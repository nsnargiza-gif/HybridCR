// HybridCR core self-tests — extracts the numerical core from the shipped HTML
// and checks closed-form reference values and edge cases. No duplication of source.
const fs=require('fs'), vm=require('vm'), path=require('path');
const HTML=path.join(__dirname,'HybridCR.html');
const html=fs.readFileSync(HTML,'utf8');
const js=[...html.matchAll(/<script[^>]*>([\s\S]*?)<\/script>/g)].map(m=>m[1]).join('\n');
function grab(src,name){
  const m=src.match(new RegExp('function\\s+'+name+'\\s*\\([^)]*\\)\\s*\\{'));
  if(!m) throw new Error('function not found: '+name);
  let depth=0,j=m.index+m[0].length-1;
  for(;j<src.length;j++){if(src[j]==='{')depth++;else if(src[j]==='}'){depth--;if(depth===0)break;}}
  return src.slice(m.index,j+1);
}
const NAMES=['erf','logGamma','gammaP','chi2cdf','chi2ppf','chi2pdf','ncx2cdf','ncx2pdf',
 'mulberry32','mleExp','momStartWeib','neg','det4','inv4','solve4','newtonWeib',
 'weibOperInfo','operVar','opBlocks6','eig2','infoQ6','lambda6','claytonBias',
 'aalenJohansen','coxSnellResiduals','simExp','simWeib','parseCSV','expSample'];
const code=NAMES.map(n=>grab(js,n)).join('\n');
const ctx={Math,console,isFinite,parseFloat,parseInt,isNaN};
vm.createContext(ctx); vm.runInContext(code,ctx);

let pass=0,fail=0; const close=(a,b,t=1e-3)=>Math.abs(a-b)<=t;
function check(name,cond,got){ if(cond){pass++; console.log('  PASS',name,got!==undefined?('= '+got):'');}
  else{fail++; console.log('  FAIL',name,'got',got);} }

console.log('== closed-form reference values ==');
check('chi2ppf(0.95,1)=3.8415', close(ctx.chi2ppf(0.95,1),3.8415,2e-3), ctx.chi2ppf(0.95,1).toFixed(4));
check('chi2ppf(0.95,4)=9.4877', close(ctx.chi2ppf(0.95,4),9.4877,3e-3), ctx.chi2ppf(0.95,4).toFixed(4));
check('ncx2cdf central==chi2cdf', close(ctx.ncx2cdf(9.4877,4,0),ctx.chi2cdf(9.4877,4),1e-6), ctx.ncx2cdf(9.4877,4,0).toFixed(5));
check('operVar(1,0.5,0.75)=[2,1]', close(ctx.operVar(1,0.5,0.75)[0],2.0)&&close(ctx.operVar(1,0.5,0.75)[1],1.0), JSON.stringify(ctx.operVar(1,0.5,0.75)));
check('erf(0)=0', close(ctx.erf(0),0), ctx.erf(0));
check('erf(1)=0.8427', close(ctx.erf(1),0.8427,1e-3), ctx.erf(1).toFixed(4));

console.log('== exact non-centrality (Module 6) ==');
const L=ctx.lambda6(200,0.70,0.15,4);
check('lambda6 balanced ~1.63', close(L.lam,1.63,0.08), L.lam.toFixed(3));
check('lambda6 lmin<=lam<=lmax', L.lmin<=L.lam+1e-9 && L.lam<=L.lmax+1e-9, `[${L.lmin.toFixed(2)},${L.lam.toFixed(2)},${L.lmax.toFixed(2)}]`);
check('eig2 identity -> [1,1]', close(ctx.eig2([[1,0],[0,1]])[0],1)&&close(ctx.eig2([[1,0],[0,1]])[1],1), JSON.stringify(ctx.eig2([[1,0],[0,1]])));

console.log('== operational information (Weibull) ==');
const oi=ctx.weibOperInfo(1.0,1.3,0.5,0.9,0.75);
check('weibOperInfo blocks finite & PD', isFinite(oi[0].I[0][0])&&oi[0].I[0][0]>0&&oi[1].I[0][0]>0, oi[0].I[0][0].toFixed(4));
check('op info grows as beta->1 (more info)', ctx.weibOperInfo(1,1.3,0.5,0.9,0.95)[0].I[0][0] > oi[0].I[0][0], 'ok');

console.log('== MLE recovers truth (seeded sim) ==');
const rng=ctx.mulberry32(42); const sim=ctx.simExp(5000,1.0,0.5,0.75,rng);
const f=ctx.mleExp(sim);
check('mleExp a1 ~1.0 (n=5000)', close(f.est[0],1.0,0.08), f.est[0].toFixed(3));
check('mleExp a2 ~0.5 (n=5000)', close(f.est[1],0.5,0.06), f.est[1].toFixed(3));

console.log('== edge cases ==');
// all censored: beta tiny -> few/no events; build manual all-censored sim
const allC={Z:[1,2,3,4,5],cause:[1,2,1,2,1],obs:[false,false,false,false,false],T0:5};
const fe=ctx.mleExp(allC);
check('all-censored: D1=D2=0, est=[0,0]', fe.D1===0&&fe.D2===0&&fe.est[0]===0, JSON.stringify(fe.est));
// parseCSV malformed + valid mix
const csv='time,cause,status\n10,1,1\n-5,1,1\nabc,2,0\n20,2,0\n30,3,1';
const pc=ctx.parseCSV(csv);
check('parseCSV keeps only valid rows (2 of 5)', pc.Z.length===2, 'rows='+pc.Z.length);
check('parseCSV drops invalid cause/t<=0/non-numeric (3 dropped)', pc.dropped===3, 'dropped='+pc.dropped);
check('parseCSV cause strictly in {1,2}', pc.cause.every(c=>c===1||c===2), JSON.stringify(pc.cause));
check('parseCSV status 0->censored', pc.obs[1]===false, 'obs='+JSON.stringify(pc.obs));
// claytonBias deterministic with seed, tau=0 -> ~unbiased
const b0=ctx.claytonBias(1.0,0.5,0.75,0.0,2000,200,7);
check('claytonBias tau=0 ~ unbiased (|bias|<0.03)', Math.abs(b0.bias)<0.03, b0.bias.toFixed(4));
const bT=ctx.claytonBias(1.0,0.5,0.75,0.4,2000,200,7);
check('claytonBias increases with tau', Math.abs(bT.bias)>Math.abs(b0.bias), `tau0=${b0.bias.toFixed(4)} tau.4=${bT.bias.toFixed(4)}`);

console.log(`\nRESULT: ${pass} passed, ${fail} failed`);
// write reference values
const ref={chi2ppf_0_95:{df1:ctx.chi2ppf(0.95,1),df4:ctx.chi2ppf(0.95,4)},
 operVar_1_0_5_0_75:ctx.operVar(1,0.5,0.75),
 lambda6_200_070_015_4:{lam:L.lam,lmin:L.lmin,lmax:L.lmax},
 weibOperInfo_block1_00:oi[0].I[0][0],
 mleExp_n5000:f.est};
fs.writeFileSync(path.join(__dirname,'reference_values.json'),JSON.stringify(ref,null,2));
process.exit(fail===0?0:1);
