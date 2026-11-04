export const distributionPolicy='regional-distribution-convergence.v3';
const limits={linf_relative:.01,l1_relative:.005,a2_relative:.005,power_relative:.005,wall_distance_m:50e-6,wall_length_difference_m:50e-6};
const demand=(x)=>{if(!x)throw new Error('분포 수렴 v3 근거 누락 또는 기준 미통과');};
export function validateDistribution(d){
 demand(d?.policy===distributionPolicy&&d.passed===true);
 const pairs=[['cubic_640','cubic_1280'],['cubic_640','quintic_640'],['cubic_1280','quintic_640'],['cubic_320','wall_refinement']];
 for(const name of ['baseline','selected']){
  const s=d[name];demand(s?.passed===true&&JSON.stringify(s.checks?.map(c=>c.pair))===JSON.stringify(pairs));
  for(const c of s.checks){
   demand(c.passed===true&&c.failures?.length===0);
   for(const [k,t] of Object.entries(limits))demand(c.limits?.[k]===t&&Number.isFinite(c[k])&&c[k]>=0&&c[k]<=t);
   demand(c.reconstruction_power_relative?.length===2&&c.reconstruction_power_relative.every(v=>Number.isFinite(v)&&v>=0&&v<=.001));
   demand(Number.isFinite(c.native_cell_diagnostic?.a2_relative)&&c.native_cell_diagnostic.a2_relative>=0&&c.native_cell_diagnostic.a2_relative<=.005);
  }
 }
}
export function applyReassessment(data,rows){
 const result=structuredClone(data),seen=new Set();
 for(const r of rows){
  demand(!seen.has(r.identity)&&typeof r.passed==='boolean'&&Number.isFinite(r.max_linf_relative)&&r.max_linf_relative>=0);seen.add(r.identity);
  const series=result.series.find(s=>s.id===r.series);demand(series&&series.stages.some(s=>s.id===r.identity));
  for(const row of [...series.stages,series.final].filter(Boolean))if(row.id===r.identity){
   row.distribution_reassessment=structuredClone(r);
   if(row===series.final&&!r.passed)row.validation='historical';
  }
 }
 return result;
}
