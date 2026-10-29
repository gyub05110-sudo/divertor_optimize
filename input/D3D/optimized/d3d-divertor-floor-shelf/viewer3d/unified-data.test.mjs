import test from 'node:test';
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {deviceOf} from './device.mjs';
import {densifyWithHeat,heatSpan} from './geometry.mjs';
const read=async path=>JSON.parse(await readFile(new URL(path,import.meta.url),'utf8'));
const catalog=await read('./catalog.json');
test('current CPD selector contains only freshly optimized unit lower-outboard branches',async()=>{
 const entries=catalog.cases.filter(c=>deviceOf(c)==='cpd2');
 assert.equal(entries.length,4);
 for(const entry of entries){
  const bundle=await read(entry.bundle_url);
  assert.equal(bundle.branch_mode,'lower_outboard');
  assert.deepEqual(bundle.branch_power_fractions,[1]);
  assert.equal(bundle.angle_floor_deg,2);
  assert.equal(bundle.geometry_reoptimized,true);
  for(const model of Object.values(bundle.models))
   assert.equal(new Set(model.sol.map(path=>path.direction)).size,1);
 }
});
function independentAreaIntegral(wall,samples){
 const arc=[0];for(let i=1;i<wall.length;i++)arc.push(arc.at(-1)+Math.hypot(wall[i][0]-wall[i-1][0],wall[i][1]-wall[i-1][1]));
 const radius=s=>{let i=0;while(i<arc.length-2&&arc[i+1]<=s)i++;return wall[i][0]+(s-arc[i])/(arc[i+1]-arc[i])*(wall[i+1][0]-wall[i][0]);};
 let total=0;for(let i=1;i<samples.length;i++){
  const a=samples[i-1],b=samples[i],ds=b[0]-a[0];
  if(a[2]!==b[2]||ds<=0||ds>.02)continue;
  total+=2*Math.PI*ds/6*(a[1]*radius(a[0])+4*(a[1]+b[1])/2*radius((a[0]+b[0])/2)+b[1]*radius(b[0]));
 }
 return total;
}
for(const entry of catalog.cases.filter(c=>c.id.includes('-unified-'))){
 test(`${entry.id}: exported toroidal power independently matches every width receipt`,async()=>{
  const bundle=await read(entry.bundle_url);
  assert.deepEqual(bundle.widths,[1,2,3,5,10]);
  assert.equal(bundle.source_model,'unified_radial_exponential_v1');
  assert.equal(bundle.power_normalization,'P_in');
  assert.deepEqual(bundle.branch_power_fractions,deviceOf(bundle)==='cpd2'&&bundle.branch_mode!=='lower_outboard'?[.5,.5]:[1]);
  for(const descriptor of bundle.stages){
   const model=await read(descriptor.detail_url);
   assert.equal(model.validation.numerically_verified,true);
   if(deviceOf(bundle)==='d3d')assert.equal(model.validation.angle_scope.floor.constraint_applied,false);
   for(const [index,width] of bundle.widths.entries())
    assert.ok(Math.abs(independentAreaIntegral(model.wall,model.heat[width])-model.validation.rendered_power[index])<1e-9);
  }
 });
}
for(const entry of catalog.cases.filter(c=>deviceOf(c)==='cpd2')){
 test(`${entry.id}: every stage retains its recorded peak and first-hit support`,async()=>{
  const bundle=await read(entry.bundle_url);
  const common=bundle.source_model==='unified_radial_exponential_v1';
  assert.equal(bundle.total_power_W,common?1:100000000);
  assert.equal(bundle.angle_floor_deg,common?2.:2.02);
  assert.equal(bundle.power_normalization,common?'P_in':'P_SOL');
  assert.equal(bundle.stages.length,common?2:7);
  if(common){
   assert.deepEqual(bundle.widths,[1,2,3,5,10]);
   assert.deepEqual(bundle.branch_power_fractions,bundle.branch_mode==='lower_outboard'?[1]:[.5,.5]);
  }
  for(const descriptor of bundle.stages){
   const model=await read(descriptor.detail_url),key=bundle.widths[0];
   if(common){
    assert.ok(model.validation.sample_count>=1280);
    assert.equal(model.validation.numerically_verified,true);
    assert.equal(model.validation.exported_integral_matches_audit,true);
    assert.ok(Math.max(...model.validation.uncovered_source_power)<=.001);
    for(const group of Object.values(model.validation.receiver_groups))
     assert.ok(Math.max(...group.integral_relative_residual.map(Math.abs))<=.01);
    assert.equal(model.validation.status,descriptor.constraint_feasible?'complete_feasible':'constraint_not_satisfied');
   }else{
    assert.equal(model.validation.sample_count,138);
    assert.equal(model.validation.status,'complete_feasible');
    assert.equal(model.validation.continuous_peak_certified,false);
   }
   if(descriptor.kind!=='baseline'||!common)assert.ok(model.validation.minimum_target_angle_deg>=bundle.angle_floor_deg);
   const samples=model.heat[key],dense=densifyWithHeat(model.wall,samples);
   const recorded=Math.max(...samples.map(s=>s[1]));
   let rendered=0;
   for(let i=1;i<dense.length;i++){
    const span=heatSpan(samples,dense[i-1][2],dense[i][2]);
    if(span)rendered=Math.max(rendered,...span);
   }
   assert.ok(Math.abs(recorded-model.summary[key].all_wall_peak_q_over_P)<1e-10);
   assert.ok(Math.abs(rendered-recorded)<1e-10,'a narrow peak must survive mesh tessellation');
   assert.deepEqual([...new Set(model.targets.map(t=>t.role))].sort(),bundle.branch_mode==='lower_outboard'?['outer_target']:['inner_target','outer_target']);
   if(!common)assert.equal(model.display_sol_omitted,0);
   assert.ok(model.sol.length>0);
  }
 });
}
