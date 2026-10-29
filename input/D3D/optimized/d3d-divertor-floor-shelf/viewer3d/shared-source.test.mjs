import test from 'node:test';
import assert from 'node:assert/strict';
import {displaySettings,sourceDescription,stageDescription,lowerBranchOf} from './device.mjs';
test('CPD lower-outboard uses one modeled branch and explicit outer target scope',()=>{
 const bundle={device_id:'cpd2',branch_mode:'lower_outboard',source_model:'unified_radial_exponential_v1',
  branch_power_fractions:[1],source_interval_mm:[0,50],angle_floor_deg:2.,power_normalization:'P_in'};
 assert.equal(lowerBranchOf(bundle),true);
 assert.ok(sourceDescription(bundle).includes('하부 외측 단일 가지'));
 assert.ok(sourceDescription(bundle).includes('100%'));
 assert.ok(displaySettings(bundle).angleText.includes('외측 타깃'));
});
test('unified D3D and CPD use explicit same source units with device branch fractions',()=>{
 for(const device_id of ['d3d','cpd2']){
  const b={device_id,source_model:'unified_radial_exponential_v1',power_normalization:'P_in',widths:[1,2,3,5,10],
   branch_power_fractions:device_id==='cpd2'?[.5,.5]:[1],source_interval_mm:[0,50],angle_floor_deg:2.02,stages:[{},{}]};
  assert.equal(displaySettings(b).unit,'q/P_in (m⁻²)');
  assert.ok(sourceDescription(b).includes('상류 반경'));
  assert.ok(stageDescription(b,{kind:'baseline',label:'기준',score:1,verified:true,constraint_feasible:false},0).includes('제약 미달'));
 }
});
