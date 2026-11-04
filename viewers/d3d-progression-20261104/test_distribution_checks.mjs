import {test} from 'node:test';
import assert from 'node:assert/strict';
import {validateDistribution,applyReassessment} from './dist/distribution-checks.mjs';
const limits={linf_relative:.01,l1_relative:.005,a2_relative:.005,power_relative:.005,wall_distance_m:50e-6,wall_length_difference_m:50e-6};
function evidence(){const checks=[['cubic_640','cubic_1280'],['cubic_640','quintic_640'],['cubic_1280','quintic_640'],['cubic_320','wall_refinement']].map(pair=>({pair,passed:true,failures:[],limits,...Object.fromEntries(Object.keys(limits).map(k=>[k,0])),reconstruction_power_relative:[0,0],native_cell_diagnostic:{a2_relative:0}}));return {policy:'regional-distribution-convergence.v3',passed:true,baseline:{passed:true,checks},selected:{passed:true,checks:structuredClone(checks)}};}
test('v3 distribution rejects hidden failure, missing pair and relaxed thresholds',()=>{
 assert.doesNotThrow(()=>validateDistribution(evidence()));
 for(const mutate of [e=>e.selected.checks[0].linf_relative=.02,e=>e.baseline.checks.pop(),e=>e.selected.checks[0].limits={...limits,linf_relative:1},e=>e.selected.checks[0].l1_relative=NaN,e=>e.selected.checks[0].reconstruction_power_relative=[.01,0]]){const e=evidence();mutate(e);assert.throws(()=>validateDistribution(e));}
});
test('failed reassessment preserves old final but marks it as old-policy only',()=>{
 const data={series:[{id:'campaign-OMP-20',final:{id:'x',validation:'current'},stages:[{id:'x'}]}]};
 const r=applyReassessment(data,[{series:'campaign-OMP-20',identity:'x',passed:false,max_linf_relative:.02}]);
 assert.equal(r.series[0].final.validation,'historical');assert.equal(r.series[0].final.id,'x');assert.equal(data.series[0].final.validation,'current');
});
