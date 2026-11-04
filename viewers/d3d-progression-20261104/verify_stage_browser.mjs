import assert from 'node:assert/strict';
import {readFileSync,writeFileSync} from 'node:fs';
import {nearestHeatSegment,formatHeatQ} from './dist/heat-probe.mjs';
const read=p=>JSON.parse(readFileSync(new URL(p,import.meta.url),'utf8'));
const evidence=read('./stage-browser-verification.json'),data=read('./dist/data.json');
const stages=read('./dist/stage-profiles.json'),saved=read('./dist/full-saved-profiles.json');
assert.equal(evidence.cases.length,12);assert.equal(evidence.console.length,0);
for(const item of evidence.cases){assert.equal(item.error,null);assert.ok(item.note.includes('표시 셀 평균'));}
for(const item of evidence.observations){
  const series=data.series.find(s=>s.id===item.series),row=item.final?series.final:series.stages[item.stage-1];
  const relative=item.final?saved.profiles[row.profile].profile:stages.cases[row.id].profile;
  const p=read('./dist/'+relative),hit=nearestHeatSegment(item.rz,p.segments,2e-5);
  assert.ok(hit);assert.equal(formatHeatQ(hit.q),item.q,JSON.stringify(item));
  item.recomputed_q=hit.q;item.coordinate_to_cell_distance_m=hit.distance;
}
evidence.numeric_observations_passed=true;
writeFileSync(new URL('./stage-browser-verification.json',import.meta.url),JSON.stringify(evidence,null,2));
console.log({cases:evidence.cases.length,actual_ui_probe_values_checked:evidence.observations.length,console_errors:0,passed:true});
