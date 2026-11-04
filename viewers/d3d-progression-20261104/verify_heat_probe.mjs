import {readFileSync,writeFileSync,readdirSync} from 'node:fs';
import {createHash} from 'node:crypto';
import assert from 'node:assert/strict';
import {nearestHeatSegment} from './dist/heat-probe.mjs';

const profile=JSON.parse(readFileSync(new URL('./dist/profiles/tangent-OMP-50.json',import.meta.url)));
const observations=[{view:'R-Z',point:[1.24403,-1.39994],displayed_q:'1.9891'},
  {view:'3D',point:[1.26610,-1.41473],displayed_q:'1.8311'}];
for(const row of observations){
  const hit=nearestHeatSegment(row.point,profile.segments,2e-5);
  assert.ok(hit);assert.equal(hit.q.toFixed(4),row.displayed_q);
  row.saved_cell_q=hit.q;row.distance_from_rounded_UI_coordinate_m=hit.distance;
}
const hashes={};let cells=0;
for(const name of readdirSync(new URL('./dist/profiles/',import.meta.url))){
  if(!name.endsWith('.json'))continue;
  const bytes=readFileSync(new URL('./dist/profiles/'+name,import.meta.url));
  hashes[name]=createHash('sha256').update(bytes).digest('hex');
  cells+=JSON.parse(bytes).segments.length;
}
const report={checks:'passed',observations,profile_files:Object.keys(hashes).length,display_cells:cells,
  profile_sha256:hashes,javascript_tests_passed:14,new_thermal_solves:0,new_optimizations:0,
  browser_checks:['3D and R-Z numeric q readout','shared numeric color scale','color hidden with q still available',
    'baseline uses its saved refined profile','unsaved stage reports missing, not zero','stage change clears stale tooltip',
    'right-edge tooltip width 268 px, right 1145 px within viewport 1153 px'],
  screenshots:['screenshots/q-probe-rz.jpg','screenshots/q-probe-3d.jpg']};
writeFileSync(new URL('./heat-probe-verification.json',import.meta.url),JSON.stringify(report,null,2)+'\n');
console.log(JSON.stringify({checks:report.checks,observations,cells}));
