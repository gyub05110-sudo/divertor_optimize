import {readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import assert from 'node:assert/strict';
import {applyQualifiedFinals,validateQualifiedProfile} from './dist/qualified-finals.mjs';
import {applyReassessment} from './dist/distribution-checks.mjs';
import {nearestHeatSegment,profileForView} from './dist/heat-probe.mjs';
const read=p=>JSON.parse(readFileSync(new URL(p,import.meta.url),'utf8'));
const hash=p=>createHash('sha256').update(readFileSync(new URL(p,import.meta.url))).digest('hex');
const original=read('./dist/data.json'),stages=read('./dist/stage-profiles.json'),dataHash=hash('./dist/data.json');
let data=structuredClone(original),saved=read('./dist/full-saved-profiles.json');
let probes=0;
for(const name of ['qualified-finals.json','distribution-finals.json']){
 const index=read('./dist/'+name),evidence={};
 for(const e of index.entries){
  assert.equal(hash('./dist/'+e.evidence.path),e.evidence.sha256);evidence[e.final.id]=read('./dist/'+e.evidence.path);
 }
 const result=applyQualifiedFinals(index,data,saved,dataHash,evidence);data=result.data;saved=result.savedProfiles;
 if(index.reassessment)data=applyReassessment(data,index.reassessment);
 for(const e of index.entries){
  assert.equal(hash('./dist/'+e.profile.profile),e.profile.sha256);
  const p=validateQualifiedProfile(read('./dist/'+e.profile.profile),e.profile);
  for(let i=0;i<p.segments.length;i+=37){const [r,z,rr,zz,q]=p.segments[i];assert.equal(nearestHeatSegment([(r+rr)/2,(z+zz)/2],p.segments,2e-5).q,q);probes++;}
 }
}
for(const s of original.series){
 const current=data.series.find(t=>t.id===s.id);
 assert.deepEqual(current.stages.map(({distribution_reassessment,...r})=>r),s.stages);
 if(!s.id.startsWith('campaign-OMP-'))assert.deepEqual(current.final,s.final);
 for(const r of current.stages)assert.equal(profileForView(data.baseline,s.source,r,'candidate',stages),stages.cases[r.id].profile);
}
const result={passed:true,probe_count:probes,original_stages_preserved:true,non_omp_finals_preserved:true,
 finals:data.series.filter(s=>s.run==='campaign'&&s.source==='OMP').map(s=>({radius:s.rmin,identity:s.final?.id??null,validation:s.final?.validation??null,v3:s.final?.distribution_reassessment?.passed??null,peak:s.final?.peak??null}))};
writeFileSync(new URL('./distribution-viewer-verification.json',import.meta.url),JSON.stringify(result,null,2));console.log(result);
