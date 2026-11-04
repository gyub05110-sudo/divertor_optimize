import assert from 'node:assert/strict';
import {readFileSync,writeFileSync} from 'node:fs';
import {createHash} from 'node:crypto';
import * as THREE from './dist/vendor/three.module.js';
import {nearestHeatSegment,profileForView,segmentTriangleRZ} from './dist/heat-probe.mjs';
import {applyQualifiedFinals,validateQualifiedProfile} from './dist/qualified-finals.mjs';
import {revolve,heatRGB,legLayers} from './dist/view-math.mjs';
const read=p=>JSON.parse(readFileSync(new URL(p,import.meta.url),'utf8'));
const hash=p=>createHash('sha256').update(readFileSync(new URL(p,import.meta.url))).digest('hex');
const original=read('./dist/data.json'),saved=read('./dist/full-saved-profiles.json'),index=read('./dist/qualified-finals.json');
const stages=read('./dist/stage-profiles.json'),legs=read('./dist/legs.json'),evidence={};
for(const e of index.entries){assert.equal(hash('./dist/'+e.evidence.path),e.evidence.sha256);evidence[e.final.id]=read('./dist/'+e.evidence.path);}
const applied=applyQualifiedFinals(index,original,saved,hash('./dist/data.json'),evidence);let probes=0,rays=0;
for(const e of index.entries){
  const series=applied.data.series.find(s=>s.id===e.series),row=series.final;
  const rel=profileForView(applied.data.baseline,'OMP',row,'candidate',stages),meta=applied.savedProfiles.profiles[rel];
  assert.equal(hash('./dist/'+rel),meta.sha256);
  const p=validateQualifiedProfile(read('./dist/'+rel),meta);assert.equal(p.donors,1280);
  const visibleLegs=legLayers(legs,row,'candidate');assert.equal(visibleLegs.length,1);
  assert.deepEqual(visibleLegs[0].legs,legs.cases[row.id]);
  for(let i=0;i<p.segments.length;i+=Math.max(1,Math.floor(p.segments.length/60))){
    const [ar,az,br,bz,q]=p.segments[i],hit=nearestHeatSegment([(ar+br)/2,(az+bz)/2],p.segments,2e-5);
    assert.equal(hit.q,q);assert.ok(heatRGB(q,applied.data.baseline.refined.OMP).every(Number.isFinite));probes++;
  }
  const s=p.segments.reduce((a,b)=>a[4]>b[4]?a:b),a=s.slice(0,2),b=s.slice(2,4);
  const positions=[revolve(a,0),revolve(b,0),revolve(b,.5),revolve(a,0),revolve(b,.5),revolve(a,.5)];
  const geometry=new THREE.BufferGeometry();geometry.setAttribute('position',new THREE.Float32BufferAttribute(positions.flat(),3));
  geometry.setAttribute('uv',new THREE.Float32BufferAttribute(segmentTriangleRZ(a,b),2));
  const mesh=new THREE.Mesh(geometry,new THREE.MeshBasicMaterial({side:THREE.DoubleSide}));
  const target=new THREE.Vector3().fromArray(positions[0]).add(new THREE.Vector3().fromArray(positions[1])).add(new THREE.Vector3().fromArray(positions[2])).multiplyScalar(1/3);
  const origin=target.clone().multiplyScalar(2),hit=new THREE.Raycaster(origin,target.clone().sub(origin).normalize()).intersectObject(mesh)[0];
  assert.ok(hit?.uv);assert.equal(nearestHeatSegment([hit.uv.x,hit.uv.y],p.segments,2e-5).q,s[4]);rays++;
  geometry.dispose();mesh.material.dispose();
  assert.equal(profileForView(applied.data.baseline,'OMP',row,'baseline',stages),original.baseline.profiles.OMP);
}
for(const s of original.series){
  const after=applied.data.series.find(t=>t.id===s.id);assert.deepEqual(after.stages,s.stages);
  if(!index.entries.some(e=>e.series===s.id))assert.deepEqual(after.final,s.final);
  for(const row of s.stages)assert.equal(profileForView(original.baseline,s.source,row,'candidate',stages),stages.cases[row.id].profile);
}
const report={passed:true,new_finals:index.entries.length,midpoint_q_probes:probes,three_ray_probes:rays,old_stages_and_finals_preserved:true,profile_and_leg_identity_bound:true};
writeFileSync(new URL('./qualified-probe-verification.json',import.meta.url),JSON.stringify(report,null,2));console.log(report);
