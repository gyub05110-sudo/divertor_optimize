import assert from 'node:assert/strict';
import {readFileSync,writeFileSync} from 'node:fs';
import * as THREE from './dist/vendor/three.module.js';
import {nearestHeatSegment,profileForView,validateStageIndex,segmentTriangleRZ} from './dist/heat-probe.mjs';
import {selectedRecord,revolve,heatRGB} from './dist/view-math.mjs';
import {createHash} from 'node:crypto';

const read=p=>JSON.parse(readFileSync(new URL(p,import.meta.url),'utf8'));
const data=read('./dist/data.json'),index=read('./dist/stage-profiles.json');
validateStageIndex(index,data,createHash('sha256').update(readFileSync(new URL('./dist/data.json',import.meta.url))).digest('hex'));
let routed=0,probes=0,rays=0;
for(const series of data.series){
  for(let n=0;n<=series.stages.length;n++){
    const row=selectedRecord(series,n,false),rel=profileForView(data.baseline,series.source,row,'candidate',index);
    assert.ok(rel);
    if(!row)continue;
    assert.equal(rel,index.cases[row.id].profile);routed++;
    const p=read('./dist/'+rel);assert.equal(p.identity,row.id);assert.equal(p.source,series.source);assert.equal(p.donors,320);
    for(let j=0;j<p.segments.length;j+=Math.max(1,Math.floor(p.segments.length/48))){
      const [r,z,r2,z2,q]=p.segments[j];
      const hit=nearestHeatSegment([(r+r2)/2,(z+z2)/2],p.segments,2e-5);
      assert.ok(hit);assert.equal(hit.q,q);assert.ok(heatRGB(q,data.baseline.refined[series.source]).every(Number.isFinite));probes++;
    }
    const s=p.segments.reduce((a,b)=>a[4]>b[4]?a:b),a=s.slice(0,2),b=s.slice(2,4);
    const positions=[revolve(a,0),revolve(b,0),revolve(b,.5),revolve(a,0),revolve(b,.5),revolve(a,.5)];
    const g=new THREE.BufferGeometry();
    g.setAttribute('position',new THREE.Float32BufferAttribute(positions.flat(),3));
    g.setAttribute('uv',new THREE.Float32BufferAttribute(segmentTriangleRZ(a,b),2));
    const mesh=new THREE.Mesh(g,new THREE.MeshBasicMaterial({side:THREE.DoubleSide}));
    const target=new THREE.Vector3().fromArray(positions[0]).add(new THREE.Vector3().fromArray(positions[1])).add(new THREE.Vector3().fromArray(positions[2])).multiplyScalar(1/3);
    const origin=target.clone().multiplyScalar(2),ray=new THREE.Raycaster(origin,target.clone().sub(origin).normalize());
    const hit=ray.intersectObject(mesh)[0];assert.ok(hit?.uv);
    const qhit=nearestHeatSegment([hit.uv.x,hit.uv.y],p.segments,2e-5);assert.ok(qhit);assert.equal(qhit.q,s[4]);rays++;
    g.dispose();mesh.material.dispose();
    assert.equal(profileForView(data.baseline,series.source,row,'baseline',index),data.baseline.profiles[series.source]);
  }
  if(series.final)assert.equal(profileForView(data.baseline,series.source,selectedRecord(series,0,true),'overlay',index),series.final.profile);
}
const report={passed:true,routed_stage_occurrences:routed,midpoint_probes:probes,three_ray_probes:rays,all_sources_radii_runs:true};
writeFileSync(new URL('./stage-probe-verification.json',import.meta.url),JSON.stringify(report,null,2));
console.log(report);
