import {test} from 'node:test';
import assert from 'node:assert/strict';
import * as probe from './dist/heat-probe.mjs';
import * as THREE from './dist/vendor/three.module.js';
import {revolve} from './dist/view-math.mjs';
import {readFileSync,readdirSync} from 'node:fs';

test('nearest segment projects to the segment, not the closest endpoint',()=>{
  assert.equal(typeof probe.nearestHeatSegment,'function');
  const hit=probe.nearestHeatSegment([1.5,-1.3],[[1,-1.3,2,-1.3,2.4]],.001);
  assert.deepEqual(hit.point,[1.5,-1.3]);assert.equal(hit.q,2.4);assert.equal(hit.t,.5);
});
test('zero is a stored value and uncovered positions stay missing',()=>{
  assert.equal(typeof probe.nearestHeatSegment,'function');
  assert.equal(probe.nearestHeatSegment([1.2,-1.3],[[1,-1.3,2,-1.3,0]],1e-5).q,0);
  assert.equal(probe.nearestHeatSegment([1.2,-1.2],[[1,-1.3,2,-1.3,2]],1e-5),null);
  assert.equal(probe.nearestHeatSegment([1.2,-1.3],[],1),null);
});
test('profile selection follows visible geometry without fabricating a stage profile',()=>{
  assert.equal(typeof probe.profileForView,'function');
  const base={profiles:{IMP:'base.json'}};
  assert.equal(probe.profileForView(base,'IMP',{profile:'selected.json'},'candidate'),'selected.json');
  assert.equal(probe.profileForView(base,'IMP',{profile:null},'overlay'),null);
  assert.equal(probe.profileForView(base,'IMP',{profile:'selected.json'},'baseline'),'base.json');
  assert.equal(probe.profileForView(base,'IMP',null,'candidate'),'base.json');
});
test('faceted toroidal triangles carry exact poloidal R/Z for picking',()=>{
  assert.equal(typeof probe.segmentTriangleRZ,'function');
  const uv=probe.segmentTriangleRZ([1.2,-1.3],[1.4,-1.4]);
  assert.deepEqual(uv,[1.2,-1.3,1.4,-1.4,1.4,-1.4,1.2,-1.3,1.4,-1.4,1.2,-1.3]);
});
test('small nonzero q is not rounded to a misleading zero',()=>{
  assert.equal(typeof probe.formatHeatQ,'function');
  assert.equal(probe.formatHeatQ(0),'0.0000');
  assert.equal(probe.formatHeatQ(.0000008),'8.000e-7');
  assert.equal(probe.formatHeatQ(null),'기록 없음');
});
test('a real ray intersection recovers R despite the toroidal chord radius bias',()=>{
  const a=[1.4,-1.3],b=[1.4,-1.4];
  const p=[revolve(a,0),revolve(b,0),revolve(b,.5),revolve(a,0),revolve(b,.5),revolve(a,.5)];
  const geometry=new THREE.BufferGeometry();
  geometry.setAttribute('position',new THREE.Float32BufferAttribute(p.flat(),3));
  geometry.setAttribute('uv',new THREE.Float32BufferAttribute(probe.segmentTriangleRZ(a,b),2));
  const mesh=new THREE.Mesh(geometry,new THREE.MeshBasicMaterial({side:THREE.DoubleSide}));
  const target=new THREE.Vector3().fromArray(p[0]).add(new THREE.Vector3().fromArray(p[1])).add(new THREE.Vector3().fromArray(p[2])).multiplyScalar(1/3);
  const origin=target.clone().multiplyScalar(2),ray=new THREE.Raycaster(origin,target.clone().sub(origin).normalize());
  const hit=ray.intersectObject(mesh)[0];assert.ok(hit);
  assert.ok(Math.abs(hit.uv.x-1.4)<1e-6);
  assert.ok(Math.abs(Math.hypot(hit.point.x,hit.point.z)-1.4)>.01);
  geometry.dispose();mesh.material.dispose();
});
test('all saved heat profiles return their own cell q at each segment midpoint',()=>{
  for(const name of readdirSync(new URL('./dist/profiles/',import.meta.url))){
    if(!name.endsWith('.json'))continue;
    const profile=JSON.parse(readFileSync(new URL('./dist/profiles/'+name,import.meta.url),'utf8'));
    for(const [r,z,r2,z2,q]of profile.segments){
      const hit=probe.nearestHeatSegment([(r+r2)/2,(z+z2)/2],profile.segments,2e-5);
      assert.ok(hit, name);assert.equal(hit.q,q,name);
    }
  }
});
