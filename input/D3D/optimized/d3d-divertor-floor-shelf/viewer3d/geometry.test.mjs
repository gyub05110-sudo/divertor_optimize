import test from 'node:test';
import assert from 'node:assert/strict';
import {revolve, heatAt, heatSpan, heatAtTriangle} from './geometry.mjs';
test('RZ rotation preserves radius and puts Z on vertical axis',()=>{
 const p=revolve(1.7,-1.3,Math.PI/2);
 assert.ok(Math.abs(p[0])<1e-12); assert.equal(p[1],-1.3);assert.ok(Math.abs(p[2]-1.7)<1e-12);
});
test('heat interpolates only inside a connected sampled interval',()=>{
 const samples=[[.1,4],[.11,2],[.2,0]];
 assert.ok(Math.abs(heatAt(samples,.105)-3)<1e-10);
 assert.equal(heatAt(samples,.09),null);assert.equal(heatAt(samples,.15),null);
 assert.equal(heatAt(samples,.2),0);
});
test('do not bridge a gap between differently tagged surfaces',()=>{
 assert.equal(heatAt([[.1,4,'floor'],[.11,2,'fixed']],.105),null);
});
test('hover interpolates rendered heat instead of projecting a chord onto the wall',()=>{
 const vertices=[revolve(3,-2,0),revolve(3.001,-2.001,0),revolve(3.001,-2.001,.071)];
 const point=vertices[0].map((_,i)=>vertices.reduce((sum,v)=>sum+v[i]/3,0));
 const hit=heatAtTriangle(vertices,[2,8,8],point);
 assert.ok(Math.abs(hit.q-6)<1e-8);
 assert.ok(Math.abs(hit.R-(3+3.001+3.001)/3)<1e-8);
 assert.ok(Math.hypot(point[0],point[2])<hit.R);
 assert.equal(heatAtTriangle(vertices,[NaN,NaN,NaN],point).q,null);
});
test('heat spans never bridge distinct receivers or unsupported sample gaps',()=>{
 const same=[[0,2,'a'],[.01,4,'a']],different=[[0,2,'a'],[.01,4,'b']],gap=[[0,2,'a'],[.03,4,'a']];
 assert.deepEqual(heatSpan(same,0,.01),[2,4]);
 assert.equal(heatSpan(different,0,.01),null);
 assert.equal(heatSpan(gap,0,.03),null);
});
