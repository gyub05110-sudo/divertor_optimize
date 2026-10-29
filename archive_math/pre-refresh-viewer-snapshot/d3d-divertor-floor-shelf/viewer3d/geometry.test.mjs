import test from 'node:test';
import assert from 'node:assert/strict';
import {revolve, heatAt} from './geometry.mjs';
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
