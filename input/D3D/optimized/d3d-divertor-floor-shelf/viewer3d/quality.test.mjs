import test from 'node:test';
import assert from 'node:assert/strict';
import {statusLabel,patchAt,compatibleShape} from './quality.mjs';
test('missing coverage states remain distinct from zero',()=>{
 assert.notEqual(statusLabel('blocked_first_hit'),statusLabel('trace_unresolved'));
 assert.ok(statusLabel('outside_source_scope').includes('범위'));
});
test('only exact compatible shapes accept a numerical audit',()=>{
 const a=[[1,0],[2,0]];
 assert.equal(compatibleShape(a,a),true);
 assert.equal(compatibleShape(a,[[1,0],[2,.001]]),false);
});
test('probe reads a certified patch and retains null for unresolved heat',()=>{
 const p={a:[1,0],b:[2,0],status:'trace_unresolved',q:{3:null}};
 assert.equal(patchAt([p],1.5,0).q[3],null);
 assert.equal(patchAt([p],1.5,.1),null);
});
