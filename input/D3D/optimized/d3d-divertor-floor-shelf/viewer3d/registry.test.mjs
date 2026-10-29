import test from 'node:test';import assert from 'node:assert/strict';
import {casesForFile} from './registry.mjs';
test('match source SHA rather than same filename; preserve distinct time slices',()=>{
 const cases=[{id:'a',shot_sha256:'aaa',time_ms:1400},{id:'b',shot_sha256:'aaa',time_ms:1500},{id:'c',shot_sha256:'bbb',time_ms:1400}];
 assert.deepEqual(casesForFile(cases,'aaa').map(x=>x.id),['a','b']);
 assert.deepEqual(casesForFile(cases,'unknown'),[]);
});
