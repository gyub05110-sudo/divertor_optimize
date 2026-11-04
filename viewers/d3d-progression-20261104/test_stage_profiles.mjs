import {test} from 'node:test';
import assert from 'node:assert/strict';
import * as heat from './dist/heat-probe.mjs';
const {profileForView,validateStageIndex}=heat;
test('stage q is looked up by exact identity and source; refined q wins',()=>{
  const base={profiles:{IMP:'base.json'}};
  const stages={cases:{a:{source:'IMP',profile:'stage-a.json'},b:{source:'OMP',profile:'stage-b.json'}}};
  assert.equal(profileForView(base,'IMP',{id:'a',profile:null},'candidate',stages),'stage-a.json');
  assert.equal(profileForView(base,'IMP',{id:'a',profile:'fine.json'},'overlay',stages),'fine.json');
  assert.equal(profileForView(base,'IMP',{id:'b',profile:null},'candidate',stages),null);
  assert.equal(profileForView(base,'IMP',{id:'missing',profile:null},'candidate',stages),null);
  assert.equal(profileForView(base,'IMP',{id:'a',profile:null},'baseline',stages),'base.json');
});
test('incomplete, wrong-source and stale stage indexes fail closed',()=>{
  const data={series:[{source:'IMP',stages:[{id:'a'}]}]};
  const index={schema:'d3d-stage-profiles.v1',complete:true,data_sha256:'hash',cases:{a:{source:'IMP',donors:320,profile:'stage-profiles/a.json',sha256:'0'.repeat(64)}}};
  assert.equal(validateStageIndex(index,data,'hash'),index);
  assert.throws(()=>validateStageIndex({...index,complete:false},data,'hash'));
  assert.throws(()=>validateStageIndex(index,data,'other-hash'));
  assert.throws(()=>validateStageIndex({...index,cases:{}},data,'hash'));
  assert.throws(()=>validateStageIndex({...index,cases:{a:{...index.cases.a,source:'OMP'}}},data,'hash'));
  assert.throws(()=>validateStageIndex({...index,cases:{a:{...index.cases.a,sha256:null}}},data,'hash'));
});
test('full-wall saved index must cover every baseline and final with matching identity and hash',()=>{
  assert.equal(typeof heat.validateSavedIndex,'function');
  const data={baseline:{profiles:{IMP:'profiles/base.json'}},series:[{source:'IMP',final:{id:'a',profile:'profiles/final.json'}}]};
  const entries={'profiles/base.json':{source:'IMP',identity:null,profile:'full-saved-profiles/base.json',sha256:'0'.repeat(64)},'profiles/final.json':{source:'IMP',identity:'a',profile:'full-saved-profiles/final.json',sha256:'1'.repeat(64)}};
  const index={data_sha256:'hash',profiles:entries};
  assert.equal(heat.validateSavedIndex(index,data,'hash'),index);
  assert.throws(()=>heat.validateSavedIndex({...index,profiles:{}},data,'hash'));
  assert.throws(()=>heat.validateSavedIndex({...index,profiles:{...entries,'profiles/final.json':{...entries['profiles/final.json'],identity:'b'}}},data,'hash'));
  assert.throws(()=>heat.validateSavedIndex({...index,profiles:{...entries,'profiles/base.json':{...entries['profiles/base.json'],sha256:null}}},data,'hash'));
});
