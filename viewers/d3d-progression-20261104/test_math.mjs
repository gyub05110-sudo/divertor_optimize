import {test} from 'node:test';
import assert from 'node:assert/strict';
import {revolve,segmentsForWall,selectedRecord,heatRGB,minRadius,maxTangentError,legLayers} from './dist/view-math.mjs';
test('revolution keeps R and Z and closes after 2 pi',()=>{
  for(const angle of [0,.7,Math.PI,2*Math.PI]){
    const [x,y,z]=revolve([1.4,-1.3],angle);
    assert.ok(Math.abs(Math.hypot(x,z)-1.4)<1e-12);assert.equal(y,-1.3);
  }
});
test('cut geometry never creates a bridge between nonadjacent points',()=>{
  const wall=[[1,-1.2],[1,-.5],[2,-.5],[2,-1.2],[2.1,-1.3]];
  assert.deepEqual(segmentsForWall(wall,true),[[wall[3],wall[4]]]);
});
test('empty, final and baseline selections cannot borrow a previous series record',()=>{
  assert.equal(selectedRecord({stages:[],final:null},2,false),null);
  assert.equal(selectedRecord({stages:[{id:'a'}],final:{id:'f'}},1,true).id,'f');
  assert.equal(selectedRecord({stages:[{id:'a'}]},0,false),null);
});
test('color scale clamps instead of changing with each candidate',()=>{
  assert.deepEqual(heatRGB(9,3),heatRGB(3,3));
  assert.deepEqual(heatRGB(-1,3),heatRGB(0,3));
  assert.ok(heatRGB(1,3).every(Number.isFinite));
});
test('IMP array and OMP scalar tangent errors are both supported',()=>{
  assert.equal(maxTangentError({tangent_angle_errors_deg:[0,1e-14,0]}),1e-14);
  assert.equal(maxTangentError({tangent_angle_error_deg:3e-15}),3e-15);
});
test('straight surfaces represented by null do not imply a zero radius',()=>{
  assert.equal(minRadius({actual_minimum_radius_m:{entry:null,floor:.0201,shelf:null}}),.0201);
  assert.equal(minRadius({actual_minimum_radius_m:{shelf:null}}),Infinity);
});
test('leg layers follow the displayed wall and never borrow an unknown candidate strike',()=>{
  const legs={cases:{baseline:{inner:'base'},abc:{inner:'selected'}}},record={id:'abc'};
  assert.deepEqual(legLayers(legs,record,'overlay').map(x=>[x.legs.inner,x.ghost]),[['base',true],['selected',false]]);
  assert.deepEqual(legLayers(legs,record,'baseline').map(x=>x.legs.inner),['base']);
  assert.deepEqual(legLayers(legs,record,'candidate').map(x=>x.legs.inner),['selected']);
  assert.deepEqual(legLayers(legs,{id:'unknown'},'candidate'),[]);
  assert.deepEqual(legLayers(legs,null,'candidate').map(x=>x.legs.inner),['base']);
  assert.deepEqual(legLayers(null,record,'overlay'),[]);
});
