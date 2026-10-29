import test from 'node:test';
import assert from 'node:assert/strict';
import {deviceOf,casesForDevice,chooseInitialCase,displaySettings,transportLabel,stageDescription} from './device.mjs';
import {densifyWithHeat,heatAt} from './geometry.mjs';

test('legacy D3D and explicit CPD cases are filtered independently',()=>{
 const cases=[{id:'a'},{id:'b',device_id:'cpd2'}];
 assert.equal(deviceOf(cases[0]),'d3d');
 assert.deepEqual(casesForDevice(cases,'cpd2'),[cases[1]]);
 assert.equal(chooseInitialCase(cases,new URLSearchParams('device=cpd2&case=a')).id,'b');
 assert.equal(chooseInitialCase(cases,new URLSearchParams('case=b')).id,'b');
 assert.equal(chooseInitialCase(cases,new URLSearchParams('device=missing')).id,'a');
});

test('CPD has actual width, correct total-power units and visible small gain',()=>{
 const meta=displaySettings({device_id:'cpd2',widths:[.767808717],power_normalization:'P_SOL',
  angle_floor_deg:2.02,section_bounds:[2.7,4.05,-2.85,-1.6]});
 assert.equal(meta.unit,'q/P_SOL (m⁻²)');
 assert.equal(meta.peakDigits,6);
 assert.equal(meta.deviceLabel,'CPD-2');
 assert.deepEqual(meta.sectionBounds,[2.7,4.05,-2.85,-1.6]);
});

test('heat anchors survive tessellation when footprint is narrower than5mm',()=>{
 const points=[[2,0],[2,.01]],samples=[[.00401,2,'target'],[.00402,10,'target'],[.00403,3,'target']];
 const dense=densifyWithHeat(points,samples,.005);
 assert.ok(dense.some(p=>Math.abs(p[2]-.00402)<1e-12));
 assert.ok(dense.some(p=>heatAt(samples,p[2])===10));
 assert.equal(heatAt(samples,.003),null);
});

test('D3D diffusion assumption and exact search budget are visible',()=>{
 const bundle={time_ms:1400,sigma_over_lambda_q:.5,evaluation_budget:7233,stages:[{},{}]};
 assert.equal(transportLabel(bundle),'1차원 횡방향 확산 · σ/λq=0.5 · 후보 7,233개');
 assert.ok(displaySettings(bundle).tag.includes('σ/λq=0.5'));
 assert.ok(stageDescription(bundle,{label:'기존 형상',kind:'recalculated_original',score:.8},0).includes('고정 형상의 확산 재계산'));
 assert.equal(transportLabel({time_ms:1400}),'횡방향 확산 없음');
});

test('source-infeasible final is distinct from an unverified geometry',()=>{
 const b={stages:[{}],selection_status:'source_accounting_infeasible'};
 const s={kind:'verified_final',verified:false,label:'잠정',score:1.1};
 const text=stageDescription(b,s,0);
 assert.ok(text.includes('기하/각도'));
 assert.ok(text.includes('소스 파워 제약 미달'));
 assert.ok(!text.includes('추가 재검증 미실시'));
});
