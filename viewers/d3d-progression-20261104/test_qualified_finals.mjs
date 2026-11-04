import {test} from 'node:test';
import assert from 'node:assert/strict';
import {applyQualifiedFinals,validateQualifiedProfile} from './dist/qualified-finals.mjs';

function fixture(){
  const id='a'.repeat(64),fp='f'.repeat(64),candidate={source:'OMP',radius_mm:20,vector:[0]};
  const stage={id,candidate,wall:[[1,-1],[2,-1]],active_curves:{floor:[[1,-1],[2,-1]]},origin:'archive/log.jsonl',peak:1.4,baseline:2.7};
  const data={reference:{power_W:1e6},baseline:{refined:{OMP:2.7}},series:[{id:'campaign-OMP-20',run:'campaign',source:'OMP',rmin:20,stages:[stage],final:null}]};
  const levels=[['cubic',320],['cubic',640],['cubic',1280],['quintic',640]].map(([method,tubes])=>({method,tubes,
    rows:{baseline:{audit_passed:true,feasible:true,peak_MW_m2:2.7},selected:{audit_passed:true,feasible:true,peak_MW_m2:1.8}},
    geometry_audit:{passed_independent:true,geometry_screen:{passed:true,policy:'full-wall-chord-1e-7.v2'}}}));
  const evidence={status:'qualified',fingerprint:fp,reference:data.reference,source:'OMP',radius_mm:20,candidate,
    origin:{stage:1,origin:stage.origin},artifact_sha256:{'selected_refined.npz':fp},
    qualification_policy:'regional-multiresolution-shared-geometry.v2',qualification_utc:'2026-11-04T00:00:00+00:00',
    selected:{peak_MW_m2:1.8},baseline:{peak_MW_m2:2.7},cpu_gpu_peak_relative_difference:0,
    cpu_gpu_distribution:{passed:true},prepublication_preservation:{inputs_unchanged:true,code_unchanged:true},
    validation:{complete:true,fingerprint:fp,identity:id,candidate,levels,verdict:{resolution_screen_passed:true,improvement_confirmed:true},
      wall_only_refinement:{geometry_audit:{passed_independent:true},rows:{baseline:{audit_passed:true},selected:{audit_passed:true}}}}};
  const final={...stage,peak:1.8,baseline:2.7,ratio:1.8/2.7,validation:'current',additional_validation:true,
    constraints:structuredClone(evidence.validation.levels[2].geometry_audit),
    archived_stage:1,archived_origin:stage.origin,
    profile:`qualified-profiles/${id}.json`,qualified_utc:evidence.qualification_utc};
  const entry={series:'campaign-OMP-20',final,evidence:{path:`qualified-evidence/${id}.json`,sha256:fp},
    profile:{identity:id,source:'OMP',donors:1280,profile:final.profile,sha256:fp,fingerprint:fp,npz_sha256:fp,peak:1.8,qualified_utc:evidence.qualification_utc}};
  const index={schema:'d3d-qualified-finals.v1',complete:true,data_sha256:'hash',fingerprint:fp,reference:data.reference,entries:[entry]};
  return {data,index,saved:{profiles:{}},evidence:{[id]:evidence}};
}

test('adds separately qualified final without mutating original stages or saved indices',()=>{
  const f=fixture(),before=structuredClone(f.data);const r=applyQualifiedFinals(f.index,f.data,f.saved,'hash',f.evidence);
  assert.equal(r.data.series[0].final.peak,1.8);assert.deepEqual(f.data,before);assert.deepEqual(f.saved,{profiles:{}});
  assert.deepEqual(r.data.series[0].stages,before.series[0].stages);
  assert.equal(r.savedProfiles.profiles[r.data.series[0].final.profile].donors,1280);
});
test('rejects non-improvement, missing independent audit, wrong wall, identity and hash',()=>{
  for(const mutate of [f=>f.index.data_sha256='stale',
    f=>f.evidence['a'.repeat(64)].validation.verdict.improvement_confirmed=false,
    f=>f.evidence['a'.repeat(64)].validation.levels[2].geometry_audit.passed_independent=false,
    f=>f.evidence['a'.repeat(64)].cpu_gpu_distribution.passed=false,
    f=>f.index.entries[0].final.wall=[[0,0],[1,1]],
    f=>f.index.entries[0].final.active_curves={floor:[[0,0],[1,1]]},
    f=>f.index.entries[0].final.archived_stage=19,
    f=>f.index.entries[0].final.constraints={},
    f=>f.index.entries[0].profile.npz_sha256='b'.repeat(64),
    f=>f.index.entries[0].final.peak=3,
    f=>f.index.entries[0].profile.identity='b'.repeat(64),
    f=>f.index.entries[0].profile.sha256=null]){
    const f=fixture();mutate(f);assert.throws(()=>applyQualifiedFinals(f.index,f.data,f.saved,'hash',f.evidence));
    assert.equal(f.data.series[0].final,null);
  }
});
test('loaded finest profile must belong to the same committed qualification',()=>{
  const f=fixture(),entry=f.index.entries[0].profile,p={...entry};
  assert.equal(validateQualifiedProfile(p,entry),p);
  for(const key of ['fingerprint','npz_sha256','peak','qualified_utc'])assert.throws(()=>validateQualifiedProfile({...p,[key]:null},entry));
});
