// Additive final overlay: leave original campaign, stage and leg hashes intact.
import {validateDistribution,distributionPolicy} from './distribution-checks.mjs';
const same=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
const sha=x=>typeof x==='string'&&/^[a-f0-9]{64}$/.test(x);
const near=(a,b)=>Number.isFinite(a)&&Number.isFinite(b)&&Math.abs(a-b)<=1e-10*Math.max(1,Math.abs(b));
const demand=(yes,message)=>{if(!yes)throw new Error('추가 정밀 결과 확인 실패: '+message);};

export function applyQualifiedFinals(index,data,savedProfiles,dataHash,evidence){
  demand(index?.schema==='d3d-qualified-finals.v1'&&index.complete===true&&index.data_sha256===dataHash&&sha(index.fingerprint),'색인/형상 버전');
  demand(same(index.reference,data.reference)&&Array.isArray(index.entries),'고정 입력');
  const pending=[],seen=new Set();
  for(const entry of index.entries){
    const series=data.series.find(s=>s.id===entry.series),f=entry.final,id=f?.id,m=evidence[id];
    demand(series?.run==='campaign'&&series.source==='OMP'&&[20,50].includes(series.rmin)&&!seen.has(series.id),'실행/영역/반경');
    seen.add(series.id);
    const stage=series.stages.find(s=>s.id===id);
    demand(stage&&same(stage.candidate,f.candidate)&&same(stage.wall,f.wall),'저장 후보/벽 좌표');
    demand(same(stage.active_curves,f.active_curves),'선택 곡선');
    demand(m?.status==='qualified'&&m.fingerprint===index.fingerprint&&same(m.reference,data.reference)&&m.source===series.source&&m.radius_mm===series.rmin,'검증 출처/입력');
    const v3=m.qualification_policy===distributionPolicy;
    demand((v3||m.qualification_policy==='regional-multiresolution-shared-geometry.v2')&&same(m.candidate,stage.candidate),'검증 정책/후보');
    if(v3)validateDistribution(m.validation?.distribution);
    demand(f.archived_stage===series.stages.indexOf(stage)+1&&f.archived_stage===m.origin?.stage&&f.archived_origin===stage.origin&&f.archived_origin===m.origin.origin,'원래 단계/출처');
    const v=m.validation,levels=v?.levels??[],wall=v?.wall_only_refinement;
    demand(v?.complete===true&&v.fingerprint===index.fingerprint&&v.identity===id&&same(v.candidate,stage.candidate),'완료된 검증');
    demand(v.verdict?.resolution_screen_passed===true&&v.verdict.improvement_confirmed===true,'수렴 및 개선');
    demand(same(levels.map(l=>[l.method,l.tubes]),[['cubic',320],['cubic',640],['cubic',1280],['quintic',640]]),'해상도');
    for(const l of levels){
      demand(l.geometry_audit?.passed_independent===true&&l.geometry_audit.geometry_screen?.passed===true&&l.geometry_audit.geometry_screen.policy==='full-wall-chord-1e-7.v2','독립 형상 검사');
      for(const key of ['baseline','selected'])demand(l.rows?.[key]?.audit_passed===true&&l.rows[key].feasible===true,'독립 열수송 검사');
    }
    demand(wall?.geometry_audit?.passed_independent===true&&wall.rows?.baseline?.audit_passed===true&&wall.rows?.selected?.audit_passed===true,'벽 격자 검사');
    demand(m.cpu_gpu_distribution?.passed===true&&Number.isFinite(m.cpu_gpu_peak_relative_difference)&&m.cpu_gpu_peak_relative_difference<=1e-10&&m.cpu_gpu_peak_relative_difference>=0,'CPU/GPU 비교');
    demand(m.prepublication_preservation?.inputs_unchanged===true&&m.prepublication_preservation.code_unchanged===true,'입력/코드 보존');
    demand(f.validation==='current'&&f.additional_validation===true&&f.qualified_utc===m.qualification_utc,'추가 검증 구분');
    demand(f.peak>0&&f.peak<f.baseline&&near(f.peak,m.selected.peak_MW_m2)&&near(f.baseline,m.baseline.peak_MW_m2)&&near(f.ratio,f.peak/f.baseline),'peak/감소율');
    demand(near(f.peak,levels[2].rows.selected.peak_MW_m2)&&near(f.baseline,levels[2].rows.baseline.peak_MW_m2)&&near(f.baseline,data.baseline.refined.OMP),'동일 해상도 기준');
    demand(same(f.constraints,levels[2].geometry_audit),'표시된 형상 제약 수치');
    demand(entry.profile?.identity===id&&entry.profile.source==='OMP'&&entry.profile.donors===1280&&sha(entry.profile.sha256)&&entry.profile.profile===`${v3?'distribution':'qualified'}-profiles/${id}.json`&&f.profile===entry.profile.profile,'열분포 후보/해상도/해시');
    demand(entry.profile.fingerprint===index.fingerprint&&entry.profile.npz_sha256===m.artifact_sha256?.['selected_refined.npz']&&sha(entry.profile.npz_sha256)&&near(entry.profile.peak,f.peak)&&entry.profile.qualified_utc===m.qualification_utc,'열분포 검증 연결');
    demand(entry.evidence?.path===`${v3?'distribution':'qualified'}-evidence/${id}.json`&&sha(entry.evidence.sha256),'검증 파일/해시');
    pending.push(entry);
  }
  const result=structuredClone(data),profiles={...savedProfiles,profiles:{...savedProfiles.profiles}};
  for(const entry of pending){
    result.series.find(s=>s.id===entry.series).final=structuredClone(entry.final);
    profiles.profiles[entry.final.profile]=structuredClone(entry.profile);
  }
  return {data:result,savedProfiles:profiles};
}

export function validateQualifiedProfile(profile,entry){
  demand(profile.identity===entry.identity&&profile.source===entry.source&&profile.donors===1280&&
    profile.fingerprint===entry.fingerprint&&profile.npz_sha256===entry.npz_sha256&&
    near(profile.peak,entry.peak)&&profile.qualified_utc===entry.qualified_utc,'로드된 열분포 검증 연결');
  return profile;
}
