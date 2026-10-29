export const deviceOf=entry=>entry.device_id||'d3d';
export const lowerBranchOf=bundle=>deviceOf(bundle)!=='cpd2'||bundle.branch_mode==='lower_outboard';
export const casesForDevice=(cases,device)=>cases.filter(c=>deviceOf(c)===device);
export function transportLabel(bundle){
 if(bundle.sigma_over_lambda_q==null)return '횡방향 확산 없음';
 return `1차원 횡방향 확산 · σ/λq=${bundle.sigma_over_lambda_q}`+
  (bundle.evaluation_budget?` · 후보 ${bundle.evaluation_budget.toLocaleString('en-US')}개`:'');
}
export function chooseInitialCase(cases,params){
 const requested=cases.find(c=>c.id===params.get('case'));
 const device=params.get('device')||deviceOf(requested||cases[0]);
 const matches=casesForDevice(cases,device);
 return (requested&&deviceOf(requested)===device?requested:matches[0])||cases[0];
}
export function displaySettings(bundle){
 const cpd=deviceOf(bundle)==='cpd2';
 return {deviceLabel:cpd?'CPD-2':'DIII-D',unit:bundle.power_normalization==='P_in'?'q/P_in (m⁻²)':cpd?'q/P_SOL (m⁻²)':'q/P_outer (m⁻²)',
  peakDigits:cpd?6:3,sectionBounds:bundle.section_bounds||[1.07,1.81,-1.46,-1.17],
  angleText:cpd?`${bundle.angle_floor_deg}° 하한은 지정 ${lowerBranchOf(bundle)?'외측 ':''}타깃의 양의 직접 소스 파워 첫 충돌점에만 적용합니다.`:
   '1.5° 하한은 shelf에만 적용합니다. 바닥 각도는 진단값입니다.',
  tag:cpd?'CPD-2 고정 평형 · 축대칭 모델':`EFIT ${bundle.time_ms} ms · 축대칭 모델`+(bundle.sigma_over_lambda_q!=null?` · σ/λq=${bundle.sigma_over_lambda_q}`:''),
  fileAccept:cpd?'.geqdsk,.eqdsk,.gfile,.g':'.parquet'};
}
export function cameraSettings(bundle){
 const cpd=deviceOf(bundle)==='cpd2';
 if(!cpd)return {radius:2.2,z:-.2,focus:[.95,-1.30,0],section:[1.25,-.35,0],
  overview:[5,3.2,-6],focusPosition:[3.2,0,-3],sectionPosition:[1.25,-.35,-6],top:[.001,7,.001],maxDistance:18};
 const wall=bundle.models.baseline.wall,radius=Math.max(...wall.map(p=>p[0])),z=(Math.min(...wall.map(p=>p[1]))+Math.max(...wall.map(p=>p[1])))/2;
 const [r0,r1,z0,z1]=displaySettings(bundle).sectionBounds,focus=[(r0+r1)*.38,(z0+z1)/2,0];
 return {radius,z,focus,section:[radius*.55,z,0],overview:[radius*2.2,z+radius*1.3,-radius*2.6],
  focusPosition:[focus[0]+radius*.8,focus[1]+radius*.5,-radius*1.1],
  sectionPosition:[radius*.55,z,-radius*3],top:[.001,z+radius*3.4,.001],maxDistance:radius*8};
}
export function formatChange(change){return Math.abs(change).toFixed(change!==0&&Math.abs(change)<.1?6:1)+'%';}
export function stageDescription(bundle,stage,index){
 if(bundle.source_model==='unified_radial_exponential_v1')return `${index+1} / ${bundle.stages.length} · ${stage.label} · 다섯 폭 최대 비율 ${stage.score.toFixed(6)} · ${stage.verified?'피크·적분 세분화 검증 통과':'수치검증 미통과'} · ${stage.constraint_feasible?'타깃 제약 충족':'기준 타깃 제약 미달'}`;
 if(deviceOf(bundle)==='d3d'&&stage.kind==='recalculated_original')return `${index+1} / ${bundle.stages.length} · ${stage.label} · 고정 형상의 확산 재계산 · 다섯 폭 worst ratio ${stage.score.toFixed(5)}`;
 if(deviceOf(bundle)!=='cpd2')return `${index+1} / ${bundle.stages.length} · ${stage.label} · 검색 당시 다섯 폭 worst ratio ${stage.score.toFixed(5)}`+
  (stage.verified?' · 세분화/보간 재검증 통과':stage.kind==='verified_final'&&bundle.selection_status==='source_accounting_infeasible'?' · 기하/각도 세분화 검증 통과 · 소스 파워 제약 미달':stage.refinement_rejected?' · 기본 해상도 후보: 세분화 재검증 탈락':stage.kind==='baseline'?' · 기준':' · 기본 해상도 후보: 추가 재검증 미실시')+
  (stage.kind==='baseline'&&stage.constraint_feasible===false?' · 기준 shelf의 1.5° 조건 미달':'');
 return `${index+1} / ${bundle.stages.length} · ${stage.label} · 실행 시작 대비 q 비율 ${stage.score.toFixed(8)} · `+
  (stage.kind==='baseline'?'이미 최적화된 실행 시작 형상':'정밀 유한 표본 제약 충족')+' · 연속 최대값 인증 미실시';
}
export function sourceDescription(bundle){
 if(bundle.source_model==='unified_radial_exponential_v1')return `공통 상류 반경 지수형 열원 · λq=1·2·3·5·10 mm · q/P_in 정규화 · ${lowerBranchOf(bundle)?'하부 외측 단일 가지':'양방향 가지'} · 가지 파워 ${bundle.branch_power_fractions.map(f=>100*f+'%').join(' / ')} · 고정 열원 범위 0–${bundle.source_interval_mm[1].toFixed(3)} mm`;
 return deviceOf(bundle)==='cpd2'?'값은 q/P_SOL로 정규화하며 P_SOL은 입력 전체 파워 100 MW입니다. 기준은 이번 최적화 실행의 시작 형상입니다.':'값은 q/P_outer로 정규화되었습니다.';
}
