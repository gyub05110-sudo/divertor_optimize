export const statusColors={blocked_first_hit:'#446b70',outside_source_scope:'#c0d3e2',trace_unresolved:'#ed65b3'};
export function statusLabel(status,reason){
 if(status==='outside_source_scope'&&reason==='separatrix_side')return 'separatrix 쪽 자속 · 선택한 외측 SOL 열원에 대응하지 않음';
 if(status==='outside_source_scope'&&reason==='beyond_source_end')return '상류 열원 범위의 바깥 끝을 넘어감';
 return ({direct_hit:'직접 첫 충돌 · 열유속 계산됨',blocked_first_hit:'같은 경로의 다른 표면에 먼저 충돌',outside_source_scope:'선택한 외측 SOL 열원 자속 범위 밖',trace_unresolved:'경로 연결 미확인 · 열유속 미정'})[status]||'추가 검증 미실시';
}
export function compatibleShape(a,b){return !!b&&a.length===b.length&&a.every((p,i)=>p.every((x,j)=>Math.abs(x-b[i][j])<1e-12));}
export function patchAt(patches,r,z){
 let result=null,best=1e-8;
 for(const p of patches||[]){const dx=p.b[0]-p.a[0],dz=p.b[1]-p.a[1],t=Math.max(0,Math.min(1,((r-p.a[0])*dx+(z-p.a[1])*dz)/(dx*dx+dz*dz))),d=(r-p.a[0]-t*dx)**2+(z-p.a[1]-t*dz)**2;
  if(d<best){best=d;result=p;}}
 return result;
}
const cache=new Map();
export async function getQuality(id){
 if(!cache.has(id))cache.set(id,(async()=>{const c=await fetch('./quality_catalog.json',{cache:'no-store'});if(!c.ok)return null;const catalog=await c.json(),url=catalog.cases[id];if(!url)return null;const r=await fetch('./'+url,{cache:'no-store'});if(!r.ok)throw Error('추가 수치검증 파일 로드 실패');return r.json();})());
 return cache.get(id);
}
export function attachQuality(model,audit){
 if(!audit||!compatibleShape(model.wall,audit.reference_wall))return model;
 const arc=point=>{let total=0,best=Infinity,result=0;for(let i=1;i<model.wall.length;i++){const a=model.wall[i-1],b=model.wall[i],dx=b[0]-a[0],dz=b[1]-a[1],len=Math.hypot(dx,dz),t=Math.max(0,Math.min(1,((point[0]-a[0])*dx+(point[1]-a[1])*dz)/(len*len))),distance=Math.hypot(point[0]-a[0]-t*dx,point[1]-a[1]-t*dz);if(distance<best){best=distance;result=total+t*len;}total+=len;}return result;};
 const auditRanges=['floor','shelf'].map(name=>{const a=arc(model[name][0]),b=arc(model[name].at(-1));return [Math.min(a,b),Math.max(a,b)];});
 return {...model,heat:audit.heat,summary:audit.summary,audit,auditRanges};
}
export function showQuality(report,model,stage,bundle){
 let panel=document.getElementById('numericalQuality');
 if(!panel){panel=document.createElement('section');panel.id='numericalQuality';panel.className='case-panel';document.querySelector('.metrics').before(panel);}
 panel.replaceChildren();
 const heading=document.createElement('h2');heading.textContent='열유속 계산 범위와 수치검증';panel.append(heading);
 const line=document.createElement('p');line.className='hint';
 if(bundle?.source_model==='unified_radial_exponential_v1'){
  const v=model.validation,residual=Math.max(...(v?.integral_relative_residual||[]).map(Math.abs));
  line.textContent=`${stage.label} · ${v?.sample_count??0}개 새 첫 충돌 표본 · ${v?.numerically_verified?'피크·면적 적분 세분화 기준 통과':'수치검증 미통과'} · ${v?.status==='complete_feasible'?'타깃 제약 충족':v?.status==='constraint_not_satisfied'?'타깃 제약 미달':'타깃 제약과 수치검증은 별도 판정'}`;panel.append(line);
  const note=document.createElement('p');note.className='hint';note.textContent=`${bundle.model_note}. 표면 열유속 적분과 연결된 열원 구간 파워의 최대 차이 ${(residual*100).toFixed(4)}%. 미표시 구간 파워 ${((Math.max(...(v?.uncovered_source_power||[])))*100).toFixed(4)}%. 연결된 구간만 보간하며 회색은 미확인 구간입니다. 반사 경계의 상류 등가 확산 근사이며 실제 SOL 체적 수송 계산은 아닙니다.`;panel.append(note);
  const link=document.createElement('a');link.href=`./cases/${bundle.id}/case.json`;link.textContent='형상·폭별 열부하·표면별 적분 검증 JSON';panel.append(link);return;
 }
 if(bundle?.device_id!=='cpd2'&&stage.kind==='recalculated_original'&&model.audit){
  line.textContent=`${stage.label} · σ/λq=${bundle.sigma_over_lambda_q} · 저장된 직접 첫 충돌 연결점의 파워 보존 재가중`;panel.append(line);
  const note=document.createElement('p');note.className='hint';note.textContent='평형과 기존 최적 형상을 고정한 확산 재계산입니다. 재최적화 결과는 마지막 단계에서 확인하세요. 반사 경계의 상류 1차원 확산 근사이며 실제 SOL 체적 수송 계산은 아닙니다. 차폐·소스 범위 밖 타깃은 미조사 상태를 유지합니다.';panel.append(note);return;
 }
 if(bundle?.device_id==='cpd2'){
  const v=model.validation,ledger=v?.power_ledgers?.['1'];
  line.textContent=`${stage.label} · ${v?.sample_count??0}개 정밀 첫 충돌 표본 · ${v?.status==='complete_feasible'?'표본 제약 충족':'미확정 / 제약 미달'} · 최소 타깃 입사각 ${v?.minimum_target_angle_deg?.toFixed(6)??'미조사'}°`;panel.append(line);
  const note=document.createElement('p');note.className='hint';note.textContent=`${bundle.model_note}. ${bundle.sigma_over_lambda_q?`확산 σ/λq=${bundle.sigma_over_lambda_q}. `:''}표본 사이에서 같은 면·같은 방향의 가까운 점만 보간합니다. 회색은 데이터가 없는 영역이며 0을 뜻하지 않습니다. 연속 열유속 최대값 또는 전체 타깃 면적 적분 인증은 아닙니다.`;panel.append(note);
  if(ledger){const p=document.createElement('p');p.className='hint';p.textContent=`입력 파워 ${(ledger.source_power_W/1e6).toFixed(3)} MW · 벽에 배정된 파워 ${(ledger.surface_power_estimate_W.reduce((a,b)=>a+b,0)/1e6).toFixed(3)} MW · 미확정 ${(ledger.unresolved_trace_weight_W/1e6).toFixed(6)} MW · 각도 하한 ${bundle.angle_floor_deg}°`;panel.append(p);}
  const link=document.createElement('a');link.href=`./cases/${bundle.id}/case.json`;link.download=`${bundle.id}.json`;link.textContent='CPD 형상·열부하 결과 JSON';panel.append(link);return;
 }
 if(!report||!model.audit){line.textContent='이 단계는 추가 표면 검증 미실시입니다. 회색 타깃은 표본/보간 범위 밖이며 열유속 0을 뜻하지 않습니다.';panel.append(line);return;}
 if(bundle?.sigma_over_lambda_q!=null){const assumption=document.createElement('p');assumption.className='hint';assumption.textContent=`1차원 횡방향 확산 σ/λq=${bundle.sigma_over_lambda_q} · 조건별 후보 ${bundle.evaluation_budget.toLocaleString('en-US')}개 · 반사 경계와 입력 파워 보존 · 기준과 후보에 같은 확산 소스를 적용합니다.`;panel.append(assumption);}
 line.textContent=`기준·최종 형상 추가 검증: ${report.verified?'설정한 수치 기준 통과':'미통과 항목 있음 — 결과는 잠정값'} · 선택 단계: ${stage.label} · 추가 검증 후 다섯 폭 worst ratio ${Math.max(...report.final_ratios).toFixed(5)}`;panel.append(line);
 if(report.source_accounting_passed===false){const warning=document.createElement('p');warning.className='hint warning';const missing=Math.max(...report.levels.at(-1).models.candidate.unmapped_power);warning.textContent=`소스 파워 제약 미달: 미포착 파워 ${(100*missing).toFixed(5)}% > 허용치 0.1%. 각도 위반과 별도 판정합니다. 소스 범위와 파워를 유지한 잠정 후보이며 실행 가능한 최적해가 아닙니다.`;panel.append(warning);}
 if(report.baseline_constraint_feasible===false){const warning=document.createElement('p');warning.className='hint warning';warning.textContent='기준 형상은 shelf의 1.5° 조건에 미달하는 비교용 참조입니다. 최종 후보에 각도 조건을 적용합니다. 수치검증 통과는 기준 대비 열유속 감소를 뜻하지 않습니다.';panel.append(warning);}
 const legend=document.createElement('p');legend.className='hint';
 for(const status of ['blocked_first_hit','outside_source_scope','trace_unresolved']){
  const s=document.createElement('span');s.style.borderLeft=`12px solid ${statusColors[status]}`;s.style.paddingLeft='6px';s.style.marginRight='12px';s.textContent=statusLabel(status);legend.append(s,document.createElement('br'));
 }
 const s=document.createElement('span');s.textContent='열유속 색상: 직접 연결된 타깃 · 연한 회색 반투명: 주변 구조';legend.append(s);panel.append(legend);
 const selected=model.audit===report.models.baseline?'baseline':'candidate';
 const residual=report.surfaces[selected].quadrature_levels.at(-1).relative_residuals;
 const text=document.createElement('p');text.className='hint';text.textContent=`최종 상류 경로 ${report.levels.at(-1).source_cells.toLocaleString()}개 · 타깃 면적 적분 최대 오차 ${(100*Math.max(...Object.values(residual).flat().map(Math.abs))).toFixed(3)}% · 연결 미확인 ${(1000*report.surfaces[selected].status_length_m.trace_unresolved).toFixed(3)} mm`;
 panel.append(text);
 if(report.source_domain_audit){const audit=report.source_domain_audit,note=document.createElement('p');note.className='hint';note.textContent=`상류 벽까지의 유효 범위: ${audit.physical_source_max_mm.toFixed(2)} mm. ${audit.physical_source_max_mm<=100?'기존 100 mm는 이미 이 범위를 덮습니다.':'100 mm를 넘어 벽까지의 범위는 별도 검사가 필요합니다.'} 상류 100 mm 끝에서 잘린 타깃 구간은 ${(audit.source_cutoff_flagged_target_length_m*1000).toFixed(3)} mm입니다. 벽 밖으로 상류를 확장하지 않습니다.`;panel.append(note);}
 const details=document.createElement('details'),summary=document.createElement('summary');summary.textContent='해상도별 최대 열유속과 검증 기준';details.append(summary);
 const table=document.createElement('table');table.style.width='100%';table.style.fontSize='11px';table.style.tableLayout='fixed';
 const header=document.createElement('tr');for(const name of ['추적','경로 수','기준 q/P (3 mm)','최종 q/P (3 mm)','최종/기준']){const cell=document.createElement('th');cell.textContent=name;header.append(cell);}table.append(header);
 for(const l of report.levels){const row=document.createElement('tr');for(const value of [l.kind==='uniform'?'균일':'추가 세분화',l.source_cells,l.models.baseline.peaks[2].toFixed(5),l.models.candidate.peaks[2].toFixed(5),l.models.candidate.ratios[2].toFixed(5)]){const cell=document.createElement('td');cell.textContent=value;row.append(cell);}table.append(row);}details.append(table);
 const note=document.createElement('p');note.className='hint';note.textContent=`최대값·개선비 변화 ≤1%, 피크 위치 변화 ≤2 mm, 타깃 면적 적분 오차 ≤1% (작은 성분은 입력 파워 10⁻⁵ 기준). 표면 연결 검사는 0.5 mm 이하 간격에서 수행합니다. 상류 파워 합 검사와 타깃 열유속 면적 적분을 별도로 확인합니다. ${report.forward_converged?'추적 수렴 통과':'추적 수렴 미통과'} / ${report.surface_integral_passed?'면적 적분 통과':'면적 적분 미통과'}. 고정된 다각형 형상과 외측 하부 SOL 모델에 대한 검증입니다.`;details.append(note);
 const link=document.createElement('a');link.href='./quality/'+report.case_id+'.json';link.download='heat-quality.json';link.textContent='상세 검증 결과 JSON';details.append(link);panel.append(details);
}
