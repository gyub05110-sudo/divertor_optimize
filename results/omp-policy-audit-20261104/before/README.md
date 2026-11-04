# D3D 형상 개선 과정 뷰어

이 뷰어는 D3D IMP/OMP 독립 1 MW 연구의 저장된 최고 기록 갱신 형상을 비교합니다. Rmin은 20·50·100 mm의 실제 연구 기록 중에서 선택합니다. 형상 사이를 보간해 새 최적해를 만들지 않습니다.

## 최신 판정: 분포 수렴 v3 (2026-11-04)

최대점 위치 2 mm 조건을 필수 합격조건에서 제외하고, 같은 물리적 벽 위치의 q 분포·전력·확산 면적 수렴으로 대체했습니다. 봉우리 순위 교체만으로 탈락시키지 않습니다. 저장된 연속 q 샘플의 재구성으로 Linf ≤1%, 면적 가중 L1 ≤0.5%, A2 변화 ≤0.5%를 검사하고 기존 형상·전력·peak·독립 first-hit/열적 검사를 유지합니다. 허용오차는 연구용 수치 설정이며 연속 오차의 수학적 인증이나 제작 표준이 아닙니다.

- **OMP 20 mm:** 보관 후보 19개 전체(기존 탈락15개, 이전 최종1개, 추가 정밀평가3개)가 v3의 위치별 q 차이 기준을 초과했습니다. 최대 차이는 전역 peak의 약2.64–2.77%입니다. 이전 4단계 형상과 1.444160060 MW/m² 결과는 그대로 보존하며 화면에 **이전 기준 통과 · v3 미통과**로 표시합니다. 기존 peak 감소가 사라졌다는 뜻은 아닙니다.
- **OMP 50 mm:** 34단계가 v3 통과. peak1.786940489 MW/m², 감소34.562428%, 최대 위치별 q 차이0.404226%, L1 0.202975%, A2 변화0.011253%. 원래1280 순위에서 최저 peak 후보가 통과하므로 나머지33개 높은 peak 후보의 전체 검증은 반복하지 않았습니다.
- 고정 기준도 v3 통과. IMP와100 mm는 이번 재검증 대상이 아닙니다. 기존237단계·원래6개 정밀 결과·이전 추가2개 결과·기준/leg/열분포 파일은 보존했습니다.

새 `distribution-finals.json`을 원래 `qualified-finals.json` 다음에 적용합니다. 기존 파일을 덮어쓰지 않고 별도 `distribution-profiles/`, `distribution-evidence/`를 사용합니다. 동일한50 mm 후보의 CPU/GPU 검사는 정확히 같은 NPZ SHA-256에 연결된 기존 검사를 재사용합니다. 초기3개 OMP20 추가 검증은 GPU를 사용해 약565초, 최종 재판정·출판은 약28초 소요됐습니다. 새 형상 최적화 탐색은 없습니다.

최종 study/검증/보관본: `D:/DivertorOptimization/studies/omp-distribution-final-20261104`.
상세 보고서: 같은 경로의 `report/REPORT.md`. 이전 연구는 변경하지 않았으며, 첫 구현의 평균 구간 차이 문제와 수정된 비교법도 각각 별도 코드·study로 남겼습니다. 뷰어 변경 전 파일은 `D:/DivertorOptimization/studies/omp-distribution-viewer-backup-20261104`에 보관했습니다.

검수: Python26개·JavaScript22개 테스트,867개 표시 셀 조회, 보호 파일/결과 체크섬, 실제 브라우저의 상태·q/좌표·중간 단계·IMP/100 mm 회귀 검사. 검수 파일은 `distribution-study-verification.json`, `distribution-viewer-verification.json`, `distribution-browser-observations.json`입니다. 아래 추가 정밀검증 설명은 **v3 이전 정책의 역사적 결과**를 포함합니다.

## 실행

- `start-viewer.cmd`를 더블클릭하거나 `launch.ps1`을 실행합니다.
- 주소: <http://127.0.0.1:8780/>
- Python: `C:/Users/FPL-user/miniforge3/python.exe`.
- 서버는 이 컴퓨터의 loopback 주소에서만 실행됩니다. 재부팅 후에는 실행 파일을 다시 열어 주세요.
- `dist/index.html`을 파일로 직접 열면 JSON 로딩이 제한되므로 위 실행 방법을 사용합니다.

모든 새 코드·뷰어 데이터·스크린샷은 이 D 드라이브 폴더에 있습니다. GitHub push나 외부 호스팅은 하지 않았습니다.

## 보는 순서

1. **IMP / OMP → Rmin → 탐색 실행 기록**을 고릅니다. 과거 접선 비교 연구와 72시간 캠페인은 서로 다른 시계열입니다.
2. 아래 **단계 선택·슬라이더·이전/다음·자동 재생**으로 기준과 기록 갱신 과정을 확인합니다. 차트의 점을 눌러 이동할 수도 있습니다.
3. 오른쪽 R–Z 단면에서 기준 점선과 선택 실선을 비교합니다. **접합부 확대**는 축척을 왜곡하지 않습니다.
4. **최종 검증 통과 / 추가 정밀검증 통과 / 과거 정밀 결과**는 별도로 선택합니다. 동일 해상도의 기준 peak와 비교하며, 탐색 단계의 연속 기록으로 붙이지 않습니다. 추가 정밀검증 결과는 원래 저장 단계 번호도 표시하며 단계 선택으로 돌아가면 해당 단계의 320 경로 분포를 보여줍니다.
5. 3D는 드래그 회전, 휠 확대, 오른쪽 드래그 이동을 지원합니다. 회전 범위, 전체 벽, LCFS, 기준/선택/겹쳐보기를 조절할 수 있습니다.
6. **q 크기에 따라 색상 표시**를 끄면 형상의 굴곡을 단색 재질로 확인할 수 있습니다. 기준 형상·모든 중간 단계·정밀 결과에 각 형상의 열분포를 연결했습니다.
7. **3D 이미지 저장**은 현재 보기를 PNG 미리보기로 엽니다. PNG 링크 또는 그림의 오른쪽 클릭 메뉴로 저장합니다. Codex 내장 브라우저에서는 다운로드 완료 이벤트가 확인되지 않았습니다. PNG 생성·미리보기·파일 내용은 확인했으며, 직접 다운로드가 동작하지 않으면 일반 브라우저에서 같은 로컬 주소를 열어 주세요.
8. **Divertor leg 표시**는 기본으로 켜져 있습니다. Inner leg는 파랑, Outer leg는 자주색입니다. ×는 X-point, 채운 점은 현재 벽의 도착점, 빈 원은 겹쳐보기의 기준 도착점입니다. 기준/선택/겹쳐보기와 단계 변경에 맞춰 leg도 바뀝니다.

## 위치별 q 조회와 색상

- 3D 표면이나 R–Z 단면의 벽에 마우스를 올리면 q [MW/m²]와 R/Z [m]가 나타납니다. 터치·클릭으로도 확인할 수 있습니다. 회전 드래그, 화면/형상 전환, 포인터 이탈, Escape에서는 이전 툴팁을 지웁니다.
- **q 크기에 따라 색상 표시**는 기본으로 켜져 있습니다. 정밀 결과가 있는 실행은 그 결과를 먼저 보여줍니다. 기록 갱신 과정은 기존 단계 선택과 슬라이더로 볼 수 있습니다.
- 3D·2D·툴팁·수치 범례가 같은 색상 함수를 사용합니다. 척도는 0에서 해당 열원의 기준 정밀 peak까지로 고정하여 후보마다 색상을 재정규화하지 않습니다. 색상을 꺼도 q 조회는 가능합니다.
- q는 저장된 연속 셀 최대 4개의 전력/면적을 합산한 **표시 셀 평균**입니다. 상단의 정밀 peak는 원래 샘플의 최대값이므로 툴팁의 최대 평균과 약간 다를 수 있습니다. 표시 셀 사이를 보간해 q값을 만들지 않습니다.
- 모든 중간 단계는 원본 후보 벡터·고정 입력·cubic 320 경로로 재평가했습니다. 원래 기록된 peak와 상대오차 1e-7 이내로 맞고 전력이 보존된 분포만 사용합니다. 화면에도 **320 경로 · 해당 단계 재평가**로 표시합니다. 최종 결과의 q를 중간 단계에 복사하지 않습니다.
- 전체 벽의 분포를 보존하여 **장치 전체 벽 표시**에서도 색상과 q 조회가 가능합니다. q=0과 파일/위치 대응 실패는 구분합니다. q는 각 셀 평균이며, 격자 아래 수준의 연속적인 점 측정값은 아닙니다.
- 기준 형상만 표시하거나 기준 단계를 선택하면 저장된 1280 경로의 기준 분포를 사용합니다. 기준 단계의 peak도 이 정밀 기준값으로 표시하며 원본 320 경로의 단계별 차트 데이터는 유지합니다. 겹쳐보기에서 기준 점선을 가리키면 기준 형상 보기로 전환하도록 안내합니다.
- 3D 조회는 실제 보이는 벽의 ray intersection을 사용합니다. 다각형 토러스의 `hypot(x,z)` 반경 오차를 피하도록 원래 R/Z를 꼭짓점 속성으로 보간합니다. LCFS·leg·기준 유령선은 벽으로 선택하지 않습니다. 저장 분포와의 위치 대응 허용오차는 기존 형상 비교와 같은 20 µm입니다.

검증: `node --test test_math.mjs test_heat_probe.mjs test_stage_profiles.mjs`의 17개 테스트와 Python 10개 테스트를 통과했습니다. 237단계 전체의 프로파일 연결·Three.js 표면 선택, 11,613개 중점 조회, 720,584개 표시 셀의 전력·면적·좌표를 대조했습니다. 실제 브라우저의 12개 열원/반경/실행 조합과 7개 q·좌표 조회도 확인했습니다. 화면은 `screenshots/all-stage-*.jpg`, 검수 기록은 `stage-profiles-verification.json`, `stage-probe-verification.json`, `stage-browser-verification.json`에 보존합니다.

237단계에는 동일 후보가 두 번 등장하므로 고유 형상 236개를 약 521초에 GPU로 재평가했습니다. IMP/OMP 각 대표 후보의 CPU/GPU q 차이는 최대 1.34e-15 MW/m², 모든 중간 단계의 원본 peak 상대오차는 최대 2.23e-16입니다. 전체 표시 셀을 회전면 면적으로 다시 적분한 전력 오차는 최대 1.07e-13입니다. 이 검사는 시각화용 단일 격자 재현 검사이며, 중간 단계의 다중 격자 정밀검증 통과를 뜻하지 않습니다.

원본 `data.json`, `legs.json`, 이전 열분포, 고정 입력은 해시로 보존을 확인했습니다. 새 원본 계산 배열은 `stage-evaluations/<candidate-id>/selected.npz`, 화면 분포는 `dist/stage-profiles/`, 기존 정밀 분포의 전체 벽 버전은 `dist/full-saved-profiles/`에 있습니다. 두 색인은 후보 ID·열원·해상도·해시가 맞아야 로딩되며, 누락 시 다른 형상의 분포로 대체하지 않습니다. 단계별 재평가 스크립트는 기존 파일과 실행 조건이 일치하는 경우에만 재개합니다.

## Leg 표시의 의미

고정 D3D 자기평형의 ψ 분리면 가지를 X-point부터 추적하고, 각 표시 벽과 처음 만나는 지점에서 자릅니다. 따라서 아래로 내려간 OMP 타깃에서는 원래 기준 도착점보다 leg가 연장됩니다. 두 가지를 기하학적으로 한 번씩 복원해 239개 벽에 대한 478개 leg를 저장했습니다. 새 최적화와 열수송 평가는 없습니다.

3D의 반투명 면은 분리면 가지를 회전한 축대칭 면이며, 색 선은 토로이달 절단면의 R–Z 가지입니다. 나선형 입자 궤적을 새로 계산한 것은 아닙니다. `poloidal_length_m`는 단면 호 길이이며 상류에서의 평행 연결 길이 L_parallel과 다릅니다. 표시용 꺾은선 단순화 허용오차는 1 µm입니다. 원래 형상/peak JSON은 그대로 보존하고 `legs.json`의 SHA-256 연결을 검사하므로, 다른 버전의 형상에 잘못된 leg를 붙이지 않습니다.

기준 도착점은 보존된 baseline review와 2.3×10⁻¹⁶ m 이내에서 일치했습니다. 이 수치는 기록 간 일치도이며, 물리적 위치 정확도를 의미하지 않습니다. 모든 leg의 좌표 유한성, 첫점 X-point, 끝점의 해당 벽 일치, 후보 ID 대응을 검사했습니다. 근거는 `legs-audit.json`과 `legs-verification.json`입니다.

## 최종 파일의 합격 범위

v3 도입 전 재검증을 통과한 개선 결과는 총 5개였습니다. 기존 3개는 캠페인 IMP 20/50 mm와 접선·곡률 비교 OMP 50 mm이며 감소율은 6.8432%, 2.6713%, 26.9121%입니다. 2026-11-04 추가 검증에서 캠페인 OMP 20/50 mm 2개를 더 연결했습니다. 나머지 정밀 결과 3개(접선 비교 OMP 20 mm, 초기 IMP 20/50 mm)는 과거 검증으로 분리합니다.

기존 3개는 설정된 수치 허용오차 안에서 접선·곡률·길이·개구부·전체 벽 충돌, 전력 보존, first-hit 독립 검사, 전도/sheath 잔차, cubic 320/640/1280·quintic 640·별도 벽 해상도 검사를 통과했습니다. 이 3개의 2026-11-03 재검증은 CPU 평가이고, 원 실행의 같은 후보 manifest에서 CPU/GPU peak 상대차는 모두 0입니다. 이 3개의 `promoted_to_qualified_store`는 false입니다. 이 화면의 '최종'은 재검증된 저장 후보를 뜻하며 전역 최적성 증명이 아닙니다.

기존 2026-11-03 재검증의 전체 코드 보존 검사에는 예외가 있습니다. 원 재검증 `summary.json`은 실행 중 미사용 모듈 두 개(`continuous_path_bounds.py`, `continuous_transport_certificate.py`)가 바뀌어 `code_preservation.unchanged=false`입니다. 별도 `code_scope_audit.json`은 수치 평가에서 사용한 코드의 해시 보존을 확인하고 accepted=true로 승인했습니다. 의존성은 사후 import/정적 분석으로 재구성했으며, 원 실행 import telemetry는 없습니다. 원 보고서나 수치 결과를 수정하지 않았습니다. 이 예외는 아래 2026-11-04 추가 검증에 해당하지 않습니다.

원 근거: `D:/DivertorOptimization/studies/candidate-requalification-20261103/summary.json`, 각 후보의 `validation.json`, `code_scope_audit.json`. 실제 summary SHA-256은 감사에 기록된 `76a40ff2e4364bd83707836b6db165325abe04a8a03055f99b0f6936effbf7fc`와 일치합니다.

### 2026-11-04 OMP 추가 정밀검증

| Rmin | 선택된 원래 단계 | 1280 peak MW/m² | 기준 대비 감소 | 전체 검증 시도 |
|---:|---:|---:|---:|---:|
| 20 mm | 4 | 1.444160060 | 47.115011% | 16 |
| 50 mm | 34 | 1.786940489 | 34.562428% | 1 |

고정 D3D 기준과 독립 1 MW 입력을 유지하고 19 + 34개 보관 후보 모두를 cubic 1280에서 평가했습니다. 낮은 peak 순으로 전체 검증하여 각 반경의 첫 통과 후보를 선택했습니다. 53개 모두 전체 검증을 통과했다는 뜻은 아닙니다. 전체 검증 17건 중 2건 통과, 15건 탈락이며 나머지 36개는 선택 후보보다 1280 peak가 낮지 않아 전체 검증을 수행하지 않았습니다.

선택된 두 후보는 기존 형상·독립 first-hit/열적 검사·전력 보존·320/640/1280·quintic·별도 벽격자·개선 확인 기준을 모두 통과했습니다. CPU/GPU q 최대 절대차는 1.34e-15 MW/m² 미만이며 셀 전력 차이는 0입니다. 실행 코드와 보호 입력의 해시 보존도 통과했고, 체크섬이 있는 불변 QualifiedStore에 영구 보존했습니다.

20 mm에서 peak가 더 낮은 15개는 최대점 위치 수렴 기준(2 mm)에 걸렸습니다. 거의 같은 높이의 두 봉우리 중 최대값 위치가 바뀌는 현상으로, 전체 분포가 그 거리만큼 이동했다는 뜻은 아닙니다. 이번에 기준을 완화하지 않았습니다. 이 결과는 **해당 보관 후보 목록에서 현재 모든 검증을 통과한 최저 1280 peak 후보**이며, 연속 형상 공간의 전역 최적성 보장이 아닙니다. 100 mm는 이번 추가 검증 대상이 아닙니다.

원본 `data.json`과 237개 단계·기존 6개 정밀 결과는 그대로 보존합니다. `qualified-finals.json`이 해시로 연결된 2개 최종 결과를 추가하며, `qualified-profiles/`와 `qualified-evidence/`에 해당 1280 분포와 검증 근거를 보존합니다. 같은 후보의 중간 단계는 원래 320 분포를 사용합니다.

결과·제외 사유·그림: `D:/DivertorOptimization/studies/omp-stage-qualification-20261104/report/REPORT.md`. 원본 배열·순위·검증·불변 보관본은 같은 study에 있습니다. 기계 판독 검수는 `qualified-final-verification.json`, `qualified-probe-verification.json`, `qualified-browser-verification.json`이며 실제 화면은 `screenshots/additional-omp*.jpg`에 있습니다.

## 포함된 기록

총 12개 실행 조건, 형상 검사를 통과한 최고 기록 갱신 단계 237개와 별도 정밀 결과 8개(기존 6개 + 추가 2개)입니다. 다른 실행에서 같은 후보가 사용될 수 있으므로 237개를 서로 다른 물리 형상의 수로 해석하면 안 됩니다.

| 실행 | 열원 | Rmin mm | 표시 단계 | 정밀 peak MW/m² | 정밀 결과 상태 |
|---|---|---:|---:|---:|---|
| 72시간 캠페인 | IMP | 20 | 92 | 2.375184 | 현재 재검증 통과 |
| 72시간 캠페인 | IMP | 50 | 61 | 2.481554 | 현재 재검증 통과 |
| 72시간 캠페인 | IMP | 100 | 0 | — | 표시할 개선 결과 없음 |
| 72시간 캠페인 | OMP | 20 | 19 | 1.444160 | 이전 기준 통과·v3 미통과·원래 4단계 |
| 72시간 캠페인 | OMP | 50 | 34 | 1.786940 | 분포 수렴 v3 통과·원래 34단계 |
| 72시간 캠페인 | OMP | 100 | 0 | — | 표시할 개선 결과 없음 |
| 접선·곡률 비교 | OMP | 20 | 16 | 1.446258 | 당시 정밀검증·현재 형상 통과 |
| 접선·곡률 비교 | OMP | 50 | 11 | 1.995853 | 현재 재검증 통과 |
| 접선·곡률 비교 | OMP | 100 | 0 | — | 표시할 개선 결과 없음 |
| IMP 초기 탐색 | IMP | 20 | 1 | 2.404015 | 당시 정밀검증·현재 형상 통과 |
| IMP 초기 탐색 | IMP | 50 | 3 | 2.530360 | 당시 정밀검증·현재 형상 통과 |
| IMP 초기 탐색 | IMP | 100 | 0 | — | 표시할 개선 결과 없음 |

OMP 캠페인의 원래 기록 갱신 후보 중 20 mm에서 34개, 50 mm에서 102개는 현재 전체 벽 형상 검사에서 제외됐습니다. 이 캠페인의 탐색 중간 기록과 과거 연구의 검증된 정밀 결과를 혼동하지 않도록 실행 선택과 검증 표시를 분리했습니다.

## 데이터 선정과 해석

- 연구 목표는 고정 D3D 기준/입력에서 peak를 낮추고 퍼뜨리는 빠른 형상 최적화입니다. 실험 측정값 일치는 필터로 사용하지 않습니다.
- 원래 `feasible=true`이며 같은 실행의 peak/baseline 비율이 이전 최고 기록보다 `1e-9` 이상 낮아진 후보를 추출합니다. 이후 중복 ID와 현재 형상 검사 실패를 제외합니다.
- 현재 형상 검사: 접선 연결, 곡률, OMP 개구부 조건, 길이 비율 ≤ 1.2, 전체 벽 자기 교차·LCFS·용기 교차. 전체 벽 정책은 `full-wall-chord-1e-7.v2`입니다.
- 이 방식은 **원래 기록 갱신 이력**을 보여줍니다. 탈락한 원래 최고값을 넘지 못한 후속 유효 후보는 표시 대상이 아닙니다. 따라서 빈 기록은 개선 형상의 부재나 불가능성 증명이 아닙니다.
- 탐색은 cubic 320 경로의 같은 기준값을 사용합니다: IMP 2.5495884235, OMP 2.7306578740 MW/m². 정밀 결과는 cubic 1280 경로: IMP 2.5496630284, OMP 2.7307561171 MW/m².
- 중간 형상에는 원래 peak만 기록돼 있어 시각화용 위치별 q를 재평가해 별도 보존했습니다. 이 단일 해상도 재평가는 전체 정밀검증을 대체하지 않습니다. 이후 별도로 수행한 OMP 20/50 mm 전체 검증 결과만 추가 정밀검증으로 승격했습니다.
- 열색상은 저장된 `cell_power × power_W / 10^6 / cell_area`를 사용합니다. 화면 경량화를 위해 같은 부모 구간의 연속 셀 최대 4개를 면적 가중 평균합니다. 전력과 면적을 함께 합산하며 비순서 donor 샘플을 연결하지 않습니다. peak 값과 ◇ 위치는 원래 정밀 샘플입니다.
- 2D와 3D는 같은 색상 스케일을 씁니다. 후보 전환 때마다 스케일을 재정규화하지 않습니다. 미세한 셀 내부 변화는 표시용 평균에서 생략됩니다.
- 3D는 R–Z 벽을 토로이달 방향으로 회전한 축대칭 표시입니다. 새 3D 수송 시뮬레이션이나 새 최적화는 실행하지 않았습니다.

## 출처와 보존

- 72시간 기록: `D:/DivertorOptimization/studies/imp-omp-global-20261029`
- 초기 IMP 기록: `D:/DivertorOptimization/studies/imp-global-pilot-20261029-v3`
- OMP 접선 기록: `D:/DivertorOptimization/worktrees/d3d-report-correction-20261104/results/omp-tangent-curvature-20261029`
- 현재 재검증: `D:/DivertorOptimization/studies/candidate-requalification-20261103`
- 추가 OMP 정밀검증: `D:/DivertorOptimization/studies/omp-stage-qualification-20261104`
- 기준 벽·형상 복원 코드: D3D worktree의 `verification/fixtures/campaign_geometry_20261029.json` 및 `d3d_optimization`.
- 기존 사이트에서는 Three.js 0.180.0, OrbitControls, 라이선스와 화면 구성 방식을 재사용했습니다. 기존 사이트의 case/catalog/형상/열분포 데이터는 사용하지 않았습니다.
- `data-audit.json`에 모든 입력 파일 SHA-256, 복원 코드 해시, 필터 개수, 원본 파일의 내보내기 전후 해시 일치를 기록했습니다. CPD_5 저장소와 D3D worktree는 수정하지 않았습니다.

## 검증과 재생성

`export_data.py`는 약 4 GB의 원본 로그를 스트리밍하여 정적 JSON을 생성합니다. Python scientific 환경은 설치된 miniforge를 사용합니다. 뷰어 사용 자체에는 Python HTTP 서버만 필요합니다.

```powershell
& 'C:/Users/FPL-user/miniforge3/python.exe' export_data.py
& 'C:/Users/FPL-user/miniforge3/python.exe' export_legs.py
& 'C:/Users/FPL-user/miniforge3/python.exe' export_stage_profiles.py
& 'C:/Users/FPL-user/miniforge3/python.exe' export_full_saved_profiles.py
& 'C:/Users/FPL-user/miniforge3/python.exe' -m unittest test_export.py test_legs.py
node --test test_math.mjs test_heat_probe.mjs test_stage_profiles.mjs
& 'C:/Users/FPL-user/miniforge3/python.exe' verify_data.py
& 'C:/Users/FPL-user/miniforge3/python.exe' verify_legs.py
& 'C:/Users/FPL-user/miniforge3/python.exe' verify_stage_profiles.py
node verify_stage_probe.mjs
node verify_stage_browser.mjs
```

명령은 이 폴더에서 실행합니다. 마지막 검증 명령은 서버가 실행 중이어야 합니다. `verification.json`은 후보 ID·원본 237행·동일 격자 비율·단조 기록·형상 정책·정밀 NPZ peak·총 전력·HTTP MIME 검증 결과입니다. `browser-verification.json`과 `screenshots/`는 실제 브라우저 확인 기록입니다.

리뷰에서 발견한 OMP 접선 오차의 scalar/array 차이와 직선 곡면의 null 반경 표시 오류를 수정하고 회귀 테스트를 추가했습니다. Windows의 기본 HTTP 서버가 `.mjs`를 text/plain으로 응답하는 문제는 `serve.py`에서 명시적인 JavaScript MIME으로 해결했습니다.

추가 OMP 결과의 검사(새 열수송 계산 없이 저장 결과를 검사):

```powershell
& 'C:/Users/FPL-user/miniforge3/python.exe' -X utf8 -m unittest test_export.py test_legs.py test_stage_profiles.py test_qualification_selection.py
node --test test_math.mjs test_heat_probe.mjs test_stage_profiles.mjs test_qualified_finals.mjs
& 'C:/Users/FPL-user/miniforge3/python.exe' -X utf8 verify_qualified_finals.py
node verify_qualified_probe.mjs
node verify_qualified_browser.mjs
```

마지막 명령은 실제 브라우저에서 저장한 관측값을 원본 표시 셀과 대조합니다. 브라우저를 새로 조작하는 자동화 명령은 아닙니다. 완료된 study의 실행 코드·입력·핵심 결과는 수정하지 마세요.
