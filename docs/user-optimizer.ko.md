# 평형 파일에서 divertor 피크를 최소화하는 명령

저장소 루트에서 실행한다. Python·NumPy·SciPy·Shapely·PyArrow 버전은 [재현 환경](../requirements-repro.txt)에 기록돼 있다. CPU가 기본이다. 호환되는 CuPy/CUDA가 있으면 `--backend gpu`로 후보의 열유속 합산을 가속할 수 있다. PowerShell에서 기존 로컬 PyArrow 복사본을 사용할 때만 `$env:PYTHONPATH = (Resolve-Path '.local_dependencies/d3d-fixed-equilibrium').Path`를 먼저 실행한다.

```powershell
python optimize_divertor.py --equilibrium input/D3D/d3d/d3d_shot_1a3155370c.parquet --time-ms 1400 --budget-seconds 600 --method screened --family hotspot --output D:/DivertorOptimization/studies/my-first-design
```

출력 폴더는 매번 새 경로여야 한다. 파일 위치·파일명은 자유다. 예산은 Python 스크립트 시작부터 입력 확인, 원벽 계산, 후보 생성과 정밀 평가까지의 경과 초다. 진행 중인 계산 하나가 예산을 넘길 수 있으며 해당 평가의 후보는 선택 대상에서 제외한다. 내보내기 시간은 따로 기록한다. `events.jsonl`은 실행 중에도 읽을 수 있고, 정상 종료 시 `result.json`이 생긴다.

## 입력과 물리 범위

- 현재 장치는 DIII-D, 열원은 외측 중간면(OMP), 설계 구역은 기존 기준벽의 하부 4개 면이다. 단일 방전의 공개 DIII-D EFIT **parquet**와 그 파일에 한 번만 존재하는 `--time-ms`가 필요하다. 지정한 시각의 R–Z 자속 격자·LCFS·자기축, 같은 시각을 포함하는 코일전류 시계열을 사용한다. G-EQDSK는 지원 범위 밖이다. 원래 분리면의 하부 두 가지와 최초 도달면을 파일마다 다시 계산한다.
- 기본 물질벽·진공용기는 [DIII-D 기하](../d3d_optimization/data/DIIID_geom.json), 설계 가능 면·고정 끝점은 [영역 정책](../baseline_review/optimization_region_policy.json)과 [기준 형상 정의](../baseline_review/baseline_review.json)를 사용한다. `--geometry`와 `--design-policy`로 같은 정의를 다른 위치에서 제공할 수 있다. 임의 장치나 당시 실제 타일을 복원하는 입력 규약은 포함하지 않는다.
- 열입력 기본값은 OMP 1 MW, 전자 밀도 10¹⁹ m⁻³, χ⊥ 0.3 m²/s, λq 3 mm, sheath 계수 7, 전도계수 κ₀ 2000이다. 각각 `--power-mw`, `--density-m3`, `--chi-m2-s`, `--lambda-q-mm`, `--gamma`, `--kappa0`로 지정한다. 주어진 입력은 실행 중 고정한다.
- `--radius-mm`는 20, 50, 100 중 하나이고 `--seed`는 음이 아닌 정수다. 고정 끝점, 바닥 G1 접선, 개구 간격, 곡률·길이, LCFS·용기 교차 및 전력·입사각 검사는 공통이다. 원래 직선 바닥과 인접 면의 각을 보존해 비교 원벽으로 사용한다.
- 각 후보에서 첫 충돌, 토로이달 평행 경로 길이, 경로별 2점 모형과 열폭을 갱신한다. 기본값 `--method screened`는 320개 자속관으로 선별하고 1280개와 정밀 벽 격자로 피크를 확인한다. `--method direct`는 `--family full`로 정밀 차등진화를 수행한다. 두 방법 모두 제한시간 안에 완료된 적격 후보의 전 벽 법선 열유속 최댓값을 최소화한다.

## 탐색할 면 선택

| `--family` | `screened`에서 움직이는 면 |
|---|---|
| `full` (기본값) | OMP 입구·선반·연결부·바닥, 총 16계수 |
| `floor-only` | 바닥 4계수 |
| `floor-entry` | 바닥 4계수와 OMP 입구 4계수 |
| `floor-shelf` | 바닥 4계수와 선반 4계수 |
| `floor-connector` | 바닥 4계수와 OMP 연결부 4계수 |
| `hotspot` | 원벽 열점으로 위 다섯 가족 중 하나를 선택 |

원벽의 정밀 계산이 예산 내에 끝나면, 선별 방식의 모든 가족에서 후보를 만들기 전에 **원벽**의 면별 법선 피크(MW/m²)와 총 침적 전력(W)을 기록한다. 전 벽 최고 피크와의 차가 그 최고값의 10⁻¹⁰ 이내인 면들을 동률로 본다. 동률에 고정벽이나 IMP 면이 포함되면 판단을 유보하고 자동 모드는 `full`을 쓴다. 나머지는 침적 전력이 가장 큰 동률 면을 택하고, 전력도 같으면 바닥·OMP 연결부·선반·OMP 입구 순으로 택한다. 바닥을 고르면 `floor-only`, 다른 면을 고르면 바닥과 그 면만 움직인다. 모든 제안에서 바닥 G1 접선 투영을 유지하고, 선택에서 제외된 면의 계수는 0으로 고정한다.

실행 중 선택한 가족을 유지한다. `result.json`의 `surface_policy`에 예상 가족(`predicted_family`), 요청한 가족(`requested_family`), 실제 가족(`selected_family`), 동률 면·측정값·유보 사유가 있다. `predicted_floor_sufficient`는 바닥 예측이면 `true`, 다른 조정 면이면 `false`, 유보면 `null`이다. 이는 주어진 수송 모형과 예산에서 어떤 면을 탐색할지에 대한 사전 예측이다. 실제 개선율은 정밀 평가를 마친 후보에서 확인한다. `events.jsonl`의 `surface_policy_decision`에 같은 결정과 시각이 남는다. 조기 종료로 결정이 없거나 직접 방식이면 `surface_policy`는 `null`이다.

## 출력 읽기

| 상태 | 뜻 |
|---|---|
| `qualified` | 시간 안에 정밀 평가와 제약을 통과한 후보 중 최소 피크. 원벽과의 개선율은 별도로 확인한다. |
| `geometry_excluded` | 입력 직후 고정 물질 구간 또는 조정 가능 면의 고정 끝점과 닫힌 LCFS 다각형의 교차·접촉을 확인했다. 모든 허용 OMP 후보가 이 교차를 공유하므로 열계산 전에 종료하며 `baseline`은 `null`이다. `fixed_wall`에 원래 선분 번호·교차 좌표가 남는다. |
| `constraint_excluded` | 열계산 가능한 원벽을 저장했다. 고정 끝점의 접선과 양의 폭을 가진 주기 Gaussian의 전 벽 지지 조건에서 허용한 어떤 선반도 입사각 하한을 충족할 수 없다. |
| `budget_exhausted` | 입력·원벽·탐색 준비 중 예산이 소진됐다. 늦게 측정된 원벽에는 `within_budget=false`를 기록한다. |
| `no_design` | 탐색했으나 시간 내 정밀 적격 후보가 없었다. 적격 형상의 부재를 입증한 것은 아니다. |
| `input_error` | 지원하지 않는 파일 구조·시각·물질벽·해석할 수 없는 자기 위상이나 잘못된 입력. 사유는 `reason`, 종료 코드는 2다. |
| `calculation_error`, `partial_result` | 수치 계산 또는 실행 환경의 실패. 후자는 앞서 완료한 적격 후보도 보존한다. 종료 코드는 1이다. |

원벽의 프로파일과 대체 설계의 적격은 별도로 기록한다. 고정 접선과 국소 자기피치가 주는 상한은 `fixed_shelf`, 열을 직접 받은 자속관 및 재분배된 셀의 선반 위반 전력은 각 형상의 `shelf_exposure`에 있다. `result.json`에는 파일·벽·정책·코드 해시, 자기축·X점·strike, SI 입력, 각 면의 첫 충돌 전력과 침적 전력, 피크와 위치, 전력수지 및 시간 기록이 있다.

- `baseline/`, `selected/`의 `geometry.csv`는 닫힌 벽의 순서 있는 R,Z 좌표와 각 선분의 이름이다. `profile.csv.gz`에는 같은 순서의 세분 셀 양 끝·호길이·회전면적·침적 W와 **셀 평균** MW/m²가 있다. 형상에 따라 셀 개수와 실제 호길이는 달라질 수 있다.
- `boundary.npz`의 `points`와 `q_MW_m2`는 정규 표본, 첫 충돌 근방 표본, 정련한 극댓값을 포함한다. 각 표본의 `sample_parent`는 0부터 시작하는 물질벽 선분 번호이며, `geometry.csv`에서 같은 `vertex` 행의 `outgoing_region`을 읽으면 면 이름을 확인할 수 있다. 보고한 피크는 이 분포의 최댓값이다. A₂는 셀 적분 전력으로 계산한 (ΣPᵢ)²/Σ(Pᵢ²/Aᵢ)다.
- `donors.csv.gz`는 각 열원의 상류 위치·전력·첫 도달면과 좌표, 입사각, 평행 경로 길이, 물리 폭을 담는다. 선반 판정에는 선분(chord)의 법선과 별도로 기록한 해석 spline 법선을 사용한다. 기준/후보의 SHA-256으로 파일을 검사할 수 있다.

모델의 범위는 고정 평형에서 지정된 열원과 수송식으로 예측한 정상상태 표면 열유속이다. 입력 파일과 시드가 같아도 다른 장비에서는 시한 안에 완료한 평가 횟수가 달라질 수 있다. 고정 형상을 비교할 때는 후보 ID·계수·입력과 수치 설정을 함께 사용한다.
