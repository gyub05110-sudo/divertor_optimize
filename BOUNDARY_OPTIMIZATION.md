# 고정 기준 형상 + 지정 입력으로 divertor 최적화

이 도구의 목적은 **주어진 기준 형상과 사용자가 정한 가상 입력에서 빠르게 divertor 형상을 비교·최적화**하는 것입니다. 실제 측정값과의 일치는 목표·선행조건·합격 기준이 아닙니다. 원래 모델 및 과거 결과는 보존하고, 새 모델을 별도 명령으로 실행합니다.

## 실행

저장소 루트에서 PowerShell 또는 터미널로 실행합니다. `--output`에는 아직 없는 새 연구 폴더를 지정합니다.

```powershell
python optimize_boundary_transport.py --source OMP --output input/D3D/optimized/my-omp-study --power-mw 1 --density 1e19 --chi 0.3 --lambda-mm 3 --gamma 7 --kappa0 2000
python optimize_boundary_transport.py --mode verify --source OMP --output input/D3D/optimized/my-omp-study
python verification/plot_boundary_transport.py --folder input/D3D/optimized/my-omp-study/OMP
```

빠른 동작 확인은 `--smoke`, 기준 형상만 계산하려면 `--mode preflight`를 사용합니다. Smoke 결과는 개선 확정 판정에 사용할 수 없습니다. 입력을 바꾸는 실행마다 새 `--output`을 지정합니다. `--source IMP`도 지원하며 `both`는 두 개의 독립 연구를 실행합니다.

| 입력 | 단위 / 의미 |
|---|---|
| `--power-mw` | 선택한 소스 계산 구간에 부여한 MW. 전체 장치 P 또는 IMP/OMP 동시 분할이 아님 |
| `--density` | upstream 대표 전자밀도, m⁻³ |
| `--chi` | 유효 전자 횡방향 열 확산계수, m²/s |
| `--lambda-mm` | 소스 exponential 폭, mm. 과거의 다섯 폭에 한정되지 않음 |
| `--gamma` | sheath 열전달 계수 |
| `--kappa0` | q∥=−κ₀ T[eV]^(5/2)∂∥T[eV]에서의 전도 계수 |

Te는 P, n 및 전도/sheath closure에서 계산되는 출력입니다. 같은 closure에 별도의 Te를 독립적으로 강제하지 않습니다. 소스의 유한 계산 구간과 정규화는 `manifest.json`의 `metadata.support_mm`에 기록합니다.

각 실행 중 입력은 고정합니다. 후보별로 형상에 따른 첫 도착점, 연결 길이, 열 상태와 폭을 다시 계산합니다. 전력 보존을 맞추기 위해 기준 형상이나 입력을 바꾸지 않습니다.

## 계산 예산과 해상도

기본은 donor 320개, trace step 0.5 mm, 형상 노드 161개, receiving-wall step 1 mm, 열적 quadrature 64입니다. `--tubes`, `--step-mm`, `--nodes`, `--wall-step-mm`, `--quadrature`로 조절합니다.

탐색 예산은 `--sobol-power`(초기 표본 수=2^값), `--population`, `--generations`로 조절합니다. 기본은 64개 초기 표본 + 24개 population × 6세대이며 중복 제안은 다시 계산하지 않습니다. 24개 spline 계수의 범위는 기본 ±100 mm이며 실제 허용 여부는 기존 형상 제약 및 추가 연속 길이 검사로 결정합니다. `--seed`로 재현합니다. 과거 적합 형상을 warm start로 사용할 때도 같은 새 모델과 입력으로 다시 평가합니다.

최적화를 실행하면 우선 탐색 후보를 저장합니다. `--mode verify`는 기준/선택 형상을 320/640/1280 donor와 quintic 보간에서 비교하고 receiving-wall만 세분화한 별도 검증도 합니다. `validation.json`의 `verdict.improvement_confirmed`를 확인해야 수치 검증을 통과한 개선이라고 말할 수 있습니다. 이 과정은 전역 최적해를 증명하지 않습니다.

## 추가한 모델

원래의 field-line first-hit 및 단방향 two-point thermal 폭 계산을 사용합니다. 그 후 첫 도착점을 중심으로, 하나의 닫힌 연결 벽의 누적 회전면적 A=∫2πR ds 위에 주기 Gaussian footprint를 둡니다.

- σphys = σflux |∂Rψsource| / |∇ψhit|
- σA = 2πRhit σphys
- 셀 전력 = donor 전력 × 주기 Gaussian의 해당 면적 구간 CDF
- q(A) = 모든 donor의 전력 밀도 합

이 정의로 모든 도착 donor의 전력이 보존되며, 살아남은 일부 셀에 전력을 재정규화하지 않습니다. 미도착 donor는 미도착 장부에 남습니다. 오목한 벽을 가로지르는 직선으로 배분하지 않고 연결된 경계를 따릅니다.

**새 footprint 가정**이며 기존 kernel과 동등하지 않습니다. σphys를 유효 벽방향 폭으로 사용하고 1/sin 입사각 폭 투영은 적용하지 않습니다. 국소 자기장/입사각은 first-hit과 기하·입사각 제약에 남습니다. 벽 위 redistribution을 실제 plasma cross-field 궤적이나 완전한 수송 방정식 해로 해석하지 않습니다.

폭이 매우 클 때 전체 연결 벽에 균일한 열유속으로 접근합니다. `broad_kernel_power_fraction`, `power_beyond_0p1m_wall_distance`, donor→region 전력 배열로 이 비국소 효과를 확인할 수 있습니다. 평균 벽 이동 거리는 셀 중심 근사이며 `wall_distance_estimate_error_bound_m`에 오차 상한을 기록합니다.

## 지원 범위와 실패의 의미

양의 유한 P, n, χ, λ, γ, κ₀를 입력받습니다. float 표현 범위와 분해능, 유한 자기장 기울기, 유효한 하나의 단순 닫힌 양의 R 벽이 필요합니다. 매우 작은 폭, underflow/overflow, 0 기울기 같은 수치 불능은 입력을 몰래 바꾸지 않고 실패로 보고합니다. 입력이 매우 좁은 소스 또는 예리한 분포를 만들면 donor/trace 해상도를 높여 검증해야 합니다.

모든 입력에서 개선이나 적합 최적해가 존재한다는 보장은 아닙니다. 탐색에서 더 좋은 후보를 찾지 못하면 적합한 기준 형상을 유지합니다. 미배분이나 수치 실패는 모델/계산 조건의 문제로 기록하며 실제 기준 형상이 잘못됐다는 뜻으로 해석하지 않습니다.

현재는 고정 자기장과 단방향 전도/sheath closure입니다. 입자·운동량·중성입자·복사·재결합·시간 변화·횡수송의 온도 피드백은 계산하지 않습니다. 서로 분리된 벽 사이의 전달도 지원하지 않습니다. 다른 종류의 기준 형상/자기장 파일을 자유롭게 교체하는 일반 geometry importer는 이번 변경에 포함되지 않았으며 현재 저장소의 D3D 형상 구성과 설계 표면을 사용합니다.

## 결과와 재현

- `manifest.json`: 입력, 기준 형상 평가, 선택 spline 계수, 실행 설정, 시간
- `search/evaluations.jsonl`: 모든 고유 후보의 제약·목표값·장부
- `validation.json`: 독립 검증, 격자/보간 비교, 최종 판정
- `baseline_refined.npz`, `selected_refined.npz`: 벽, 열유속, 셀 전력, 도착점, σ, 열 상태, donor→region 전력
- `REPORT.ko.md`, `comparison.png`, `arrival_map.png`: 결과와 그림
- `provenance.json`: 실행 전에 기록한 원래 코드·새 모델·검증 코드·입력 파일 SHA-256

새 검증은 코드/입력 hash가 달라졌으면 중단합니다. 원래 Python 169개 파일의 hash는 `docs/superpowers/plans/2026-10-28-boundary-preservation.json`으로 검사합니다. 신규 시험은 `python -m pytest verification/test_boundary_transport.py verification/test_boundary_evaluator.py verification/test_boundary_diagnostic.py verification/test_boundary_study.py -q`로 실행합니다.
