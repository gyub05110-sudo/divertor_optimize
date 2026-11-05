# Divertor geometry optimization

주어진 DIII-D 기준 형상과 자기 평형, 입력 전력·밀도·횡방향 확산계수에서 **도착 열유속을 넓게 분산하여 peak를 낮추는 형상**을 빠르게 탐색합니다. 실제 측정값과의 일치는 이 프로젝트의 최적화 목표나 실행 조건이 아닙니다.

탐색 정확도는 같은 수송식의 정밀 계산에서 얻는 최저 피크와 후보의 차이로, 속도는 그 해를 얻기까지의 준비·탐색·재평가 시간으로 평가합니다. FIREFLY는 고정한 후보를 별도 수송식으로 다시 계산해 개선 방향과 적용 범위를 확인하는 보강 자료입니다. [NF 방법 검증 결과](results/nf-method-20261105/REPORT.ko.md) · [설계와 재현 방법](docs/nf-method-design.ko.md)

원래 캠페인과 후속 시험의 모든 형상은 제공된 DIII-D 벽을 기준으로 정의됩니다. 일부 후속 캠페인은 그 기준에서 이미 탐색한 형상도 초기 후보에 포함했습니다. 28점 검증은 과거 형상 사이에서 정의한 유한 설계 영역의 시험입니다. [기준과 초기 후보의 출처](results/nf-method-20261105/REPORT.ko.md) · [최초 도달면 검증](results/preimpact-firstwall-20261105/REPORT.ko.md) · [FIREFLY 보고서·연구 원자료와 복원 방법](results/research-archive-20261105/README.md)

같은 DIII-D 벽에서 과거 후보를 읽지 않고 생성한 72개 형상으로 200초 예산을 비교했습니다. 여섯 비용 재생 시험에서 선별과 직접 탐색은 같은 형상을 골랐고, 추가 OMP 실시간 시험의 기준 대비 피크 감소율은 각각 2.10%와 2.00%였습니다. [동일 예산 시험과 전체 원자료](results/d3d-base-budget-20261105/REPORT.ko.md)

최종 후보에는 기준벽보다 피크가 낮은지, 인접 후보까지 구별할 수 있는지, 국소 열유속이 수렴했는지를 따로 판정합니다. OMP20은 4개 수치 설정과 두 횡확산계수에서 기준벽의 피크를 45.64–47.26% 낮췄지만, 1% 국소 분포 기준을 넘었습니다. [최종 후보 재검증·추가 평형에서의 선별 결과](results/nf-finalist-robustness-20261106/REPORT.ko.md)

이 저장소는 `freegs_divertor_no-transport` 작업 폴더의 2026-10-29 코드와 작업 기록을 보존합니다. 폴더 이름 및 과거 문서와 달리, 최신 실행 경로에는 연결 길이에 의존하는 수송 폭과 전력 보존형 경계 분배 모델이 포함됩니다. 과거 직접 매핑 모델도 그대로 남아 있습니다.

## 2026-11-04 OMP 판정 재검수와 공개 뷰어

[공개 형상·열분포 뷰어](https://d3d-divertor-shape-explorer.sorainpublic.chatgpt.site)는 로컬 서버 없이 열립니다. 현재 뷰어와 데이터 전체를 [`viewers/d3d-progression-20261104/`](viewers/d3d-progression-20261104/)에 보존했습니다. 저장소 루트에서 `python -m http.server 8781 --bind 127.0.0.1 --directory viewers/d3d-progression-20261104/dist`로도 열 수 있습니다.

OMP의 argmax 위치 2 mm를 필수 조건에서 제외한 방향은 유지하되, v3 1%/.5%는 잠정 수치 허용오차로 구분했습니다. 20 mm 18단계는 계산상 peak 47.26% 감소·A2 53.23% 증가이지만 국소 q 차이 약 2.693%가 남아 v3 미통과입니다. 50 mm 34단계는 peak 1.786940489 MW/m², 34.56% 감소로 v3 통과입니다. 기존 원 캠페인 최종 후보와 이번 재검증 대상인 보관 중간 후보는 서로 다릅니다. 이번 감사는 합격 여부나 최종 형상을 바꾸지 않았습니다.

[작업·배포 보고서](docs/publication/d3d-progression-20261104/REPORT.ko.md) · [수치 재검수 상세](results/omp-policy-audit-20261104/REPORT.md). 원시 NPZ·계산 캐시는 기존 D 드라이브 study에 보존했으며, 웹 뷰어 실행에는 필요하지 않습니다. 아래 날짜별 내용은 당시 결과 기록입니다.

## 2026-11-03 우열 판정 조기 종료

최종 후보 비교에 공유 시간 예산과 조기 종료를 추가했습니다. 같은 320 donor 입력의 새 대조 실험에서 기존 방식은 335.99초에 4/5개 비교를 판별했고, 새 방식은 173.65초에 5/5개를 판별했습니다. replay는 10회에서 6회로 줄었습니다. 동일한 구간 폭이 아니라 엄격한 상·하한 분리로 후보 우열을 판별하는 기준의 비교이며, 전체 형상 탐색 가속률을 뜻하지 않습니다. 기존 탐색 루프는 유지합니다. [실험 범위와 비용 보고서](results/adaptive-peak-comparison-20261103/REPORT.ko.md).

## 2026-11-03 후보 재검증과 선택적 peak 검증

현재 형상 정책으로 기준 형상과 저장 후보를 같은 조건에서 재검증했습니다. IMP 20/50 mm는 각각 2.375184/2.481554 MW/m²(기준 대비 6.84%/2.67% 감소), 과거 OMP 50 mm는 1.995853 MW/m²(26.91% 감소)로 해상도·보간·벽 격자 검사를 통과했습니다. 기존 캠페인의 OMP 50 mm는 3.003868 MW/m²로 기준보다 높습니다. 아래 과거 OMP 20 mm 결과는 이후 강화된 전체 벽 검사에서 탈락했으므로 현재 유효 최선값으로 사용하지 않습니다.

빠른 탐색 경로는 유지합니다. 새 `verification.run_selective_peak_comparison` 및 `verification.run_selective_transport_study`는 소수 최종 후보용 별도 실행이며 기본 탐색에 연결하지 않았습니다. 고정된 수치 donor 분포의 공간 peak 판정과 연속 경로 오차를 포함하는 구간은 서로 다른 보증입니다. 연속 경로 구간은 줄었지만 아직 넓고, 전역 최적성은 미증명입니다. 실행 예시, 실제 비용과 한계는 [이번 보고서](results/selective-candidate-validation-20261103/REPORT.ko.md)에 기록합니다.

## IMP 접선 제약과 장시간 전역 탐색

`run_global_campaign.py`는 IMP 경사면–수직면 접합부를 공통 접선으로 연결하고, 외측 고정면과의 두 접선도 유지합니다. 공통 접선 방향과 IMP 계수 8개를 함께 탐색합니다. OMP는 이전 접선·개구부 제약을 유지합니다. 두 영역은 각각 독립적인 1 MW 실험이며, 각 영역의 최소 곡률반경 20·50·100 mm를 비교합니다. 원본 기준 형상은 그대로입니다.

```powershell
python run_global_campaign.py --output D:/DivertorOptimization/studies/imp-omp-global-20261029 --hours 72
python summarize_global_campaign.py --output D:/DivertorOptimization/studies/imp-omp-global-20261029
```

첫 명령은 코드·입력·PyArrow 복사본과 해시를 저장한 후 숨겨진 IMP/OMP 프로세스를 시작합니다. 두 프로세스는 RTX 3090의 FP64 열유속 분포 합산을 사용하고, 형상·첫 충돌·열적분·전력 CDF 계산은 CPU에서 수행합니다. CPU/GPU 분포 배열 일치를 확인한 뒤 진행합니다. 각 조건은 64개 집단의 differential evolution, 독립 재시작, 전역 범위 초기값, 양방향 변이와 국소 개선을 사용합니다. 전역 최적성을 증명하지 않습니다.

체크포인트는 약 2초마다 및 구간 전환 시 저장합니다. 프로세스 중단 후 `--mode launch --source IMP` 또는 `OMP`로 같은 출력 폴더를 지정하면 저장 상태에서 재개하며 종료 시한은 연장되지 않습니다. 같은 영역의 중복 쓰기는 OS 잠금으로 차단합니다. 원래 실행 파일·수치 라이브러리 버전·입력 해시가 달라지면 재개를 거부합니다. `STOP` 파일을 출력 폴더에 만들면 탐색을 마치고 정밀 검증으로 넘어갑니다. 정상 실행 중에는 중복 재시작하지 마세요.

72시간 중 마지막 1시간은 정밀 검증에 예약합니다. 640개 격자에서 후보를 다시 고른 뒤 320/640/1280개 격자, cubic/quintic 보간 및 벽 격자 세분화 검사를 통과한 후보만 최종 채택합니다. 곡률 조건별 후보를 번갈아 검증하고 이미 검증된 최선 후보는 재개해도 보존합니다. 조건을 만족하는 후보를 못 찾은 경우 제약을 완화하지 않고 미발견으로 기록합니다. 프로세스가 실행되는 동안 자동 시스템 절전을 억제하고 종료 시 해제합니다.

전체 결과·로그·체크포인트는 지정한 D 드라이브 폴더에 저장합니다. `summary.json/csv`, `REPORT.ko.md`, `qualified_geometry.png`, `search_progress.png`는 두 번째 명령으로 생성합니다. 결과 표의 정밀 peak와 탐색 중 거친 격자의 peak는 구분합니다.

2026-10-29 15:10 KST에 `D:/DivertorOptimization/studies/imp-omp-global-20261029`에서 본 탐색을 시작했습니다. 절대 종료 시한은 2026-11-01 15:10 KST입니다. 실제 중단/재개 시험에서 탐색 기록과 종료 시한 보존을 확인했고, 전체 테스트 197 passed / 1 skipped 및 원본 보존 검사를 통과했습니다. 1시간 간격 점검을 설정했습니다. 초기 후보의 CPU/GPU 분포 차이는 peak 기준 약 4.4e-16~6.1e-16이었고 초기 전체 후보 평가 가속은 약 1.35~1.44배였습니다.

[짧은 IMP 초기 탐색 보고서](results/imp-tangent-pilot-20261029/REPORT.ko.md): 기준 2.549663 MW/m²에서 20 mm는 2.404015 MW/m²(5.712% 감소), 50 mm는 2.530360 MW/m²(0.757% 감소)로 정밀 검증을 통과했습니다. 100 mm는 아직 유효 후보를 찾지 못했습니다. 짧은 시험은 정해진 18분 시한에 종료됐으며 일부 OMP 후속 검증은 시간 내 완료하지 못했습니다. 이 결과를 72시간 탐색의 최종 결과로 해석하지 않습니다. 본 탐색은 초기 후보를 재평가하여 활용하고, 새 결과는 D 드라이브 실행 폴더에 계속 저장합니다.

## 직전 OMP 제약 최적화 결과

OMP 왼쪽 접합부의 접선 연결, 최소 곡률반경, 오른쪽 개구부 간격을 함께 적용한 재최적화입니다. 원본 기준 형상과 IMP 영역은 보존했습니다. 중앙 바닥이 아래로 내려가는 후보도 탐색했습니다.

| 최소 곡률반경 제약 | peak [MW/m²] | 기준 대비 | 오른쪽 최소 간격 |
|---|---:|---:|---:|
| 20 mm | 1.446258 | 47.038% 감소 | 42.577 mm |
| 50 mm | 1.995853 | 26.912% 감소 | 44.781 mm |
| 100 mm | 3.111172 | 13.931% 증가 | 41.239 mm |

기준 peak는 2.730756 MW/m², 요구 간격은 33.354336 mm입니다. 세 후보 모두 형상·전력·격자·보간 검증을 통과했으며, 20·50 mm에서 개선을 확인했습니다. 접선 연속(G1)을 적용했고, 접합부 곡률 연속(G2)은 강제하지 않았습니다. 기록된 탐색 예산에서 찾은 결과이며 전역 최적을 증명하지 않습니다.

![접선·곡률·간격 제약 비교](results/omp-tangent-curvature-20261029/geometry_constraints.png)

- [재최적화 보고서](results/omp-tangent-curvature-20261029/REPORT.ko.md)
- [비교 데이터](results/omp-tangent-curvature-20261029/summary.csv)
- [검증 기록](results/omp-tangent-curvature-20261029/completion_checks.json)
- [설계·구현 계획](docs/superpowers/plans/2026-10-29-omp-tangent-curvature.md)

동일한 고정 입력으로 새 비교를 실행하려면 다음 순서를 사용합니다. 후보 선별·검증 스크립트는 현재 CuPy/CUDA 환경을 사용합니다. 전체 원본 실행 데이터와 코드/입력 스냅샷은 `D:/DivertorOptimization/studies/omp-tangent-curvature-20261029-v2`에 보존했습니다. 저장소 결과 폴더에서는 중복 입력 스냅샷과 테스트 사본을 제외했습니다. 테스트 원본은 `verification/`에 있습니다.

```powershell
python run_tangent_curvature_study.py --output D:/DivertorOptimization/studies/new-constrained-run
python qualify_tangent_candidates.py --output D:/DivertorOptimization/studies/new-constrained-run
python verification/validate_tangent_curvature.py --output D:/DivertorOptimization/studies/new-constrained-run
python visualize_tangent_curvature.py --output D:/DivertorOptimization/studies/new-constrained-run
```

## 이전 IMP·OMP 독립 최적화 — 추가 접선·곡률·간격 제약 전

각 열원에 **1 MW씩 별도로** 부여한 독립 실험입니다. IMP 실행은 IMP의 두 표면만, OMP 실행은 OMP의 네 표면만 바꾸고 반대 영역은 고정했습니다.

| 영역 | 기준 peak [MW/m²] | 검증 후 채택 peak [MW/m²] | 감소율 |
|---|---:|---:|---:|
| IMP | 2.549663 | 2.009104 | 21.201% |
| OMP | 2.730756 | 1.357740 | 50.280% |

고정 입력: n=1e19 m⁻³, χ⊥=0.3 m²/s, λq=3 mm, γ=7, κ₀=2000. 두 후보 모두 전력 보존, 독립 first-hit/열 상태 계산, 격자·보간 검사 및 고정 영역 불변 검사를 통과했습니다. 원래 탐색 최저 peak 후보 중 검증 기준에 미달한 결과도 함께 보존했습니다.

![최종 열유속 비교](results/regional-20261028/heat_flux_comparison.png)

- [최적화 보고서](results/regional-20261028/REPORT.ko.md)
- [IMP 최종 좌표와 계수](results/regional-20261028/final_geometry/IMP/)
- [OMP 최종 좌표와 계수](results/regional-20261028/final_geometry/OMP/)
- [요약 데이터](results/regional-20261028/summary.csv)
- [코너·협소 구간 분석](results/omp-corner-gap-20261029/REPORT.ko.md)
- [하향 변형 진단 데이터](results/omp-corner-gap-20261029/downward_probe.json)
- [전체 작업 기록과 다음 과제](docs/WORKLOG.ko.md)

## 설치 및 실행

검증한 환경은 Windows, Python 3.13.12입니다. 프로젝트를 D 드라이브에 clone하고 결과도 D에 저장하는 구성을 권장합니다.

```powershell
python -m venv D:/DivertorOptimization/venv
D:/DivertorOptimization/venv/Scripts/Activate.ps1
python -m pip install -r requirements-repro.txt
python -m pip install --no-deps --target .local_dependencies/d3d-fixed-equilibrium pyarrow==25.0.1
$env:MKL_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
python tools/verify_publication.py
```

두 번째 설치는 보존된 PyArrow 소스 파일의 해시 검사를 위해 필요합니다. `.local_dependencies`는 Git에 포함하지 않습니다. 다른 PyArrow 버전으로 대체하여 해시 검사를 무시하지 마세요.

CPU로 새 연구를 실행하려면 **아직 존재하지 않는 출력 폴더**를 사용합니다.

```powershell
python optimize_regional_transport.py --output D:/DivertorOptimization/studies/new-run --source both --backend cpu --generations 10
python optimize_regional_transport.py --output D:/DivertorOptimization/studies/new-run --source both --mode verify
python visualize_regional_transport.py --output D:/DivertorOptimization/studies/new-run
```

`--power-mw`, `--density`, `--chi`, `--lambda-mm`로 기준 입력을 지정합니다. `--source IMP` 또는 `OMP`로 한 영역만 실행할 수 있습니다. 최저 peak 후보가 정밀 검증을 통과하지 않으면 `qualify_regional_candidates.py`로 다른 후보를 선별하고 **선별된 후보에 전체 검증을 다시 수행**합니다. 선별만으로 개선 확인을 선언하지 않습니다.

GPU는 선택 기능입니다. RTX 3090에서 FP64 Gaussian 합산을 가속했으며, 같은 OMP 후보 3개 평가의 중앙시간 기준 약 1.46배 빨랐습니다. 전체 최적화 가속비는 아닙니다. 형상 교차·열 상태·CDF 전력·최종 검증은 CPU에서 수행합니다. [GPU 설치 및 재현 범위](docs/REPRODUCIBILITY.md)를 참고하세요.

```powershell
# CuPy/CUDA 환경을 설정한 경우
python optimize_regional_transport.py --output D:/DivertorOptimization/studies/new-gpu-run --source both --backend gpu
python -m pytest -q
```

GPU가 없는 환경에서는 GPU 전용 테스트 두 개를 제외할 수 있습니다.

```powershell
python -m pytest -q -k "not test_gpu_density_matches_cpu_narrow_broad_and_seam and not test_gpu_deposition_matches_preserved_cpu"
```

## 모델 범위와 현재 한계

- 자기 평형은 고정하고 타겟 형상을 바꿉니다. 기준 형상을 비교용으로 보존합니다.
- 연결 길이는 토로이달 감김을 포함합니다. 길이에 비례한 임의의 전력 손실을 주는 모델은 아닙니다.
- 열전도/sheath 기반 열 상태와 확산 폭을 이용하고, 연결된 벽의 면적 좌표에 Gaussian 전력을 분배합니다. 전체 2D 플라스마·중성입자 수송 방정식이나 방사 손실을 풀지는 않습니다.
- 최적화 목적은 전체 벽의 peak 감소입니다. 곡률반경·개구부 최소 간격 제약은 **아직 추가하지 않았습니다**. 기존 접합부의 코너도 남습니다.
- 두 결과는 독립 최적화입니다. 합친 형상에서 두 열원을 동시에 켠 최적해나 전역 최적해를 뜻하지 않습니다.

## 폴더 안내

| 경로 | 내용 |
|---|---|
| `d3d_optimization/` | 보존한 직접 매핑 및 추가 수송·GPU·영역 제한 모듈 |
| `verification/` | 회귀 테스트, 독립 물리량 대조 및 수치 검증 |
| `input/` | 실제 실행 입력과 과거 연구 결과 |
| `baseline_review/` | 기준 형상·자기장선 검토 자료 |
| `results/` | 2026-10-28 최종 연구와 2026-10-29 국소 형상 분석 |
| `docs/superpowers/plans/` | 설계 계획·검토 기록·원본 보존 해시 |
| `docs/publication/` | 이번 GitHub 배포의 파일 복사·검증 기록 |

과거 README는 [docs/README_LEGACY.md](docs/README_LEGACY.md)에 보존했습니다. 과거 문서의 '수송 없음' 설명은 당시 실행에 해당합니다. 기록에 남은 로컬 절대경로와 중복 snapshot 제외 범위는 [재현 안내](docs/REPRODUCIBILITY.md)에 명시했습니다.
