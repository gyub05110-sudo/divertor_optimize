# NF 방법 논문용 검증 설계
2026-11-05. 원격 기준 dc8ab43. 제공된 기준 벽·고정 평형·입력별 독립 열원과 공학 제약을 사용한다. 목적량은 전 벽의 최대 열유속이며 같은 전력에서 A₂를 함께 보고한다. 실험 보정과 FIREFLY 점수 학습은 포함하지 않는다.

1. **탐색 자체의 시험.** 공개된 과거 형상 세 개를 출발점으로 열원별 2차원 유한 설계공간(삼각 격자 한 변 6등분, 28점)을 사전에 고정한다. G1 등 원 제약과 전체 벽 충돌 검사를 유지한다. 최소 곡률반경은 양쪽 모두 20 mm. 320경로 GPU로 모든 점을 평가하고 유효 상위 3개만 같은 입력의 cubic1280으로 다시 평가한다. 선택한 뒤 모든 28개를 같은 cubic1280으로 별도 평가해 격자 내 최저값 대비 초과율, 적격 후보 수, 기준 대비 피크와 A₂, 준비 포함 시간과 캐시 준비 후 시간을 보고한다. 정확히 같은 정밀 평가를 모두에게 적용한다. 1% 초과율은 이 시험의 설계 판단 기준이며 수학적 오차상한이 아니다. cubic1280은 같은 모델의 수치 참조다. 전역 최적이나 완전한 다중해상도 적격 인증을 주장하지 않는다. 누락/실패 표본은 표시한다.
2. **재사용성.** 기준값 1 MW, n=10¹⁹ m⁻³, χ=0.3 m²/s, λq=3 mm와 단일 변화 χ=0.6 m²/s로 각각 재탐색한다. OMP·IMP는 개별 부하다. 고정 입력을 명시하는 원 캠페인 CLI와 스냅샷/재시작 계약을 검증한다.
3. **어느 물리 가정이 위험한가.** 마지막 횡방향 변위가 자기선 라벨에 생긴 후 타깃에 도달하는 통제된 수송 한계를 사용한다. 경사·곡면에서 단위 전력의 해석적 첫 충돌 면적분포와 현 벽면 Gaussian의 분포·순위를 대조한다. 실제 DIII-D 계산에 이 한계의 가정을 강제하지 않는다. 전력 보존만으로 통과시키지 않는다.
4. **외부 보강.** 이미 공개된 FIREFLY 결과를 양·음 사례 모두 그대로 사용한다. 같은 입력 정의/배경 차이/격자 의존성을 표기한다. 현 원고가 주장할 수 있는 성능과 추가 자원으로 검증할 부분을 구분한다. 이후 수송식을 바꿀 경우 미사용 형상·운전점에서 사전 고정 재시험한다.

## 재현

PowerShell에서 저장소 루트와 기존 Python/CuPy 의존성을 사용한다. 다음 출력 이름은 새 경로로 지정한다. 기존 결과에 덮어쓰지 않는다.

```powershell
$study = 'D:/DivertorOptimization/studies/nf-reproduction'
New-Item -ItemType Directory -Path $study
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
python -m verification.nf_search_benchmark --output "$study/grid"
python -m verification.nf_finalist_sensitivity --grid "$study/grid" --output "$study/sensitivity"
python -m verification.nf_width_ablation --grid "$study/grid" --output "$study/width-ablation"
python -m verification.nf_operating_check --output "$study/operating"
python -m verification.nf_incidence_limit --output "$study/inclination.json"
python -m verification.nf_branch_coverage --shot input/D3D/d3d/d3d_shot_1a3155370c.parquet --evidence D:/DivertorOptimization/analyses/submission-adversarial-20261105/audit.json --source-csv D:/DivertorOptimization/studies/firefly-baseline-20261105/OMP.csv --output "$study/coverage.json"
python run_global_campaign.py --mode create --source IMP --hours .45 --population 12 --projected-seeds 2 --slice-seconds 5 --search-seconds 120 --validation-reserve-seconds 120 --seed 11052026 --backend gpu --power-mw 2 --density 1.3e19 --chi .6 --output "$study/campaign-2MW"
```

같은 자기장·같은 λq로 운전점을 여럿 평가할 때는 `prepare_boundary`의 기준 객체와 새 `ReferencePoint`를 `d3d_optimization.operating_sweep.reuse_boundary`에 전달한다. 반환한 기준 객체로 기존 `evaluator`와 `evaluate_candidate`를 호출한다. 평형·벽·격자·열원 폭이 달라지면 새 기준 객체를 준비한다. `verification.nf_operating_check`가 2 MW에서 독립 재계산과 이 경로의 동등성을 확인한다. 캠페인 생성 시 지정 가능한 물리 인자는 `--power-mw`, `--density`, `--chi`, `--lambda-q-mm`, `--gamma`, `--kappa0`이며 이미 생성된 캠페인의 운전점은 고정된다.
