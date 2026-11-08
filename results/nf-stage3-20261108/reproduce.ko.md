# 3단계 자료 재현

[최초 비교계획](protocol.ko.md)은 커밋 `b75c2a9`에 있다. [입력 등록](inputs/registration.json)은 커밋 `5530789`에 있다. [실행 환경](execution-environment.json)과 [개발 예 시운전](commissioning.json)은 새 평형 열계산 전에 작성했다. 연구 방향을 다시 선택할 때 사용할 이전 체크포인트는 `codex/selective-candidate-validation`의 `8e63e34`, 두 개발 평형 비교는 `codex/d3d-peak-input-and-shape-study`의 `2c5e0d2`다.

## 원자료와 측정

공개 원파일의 고정 revision, 개별 다운로드 URL, 크기와 SHA-256은 [선택 목록](inputs/selection.json)에 있다. 같은 바이트를 되받으면 원래 익명 파일명으로 사용할 수 있다. 이전 입력 24개 원파일에서 얻은 프레임 지문과 새 입력 12개의 1400 ms 프레임을 비교한 절차는 [등록 코드](../../verification/stage3_inputs.py)에 있다. 등록 재연에는 이 코드가 참조하는 저장소의 이전 연구 증거도 필요하다. [공개 데이터셋의 고정 revision](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge/tree/1e280905b85f2a6fdde7e06fca8cf3a1edf447cb)

원실행은 Windows에서 아래 순서로 수행했다. `$study`는 저장소 밖의 증거 디렉터리다. 봉인된 `expected.json`과 캠페인 실행기는 등록 JSON의 정확한 바이트, 절대 경로, 인터프리터 및 코드에 결합된다. 원실행 디렉터리에서 중단된 작업은 마지막 명령으로 이어 간다. 다른 위치나 새 다운로드에서 만든 등록 기록은 시각·경로가 달라지므로 보존된 등록 해시를 재생성할 수 없다. 다른 컴퓨터에서 원자료를 검증하거나 탐색을 새로 실행하는 방법은 아래에 있다.

```powershell
$py = 'C:/Users/FPL-user/miniforge3/python.exe'
$study = 'D:/DivertorOptimization/studies/nf-stage3-20261108'
$env:PYTHONPATH = (Resolve-Path '.local_dependencies/d3d-fixed-equilibrium').Path
$env:DIVERTOR_GPU_RUNTIME = 'D:/DivertorOptimization/runtime/cupy'
& $py -B -m verification.stage3_inputs freeze --output "$study/inputs"
& $py -B -m verification.stage3_inputs acquire --output "$study/inputs"
& $py -B -m verification.stage3_inputs verify --output "$study/inputs"
& $py -B -m verification.stage3_campaign seal --study $study
& './results/nf-stage3-20261108/run-study.ps1' -Python $py -Study $study -GpuRuntime $env:DIVERTOR_GPU_RUNTIME
```

`run-study.ps1`은 시간 제한 탐색, 실행 증거 확인, 시간별 결과 집계, 고정 형상 재평가 순서로 실행한다. 한 단계가 실패하면 종료한다. 검산만 따로 할 때는 아래 명령을 사용한다.

```powershell
& $py -B -m verification.stage3_campaign check --study $study
& $py -B -m verification.stage3_report --study $study
& $py -B -m verification.stage3_replay --study $study
```

`runs/<입력>/<시드>/<방법>/`의 `result.json`은 기준벽과 선택 형상, `events.jsonl`은 각 제안·정밀 평가의 실제 완료 시각을 보존한다. `attempts/`에는 최초 호출 인자와 중단 기록이, `observations/`에는 작업자 생성부터 종료까지의 외부 경과가 남는다. `replays/`는 시간 제한 탐색이 끝난 뒤 고정한 형상의 결과다. 엄밀한 전체 승리가 발생하지 않은 쌍도 같은 표에 남는다.

## 원고 그림과 사후 검사

저장소의 `results/nf-stage3-20261108`에서 생성하는 그림과 표는 완료된 정밀 계산만 읽는다. 사례별 열유속 그림은 개발 단계와 같은 물질좌표 및 셀 평균 집계 함수를 사용한다. 연결부와 바닥 열점 유형별 첫 적격 파일, 시드 20711을 사용한다.

```powershell
$paper = 'results/nf-stage3-20261108'
& $py -B "$paper/plot_anytime.py" --study $study
& $py -B "$paper/regional_power.py" --study $study
& $py -B "$paper/plot_example_profiles.py" --study $study
& $py -B "$paper/input_diagnostics.py" --study $study
& $py -B "$paper/fixed_wall_diagnostics.py" --study $study
& $py -B "$paper/tabulate_ablation.py" --summary "$study/replays/summary.json"
```

## 다른 컴퓨터에서 보존 자료 확인

아래 검증에는 Python 3.11 이상 표준 라이브러리와 `evidence/`의 전체 분할 파일이 필요하다. 원래 D: 경로나 GPU가 없어도 실행할 수 있다. 실제 파일과 해시를 원자료까지 대조하려면 [선택 목록](inputs/selection.json)의 고정 URL에서 12개 parquet를 원래 파일명으로 받은 후 `restore`에 `--raw-dir`를 넘긴다. 원자료의 크기·SHA가 틀리면 복원을 중단한다. 과거 원파일과의 중복 비교 근거는 프레임 지문으로 포함한다.

```powershell
$paper = 'results/nf-stage3-20261108'
$bundle = "$paper/evidence"
$sha = (Get-Content -Raw "$paper/evidence.sha256").Trim()
& $py -B "$paper/package_evidence.py" verify --bundle $bundle --manifest-sha256 $sha
& $py -B "$paper/package_evidence.py" restore --bundle $bundle --manifest-sha256 $sha --out 'D:/reproductions/nf-stage3-original'
```

고정한 원래 경로에 모든 자료가 있고 자체 원실행에서 증거 묶음을 다시 만들 때에는 `build --study $study --bundle <새 디렉터리>`를 사용한다. 보존본 검증은 새 최적화의 벽시계 성능을 측정하는 절차와 독립적이다.

검증한 분할 파일은 19개(총 863,740,662바이트)이며 2,256개의 논리 파일을 복원한다. 12개 원본 parquet는 고정 외부 참조다. 시운전·실행환경·소스 스냅샷, 108개 슬롯, 중단 여부와 등록 제외, 각 정밀 결과의 기하·열분포·도달 경로, 99개 고정 형상 재계산 측정까지 연결된다. 최종 회귀시험은 844개 통과, 기존 조건부 시험 1개 건너뜀, 하위 시험 12개 통과였다. [manifest 해시](evidence.sha256) · [manifest 내용](evidence/evidence-manifest.json)

## 다른 환경에서 파일 하나를 다시 최적화

같은 코드의 공개 CLI에 내려받은 실제 파일 경로를 넘긴다. 다음은 한 시드에서 세 방법을 순차 실행하는 예다. 각 출력에는 새 경과와 새 결정이 기록되므로 원 측정과 구분한다. GPU 설치와 지원 파일 형식은 [사용자 안내](../../docs/user-optimizer.ko.md)에 있다.

```powershell
$eq = 'D:/data/d3d_shot_64058989a5.parquet'
$base = @('--equilibrium', $eq, '--time-ms', '1400', '--budget-seconds', '600',
          '--seed', '20711', '--backend', 'gpu', '--radius-mm', '20', '--power-mw', '1',
          '--density-m3', '1e19', '--chi-m2-s', '0.3', '--lambda-q-mm', '3', '--gamma', '7', '--kappa0', '2000')
& $py -B optimize_divertor.py @base --family floor-only --method screened --output 'D:/reproductions/case-floor'
& $py -B optimize_divertor.py @base --family full --method screened --output 'D:/reproductions/case-full'
& $py -B optimize_divertor.py @base --family full --method direct --output 'D:/reproductions/case-direct'
```
