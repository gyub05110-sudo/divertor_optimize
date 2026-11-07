# DIII-D 평형 전향 시험의 재현

검색 코드와 규칙을 `86dd454`에서 고정했다. 원래 공급된 DIII-D 벽, OMP 최적화 끝점과 운전 입력을 공통으로 사용하고 공개 EFIT 프레임만 바꿨다. 각 익명 파일의 실제 방전에서 같은 벽·운전점이 측정됐다는 전제는 없다. 새 세 파일의 식별자는 `332a6d97fc`, `e806f47778`, `51d83022ca`다. 세 프레임의 시각은 각각 1,400 ms다. 입력 파일명에 남아 있는 `1a3155370c`는 기존 로더가 읽는 경로다. 각 실행에 실제 공급한 바이트의 출처는 `preflight.json`의 원본 경로와 `frozen.json`의 파일 해시로 확인한다. [고정된 공개 자료](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge/tree/1e280905b85f2a6fdde7e06fca8cf3a1edf447cb) · [입력 선별 및 원본 보관 기록](D:/DivertorOptimization/studies/nf-prospective-inputs-20261107/reproduce.ko.md)

## 환경

아래 명령은 PowerShell에서 수행한다. 연구 코드가 고정한 Python 패키지 버전과 장비는 `development-plan.json`과 `frozen.json`에 들어 있다. 실제 비교는 단일 RTX 3090을 사용했다. 시한은 자식 프로세스 기동부터 재고 매 실행 전용 GPU 컴파일 캐시와 임시 폴더를 비운 상태에서 시작한다. 각 방법의 파일 시스템 읽기에는 OS 페이지 캐시가 영향을 줄 수 있다.

```powershell
$repo = 'D:/DivertorOptimization/worktrees/d3d-report-correction-20261104'
$py = 'C:/Users/FPL-user/miniforge3/python.exe'
$study = 'D:/DivertorOptimization/studies/nf-prospective-search-20261107'
$inputs = 'D:/DivertorOptimization/studies/nf-prospective-inputs-20261107'
Set-Location -LiteralPath $repo
$env:PYTHONPATH = "$repo/.local_dependencies/d3d-fixed-equilibrium"
$env:DIVERTOR_GPU_RUNTIME = 'D:/DivertorOptimization/runtime/cupy'
$env:TEMP = $env:TMP = 'D:/DivertorOptimization/tmp'
$env:CUPY_CACHE_DIR = 'D:/DivertorOptimization/cache/cupy'
$env:CUDA_CACHE_PATH = 'D:/DivertorOptimization/cache/cuda'
$env:MPLCONFIGDIR = 'D:/DivertorOptimization/cache/matplotlib'
$env:HF_HOME = 'D:/DivertorOptimization/cache/huggingface'
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:PYTHONIOENCODING = 'utf-8'
$env:OPENBLAS_NUM_THREADS = $env:MKL_NUM_THREADS = $env:OMP_NUM_THREADS = '1'
```

## 저장된 결과 검산

입력 검사는 이전에 보관한 18개 파일의 3,903개 프레임과 새 자료의 배열을 재검사한다. 캠페인 검사는 봉인한 코드·물리 입력·작업 명세·부모 로그·원 이벤트를 읽어 각 시점의 실제 정밀 결과를 재구성한다. 검산할 때 물리 계산을 다시 시작할 필요는 없다. 보고서를 새 폴더에 생성하여 기존 파일과 비교할 수 있다.

```powershell
& $py -B -m verification.nf_prospective_inputs --study $inputs --phase verify
& $py -B -m verification.nf_prospective_campaign --study $study --action audit
& $py -B -m verification.nf_prospective_report --study $study --output 'D:/DivertorOptimization/studies/nf-prospective-report-replay'
```

## 실제 검색과 사후 물리 계산 재현

새로 검색할 때는 이전 입력의 원본 바이트를 검증하고 비어 있는 별도 D 연구 경로를 사용한다. 2026-11-06 개발 입력 해시와 비교하는 코드는 당시에 봉인한 `nf-anytime-fullspace-20261106/frozen.json`도 읽는다. 신규 입력을 새로 취득할 때는 링크한 입력 안내의 등록·선별 순서를 사용하고, `--phase verify`가 성공하면 그 JSON 출력을 새 입력 폴더의 `verification.json`으로 UTF-8 저장한다. 각 파일의 SHA와 예약된 시드·작업 순서를 바꾼 실행에는 새 등록 문서가 필요하다.

```powershell
$repeatStudy = 'D:/DivertorOptimization/studies/nf-prospective-search-repeat'
& $py -B -m verification.nf_prospective_campaign --study $repeatStudy --action freeze-development
& $py -B -m verification.nf_prospective_campaign --study $repeatStudy --action development
& $py -B -m verification.nf_prospective_campaign --study $repeatStudy --action freeze-holdout --holdout-source $inputs
& $py -B -m verification.nf_prospective_campaign --study $repeatStudy --action holdout --holdout-source $inputs
& $py -B -m verification.nf_prospective_campaign --study $repeatStudy --action audit
& $py -B -m verification.nf_prospective_numerics --study $repeatStudy --lock-only
& $py -B -m verification.nf_prospective_numerics --study $repeatStudy --score
& $py -B -m verification.nf_prospective_report --study $repeatStudy --output "$repeatStudy/report"
```

개발 단계의 두 설정이 사전 조건에 미달하면 `freeze-holdout`에서 멈춘다. 고정 시드라도 시간 제한이 있는 탐색에서는 하드웨어 부하, JIT 컴파일과 디스크 시간에 따라 시한 안에 마친 후보가 달라질 수 있다. 원기록을 읽어 얻은 선택의 재현과 새 벽시계 실행의 재현을 구분한다.
