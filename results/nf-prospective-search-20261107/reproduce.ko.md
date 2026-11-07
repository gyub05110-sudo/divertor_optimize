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
& $py -B -m verification.nf_prospective_report --study $repeatStudy --output 'D:/DivertorOptimization/studies/nf-prospective-search-repeat-report'
```

개발 단계의 두 설정이 사전 조건에 미달하면 `freeze-holdout`에서 멈춘다. 고정 시드라도 시간 제한이 있는 탐색에서는 하드웨어 부하, JIT 컴파일과 디스크 시간에 따라 시한 안에 마친 후보가 달라질 수 있다. 원기록을 읽어 얻은 선택의 재현과 새 벽시계 실행의 재현을 구분한다.

## 원벽 실패 및 고정 선반 끝점

입사각 입력 진단은 후보나 표적 열유속 없이 실행할 수 있다. 실제 세 입력의 실행 순서는 아래와 같다. `below_floor_if_powered`는 고정 끝점 근방까지 전력이 도달할 때 1.5° 하한을 위반한다는 뜻이다. 표적의 차폐와 최초 도달면이 바뀔 수 있으므로 이 표시는 전체 형상 공간의 적격성 판정으로 사용하지 않는다.

```powershell
& $py -B -m verification.nf_shelf_preflight --input-root "$inputs/cases/332a6d97fc/input-root" "$inputs/cases/e806f47778/input-root" "$inputs/cases/51d83022ca/input-root"
```

캠페인 감사를 마친 뒤 두 실패 평형의 실제 원벽 생성자를 같은 인자로 재실행했다. `register`는 18회 전체의 작업·종료·원 이벤트 해시에서 서로 다른 두 실패 대상을 고정하고, `replay`는 기존 예외를 보존하면서 폭별 직접 매핑·각도·도달한 단계를 기록한다. `nf_shelf_corner_check`는 실제 선반 기저 벡터와 cubic/quintic 자기장을 대조한다. 사후 진단의 결과는 600초 및 1800초 선택에 편입하지 않는다.

```powershell
& $py -B -m verification.nf_prospective_baseline_diagnostic --study $study --action register
& $py -B -m verification.nf_prospective_baseline_diagnostic --study $study --action replay
& $py -B -m verification.nf_shelf_corner_check --study $study
& $py -B "$study/paper_checks.py" --study $study --summary "$repo/results/nf-prospective-search-20261107/summary.json"
```

## 비용 및 증거 보관

부모 프로세스에서 잰 시간은 개발 3,651.64초, 전향 시험 6,884.14초, 네 설정 재계산 1,081.16초, 실패 진단 재현 139.61초였다. 전향 시험의 명목 예산 합계는 18,000초다. 12회가 원벽 준비에서 일찍 끝난 시간이 실측 합계에 반영됐다. `*-parent-0000.json`에 종료 시각·실측 시간·로그 SHA-256을 보관했다. 사후 수치 점수 44건 중 8건의 원벽 계산 실패를 알리는 점수기 종료 코드는 2였고, 완료된 보고서 상태는 `closed_with_exceptions`다.

`evidence.zip`은 종료한 개발·전향 작업의 원 이벤트·반환값·실패·사후 점수·진단과 기준 코드 스냅샷, 입력 원장과 채택한 입력 파일 트리를 담고 있다. 해제 경로는 `D:/DivertorOptimization` 기준이다. 공개 원본 다섯 개는 `evidence-manifest.json`의 `public_raw_files`에 고정 URL·길이·SHA-256을 기록했다. 로컬 D 원본은 입력 연구의 `raw` 폴더에 보존했다. 캐시·잠금 파일은 보관 대상에서 제외했다. 이전 개발 입력과 그 출처는 [이전 입력 연구](../unused-equilibrium-20261106/reproduce.ko.md) 및 [이전 600초 연구](../nf-anytime-fullspace-20261106/reproduce.ko.md)에 있다. 각 ZIP 멤버의 SHA-256을 압축 해제한 바이트로 다시 검사했다. [압축 파일](evidence.zip) · [멤버와 공개 원자료 목록](evidence-manifest.json) · [압축 파일 SHA-256](evidence.zip.sha256)
