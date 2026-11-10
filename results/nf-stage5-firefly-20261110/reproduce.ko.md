# FIREFLY 고정 형상 비교 재현

공개한 [형상·시드 명단](wall-registration.json)은 4단계가 끝난 시점에 알려진 결과로 구성했다. [입력 연결 기록](bridge.json)은 자기장과 기하, 가중 열원, 하부 leg 착지를 검산한 결과다. [실행 봉인](execution.json)은 원벽과 후보의 입력 바이트, 함수 코드, 네이티브 라이브러리 환경, 시드별 실행 순서를 고정한다. [4단계 재현 안내](../nf-stage4-20261109/reproduce.ko.md)에는 이전 탐색 결과를 복원하는 방법이 있다.

## 원래 연구 컴퓨터에서 검산

Windows PowerShell에서 저장소 루트로 이동한 뒤 아래 변수를 설정한다. 연구 폴더는 D 드라이브에 있고, Linux FIREFLY는 WSL Ubuntu 안에서 11월 5일자 FLARE/FIREFLY 루트 파일시스템을 사용한다. 파일별 자기격자는 해당 평형의 원벽과 모든 고정 후보에 공통이다. 각 후보의 열 집계에는 그 후보의 물질벽을 쓴다. [실행 환경과 바이너리 지문](execution-environment.json)

```powershell
$py = 'C:/Users/FPL-user/miniforge3/python.exe'
$study = 'D:/DivertorOptimization/studies/nf-stage5-firefly-20261110-r2'
$paper = 'results/nf-stage5-firefly-20261110'
$env:PYTHONPATH = (Resolve-Path '.local_dependencies/d3d-fixed-equilibrium').Path
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:PYTHONIOENCODING = 'utf-8'
& $py -B -m verification.stage5_inputs --study $study --verify
& $py -B -m verification.stage5_bridge_report --study $study
& $py -B -m verification.stage5_dispatch check --study $study
& $py -B -m verification.stage5_metrics --study $study
& $py -B "$paper/summarize.py" --study $study --output $paper
```

입력 검산은 원래 parquet와 4단계 영수증에서 평형, 개별 열원, 형상을 다시 구성한다. 실행 검산은 276개 결과 파일 및 프로세스 영수증의 해시를 확인한다. 통계 분석은 원래 NetCDF 집계값을 읽어 벽의 고리마다 W를 복원한 다음 대조량과 신뢰구간을 다시 계산한다. `profiles/<평형>/<형상>.csv.gz`에서 각 고리의 면적과 위치, 물질면, 축소 모델 W, FIREFLY 시드별 W를 확인할 수 있다. [고리별 분석 코드](../../verification/stage5_metrics.py)

## 보존본 검증과 경로를 옮긴 재분석

Git의 [공개 보존본](evidence-compact/)은 시드별 원본 NetCDF 276개, 비교에 사용한 고리별 W, 두 실패 원본과 경계 수정 제어 결과, 코드와 영수증을 포함한다. 별도 full 보존본에는 원래 parquet와 4단계 후보 측정 배열도 들어 있다. 두 보존본은 큰 자기 격자 7개의 생성 입력과 크기·해시를 담고, 격자 원본은 연구 폴더에 보존한다. 해시와 공개된 원본 수는 [배포 기록](evidence-checksums.json)에서 확인할 수 있다.

```powershell
$compact = "$paper/evidence-compact"
& $py -I -B -S "$compact/stage5_evidence.py" verify --bundle $compact
$full = "$study/evidence-full"
$copy = 'D:/DivertorOptimization/tmp/stage5-restore-new'
$replay = 'D:/DivertorOptimization/tmp/stage5-replay-new'
& $py -I -B -S "$full/stage5_evidence.py" restore --bundle $full --out $copy
& $py -B "$copy/analysis-code/verification/stage5_portable_metrics.py" --restored $copy --out $replay
& $py -B -m verification.stage5_decision_check --study $study
```

복원과 재분석 출력은 새 디렉터리를 지정한다. 표준 라이브러리 검증기는 각 압축 조각과 논리 파일의 SHA-256을 대조한다. 경로 독립 재분석에는 NumPy와 h5py가 필요하다. [시드별 원시 전력에서 다시 계산한 45개 대조량](decision-crosscheck.json)에는 짝지은 전체 재표집과 정밀도 판정 일치 여부가 기록돼 있다.

## 중단된 원래 실행 재개

```powershell
& "$paper/run-study.ps1" -Study $study -Python $py -Workers 2
```

위 스크립트는 현재 코드, 입력, 준비 기록을 봉인과 대조하고 영수증이 있는 슬롯은 원본 바이트를 검사한다. 누락된 슬롯을 원래 시드로 계산하고 마지막에 다시 검산 및 분석한다. 각 시드의 MPI 4랭크가 합계 4,096개 이력을 추적한다. 실행 시간에는 각 프로세스의 준비·계산 시간이 포함된다. 4단계에서 사용자가 형상을 얻는 데 적용한 600초 시간 한도는 별도 기록이다.

## 입력 연결을 새로 구성하는 순서

원본 루트에서 `verification.stage5_campaign register`, `verification.stage5_inputs`, `verification.stage5_native smoke` 순으로 실행한다. 원본 `smoke`는 11월 5일자 원래 라이브러리와 기준 NetCDF의 비트 일치 검사용이다. `verification.stage5_native prepare --case <평형 ID>`를 각 평형에 적용하고 `verification.stage5_bridge_report`를 실행한다. 각 모듈은 `--study`로 연구 디렉터리를 받는다. 대표 평형은 선반 열점 `3264fdd71f`와 바닥 열점 `70481becc0`이다.

수정본의 재현에는 별도 새 연구 디렉터리를 쓴다. 원본에서 검산한 `inputs/`, `wall-registration.json`, `native/cases/`, `native/historical_smoke/`를 복사한다. [반사 경계 수정 재현 안내](native-repair.md)에 따라 수정된 라이브러리와 세 제어 계산을 만들거나 보존본의 `native-physics/`를 복원한다. [재실행 환경](execution-environment-r2.json)을 새 디렉터리의 `execution-environment.json`으로 복사한다. 입력 연결을 다시 검산해 `bridge.json`을 만들고, 첫 형상 열계산 전에 `verification.stage5_dispatch seal`로 새 실행 기록을 만든다. [발표된 재실행 봉인](execution-r2.json)과 입력·분석 코드·열계산 라이브러리 지문을 대조할 수 있다. 봉인의 시각과 커밋은 해당 재현 시점의 값이다. WSL 진입 스크립트는 수정 실행마다 전용 마운트 네임스페이스를 열며, 각 MPI 프로세스는 버전 경로의 라이브러리를 사용한다. [진입 스크립트](../../verification/stage5_native_enter.sh) · [네이티브 검산 코드](../../verification/stage5_native.py)
