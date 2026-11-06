# 재현 절차와 원기록 보존

이 안내는 저장된 판단을 다시 계산하는 감사와, 같은 설정으로 실제 탐색을 새로 실행하는 절차를 구분한다. 기준 연구는 `D:/DivertorOptimization/studies/nf-anytime-fullspace-20261106`이다. 이 보고서 폴더에는 작은 요약·그림·출처 자료를 모으고, 전체 원기록·코드 스냅샷·캐시는 D의 연구 폴더에 보존한다. 명령은 아래 저장소 루트에서 실행한다. 기존 연구의 `frozen.json`, `protocol.ko.md`, 입력, 이벤트, 봉인된 사후 자료를 수정해 재실행 조건을 맞추지 않는다.

**환경과 경로**

동결 기록의 환경은 Windows 11 10.0.26200, Python 3.13.12, NVIDIA GeForce RTX 3090 24,576 MiB, NVIDIA 드라이버 610.88이다. Python 실행 파일은 `C:/Users/FPL-user/miniforge3/python.exe`, 추가 평형 입출력 의존성은 저장소의 `.local_dependencies/d3d-fixed-equilibrium`, GPU 런타임은 `D:/DivertorOptimization/runtime/cupy`를 사용한다. C의 Python을 실행하더라도 연구 출력·임시 파일·가속기 캐시는 D에 둔다.

| 동결 패키지 | 버전 |
| --- | --- |
| NumPy / SciPy | 2.4.3 / 1.17.1 |
| PyArrow / Shapely | 25.0.1 / 2.1.2 |
| cupy-cuda12x | 14.2.0 |
| cuda-pathfinder / fastrlock | 1.8.4 / 0.8.4 |
| nvidia-cuda-nvrtc-cu12 | 12.9.86 |
| nvidia-cuda-runtime-cu12 | 12.9.79 |

전체 버전·장비·수치 튜플·운전 입력·파일별 SHA-256의 기준은 연구 폴더의 `frozen.json`이다. 이미 실행한 연구의 재개는 동결된 코드·입력·환경과 일치해야 한다. 다른 컴퓨터나 경로로 옮겨 실행하려면 원기록을 고치는 대신 검증한 입력과 별도 출력 폴더로 새 연구를 동결한다.

```powershell
$anytimeRepo = 'D:/DivertorOptimization/worktrees/d3d-report-correction-20261104'
$anytimePython = 'C:/Users/FPL-user/miniforge3/python.exe'
$anytimeStudy = 'D:/DivertorOptimization/studies/nf-anytime-fullspace-20261106'
$anytimeInputs = 'D:/DivertorOptimization/studies/unused-equilibrium-20261106'
$anytimeRuntime = 'D:/DivertorOptimization/runtime/anytime-reproduction'
Set-Location -LiteralPath $anytimeRepo
$env:PYTHONPATH = "$anytimeRepo/.local_dependencies/d3d-fixed-equilibrium"
$env:DIVERTOR_GPU_RUNTIME = 'D:/DivertorOptimization/runtime/cupy'
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:PYTHONIOENCODING = 'utf-8'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:TEMP = "$anytimeRuntime/tmp"
$env:TMP = $env:TEMP
$env:MPLCONFIGDIR = "$anytimeRuntime/matplotlib"
$env:CUPY_CACHE_DIR = "$anytimeRuntime/cupy"
$env:CUDA_CACHE_PATH = "$anytimeRuntime/cuda"
$env:HF_HOME = "$anytimeRuntime/huggingface"
foreach ($anytimeDir in @($env:TEMP, $env:MPLCONFIGDIR, $env:CUPY_CACHE_DIR, $env:CUDA_CACHE_PATH, $env:HF_HOME)) {
    $null = New-Item -ItemType Directory -Path $anytimeDir -Force
}
```

주 실행 감독자는 각 자식의 `TEMP`, `TMP`, `CUPY_CACHE_DIR`, `CUDA_CACHE_PATH` 등을 해당 `live/<case>/<seed>/<policy>/cache` 아래의 새 전용 디렉터리로 다시 지정한다. 사후 계산은 `posttrial-runtime/reference`, `posttrial-runtime/uncertainty`, `posttrial-runtime/screen-rejections`를 쓴다. GPU/CPU 계산을 다른 연구와 동시에 실행하지 않는다. OS 파일 캐시나 장비 전체의 GPU 상태를 초기화하는 절차는 없다. 여기서의 시작 조건은 새 프로세스와 빈 전용 CUDA/CuPy 디스크 캐시다.

**보존 자료와 출처**

| 위치 | 용도 |
| --- | --- |
| 주 연구의 `frozen.json`, `protocol.ko.md`, `code-snapshot/` | 주 실행의 규칙, 코드·입력 바이트와 해시 |
| `live/<case>/<seed>/<policy>/events.jsonl`, `result.json`, `process.json`, `stdout.log` | 완료·초과·실패 이벤트, 후보·초기 집단, 실제 프로세스 비용과 로그 |
| `audit.json`, `final-independent-audit.json` | 감독자 감사와 별도 구현의 원이벤트 재집계 |
| `reference/` | 패널별 후보·계수·해시, coarse/fine 전수 점수, live 후보와의 중복, 비용 |
| `uncertainty/` | 고정 후보 목록·선택 근거, 네 수치 설정의 재평가, 일치 검사와 비용 |
| `screen-rejections/` | 원기록에서 전부 추출한 coarse 열계산 후 탈락 후보의 별도 native 정밀 검사 |
| 세 사후 폴더의 `protocol.json`, `code-snapshot/`, `*.sha256` | 주 실행과 별도로 봉인한 사후 코드·규칙·산출물 |
| `amendment.json`과 이전 연구의 `frozen.json`, `interruption.json` | 12/12에서 생산 기본값 64/12로 정정한 이력 |

평형 원자료와 사례별 `input-root`는 `$anytimeInputs`의 `raw/` 및 `cases/`에 있다. [Sophelio/fusion-equilibrium-challenge](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge)의 [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) 배포본이며, 고정 revision은 `1e280905b85f2a6fdde7e06fca8cf3a1edf447cb`다. 원자료의 저작·배포 출처를 유지하고, 본 연구가 1400 ms 평형을 선택해 고정 DIII-D 벽과 가정한 운전 입력을 결합했음을 명시한다. 원자료의 측정 전력·밀도에 맞춰 운전점을 보정하지 않았다.

| 사용 사례의 고정 다운로드 URL | 바이트 수 | SHA-256 |
| --- | ---: | --- |
| [deb842cf1e.parquet](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge/resolve/1e280905b85f2a6fdde7e06fca8cf3a1edf447cb/data/diii_d_train/d3d_shot_deb842cf1e.parquet) | 6,898,095 | `24000d2a2fea8838e5226f63d7914e2403c95af5ac4c1d3c9cafd2570b29952d` |
| [fad1624516.parquet](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge/resolve/1e280905b85f2a6fdde7e06fca8cf3a1edf447cb/data/diii_d_train/d3d_shot_fad1624516.parquet) | 7,134,726 | `42725efdd7f426e212e47bc95b021666651866594b7e557209b3e177920514e6` |
| [4d14a35fcb.parquet](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge/resolve/1e280905b85f2a6fdde7e06fca8cf3a1edf447cb/data/diii_d_train/d3d_shot_4d14a35fcb.parquet) | 17,627,110 | `7d9a12aef9d23115fae9b367f128f0a1b7a20b934cee45fe5a94d453373c8b69` |

선정·제외 이력 전체와 원본 12개 파일의 URL·크기·해시는 원자료 연구의 `selection.json`, `downloads.json`, `preflight.json`, `protocol.ko.md`에 있다. 주 연구는 이 네 파일과 사례별 입력의 해시를 연결한다. 작은 보고서 묶음만으로 원자료와 전체 실행을 대체할 수 없다. 복원 시에는 위 고정 URL에서 받은 바이트를 기록과 대조하고, 기존 연구 감사에 필요한 원래 경로를 복원하거나 새 입력 연구·새 출력을 만든다.

사례별 입력 root를 새로 만들 때에는 [원자료 준비 절차](../unused-equilibrium-20261106/reproduce.ko.md)의 새 폴더 복사와 `python -m verification.nf_unused_inputs --study <새 D 원자료 폴더>`를 사용한다. 선정/제외를 재현하려면 선택된 세 파일만이 아니라 `downloads.json`에 열거한 12개 원자료를 복원하고, 새로 복사한 목록의 `local` 경로를 검증한 다운로드 위치로 맞춘다. 저장소의 원래 DIII-D 입력·벽 파일도 필요하다. 새 전처리 기록은 생성 시각·경로 때문에 기존 봉인 기록과 바이트가 다를 수 있으므로, 이렇게 재구성한 입력은 아래 `--source-study`에 지정해 **새 연구를 동결**한다. 기존 연구 감사를 위해 봉인된 해시를 고치는 용도가 아니다.

보고서 생성기는 `amendment.json`의 `prior_study`도 읽는다. 따라서 `D:/DivertorOptimization/studies/nf-anytime-20261106/frozen.json`과 `interruption.json` 및 기록된 해시 연결을 보존해야 한다. 이 이전 부분 실행은 비교 표본에 합치지 않으며, 초기 집단 설정 정정의 출처로만 사용한다.

**저장된 판단의 재현: 추가 열계산 없음**

독립 감사기는 원래 `events.jsonl`에서 200/400/600초 선택을 다시 계산하고 입력·코드 스냅샷·초기화·프로세스 시각·저장 선택의 연결을 검산한다. 새 감사 출력 이름을 사용한다. 아래 명령은 기존 연구나 봉인된 사후 자료를 덮어쓰지 않는다.

```powershell
$anytimeAuditOut = 'D:/DivertorOptimization/studies/nf-anytime-fullspace-reaudit-01.json'
& $anytimePython -m verification.nf_anytime_final_audit --study $anytimeStudy --out $anytimeAuditOut
if ($LASTEXITCODE -ne 0) { throw 'Independent audit failed; inspect its output' }
```

`--out` 파일이 이미 있으면 감사기는 덮어쓰지 않고 중단한다. 원래 감독자의 `--audit-only`는 `audit.json`을 다시 쓰므로, 사후 `live-evidence.json`이 이미 그 바이트를 참조하는 보존 연구에는 위 별도 감사 출력을 사용한다. 물리 계산을 새로 실행해 당시 벽시계 종료 시각을 재현하려고 하지 않는다. 동일 시드·장비에서도 운영체제와 실행 비용의 변동으로 시한 내 완료 후보 수와 최종 선택이 달라질 수 있다. 당시 선택의 재현은 보존한 완료 시각과 후보를 재집계하는 방식이다.

기존 연구의 세 사후 `summary.json`까지 보존되어 있다면, 아래 명령으로 추가 열계산 없이 보고서와 그림을 새 결과 폴더에 다시 만들 수 있다. 생성 시각 같은 메타데이터가 바뀔 수 있으므로 문서 파일의 바이트 동일성과 저장된 선택·수치의 재집계를 구분한다.

```powershell
$anytimeReReport = "$anytimeRepo/results/nf-anytime-fullspace-rereport-01"
if (Test-Path -LiteralPath $anytimeReReport) { throw 'Choose a fresh report directory under repository results' }
& $anytimePython -m verification.nf_anytime_report --study $anytimeStudy --out $anytimeReReport
if ($LASTEXITCODE -ne 0) { throw 'Stored-evidence report generation failed' }
```

**같은 고정 입력으로 새 27회 실행**

다음 예시는 위 원자료 연구가 보존된 같은 컴퓨터에서 새 출력으로 반복하는 절차다. `--source-study`는 세 평형의 이미 검증된 입력을 읽기 전용으로 재사용한다. 이전 `frozen.json`, `live/`, 캐시, 사후 결과는 복사하지 않는다. 복사하는 `amendment.json`은 채택한 64/12 규칙의 역사적 출처이며, 새 반복 중에 설정을 다시 변경했다는 뜻이 아니다.

```powershell
$anytimeRepeat = 'D:/DivertorOptimization/studies/nf-anytime-fullspace-repeat-01'
if (Test-Path -LiteralPath $anytimeRepeat) { throw 'Choose a fresh study directory' }
$null = New-Item -ItemType Directory -Path $anytimeRepeat
foreach ($anytimeName in @('protocol.ko.md', 'amendment.json')) {
    Copy-Item -LiteralPath "$anytimeStudy/$anytimeName" -Destination "$anytimeRepeat/$anytimeName"
}
& $anytimePython -m verification.nf_anytime_campaign --study $anytimeRepeat --source-study $anytimeInputs --freeze-only
if ($LASTEXITCODE -ne 0) { throw 'Main protocol freeze failed' }
& $anytimePython -m verification.nf_anytime_campaign --study $anytimeRepeat
if ($LASTEXITCODE -ne 0) { throw 'Main supervisor stopped; preserve and inspect its records' }
& $anytimePython -m verification.nf_anytime_final_audit --study $anytimeRepeat --out "$anytimeRepeat/final-independent-audit.json"
if ($LASTEXITCODE -ne 0) { throw 'Main evidence failed independent audit' }
```

등록된 세 평형 × 세 시드(20611,20627,20639) × 세 방법을 직렬 fresh process로 실행한다. 방법은 `coarse_adaptive`, `fine_adaptive`, `random_fine`이며 적응 방법의 설정은 population=64/projected_seeds=12다. 초기 집단의 범위·중심화 affine rank와 실제 시간 내 제안·평가 범위를 구분해 보존한다. 전체 비교에는 `--case`, `--seed`, `--policy` 필터를 사용하지 않는다. 중단 뒤 동일 연구를 이어갈 때에는 변경되지 않은 코드·입력으로 `--skip-completed`를 지정할 수 있다. 이 옵션은 감사 가능한 기존 완료·실패 기록을 보존하며, 실패를 성공할 때까지 반복하는 옵션이 아니다.

각 600초는 부모의 프로세스 시작 직전 monotonic 시각부터 임포트, 기준벽·자기선 준비, GPU 초기화, 후보 생성·탈락, 수송 계산, 첫 적격 결과의 CPU/GPU 전체 배열·피크 대조, 실행 중 로그 비용을 포함한다. 입력 취득·기하 사전 검사, 부모 해시/빈 캐시 폴더 준비, 최종 직렬화·보고서와 사후 계산 비용은 별도다. 600초 뒤 끝난 값은 보존하되 시간 내 선택에서 제외하고, 1200초 watchdog이 미완료 프로세스를 중단하면 실패·검열로 남긴다. 따라서 명목 27×600초와 실제 전체 소요시간은 같지 않을 수 있다.

**별도 패널·수치 검증과 최종 보고서**

이후 명령의 대상은 새 반복인 `$anytimeRepeat`다. 보존 연구의 현재 주 실행 기록을 수정하지 않는다. 주 실행의 규칙·코드 동결과 별도로 사후 스크립트 및 규칙을 먼저 봉인한다. `--freeze-only`는 추가 물리 계산을 하지 않는다. `uncertainty --select-only`는 완료된 주 실행에서 후보 목록과 근거를 잠그며 추가 열계산을 하지 않는다.

```powershell
& $anytimePython -m verification.nf_anytime_reference --study $anytimeRepeat --freeze-only
if ($LASTEXITCODE -ne 0) { throw 'Reference protocol freeze failed' }
& $anytimePython -m verification.nf_anytime_uncertainty --study $anytimeRepeat --freeze-only
if ($LASTEXITCODE -ne 0) { throw 'Uncertainty protocol freeze failed' }
& $anytimePython -m verification.nf_anytime_uncertainty --study $anytimeRepeat --select-only
if ($LASTEXITCODE -ne 0) { throw 'Numerical-check selection lock failed' }
& $anytimePython -m verification.nf_anytime_reference --study $anytimeRepeat
if ($LASTEXITCODE -ne 0) { throw 'Reference calculation stopped; inspect preserved evidence' }
& $anytimePython -m verification.nf_anytime_uncertainty --study $anytimeRepeat
if ($LASTEXITCODE -ne 0) { throw 'Numerical checks stopped; inspect preserved evidence' }
& $anytimePython -m verification.nf_anytime_rejections --study $anytimeRepeat --select-only
if ($LASTEXITCODE -ne 0) { throw 'Coarse rejection selection lock failed' }
& $anytimePython -m verification.nf_anytime_rejections --study $anytimeRepeat
if ($LASTEXITCODE -ne 0) { throw 'Coarse rejection checks stopped; inspect preserved evidence' }
& $anytimePython -m verification.nf_anytime_final_audit --study $anytimeRepeat --out "$anytimeRepeat/final-independent-audit-after-post.json"
if ($LASTEXITCODE -ne 0) { throw 'Post-calculation main-evidence audit failed' }
$anytimeReportOut = "$anytimeRepo/results/nf-anytime-fullspace-repeat-01"
if (Test-Path -LiteralPath $anytimeReportOut) { throw 'Choose a fresh report directory under repository results' }
& $anytimePython -m verification.nf_anytime_report --study $anytimeRepeat --out $anytimeReportOut
if ($LASTEXITCODE -ne 0) { throw 'Report generation failed' }
```

`reference`는 평형/시드마다 최대 24개 기하 적격 후보를 생성한다. 각 평형의 세 시드 패널을 해당 후보의 **추가** 열평가 전에 파일·해시로 봉인한 뒤 coarse/fine 전수 계산을 수행한다. 기준벽 준비에는 열계산이 포함된다. 같은 시드·생성기를 재사용하므로 패널 후보는 주 실행 또는 이전 연구에서 이미 제안·평가된 후보와 겹칠 수 있다. 중복 교집합을 기록하며, 이 패널을 독립 holdout 검증으로 부르지 않는다. 패널의 무작위 기대값은 공통 coarse 적격 모집단에서 같은 fine 요청 수를 비복원 균등 추출한 비교로, 실시간 `random_fine`이나 동일 시간 비교를 대신하지 않는다.

`uncertainty`는 기준벽, 시간 지점별 선택, 마지막 초과 평가 후보와 사전에 정한 관측 후보 내 근접 피크 비교쌍을 잠근다. 그 목록을 cubic640, cubic1280, quintic640, quintic1280에서 계산한다. native cubic1280의 재평가 일치는 실제로 이전 관측이 있는 후보·관측량에 한해 판정한다. 아직 관측하지 않은 후보나 기하 적격성만 확인한 경우를 열유속 재현으로 세지 않는다. 추가 결과로 주 실행의 선택이나 시간 내 성과를 바꾸지 않는다.

`screen-rejections`는 주 실행 후 적대적 검토에서 추가한 별도 검사다. 시간 내 완료된 새 coarse 침적 결과가 부적격인 경우를 전부 추출하며, 열유속 크기나 오류 유형으로 대상을 골라내지 않는다. 새 native cubic1280 계산 전에 목록을 봉인하고 해당 평형의 기준벽과 함께 계산한다. 이번 기록에는 입사각 제약으로 탈락한 세 후보/두 평형이 있었다. fine에서도 부적격인지, fine에서는 통과하는 오탈락인지, 추가 계산이 실패해 미확인인지 구분한다. 이는 원래 사전 등록한 600초 성과도, 전체 형상 공간의 오탈락률 추정도 아니다.

사후 명령은 봉인된 기존 결과를 확인하며 재사용하고, 남은 계산을 이어간다. 완료된 사후 폴더에서 같은 명령을 다시 호출하는 것은 모든 후보의 새 물리 재실행이 아니다. 새로운 반복 계산에는 위와 같이 새 연구 폴더를 쓴다. 코드·프로토콜·후보 목록·해시가 달라진 상태를 같은 연구에 섞으면 중단한다. 세 물리 사후 단계는 주 실행과 서로 동시 실행하지 않는다.

최종 전체 보고서는 `reference/summary.json`, `uncertainty/summary.json`, `screen-rejections/summary.json`이 보존된 뒤 생성한다. 보고서 명령은 주 실행과 봉인된 사후 근거를 읽어 `REPORT.ko.md`, `summary.json`, `anytime_curves.png`, `anytime_curves.pdf`를 지정한 저장소 `results` 폴더에 만든다. 독립 최종 감사기는 주 실행 원기록을 검산하고, 보고서 생성기는 사후 해시 연결과 저장된 판정을 다시 확인한다. 어느 명령도 새 물리 결과를 만들어 누락을 채우지 않는다.

비용은 주 실행의 `process.json`과 이벤트 구간, 사후 각 `runs/*.started.json`·`*.finished.json`, `preparation/`, 후보 `scores/`, 사후 `summary.json`의 `costs`에서 구분한다. 주 실행 600초 성능, 초과 및 프로세스 종료 비용, 패널 생성·전수 계산, 추가 수치 설정 검증을 각각 보고한다. 같은 고정 수송식의 유한 설정 범위는 연속체 오차상한이나 실험적 정확도 인증이 아니며, FIREFLY는 별도 물리 교차 검증 범위로 유지한다.

**압축 원기록과 검증 코드**

`evidence.zip`은 완료된 주 실행·사후 계산과 중단한 초기 설정의 원기록, 코드 스냅샷, 원자료 선정 메타데이터를 담는다. CUDA/CuPy 캐시·임시 파일·잠금 파일과 외부 parquet/사례별 입력 root는 제외한다. 파일별 바이트 수와 SHA-256은 `evidence-manifest.json`, 압축 파일 자체의 SHA-256은 `evidence.zip.sha256`에 있다. 압축 묶음은 기존 연구에 덮어쓰지 말고 새 디렉터리에서 검사한다. 외부 입력은 위 다운로드·복원 절차가 여전히 필요하다.

모든 사후 계산이 끝난 원기록에서 새 묶음을 만들려면 다음과 같이 별도 출력을 지정한다. 생성기는 압축 후 모든 항목을 다시 읽어 바이트 수와 SHA-256을 확인하며 기존 묶음을 덮어쓰지 않는다.

```powershell
$anytimeBundleOut = 'D:/DivertorOptimization/exports/anytime-evidence-repeat-01'
& $anytimePython results/nf-anytime-fullspace-20261106/build_evidence_archive.py --study $anytimeStudy --out $anytimeBundleOut
if ($LASTEXITCODE -ne 0) { throw 'Evidence bundle verification failed' }
```

새 검증 모듈의 단위 테스트는 물리 캠페인을 재실행하지 않는다. 실제 장비의 시간 비교가 진행 중일 때는 테스트를 동시에 실행하지 않는다.

```powershell
& $anytimePython -m pytest -q verification/test_nf_anytime_search.py verification/test_nf_anytime_campaign.py verification/test_nf_anytime_reference.py verification/test_nf_anytime_uncertainty.py verification/test_nf_anytime_report.py verification/test_nf_anytime_final_audit.py verification/test_nf_anytime_rejections.py verification/test_nf_anytime_bundle.py verification/test_publication_verifier.py
if ($LASTEXITCODE -ne 0) { throw 'Anytime validation tests failed' }
& $anytimePython tools/verify_publication.py
if ($LASTEXITCODE -ne 0) { throw 'Preserved-source verification failed' }
```

마지막 명령은 과거 보존 코드·연구의 지문과 저장된 판정을 확인한다. 현재 AGENTS 지침과 NPZ 출력 메타데이터는 과거 배포 뒤 의도적으로 갱신됐으므로, 그 두 파일의 원본은 `docs/publication/legacy-20261029/`에서 원래 해시로 대조한다. 알려진 현재 바이트도 별도로 제한하며 임의의 추가 변경을 허용하는 우회가 아니다. 원래 보존 manifest는 바꾸지 않는다. 이번 캠페인의 실제 실행 코드는 해당 `frozen.json`과 독립 감사로 따로 검증한다. 보존 검사는 이번 형상의 새 수치 수렴 계산을 대신하지 않는다.
