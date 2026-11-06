# 재현과 보존

실행 자료의 위치는 D:/DivertorOptimization/studies/unused-equilibrium-20261106/이다. 이 보고서 폴더에는 작은 원기록과 그림을 복사했다. 원본 parquet 12개, 이전 원자료와의 대조용 다운로드, 사례별 입력 root, 코드 스냅샷과 전체 로그는 위 D 경로에 있다. 공개 자료는 [Sophelio 데이터셋](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge)의 CC-BY-4.0 배포본이며, 파일별 고정 revision URL과 SHA-256은 downloads.json에 있다.

저장소 루트에서 현재 실행 기록을 다시 감사할 수 있다.

```powershell
$study = 'D:/DivertorOptimization/studies/unused-equilibrium-20261106'
$env:PYTHONPATH = "$PWD/.local_dependencies/d3d-fixed-equilibrium"
$env:PYTHONDONTWRITEBYTECODE = '1'
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
$env:TEMP = 'D:/DivertorOptimization/tmp'
$env:TMP = $env:TEMP
$env:MPLCONFIGDIR = 'D:/DivertorOptimization/cache/matplotlib'
$env:CUPY_CACHE_DIR = 'D:/DivertorOptimization/cache/cupy'
$env:CUDA_CACHE_PATH = 'D:/DivertorOptimization/cache/cuda'
python -m verification.nf_unused_live --study $study --audit-only
python -m verification.nf_unused_numerics_compare --study $study
python -m pytest -q verification/test_nf_unused_inputs.py verification/test_nf_unused_live.py verification/test_nf_base_live.py verification/test_nf_base_budget.py verification/test_nf_base_budget_audit.py
```

이번 실행에서 이 대상 테스트 30개가 통과했다. 원래 nf_base_live와 nf_base_budget의 해시를 확인하고 새 실행기는 별도 모듈로 추가했다. 초기 검토에서 지적된 사전 고정 해시 연결과 중간 종료 JSON 처리는 실시간 실행 전에 보완했다. 각 자식 프로세스는 D의 전용 임시 폴더와 빈 GPU 캐시를 사용한다.

동일 입력으로 새 실행을 만들 때에는 기존 원자료를 읽기 전용으로 참조하고 별도 출력 폴더를 만든다. 아래는 이 컴퓨터에 위 study 원자료가 보존돼 있다는 조건이다. 새 폴더가 이미 있으면 중단한다. 각 원자료를 다른 위치에 복원할 때에는 고정 URL에서 받은 파일의 바이트 수와 SHA-256을 downloads.json과 대조한 뒤 그 파일의 local 경로를 새 목록에 반영해야 한다. 열계산에는 원자료의 전력·밀도 측정값을 대입하지 않으며 고정 ReferencePoint를 사용한다.

```powershell
$trial = 'D:/DivertorOptimization/studies/unused-equilibrium-repeat'
if (Test-Path -LiteralPath $trial) { throw 'Choose a fresh output directory' }
$null = New-Item -ItemType Directory -Path $trial
foreach ($name in 'protocol.ko.md','selection.json','downloads.json','numerics-protocol.ko.md') {
    Copy-Item -LiteralPath "$study/$name" -Destination "$trial/$name"
}
python -m verification.nf_unused_inputs --study $trial
if ($LASTEXITCODE -ne 0) { throw 'Input preparation failed' }
python -m verification.nf_unused_live --study $trial
if ($LASTEXITCODE -ne 0) { throw 'Inspect preserved trial records' }
python -m verification.nf_unused_summary --study $trial
python -m verification.nf_unused_numerics --study $trial
python -m verification.nf_unused_numerics_compare --study $trial
python -m verification.nf_unused_figures --study $trial
```

실시간 시험은 머신의 실제 경과 시간을 사용하므로 재실행 때 같은 후보까지 정밀 평가를 마치는 횟수와 최종 선택은 달라질 수 있다. 당시 판단의 재현은 보존한 이벤트를 감사하는 첫 번째 명령으로 수행한다. --case와 --seed로 사전에 등록한 한 쌍만 실행할 수 있고, --skip-completed를 명시하면 감사 가능한 완료 쌍은 보존한다. 실행 중 결과나 코드를 바꿔 다른 실행과 혼합하는 경우에는 해시 검사에서 중단한다.

