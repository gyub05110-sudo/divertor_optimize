# 3단계 자료 재현

[최초 비교계획](protocol.ko.md)은 커밋 `b75c2a9`에 있다. [입력 등록](inputs/registration.json)은 커밋 `5530789`에 있다. [실행 환경](execution-environment.json)과 [개발 예 시운전](commissioning.json)은 새 평형 열계산 전에 작성했다. 연구 방향을 다시 선택할 때 사용할 이전 체크포인트는 `codex/selective-candidate-validation`의 `8e63e34`, 두 개발 평형 비교는 `codex/d3d-peak-input-and-shape-study`의 `2c5e0d2`다.

## 원자료와 측정

공개 원파일의 고정 revision, 개별 다운로드 URL, 크기와 SHA-256은 [선택 목록](inputs/selection.json)에 있다. 같은 바이트를 되받으면 원래 익명 파일명으로 사용할 수 있다. 이전 입력 24개 원파일에서 얻은 프레임 지문과 새 입력 12개의 1400 ms 프레임을 비교한 절차는 [등록 코드](../../verification/stage3_inputs.py)에 있다. 등록 재연에는 이 코드가 참조하는 저장소의 이전 연구 증거도 필요하다. [공개 데이터셋의 고정 revision](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge/tree/1e280905b85f2a6fdde7e06fca8cf3a1edf447cb)

원실행은 Windows에서 아래 절차로 수행한다. `$study`는 저장소 밖의 증거 디렉터리다. 일단 봉인된 `expected.json`은 디렉터리, 인터프리터, 과학 코드 및 원본 입력 경로에 결합된다. 동일 호스트에서 중단된 실행은 마지막 명령으로 이어 간다. 다른 환경의 속도 재측정은 새 연구 디렉터리에서 따로 기록한다. 600초는 벽시계 예산이므로 같은 시드에서도 제한시간까지 완료한 시도가 달라질 수 있다.

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
