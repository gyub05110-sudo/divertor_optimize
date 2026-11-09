# 4단계 기록 검산과 복원

분석 단위는 등록한 24개 원파일이다. 계획은 파일당 10개 슬롯과 첫 시드의 고정 형상 재계산 조건을 정한다. 원본 parquet의 고정 URL·바이트 수·SHA-256, 1400 ms 등록 상태 및 실험 환경을 각각 확인할 수 있다. [전향 계획](protocol.ko.md) · [선택 명단](inputs/selection.json) · [등록 명단](inputs/registration.json) · [실행 환경](execution-environment.json)

## 원래 실행 위치에서 재검산하기

원래 연구 폴더는 `D:/DivertorOptimization/studies/nf-stage4-20261109`이다. `expected.json`과 `expected.sha256`은 실행 순서, 코드·입력·설정, 인터프리터와 원래 경로를 봉인한다. 처음 두 명령은 같은 환경에서 완료된 슬롯과 제품 비교·예측 표를 검산한다. 마지막 명령은 시간 제한 탐색과 분리된 고정 형상 재계산을 검산하고 누락된 측정을 수행한다. [실행 스크립트](run-study.ps1)

```powershell
$py = 'C:/Users/FPL-user/miniforge3/python.exe'
$study = 'D:/DivertorOptimization/studies/nf-stage4-20261109'
$paper = 'results/nf-stage4-20261109'
$env:PYTHONPATH = (Resolve-Path '.local_dependencies/d3d-fixed-equilibrium').Path
$env:DIVERTOR_GPU_RUNTIME = 'D:/DivertorOptimization/runtime/cupy'
& $py -B -m verification.stage4_campaign check --study $study
& $py -B -m verification.stage4_report --study $study
& $py -B -m verification.stage4_replay --study $study
```

`registered-runs.csv`는 각 슬롯의 상태와 시각, `product-pairs.csv`는 같은 파일·시드·시간에서 자동과 전체의 차이, `policy-predictions.csv`는 첫 시드의 네 강제 탐색과 사전 예측을 연결한다. `stage4_report.json`은 기준벽·체크포인트·원기록 해시를, `stage4_replay.json`은 사후 재계산 조건의 발동 여부와 해당 비교를 담는다. `runs/<파일>/<시드>/<가족>/events.jsonl`에는 후보 생성 전에 기록한 정책과 후보별 시각이 있다. `source-code/`는 실행계획에 결합된 코드 스냅샷이다.

## 보고서 그림과 면별 진단 재생성

최종 보고서와 원래 표면 열유속 배열에서 다음 명령으로 두 PDF와 PNG 및 첫 시드의 면별 진단을 다시 만들 수 있다. `$py`, `$study`, `$paper`는 위와 같다. 진단 명령은 JSON을 표준출력으로 내보낸다.

```powershell
& $py -B "$paper/plot_policy_validation.py" --study $study --output $paper
& $py -B "$paper/plot_hotspot_profiles.py" --study $study --output $paper
& $py -B "$paper/region_mechanism.py" --study $study
```

## 완료된 기록을 보존하기

Python 3.11 이상 표준 라이브러리로 실행한다. 저장 위치는 매번 새 디렉터리를 지정한다. 원래 연구 폴더에서 모든 과학·후처리 작업이 끝나면 아래 두 명령을 실행한다. 두 보존본의 `stage4-manifest.json`에서 `inventory_sha256`가 같으면 전체 상대 경로·크기·파일 해시 목록이 일치한다. 각 명령이 출력하는 **별도의 manifest SHA-256**을 배포 기록에 보관한다. [보존 도구](package_evidence.py)

```powershell
$compact = "$paper/evidence-compact"
$full = 'D:/DivertorOptimization/exports/nf-stage4-20261109-full'
& $py -B "$paper/package_evidence.py" build --mode compact --study $study --bundle $compact
& $py -B "$paper/package_evidence.py" build --mode full --study $study --bundle $full
```

Git에서 압축을 풀지 않고 결과 표를 열 수 있게 하려면 완료된 연구 폴더의 `stage4_report.json`, `stage4_replay.json`, `registered-runs.csv`, `product-pairs.csv`, `policy-predictions.csv`를 이 안내와 같은 디렉터리에 복사한다. 보존본의 파일 목록에서 해당 상대 경로의 SHA-256을 대조한다.

| 모드 | 포함하는 바이트 |
|---|---|
| `compact` | 원파일 등록·요청·해시와 프레임 지문, 봉인 코드와 계획, 실행 영수증·표준출력·후보별 이벤트·결과 JSON, 고정 형상 재계산 기록과 비교 표. 생략된 원파일·수치 배열도 파일 목록에 크기와 SHA-256을 남긴다. |
| `full` | 같은 기록과 원본 parquet, `geometry.csv`, 압축된 벽 프로파일·도달 경로, `boundary.npz` 등 원래 연구 폴더의 전체 증거. |

두 모드 모두 `runtime/` 아래 머신별 캐시는 생략한다. 연구 폴더의 다른 파일은 전부 목록에 넣는다. 같은 내용은 분할 ZIP 안에서 한 번만 저장한다. 완결된 240개 슬롯의 결과와 24개 고정 형상 재계산 결정이 필요하며, 선언된 원파일 또는 내보낸 측정 파일이 사라지거나 바뀌면 작성이 실패한다.

## 다른 컴퓨터에서 보존본 확인하기

보존본 루트에 복사한 `package_evidence.py`를 사용한다. 이 명령은 GPU·원래 절대 경로 없이 저장된 바이트, 압축 조각, 입력·코드·슬롯·재계산 결과 사이의 해시 연결을 검사한다. `--manifest-sha256`에는 배포 기록에 보관한 해당 보존본의 값을 넣는다.

```powershell
$bundle = 'D:/received/nf-stage4-evidence'
$sha = '<해당 보존본의 발표된 manifest SHA-256>'
python -B "$bundle/package_evidence.py" verify --bundle $bundle --manifest-sha256 $sha
python -B "$bundle/package_evidence.py" restore --bundle $bundle --manifest-sha256 $sha --out 'D:/reproductions/stage4-records'
```

`full`로 복원하면 원본 입력과 모든 수치 배열을 함께 얻는다. `compact`로 복원하면 수록한 기록만 생성한다. 원래 연구 폴더나 그 전체 복사본에도 접근할 수 있으면 `restore` 명령에 `--study D:/DivertorOptimization/studies/nf-stage4-20261109`를 추가해 누락된 대용량 파일을 해시 검증 후 채울 수 있다. 보존된 측정 시각은 원실행의 값이다. 새 장비에서 동일한 벽시계 예산으로 최적화를 다시 수행하면 완료 후보 수와 선택 형상이 달라질 수 있다. [다시 최적화하는 공개 CLI](../../docs/user-optimizer.ko.md)
