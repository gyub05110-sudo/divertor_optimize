# 공통 격자 재현

보존 위치는 `D:/DivertorOptimization/studies/firefly-baseline-20261105/common_rank`다. 형상 최적화는 `D:/DivertorOptimization/worktrees/d3d-report-correction-20261104`의 커밋 `dc8ab43a9a930ca155813cb53a448bba39f6291f`에서 읽었다. 원래 모델과 이 연구 폴더의 연결은 `metrics.original` 및 `original_imp.py`에서 찾을 수 있다. 부모의 `REPRODUCE.md`, `firefly-midplane.patch`, `checksums.json`은 FIREFLY·FLARE 소스/빌드 및 앞선 기준 계산을 기록한다.

## 계산

기존 D 드라이브 환경을 사용한다. PowerShell에서 이 폴더로 이동한 뒤 실행한다. 자기 격자 두 개는 이미 보관했다. `start.sh`는 입력과 경계가 준비된 후 격자를 생성하는 래퍼다. 열수송은 CPU MPI 4랭크/계산이다. 서로 다른 형상을 병렬로 진행할 수 있으며 모든 계산에 4랭크를 유지한다. 기존 저장 파일을 건너뛰는 스크립트와 덮어쓰는 대조군 스크립트가 있으므로, 원본을 보존하려면 별도 복사본에서 경로를 조정한다.

```powershell
$enter = '/mnt/d/DivertorOptimization/studies/firefly-baseline-20261105/enter.sh'
$env:OPENBLAS_NUM_THREADS='1'
python prepare.py
wsl -d Ubuntu -u root -- bash $enter common_rank/start.sh
wsl -d Ubuntu -u root -- bash $enter common_rank/dispatch.sh
wsl -d Ubuntu -u root -- bash $enter common_rank/stress.sh
wsl -d Ubuntu -u root -- bash $enter common_rank/cold_followup.sh
python split_wall.py
wsl -d Ubuntu -u root -- bash $enter common_rank/split_control.sh
wsl -d Ubuntu -u root -- bash $enter common_rank/replay.sh
python original_imp.py
python thomson_registration.py
python flux_label_comparison.py
```

`design.md`에는 주 분석 계획이 있다. `stress.sh`, `cold_followup.sh`, `split_control.sh`의 저온·고리 해상도·동일 벽 재분할 시험은 사후 탐색이다. 첫 3개 저온 시드의 피크 부호 역전을 보고 나머지 3개 및 1 m 스텝 대조군을 추가했다. 6시드 저온 결과에 대한 구간의 해석에는 이 순차 선택을 고려한다.

## 통계 및 검산

```powershell
$env:OPENBLAS_NUM_THREADS='1'
python analyze_rank.py
python rank_figure.py
python tabulate_report.py
python assemble_report.py
python verify.py
```

CUDA의 FP64로 시드별 프로파일을 가중합했다. 모든 조합의 두 지표를 CPU에서 독립적으로 중복 계산해 수치 일치를 요구한다. `verify.py`는 이와 별도로 각 native 고리를 교차하는 출력 고리로 하나씩 나누어 부분 원뿔대 면적에 비례해 전력을 분배하고, 결과를 내보낸 공간 프로파일과 대조한다. 저장된 모든 입력·결과의 해시 목록은 `checksums.json`, 검사 요약은 `verified.json`이다. 부모 연구의 보호 대상 161개도 검사한다.

`profile_열원_형상.csv`는 고리 중심 거리와 R, Z, 면적, 원래 모델 열유속, FIREFLY 통합 열유속과 개별 시드 열유속을 담고 있다. `rank.json`에는 공간 해상도·시드별·수송 입력별 결과와 수지가 들어 있다. `results/*.nc`는 보정 전 집계값, `results/*.dat`는 자기선 결과다.

## 공개 자료 대조

`github_audit.py`, `public_id.py`, `compare_public_shots.py`, `public_benchmarks.py`, `inspect_alternate.py`가 각 JSON 기록을 작성한다. GitHub 예제 파일은 Git LFS SHA-256을 먼저 확인하고, 방전별 자기장·Thomson 시각 및 온도 등 5개 배열의 지문을 대조했다. `compare_public_shots.py`가 읽는 `public-starter`는 `https://github.com/Sophelio/fusion-equilibrium-challenge-starter.git`의 이력이 있는 클론이다.

