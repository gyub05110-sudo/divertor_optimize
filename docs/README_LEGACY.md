# 사진 상태의 D3D 직접 자기장 매핑 복원

사진의 d3d-divertor-floor-shelf / 1400 ms / case 6758 / 29단계 이력을 복원했습니다.
열모델은 외측 중간면의 지수 열원 → 하부 최초 벽 충돌이며 횡방향 수송은 없습니다.
archived math는 archive_math에 원본 그대로, 실행용 수학은 d3d_optimization에 import 이름만 바꾸어 보존했습니다.
최근 열전도·확산 코드와 계산 결과는 복사하지 않았습니다. 원래 freegs_divertor는 수정하지 않았습니다.
복원 범위는 해당 D3D 실행/뷰어/재현에 필요한 코드와 입력입니다. 과거 전체 CPD/FreeGS 폴더의 시점별 백업은 아닙니다.
원본 manifest와 탐색 결과는 수정하지 않았으며 local_manifest.json만 이 폴더의 parquet 경로를 가리킵니다.
viewer UI 코드는 현재 보존본이며, 사진의 과거 열분포/형상/검증 파일을 로드하도록 등록 목록을 제한했습니다.

## 기존 결과 보기
PowerShell에서 input/D3D/optimized/d3d-divertor-floor-shelf/viewer3d/start_viewer.ps1을 실행합니다.
로컬 주소는 http://127.0.0.1:8776/ 입니다.

## 새로운 계산
python -X utf8 optimize_d3d.py --mode preflight --output input/D3D/optimized/NEW_PREFLIGHT
위 명령은 보존한 원래 floor+shelf 실행입니다.

## IMP·OMP 독립 최적화

새 실행은 실제 IMP와 OMP에 각각 단위 열원을 설정하며 50:50 분배를 적용하지 않습니다. 횡방향 열수송은 없습니다.
기준 형상에서 선택한 6개 구간(floor, shelf, IMP 69/70 일부, OMP 58/60 일부)의 24개 B-spline 계수를 사용합니다.
IMP Z > -1.15 m와 OMP 53/54/55는 고정하며 전체 벽의 열유속을 목적함수에 포함합니다.
λq = 1, 2, 3, 5, 10 mm의 최대 q/P를 각각의 열원 기준 형상과 비교합니다. 두 결과를 합산하지 않습니다.

```powershell
python -X utf8 optimize_independent_sources.py --mode preflight --source both --output input/D3D/optimized/my-independent-preflight
python -X utf8 verification/run_continuous_search.py --mode optimize --source both --boundary-mode continuous --output input/D3D/optimized/my-independent-search
python -X utf8 optimize_independent_sources.py --mode verify --source both --output input/D3D/optimized/my-independent-search
```

기본 탐색은 기존과 같은 Sobol 1024/512/512개(20/50/100 mm), DE seed 42/93, 최대 40세대, 초기 최대 64개 개체, polish 없음입니다.
유효 형상 1000개를 확보하는 추가 Sobol 단계는 짧은 구간의 길이를 고려해 ±min(20 mm, 최소 구간 길이×0.15)로 샘플링합니다.
최대 허용 계수 범위와 형상/각도/에너지 제약은 유지합니다. 기본 DE 초기 개체 수는 확보한 feasible 후보 수에 따라 32–64개입니다.
추가 샘플링 100000개 한도에 도달하면 실패를 보고하며 로그를 보존합니다. 후보 중복은 12자리 반올림 계수로 제거합니다.
현재 기본값은 continuous 경계입니다. IMP 69번의 시작점을 ψ=ψX의 연속 strike로 맞추고 69/70 연결점의 샘플 틈을 제거합니다.
형상 끝점 변위·접선 고정은 유지하며, SOL측 separatrix 열유속 극한도 전체 벽 최대값에 포함합니다. 이 극한점에 유한한 열원 전력을 추가하지 않습니다.
기본 imp-buffer-mm=0으로 PFR 쪽 설계 영역을 연장하지 않습니다. --imp-buffer-mm는 별도 설계 선택이며 물리 검증을 대신하지 않습니다.

현재 검증은 cubic 320/640/1280/2560 및 quintic1280 경로에서 선택 후보를 다시 계산하고 전체 벽 first-hit를 독립 대조합니다.
각 폭의 1280→2560 비율 변화와 cubic→quintic 변화가 각각 0.002 이하여야 해상도 검사를 통과합니다.
개선 확인에는 감소량이 관측 변화의 두 배 + 0.0001보다 커야 하며 모든 검사에서 기준 형상보다 악화되지 않아야 합니다.
이는 경험적 수치 선별 기준이며 엄밀한 오차 상한·전역 최적해 증명이 아닙니다. 형상 제약 통과와 개선 확인은 별도 필드입니다.
위 실행 래퍼는 전체 계산 모듈 및 실제 기준벽 입력의 해시를 탐색 전에 저장합니다. 검증은 누락·불일치를 거부합니다.
보존된 legacy 실행은 --boundary-mode legacy로 재현할 수 있습니다. legacy 검증만 cubic320/640, quintic640 및 --verification-attempts 최대64 후보 방식을 사용합니다.

연결 시험은 --smoke 옵션으로 실행합니다. 48경로·작은 Sobol 표본·DE 1세대이므로 최종 최적화 성능으로 해석하지 마세요.
구현 확인 결과는 input/D3D/optimized/independent-smoke-v2-20261006 및 independent-preflight-20261006에 있습니다.
각 실행의 IMP/OMP 하위 폴더에 manifest.json, evaluations.jsonl, case_arrays.npz, viewer_result.json이 저장됩니다.
http://127.0.0.1:8777/baseline_review/의 ‘독립 최적화 결과’에서 viewer_result.json을 선택하면 후보 형상과 q/P를 표시합니다.
기준 자기장에 따른 두 leg는 후보 형상과 별개로 기존 기준 형상의 strike 좌표를 유지합니다.
연속 경계 결과의 열유속 점에는 후보의 실제 SOL 극한 충돌점이 추가됩니다. 색상 설계영역은 해당 결과에 저장된 경계를 사용합니다.

2026-10-24 문헌·계획·적대적 검토·실행 기록은 docs/superpowers/plans/2026-10-24-*.md에 있습니다.
본 탐색은 input/D3D/optimized/continuous-boundary-20261024, 기존 OMP 후보의 수정 모델 재검증은 continuous-boundary-incumbent-20261024에 보존합니다.
이 두 실행의 전체 해시는 탐색 진행 중 보강됐고 sidecar에 해당 시점이 명시됩니다. 이후 실행부터 위 래퍼를 사용합니다.

과거 시간 기록: 사진의 floor+shelf 실행 7233 unique 후보 / 600.565초(10.0분), expanded 실행 14702 unique 후보 / 734.699초(12.2분).
새 실행은 24변수이며 두 열원을 별도로 계산하므로 과거 실행 시간과 동일하지 않습니다.
