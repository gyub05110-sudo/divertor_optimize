# 기준 형상 검토 — 2026-10-06

사진의 실행은 `d3d-divertor-floor-shelf`, `d3d_shot_1a3155370c.parquet`, EFIT 1400 ms, retained case 6758이다. 새 요청의 그림은 **기준 형상**을 사용한다. 과거 최적화 형상을 새 기준으로 바꾸지 않았다.

## 복원 범위

`freegs_divertor_no-transport`에는 해당 D3D 직접 자기장 매핑을 실행하고 과거 결과를 보는 데 필요한 코드·입력·원본 결과·뷰어가 있다. 당시 전체 FreeGS/CPD 폴더의 완전한 시점 백업을 확보한 것은 아니다.

수학 커널은 `migration-20261003-152601`의 보존본이다. 실행 패키지는 import 이름만 변경했다. search/validate 래퍼는 unified 모델 변경 이전 백업에서 복원했다. 뷰어의 현재 보존 UI는 유지하되 등록 목록을 사진 실행 하나로 제한했다. 원본 manifest, 탐색 로그, 형상 및 품질 검증 결과는 수정하지 않았다. 원래 `freegs_divertor`에는 쓰지 않았다.

## 요청사항을 이해한 방식

1. 하부 X-point에서 실제 ψX 등고선을 따라 inner/outer separatrix leg를 각각 벽의 최초 충돌까지 표시한다. LCFS 배열의 가장 낮은 점에서 임의의 직선을 그리지 않는다.
2. 설계영역은 원래 floor(segment 68)와 shelf(segment 59) 전체에 사용자가 선택한 기준 SOL 최초 충돌 벽 부분을 더한다. 2026-10-06 추가 지시에 따라 IMP 쪽은 **Z<=-1.15 m**만 포함한다. 첨부한 OMP 바깥 경사 벽 **segment 53/54/55**는 최적화에서 제외하고 고정한다. shelf 인근 58/60의 선택 충돌 부분은 유지한다. 원래 segment 전체가 아니라 선택한 충돌 부분만 새 설계영역으로 본다. 후보에서는 최초 충돌·차폐·전체 벽 열부하를 다시 계산하고, 제외한 고정 벽의 열부하도 계속 포함한다.
3. 사용자가 실제 IMP·OMP 열원을 요청했다고 확인했다. `IMP: R_IMP−x`, `OMP: R_OMP+x`의 별도 상류 프로파일을 만들고 각각 독립적으로 정규화/최적화한다. 두 입력을 하나의 열원에서 방향별 50:50로 분리하지 않는다. 미래 IMP/OMP 해의 총 열부하를 합치려면 실제 입력 전력이 따로 필요하다. 사진의 과거 결과도 이미 outer-branch power=1이며 50:50 결과가 아니다.
4. 두 leg 사이 공간은 PFR(private flux region), 주변 벽 부품은 PFC(plasma-facing component)다. 직접 SOL 매핑은 이 공간의 열분포를 결정하지 못한다. 그림에서 분홍색은 '열분포 미정'이며 물리적인 q=0을 뜻하지 않는다.

## 계산으로 확인한 위치

| 항목 | R (m) | Z (m) | 원본 벽 segment |
|---|---:|---:|---:|
| X-point | 1.17659194 | -1.25239897 | — |
| inner leg 최초 충돌 | 1.08825770 | -1.28532210 | 69 |
| outer leg 최초 충돌 | 1.22136842 | -1.35869801 | 68 |
| IMP separatrix / 중간면 | 1.10136428 | -0.06456412 | — |
| OMP separatrix / 중간면 | 2.25131111 | -0.06456412 | — |

IMP SOL 최초 충돌 표본은 경사면 69의 일부와 수직 벽 70의 일부에 있다. 원래 수직 벽의 표본 포락 범위는 대략 Z=-1.2163~-0.1669 m지만, 사용자가 선택한 설계영역은 **Z=-1.15 m부터 아래쪽**으로 제한한다. 이것은 열원이나 최초 충돌의 물리적 범위를 변경한 것이 아니다.

OMP SOL은 floor 68, shelf 59와 60/58/55/54/53의 일부에 먼저 닿는다. **53/54/55의 최초 충돌과 열부하 계산은 유지하되 형상 최적화에서 제외한다.** 61~67 등의 차폐된 후속 교차점은 최초 충돌 집합에 포함하지 않는다. source 최대 범위는 IMP 약 84.054 mm, OMP 약 98.776 mm다. 전체 유효 source 범위를 사용하므로 좁은 λq에서는 매우 약한 tail의 충돌 벽도 보인다.

## 그림과 수치의 적용 범위

각 중간면에서 800개 SOL 경로를 과거 RK4 커널로 새로 계산했다. trace step은 1 mm다. 두 leg는 saddle Hessian의 국소 방향에서 ψX 등고선으로 투영한 뒤 DOP853으로 추적했고, 출발 거리/step을 절반으로 바꿔 타점 변화를 확인했다. 타점 변화는 inner 1.82e-8 m, outer 2.46e-8 m다. 이것은 주어진 cubic EFIT 보간 내의 수치 재현성이며 실험 위치 정확도 인증이 아니다. 원본 EFIT 격자와 형상 한계는 유지한다.

SOL 표본의 미해결 최초 충돌은 두 입력에서 모두 0이지만, 연속 전체 footprint 경계의 수렴을 확인했다는 의미는 아니다. 새 설계영역을 실제 코드에 연결하기 전에는 전이점과 매우 얕은 입사 구간을 추가 세분화해야 한다. 그림의 선 색은 열유속이 아니라 경로/영역 분류다. 새 IMP/OMP 열유속과 최적화는 아직 계산하지 않았다.

## 실행 및 확인

- 검토 페이지: `http://127.0.0.1:8777/baseline_review/`
- 재실행: 루트 `start_review.ps1`
- 좌표 재계산: `python -X utf8 review_baseline.py`
- 그림 재생성: `python -X utf8 plot_baseline_review.py`
- 저장된 최초 충돌에 선택 영역 적용: `python -X utf8 design_region_policy.py` (`optimization_region_policy.json` 사용). 새 trace를 재계산할 때도 같은 선택 정책을 적용한다.
- 복원 확인: `python -X utf8 verify_restoration.py`
- 원래 수학/형상 테스트: `python -X utf8 -m pytest verification/test_d3d_frozen_heat.py verification/test_d3d_multisurface.py -q` → 33 passed
- CLI 사전 점검: `restored-preflight-20261006-v2` → exit 0, 원래 baseline peak와 동일. v1은 후대 래퍼 인자 불일치로 실패한 기록이며 성공 결과로 취급하지 않는다.

검토용 3D 페이지는 별도 페이지다. 과거 열분포/29단계 뷰어는 보존했고, leg와 새 설계영역의 생산 경로 통합은 이 기준 형상 표시를 확인한 후 진행한다.
