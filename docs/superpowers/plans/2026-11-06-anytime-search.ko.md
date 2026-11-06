# 여러 시간 예산의 탐색 및 최종 판정 구현 계획

> **For agentic workers:** Use superpowers:subagent-driven-development to implement and review independent tasks. The user has approved execution of the preceding design; no repeated approval is required.

**Goal:** DIII-D 기준벽을 유지한 전체 탐색에서 시간에 따른 정밀 확인 결과를 비교하고, 최종 선택의 수치 설정 민감도를 판정한다.

**Architecture:** 과거 실행기와 물리식을 변경하지 않는 별도 검증 모듈을 추가한다. 기존 SearchCase의 적응 탐색을 저비용/정밀 설정으로 각각 실행하고, 기존 투영 무작위 생성기의 직접 정밀 탐색을 대조한다. 정밀 확인 완료 시각으로 중간 결과를 재구성한다. 결과는 모두 D에 보존한다.

**Tech Stack:** Python, NumPy/SciPy, CuPy RTX 3090, pytest, JSON/JSONL, Matplotlib.

## 승인된 범위와 이번 실행

- 주장은 지정한 가정에서 좋은 형상을 빠르고 정확히 찾는 방법의 성능이다. FIREFLY 결과·실측은 탐색 점수와 계수 조정에 사용하지 않는다.
- 첫 단계는 세 기존 외부 방전 평형에서 세 고정 시드(20611,20627,20639), 세 방법의 600초 실행이다. 200/400/600초 지점을 계산 전에 고정한다. 27개 직렬 실행으로 명목 4.5시간이다. 이전 200초 시험은 별도 보존한다. 1800초 확장과 추가 방전은 이 단계의 완료 조건에 포함하지 않으며, 미실행으로 명시한다.
- 여기 사용한 평형의 과거 결과는 이미 공개되었다. 새 후보를 제안하지만 완전히 보지 않은 평형의 전향 시험으로 부르지 않는다. OMP 1 MW, 동일 ReferencePoint, Rmin 20 mm, 기존 원벽·제약을 사용한다.
- 방법: coarse_adaptive(기존 SearchCase, cubic320, 새 적격 저비용 후보 네 개마다 그 묶음 최저 후보를 cubic1280으로 확인), fine_adaptive(동일 SearchCase와 초기화 정책, cubic1280), random_fine(기존 투영 무작위 생성기, cubic1280). 모델과 탐색기 튜닝은 하지 않는다. 각 방법에서 기준벽 유지가 가능하다.
- 적응 탐색은 동일한 형상 초기화 설정 population=12, projected_seeds=12, 과거 최적 형상/파일을 읽지 않는 staging root를 사용한다. 최적화 피드백에 따라 생성 경로는 달라질 수 있어 이후 후보 접두사 동일성을 요구하지 않는다. 무작위 방법은 같은 RNG 시드·기하 제약을 유지한다. 후보 12개 상한은 제거한다.
- 프로세스 시작부터 임포트·자기선 준비·GPU 초기화·후보 생성/탈락·검증·실행 중 로그 비용을 포함한다. 장비는 공유 실행하지 않는다. 전용 빈 D 캐시를 사용하고 방법 순서를 세 방법 Latin rotation으로 교차한다. 생성·평가 모두 최대 시각 뒤 새 작업을 시작하지 않는다.
- 초과된 평가도 값·시각·기하 적격성을 기록하지만 시한 내 선택·최초 개선 시각에 넣지 않는다. 중간 시점은 결과 확인만 하며 선택 규칙을 변경하지 않는다. 실패한 방법도 남기고 평형 교체를 금지한다.
- 별도 사후 계산에서 각 시점 최종 후보와 마지막 초과 후보를 고정하고 cubic640/quintic640/quintic1280 및 cubic1280 재계산을 대조한다. 같은 평형에서 시간 내 정밀 평가했던 별개 후보 중 native 피크가 가장 가까운 후보도 미리 고정한다. 관측된 집합 안에서 기준 대비 개선과 후보 간 순위 분리를 판정하며, 형상 공간의 모든 인접 후보나 연속체 오차 보증으로 확대하지 않는다.
- 추가 참조는 각 case/seed의 고정 투영 후보 24개를 coarse/fine로 전수 계산하는 패널이다. 이 패널에서만 전수 최저 대비 선별 손실을 주장한다. 적응 탐색 전체에서의 최저는 관측된 후보 집합의 best-known이며 전역 최적이 아니다.

## Task 1: 기록과 실제 탐색 실행기

- [x] `verification/nf_anytime_search.py`, `verification/test_nf_anytime_search.py`: 먼저 모의 시계로 200초 경계, 초과 값 보존, partial batch, 기준벽 fallback, 캐시 중복, 실패 처리, 시드 반복을 검증하는 실패 테스트를 작성한다.
- [x] 기존 SearchCase를 조합/최소 하위 클래스로 감싸며 기존 생산 파일은 수정하지 않는다. `run_arm(policy, input_root, seed, out, started_ns, budget_s=600, clock_ns=perf_counter_ns)`를 구현한다. 반환 JSON에 events, candidates, baseline, checkpoints를 포함한다.
- [x] `snapshot(record,t)`는 `end_s <= t`인 적격 정밀 결과와 시간 내 준비된 기준벽만 본다. late 값은 기록하되 선택에서 제외한다.
- [x] `python -m pytest -q verification/test_nf_anytime_search.py`의 실패/성공을 확인하고 독립 사양·품질 검토를 받는다. 별도 원래 평형의 240초 실제 연결 시험에서 GPU/CPU parity와 원 event 감사를 통과했다.

## Task 2: 실행 고정과 감독

- [ ] `verification/nf_anytime_campaign.py`: 원자료와 원벽 해시를 기존 `nf_unused_live.load_cases`로 검증한다. 코드/입력/프로토콜 해시를 동결하고 재개 시 바이트 변경을 거부한다. 새 디렉터리의 raw event와 프로세스 로그를 보존한다.
- [ ] 27개 serial fresh-process 실행을 감독한다. 자식 오류·절단 JSON은 실패로 남기며 다음 방법은 계속 실행한다. 요약은 원 event에서 독립 재집계한다.
- [ ] D 임시/캐시만 사용, 원래 프로덕션 코드 해시 불변과 cold-start 비용 경계를 확인한다.

## Task 3: 추가 정확도 및 수치 판정

- [ ] `verification/nf_anytime_reference.py`: 고정 24개 패널의 열점수 보기 전 후보/해시를 저장하고 전수 coarse/fine를 계산한다. 같은 정밀 예산의 무작위 선택 기대값을 정확히 계산한다.
- [ ] `verification/nf_anytime_uncertainty.py`: 최종 채택/초과 후보 목록을 고정하고 각 평형 기준과 함께 추가 수치 설정으로 계산한다. 비교된 설정마다 개선 부호를 표시하고 범위가 겹치면 미세 우열 미확정으로 기록한다.
- [ ] 후보 구별과 기준 대비 큰 개선을 분리하고 A2 감소/부적격/시간초과를 숨기지 않는다. 모델 가정 차이와 수치 민감도를 혼동하지 않는다.

## Task 4: 검산과 보고

- [ ] 200/400/600초 쌍별 결과, 전체 anytime 곡선, 최초 5% 개선 시각, 실제 비용, 평가 수를 보고한다. 여러 시드를 독립 방전으로 세지 않는다.
- [ ] `results/nf-anytime-20261106/`에 작은 결과·한국어 보고서·재현 명령·출처/해시를 보존한다. 큰 입력·스냅샷·로그는 `D:/DivertorOptimization/studies/nf-anytime-20261106`에 둔다.
- [ ] 원기록 독립 감사, 새/기존 관련 회귀 테스트, 그림 육안 검사 후 로컬 커밋한다. FIREFLY 추가 교차 검증·연속 횡확산·1800초·다른 열원으로 해결 범위를 확대하지 않는다.
