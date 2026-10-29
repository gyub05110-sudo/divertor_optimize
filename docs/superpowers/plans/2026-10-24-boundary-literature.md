# IMP 경계 문제: 문헌 조사 (2026-10-24)

## 조사 결과

1. Dekeyser, Reiter & Baelmans (2014), *Divertor target shape optimization in realistic edge plasma geometry*, Nuclear Fusion 54, 073022. DOI: https://doi.org/10.1088/0029-5515/54/7/073022 . 확인한 저자 제공 본문: https://www.researchgate.net/publication/262150205_Divertor_target_shape_optimization_in_realistic_edge_plasma_geometry . 타깃 끝점을 고정하여 주변 구조와 연결하고, 자기장선을 따른 형상 자유도와 separatrix 부근 해상도를 다룬다. 해당 논문의 플라즈마 유체 모델은 현재 코드의 prescribed direct mapping과 다르므로 정량 결과를 가져오지 않는다.
2. Frerichs (2025), *Particle swarm optimization of divertor targets for heat load control*, preprint: https://arxiv.org/html/2509.00206v1 . 타깃의 위치와 방향을 탐색하고 floor/baffle 열부하도 비용에 포함한다. 현재 코드에서도 선택하지 않은 벽으로 열을 보내는 것을 개선으로 취급하면 안 된다는 설계 근거다. PSO로 알고리즘을 바꾸라는 근거는 아니다.
3. *Stellarator island divertor shape optimization for reduced peak heat fluxes* (2026), preprint: https://arxiv.org/html/2602.24049v1 . 고정 자기장 속 형상 최적화, 매끄러운 연결, 자기 차폐와 다른 표면의 열부하를 다룬다. Stellarator 및 횡방향 확산 모델이므로 DIII-D 무수송 모델의 실험 검증 자료가 아니다.
4. Scarabosio et al., *Scaling of the Divertor Power Spreading (S-factor) in Open and Closed Divertor Operation in JET and ASDEX Upgrade*: https://scipub.euro-fusion.org/archives/jet-archive/scaling-of-the-divertor-power-spreading-s-factor-in-open-and-closed-divertor-operation-in-jet-and-asdex-upgrade . PFR 방향 spreading S와 upstream decay width λq는 구별된다. 현재 오류를 고치기 위해 Gaussian 확산을 추가하면 사용자가 복원한 물리 모델 자체가 바뀐다.

## 코드에 적용하는 판단 (문헌의 직접 결론과 구분)

- 샘플 첫 충돌점을 설계 경계로 사용하지 않고 ψ=ψX의 연속 strike를 계산한다.
- 수송, 열원, 24개 B-spline 자유도, Sobol+DE 방식은 유지한다.
- 끝점 접선을 고정한 채 경계를 strike에 맞추면, 그 점의 SOL측 열유속 극한이 변하지 않을 수 있다. 이 값도 목적함수에서 검사해야 한다.
- 주변 벽에서의 first hit를 검사한 후 경계 극한값을 계산한다. 가려진 과거 strike에 열을 계속 부여하지 않는다.
- 5 mm PFR 쪽 형상 지지구간은 본 작업의 설계 제안이지 문헌이 보장한 길이가 아니다. 사용자 응답 전에는 적용하지 않는다.
- 열유속 실측 자료가 없으므로 가능한 것은 수치 검증이다. 확률적 탐색 및 유한 해상도 검사는 전역 최적해 증명이 아니다.

## 확인된 원인

이전 IMP 320 경로 결과는 약 2% 개선, 640에서는 약 1%, 1280/2560에서는 0%였다. 69번 벽의 연속 strike와 샘플 기반 설계 시작점 사이에는 약 0.9823 μm의 고정 구간이 남았다. 현재 spline은 양 끝 변위와 기울기를 모두 0으로 고정한다. 이 두 조건을 각각 다뤄야 한다. 이전 생산 결과는 수정하거나 덮어쓰지 않는다.
