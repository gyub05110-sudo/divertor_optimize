# 가상 기준 운전점에서의 수송·형상 비교

기존 direct와 2026-10-27 broadening 계산을 보존하고 별도 실행기 `compare_reference_transport.py`를 추가했다. 실험값이 없어도 동일한 가상 입력·경계조건에서 후보를 비교하는 연구용 모델이다.

2026-10-28 실행 결과는 [해석과 한계](../input/D3D/optimized/reference-comparison-20261028/INTERPRETATION.md), [수치 결과와 그림](../input/D3D/optimized/reference-comparison-20261028/REPORT.md)에 있다. IMP의 기준점 피크 감소 13.31%는 수치 검증을 통과했지만, 폭 고정 대조군과 같은 형상을 선택했고 일부 운전점·분배 가정의 민감도에서 부적합했다. OMP 수송 모델은 기준 형상부터 부적합해 최적화에서 제외했다.

## 고정 기준점

| 입력 | 값 |
|---|---:|
| 독립 열원 한쪽의 입력 전력 | 1 MW |
| 상류 전자밀도 | 1e19 m⁻³ |
| 전자 횡방향 열확산계수 | 0.3 m²/s |
| 상류 exponential 폭 | 3 mm |
| 이온 | deuterium, Ti=Te 근사 |
| sheath 열전달계수 | 7 |
| 평행 열전도 계수 κ0 | 2000 W m⁻¹ eV⁻⁷ᐟ² |

IMP와 OMP는 각각 1 MW인 독립 시나리오다. 실제 장치 총전력이나 실측 상태를 뜻하지 않는다. 상류/표적 온도를 둘 다 입력하지 않으며 온도는 경로별 에너지·sheath 조건에서 계산한다. chi는 전자 열용량 1.5 e n에 대한 확산계수다.

## 무엇을 계산하는가

1. 후보 벽의 실제 첫 충돌점을 찾는다.
2. 자속관의 q_parallel/B를 보존하며 실제 3차원 경로에서 평행 열전도를 적분한다.
3. 상류 밀도, 압력 근사 및 sheath 조건으로 Tu, Tt, nt를 계산한다.
4. 경로의 온도·밀도·자기선속 좌표 변환을 이용해 횡방향 퍼짐 폭을 적분한다.
5. 같은 연속 벽 가지에 Gaussian 전력을 배분하고 면적·입사각·차폐를 평가한다. 미배분 전력은 따로 남기며 0.1%를 넘으면 부적합이다.

`direct`: 원래 직접 매핑의 3 mm 조건. `frozen_transport`: 기준 형상에서 구한 경로별 sigma를 고정. `thermal`: 후보마다 온도와 sigma를 재계산. 동일 후보의 세 모델 교차 평가로 단순 분포 평활화와 형상에 따른 수송 반응을 구분한다. 여기서는 3 mm 단일 조건을 최적화하므로 과거 다섯 폭 minimax 결과와 수치를 직접 비교하지 않는다.

## 실행

프로젝트 폴더에서 새 출력 경로를 지정한다.

```powershell
python compare_reference_transport.py --mode preflight --source both --output input/D3D/optimized/my-reference-preflight
python compare_reference_transport.py --mode optimize --source both --output input/D3D/optimized/my-reference-study
python verification/validate_reference_transport.py --output input/D3D/optimized/my-reference-study
```

`--power-mw`, `--density`, `--chi`로 별도의 기준점을 명시할 수 있다. 다른 기준점은 새 출력 폴더에서 실행한다. `--smoke`는 짧은 연결 시험이며 개선 검증에 쓰지 않는다. 출력을 덮어쓰지 않는다. 부분/완료 validation도 덮어쓰지 않는다.

실제 적용된 전체 reference, 수치 설정, seed와 후보 계수는 각 source의 `manifest.json`에 있다. 검색 전에 입력·코드 해시를 `provenance.json`에 고정한다. 코드를 변경한 뒤 과거 실행에 검증을 붙이려 하면 해시 검사에서 거부한다.

## 결과 읽기

- `REPORT.md`: 기준점, 후보/모델 교차 평가, 수치 판정, 형상·열부하·폭·온도 그림, 운전점 민감도.
- `manifest.json`: 탐색 원자료와 후보 계수. 비율은 같은 모델의 기준 형상 대비다.
- `evaluations.jsonl`: 모든 후보의 점수, 전력 분배와 수송 진단.
- `validation.json`: 320/640/1280경로, quintic640 교차 검증과 고정 후보의 전력/밀도/chi 절반·두 배 평가.
- `*.npz`, `*_refined.npz`: wall, q_over_P, 절대 q_MW_m2, 실제 침적 sigma, 경로별 pre-redistribution Tu/Tt/nt 및 길이.

frozen 모델의 `thermal` 진단은 해당 후보에서 다시 계산했을 경우의 상태다. 실제 침적에 사용한 폭은 `deposited_width_statistics` 또는 NPZ의 `sigma_m`이며, `thermal_sigma_m`과 구분한다.

수치 개선 확인은 원래 및 모든 검증 격자의 적합성, 독립 적분과 격자 수렴을 요구한다. 민감도 표는 단일 해상도 고정 후보 재평가이며 새로운 최적해나 전 격자 검증이 아니다. 실패한 행의 낮은 피크를 개선으로 해석하지 않는다.

## 연결 길이가 폭에 들어가는 방식

입력으로 고정하는 것은 전력·상류 밀도·수송계수와 경계조건의 **규칙**이다. 모든 후보의 온도나 sheath 밀도까지 같게 강제하지는 않는다. 이들은 실제 충돌점과 경로 자기장에 따라 다시 계산된다.

경로 거리 s는 3차원 자기장선 거리이고, B는 자기장 크기다. 경로별 상류 열유속 q_u는 동일한 exponential 열원에서 정한다. 전력 손실 없는 자속관 근사에서는 다음을 사용한다.

```text
q_parallel(s) / B(s) = q_u / B_u
J = integral_0^L [B(s)/B_u] ds
T_u^(7/2) - T_t^(7/2) = (7/2) q_u J / kappa0
n_t = n_u T_u / (2 T_t)
q_t = gamma e n_t T_t sqrt(2 e T_t / m_i)
```

온도는 eV, e는 전하량이다. bulk 밀도는 n(s)=n_u T_u/T(s)로 근사하며, 마지막 bulk 밀도와 별도 sonic sheath 밀도의 factor-two 차이를 경로에 걸쳐 보간하지 않는다. 전자 열용량 기준 평행 확산계수와 횡방향 퍼짐은 다음과 같다.

```text
chi_parallel(s) = kappa0 T(s)^(5/2) / [1.5 e n(s)]
g(s) = |grad psi(s)| / |partial_R psi at source separatrix|
S = integral_0^L g(s) sqrt[chi_perp / chi_parallel(s)] ds
sigma = S / sqrt(2)
```

S는 사용한 전도 broadening 근사의 폭 관례이고, sigma는 실제 Gaussian 커널의 표준편차다. 따라서 긴 경로가 직접 적분 구간을 늘리는 효과와 온도·밀도·자속 좌표 변환을 바꾸는 효과가 함께 들어간다. 길이만으로 모든 형상에 단조로운 이득을 보장하지 않는다. 이 모델에서 별도의 유동 속도나 입자 체류시간은 풀지 않는다.

## 한계

이 모델은 **one-way effective two-point closure**다. 계산된 Tt는 횡방향 전력 재분배 전의 경로 상태이며, 재분배 후 열유속으로 sheath 온도를 다시 풀지 않는다. 입자·운동량·유동·복사·중성입자·detachment를 풀지 않고, 자기장도 고정이다. 따라서 물리적으로 완결된 2D 수송 또는 Super-X 자기배위 최적화가 아니다.

낮은 전력 꼬리에서 극저온/큰 밀도가 생길 수 있다. 온도 하한·폭 상한으로 숨기지 않으며, Tt<5 eV 전력 비중, 작은 전도 온도강하 비중, free-streaming 비율, sigma/lambda를 진단한다. 이런 지표가 나쁘면 물리적 적용 범위가 좁다는 뜻이지, 수치식의 해가 실험적 detachment를 입증한다는 뜻은 아니다.

설계·리뷰·실행 기록은 `docs/superpowers/specs/2026-10-28-reference-transport-design.md`, `docs/superpowers/plans/2026-10-28-reference-transport.md`, `docs/superpowers/plans/2026-10-28-reference-reviews-execution.md`에 있다.

참고: [Nille–Sieglin–Eich conductive broadening](https://arxiv.org/html/1610.04148), [JET-ILW two-point power/momentum conventions](https://www.sciencedirect.com/science/article/am/pii/S2352179123000649). 전 경로와 가변 자기장에 대한 적용은 이 구현의 근사 확장이다.
