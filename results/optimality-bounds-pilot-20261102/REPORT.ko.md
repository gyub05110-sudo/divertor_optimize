# 최적성 하한 추적: 1변수 IMP 시범 결과

## 최적값을 미리 알아야 하는가?

아니다. `x`는 형상 계수, `f(x)`는 그 형상의 최대 열유속, `f*`는 허용 형상 전체에서 가능한 최소 최대 열유속이다. 실행 가능한 후보 `x_hat`의 성능에 대한 검증된 상한을 `U`, 모든 허용 형상에 공통으로 적용되는 검증된 하한을 `L`이라고 하면 `L <= f* <= U`이다. `U-L`은 최적값을 모르더라도 계산할 수 있다.

예를 들어 U=2.40, L=2.30 MW/m²이면 실제 최적값은 그 사이에 있다. L=2.30을 달성하는 형상이 발견됐다는 뜻은 아니다. 탐색으로 U를 낮추고, 증명으로 L을 높여 허용 오차 이내로 만나게 하는 것이 목표다. 최적값은 숫자이고 최적 형상은 그 값을 달성하는 형상이므로, 숫자만 안다고 형상이 자동으로 주어지는 것도 아니다. 일반적인 부동소수점 계산에서는 정확한 등식보다 명시한 오차 안의 최적성 인증을 목표로 삼는다.

## 이번에 구현한 범위

- 기존 탐색·최종 판정 정합성과 검증 후보 보존 코드는 이전 커밋 d5c47e1에 보존했다.
- 새 모듈은 고정한 spline 상수와 벽 좌표를 정확한 유리수로 해석한다. 구간 de Boor 계산으로 곡선의 반경과 속도를 감싸고, 회전 표면적 상한 A_upper를 계산한다.
- 열유속이 비음수이며 전력 P가 전부 이 표면에 도착한다면 `peak >= P/A >= P/A_upper`이다. 형상 제약을 일부 완화해 계산해도 실제 허용 형상 집합의 하한으로 사용할 수 있다.
- 계수 구간을 이분할하고 양쪽 자식 구간을 모두 보존한다. 전역 하한은 남아 있는 구간 하한들의 최솟값이다. 샘플링한 값으로 구간을 제거하지 않는다.
- `bounds.json` 재검증은 입력·전체 분할·하한을 재계산한다. 수치 peak를 인증된 상한으로 바꾸는 기능은 없다. `certified_gap_MW_m2=null`, `global_optimum_proven=false`로 기록한다.

이것은 **정확한 고정 spline 벽 + 전량 전달 가정에 대한 하한**이다. 운영 중인 부동소수점 벽 생성, 자기장 보간, 경로 추적, 열수송 및 연속 peak 계산의 오차 경계까지 증명한 결과는 아니다. 현재 기능은 하한 세분화이며 완성된 수송 모델 전역 최적화 인증기는 아니다.

## 실제 실행

기준: D3D 원형 유지, 독립 IMP 1 MW, Rmin=20 mm, 기존 검증 후보에서 시작. IMP 계수 9개 중 인덱스 2인 offset만 ±0.5 mm 변화시키고 나머지 계수 및 접합각을 고정했다. 기존 캠페인을 재시작하지 않고 별도 D 폴더에서 320개 경로·9개 표본만 계산했다. 이 작은 실행은 CPU로 14.6초 소요됐다.

| 항목 | 결과 |
|---|---:|
| 탐색한 표본 | 9 |
| 수치 판정을 통과한 표본 | 1 (기존 후보) |
| 수치 후보 peak | 2.375044182 MW/m² |
| 처음의 조건부 모델 하한 | 0.01348487644 MW/m² |
| 16회 이분할 후 조건부 모델 하한 | 0.01348572638 MW/m² |
| 남은 구간 | 17 |
| 인증된 상한·최적성 gap | 없음 |
| 기존 입력/스냅샷 해시 보존 확인 | 876개 |

음의 방향 표본 4개는 곡률, 양의 방향 표본 4개는 재료 길이 제약에서 탈락했다. 샘플 사이의 다른 허용 형상 존재 여부를 증명한 것은 아니다. 수치 후보의 전력 장부는 미도착 0, 미해결 0, 보존 오차 1.11e-16이다. 이 값은 수치 검증 결과이며 정확한 산술의 전력 보존 증명은 아니다.

전체 벽 면적만 이용한 하한은 매우 느슨하다. 수치 후보와 하한의 차이 2.36156 MW/m²도 **인증된 gap이 아니다**. 이번 결과는 새 최적 형상을 찾거나 기존 형상의 최적성을 입증하지 않았다. 한 변수로 제한한 결과를 IMP 전체 9변수나 OMP 16변수에 확대 해석할 수 없다.

## 다음에 강화할 부분

1. 제약 경계에서 여러 계수를 함께 움직이는 2~3변수 부분문제를 정의한다.
2. 열이 실제로 도달할 수 있는 영역, 전달 전력, 확산 폭의 구간 경계를 이용해 하한을 강화한다. 단순히 낮은 수치 표본값을 하한으로 쓰면 안 된다.
3. 동일한 모델에서 후보의 형상 허용성과 연속 peak 상한을 수치 오차까지 포함해 검증한다.
4. 이 조건이 갖춰져야 `L <= f* <= U`와 인증 gap을 출력한다. 그 후에 더 큰 변수 공간으로 확장한다.

## 저장·재현

데이터: `D:/DivertorOptimization/studies/optimality-bounds-pilot-20261102`.
`settings.json`, `bounds.json`, `observations.json`, `summary.json`, `artifact_sha256.json` 및 실행 당시 소스를 보존한다. 이후 sqrt 극소수/초대형수 처리와 자동 소스 스냅샷 저장을 강화했으며, 저장된 D3D 하한은 현재 코드로도 동일하게 재검증된다. 전체 원본 결과는 덮어쓰지 않는다.

```powershell
python -m verification.run_optimality_bounds --replay D:/DivertorOptimization/studies/optimality-bounds-pilot-20261102/bounds.json
python -m pytest verification/test_optimality_bounds.py verification/test_campaign_incumbents.py verification/test_campaign_consistency.py verification/test_global_campaign.py verification/test_tangent_geometry.py verification/test_global_certificate_geometry.py verification/test_run_provenance.py -q
```

새 테스트 10개와 기존 관련 회귀 테스트 52개, 총 62개 통과. 독립 코드 검토에서 지적된 subnormal sqrt 정체 문제를 재현 후 수정했고, 구간 전체 덮음·변조 거부·기존 실제 형상 면적 포함 여부를 확인했다.

최적성 상하한 및 공간 분할의 일반적인 설명: [Gurobi nonlinear optimization 공식 문서](https://docs.gurobi.com/projects/optimizer/en/current/features/nonlinear.html). 이 문서를 인용한 것이 현재 열수송 모델에 대한 인증을 대신하지 않는다.
