# 수송 분포 peak 구간 검증과 3계수 IMP 탐색

## 결과

IMP, Rmin=20 mm, 독립 입력 1 MW에서 세 형상 계수를 함께 변화시키는 시범 탐색과 고정 수송 분포의 연속 peak 검증을 완료했다. **125개 표본 중 19개가 수치 판정을 통과했고, 가장 좋은 형상은 기존 검증 후보 그대로였다.** 새로운 개선 형상이나 형상 전역 최적성 증명을 얻은 것은 아니다.

### 형상 탐색

- 실제 D3D 기준 형상과 기존 수송 입력을 유지했다.
- IMP offset 계수 인덱스 0, 1, 2를 각각 ±0.5 mm, 5단계로 변화시켰다. 나머지 계수와 접합각은 고정했다.
- 공통 접선·곡률·재료 길이·전체 벽 판정과 수송 평가를 사용했다.
- 탈락 106개 중 92개는 연속 형상 제약, 14개는 재료 길이 제한에서 탈락했다.
- 통과 후보의 320 경로 peak 범위는 2.375044182~2.387526606 MW/m²였다. 가장 좋은 값은 기존 후보와 동일했다.
- 이 표본 결과로 표본 사이의 더 좋은 해나 다른 변수 방향의 개선 가능성을 배제할 수 없다.

### 고정 수송 분포의 연속 peak 검증

| 경로 수 | 수치 peak [MW/m²] | 검증 구간 [MW/m²], 바깥쪽 반올림 | 검증 시간 |
|---:|---:|---:|---:|
| 320 | 2.375044182 | [2.375044, 2.376031] | 124.5초 |
| 640 | 2.375136847 | [2.375136, 2.376119] | 317.4초 |
| 1280 | 2.375134353 | [2.375134, 2.376132] | 500.1초 |

세 경우 모두 원래 정밀도로 기록한 구간 폭이 0.001 MW/m² 이하라는 종료 조건을 만족했다. 1280 경로의 정확한 저장 구간은 [2.37513433552238417998170636251989872673813, 2.37613105556017023194792747363266371842277]이다.

이 구간은 **계산된 도착 위치·폭·전력을 정확한 고정 상수로 해석한 무한 주기 Gaussian 혼합 분포의 공간 최대값**을 감싼다. 기존 수치 peak가 구간 안에 들어오는 것을 확인했다. 구간 계산에는 방향성 반올림, 주기 경계의 생략한 Gaussian 꼬리, 전체 수신 좌표의 분할이 포함된다.

640→1280 경로 수치 peak 변화는 약 0.000002494 MW/m²이다. 모든 경우 미도착·미해결 전력은 0, 전력 장부 오차는 최대 6.66e-16이었다. 이는 수치 관측이며 실제 수송 모델 전체의 오차를 엄밀하게 제한한 결과는 아니다.

**이번에는 source 경로 수만 320/640/1280으로 바꿨다.** 경로 추적 step=0.0005 m, target nodes=161, receiving wall step=0.001 m, cubic 보간은 고정했다. 기존 최종 검증의 1280 설정(step=0.000125 m, nodes=641, wall step=0.00025 m)과 다르므로, 기존 저장 peak 2.375184231을 이번 2.375134353으로 대체하지 않는다. 같은 이번 설정의 기준 peak 2.549656690 대비 감소율은 약 6.84493%이다.

### 수송 정보를 이용한 조건부 하한

1280 경로의 고정 donor 전력과 전체 면적을 유지하고, 각 확산 폭은 ±1%, 도착 위치는 해당 donor의 원래 폭의 ±1%만큼 변화한다고 **가정한** 분포 집합을 정의했다. 도착 좌표는 [0, 전체 면적] 범위로 제한했다.

이 가정된 집합의 모든 분포에 대해 peak는 다음 구간에 포함된다.

```
2.31255146056807065846737502833988603379788
  <= peak <=
2.56626605150330743137049929465377839161437  MW/m²
```

64회 분할 예산에 도달해 멈췄으며 0.001 종료 폭은 달성하지 못했다. 그럼에도 현재의 넓은 구간은 해당 가정 아래에서 유효하다.

**이 2.31255를 형상 최적화 문제의 하한으로 사용할 수는 없다.** 실제 형상 변수 구간의 모든 도착 위치·폭·전체 면적이 이 가정된 수송 범위에 포함된다는 증명이 아직 없기 때문이다. 현재 구조는 그러한 포함 관계가 별도로 증명되었을 때 수송 하한을 계산할 수 있는 기반이다. `family_covers_geometry_box_proven=false`, `certified_geometry_optimality_gap=null`로 저장했다.

## 구현과 확인

- `d3d_optimization/transport_peak_bounds.py`: 방향성 Decimal 구간, 무한 주기 꼬리 경계, 넓은 분포의 Fourier 경계, 전체 공간 peak 분할 및 재검증.
- `verification/run_transport_peak_pilot.py`: 3계수 탐색, 가장 좋은 수치 후보 저장, 경로 수 비교, 고정/조건부 분포 인증, 코드·입력 스냅샷.
- 상세 식과 증명 범위: `docs/transport_peak_enclosure.md`.
- 관련 테스트 76개 통과. 독립 검토에서 90자리 고정밀 비교 320개 및 정확한 유리수 산술 비교 1,000개도 통과했다.
- 검토 중 경계 끝점 c=A가 방향성 나눗셈 때문에 거부되는 문제를 재현하고 수정했다. 실제 실행 좌표는 modulo 구간 안에 있어 이 문제의 영향을 받지 않았다.
- 반복되는 exp(-100) 상한을 캐시했다. 고정 상수의 동일성이 확인됐으며 결과 정의는 유지된다.
- 과거 검증 후보·입력·스냅샷 900개를 실행 전후 확인했고, campaign/provenance JSON 3개의 해시도 별도 기록했다. 기존 qualified store 승격은 수행하지 않았다.

초기 320 경로 준비는 8.61초, 시범 실행 전체는 약 1,123.7초(18.7분)였다. 구간 검증은 CPU에서 실행했다. 탐색 내부에 이 비용을 추가하지 않고, 선별된 최종 후보에 적용하는 후처리로 분리하는 것이 적합하다.

## 저장과 재현

원본 데이터 및 실행 당시 소스:
`D:/DivertorOptimization/studies/transport-peak-pilot-20261102`

`settings.json`, `observations.json`, `best_numerical.json`, `resolution.json`, 세 `kernel_*.json`과 `certificate_*.json`, `conditional_family_certificate.json`, `summary.json`, 입력·런타임 해시 및 `source/`를 보존했다. 실행 후의 끝점 처리·메타데이터·캐시 보완은 저장소 코드에 반영했고, 실제 실행 소스 스냅샷은 그대로 남겼다.

```powershell
python -m verification.run_transport_peak_pilot --replay D:/DivertorOptimization/studies/transport-peak-pilot-20261102/certificate_320.json
python -m pytest verification/test_transport_peak_bounds.py verification/test_optimality_bounds.py verification/test_campaign_incumbents.py verification/test_campaign_consistency.py verification/test_global_campaign.py verification/test_tangent_geometry.py verification/test_global_certificate_geometry.py verification/test_run_provenance.py -q
```

그림은 저장된 kernel/certificate에서 생성했다. 그림의 구간 표시는 일반 반올림한 라벨이며, 엄밀한 끝점은 JSON과 위 표의 바깥쪽 반올림 값을 사용한다.

## 남은 핵심 작업

형상 변수 구간 → 첫 도착 경로/위치의 구간 → 연결 길이·열수송 폭의 구간 → 수신 면적의 구간을 오차까지 포함해 연결해야 한다. 그 다음에야 이번 수송 하한을 형상 최적화의 하한으로 사용할 수 있다. 물리 모델 전체의 후보 상한도 같은 오차 체계에서 검증해야 한다. 이번 결과는 그중 **고정된 수송 입력 이후의 연속 peak 계산**을 검증한 단계다.
