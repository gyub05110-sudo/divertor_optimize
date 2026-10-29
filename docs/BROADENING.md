# 연결 길이 의존 broadening 실행 안내

기존 실행기와 계산 모듈은 변경하지 않았다. 새 모델은 `compare_broadening.py`로 별도 실행한다. `sigma_ref=0`이면 기존 continuous-boundary 직접 매핑 모델을 그대로 호출한다.

## 모델

- 상류 지수 열원, λq=1/2/3/5/10 mm와 기존 단위 전력을 유지한다. IMP/OMP는 독립 계산이다.
- 연결 길이는 `integral |B|/Bp ds_p`, 중간면에서 실제 최초 벽 충돌까지의 3차원 길이다. 단면상의 leg 길이를 그대로 쓰지 않는다.
- 퍼짐 표준편차는 `sigma_ref * (L / L_ref)**exponent`이다. 표준값 `L_ref=10 m`, `exponent=1`이며 `0.5`도 지원한다. 폭은 upstream-equivalent 자기선속 좌표의 길이로, 물리적 확산계수나 Eich S가 아니다.
- 같은 연속·단조 자기선속 벽 가지에서 Gaussian 전력을 적분한다. 벽의 토로이달 면적과 입사 방향을 Jacobian에 한 번 반영한다.
- 고정 벽을 포함한 전체 벽을 평가한다. 같은 벽 가지의 PFR 쪽으로 침적 지지를 확장할 수 있으나 PFR 수송 방정식을 푸는 것은 아니다.
- 벽 뒤로 통과하는 재분배를 막기 위해 5 mm poloidal 상류점에서 받는 벽 구간 전체까지의 삼각형 영역을 검사한다. 이는 근사 모델의 국소성 제약이며 실제 확산 경로를 입증하지 않는다.
- 차폐·가지 바깥·미충돌 전력을 재정규화하지 않는다. 미해결 전력이 0.1%보다 크면 부적합이다. 강한 퍼짐에서는 이 근사의 적용 범위가 제한될 수 있다.
- 자기장은 고정돼 있다. plate 형상이 바뀌며, Super-X 자기배위를 새로 만드는 계산은 아니다.

## 실행

프로젝트 폴더에서 실행하며 출력 폴더는 새 이름을 사용한다.

```powershell
python compare_broadening.py --mode preflight --source both --output input/D3D/optimized/my-broadening-preflight
python compare_broadening.py --mode optimize --source both --tubes 320 --nodes 161 --step-mm 1 --wall-step-mm 2 --output input/D3D/optimized/my-broadening-study
python verification/validate_broadening.py --output input/D3D/optimized/my-broadening-study
```

기본 비교 폭은 `--sigma-mm 0 0.1 0.25 0.5` mm이다. `--exponent 0.5`, `--length-ref-m 10`, `--anchor-mm 5`를 명시적으로 바꿀 수 있다. 실제 사용 인수는 manifest에 저장한다.

탐색 기본 예산은 Sobol 128개, DE population 24, 최대 8세대이다. 검증된 과거 직접 매핑 후보를 모든 모델의 동일 초기 후보로 포함한다. 더 작은 `--bound-mm`를 요청해 과거 후보가 범위를 벗어나면 자르지 않고 제외 사유를 기록한다. 이 탐색은 제한된 비교 실험이며 전역 최적화 인증이 아니다.

`--smoke`는 48경로·8개 Sobol·DE 1세대의 연결 시험이다. 수치/물리 개선의 증거로 사용하지 않는다. 기존 출력, 완료/부분 검증 결과는 덮어쓰지 않는다.

## 결과와 검증

- `REPORT.md`: 비교 표, 형상과 열분포 그림, 적용 범위.
- `IMP/manifest.json`, `OMP/manifest.json`: 설정, 기준값, 후보 계수, 모든 후보의 모델 간 교차 평가.
- 각 모델의 `evaluations.jsonl`: 평가한 모든 후보와 부적합 사유.
- `provenance.json`: 실행 전에 고정한 코드와 실제 입력의 SHA256.
- `validation.json`: 320/640/1280경로 및 quintic640 재평가, 독립 최초 충돌·연결 길이·면적 적분·입사각·고정 벽 검사.
- `*_refined.npz`: 가장 미세한 cubic 격자에서의 벽과 열분포.

검증기는 기존 직접 모델의 0전력 separatrix 극한 q를 양의 broadening에 섞지 않는다. 대신 새 연속 밀도의 최대값을 양쪽 끝점과 보강 샘플에서 평가한다. 퍼짐이 작아 개별 source ray 피크가 나타나면 source/벽 격자 수렴 검사가 실패할 수 있으며 이를 보존한다.

수치 개선 인정 조건은 원래/모든 재평가 수준에서 적합, 최대 비율 변화≤0.002, 절대 피크 상대 변화≤0.005, 전폭 악화 없음, 감소량>관측 비율 변화의 2배+0.0001이다. 독립 밀도 적분/CDF 전력 차이는≤1e-4, 연결 길이 상대 차이는≤0.005, flux drift/sigma는≤0.05여야 한다. 가장 거친 검증 수준도 실패하면 개선 확인으로 올리지 않는다.

2.5/10 mm 상류 anchor와 다른 길이 지수의 민감도는 별도 항목이다. 수치 수렴과 수송 가정의 타당성은 다르다. 검증 완료가 실험 검증 또는 전역 최적해를 뜻하지 않는다.

2026-10-27 완료 비교는 `input/D3D/optimized/broadening-comparison-20261027/REPORT.md`에 있다. 기본 2 mm 추적 검증에서 연결 길이 오차가 기준을 초과하여, 이 실행에서는 `refinement_plan.json`의 더 세밀한 추적 간격을 모든 후보에 적용했다. 임계값은 바꾸지 않았다. 정확한 실행 명령은 `docs/superpowers/plans/2026-10-27-broadening-execution.md`에 기록했다.

해당 실행의 계산 완료 후 보고서 출력 코드만 수정했다. 당시 코드 사본과 `presentation_amendment.json`을 결과 폴더에 보존했으며 과거 provenance는 그대로다. 따라서 현재 출력 코드로 과거 실행의 엄격한 해시 검사를 통과시키려 하지 말고, 정확한 재현에는 별도 사본에서 보관된 출력 코드를 사용하거나 새 실행 폴더를 생성한다.

## 개발·검토 기록

설계: `docs/superpowers/specs/2026-10-27-connection-broadening-design.md`

구현 계획: `docs/superpowers/plans/2026-10-27-connection-broadening.md`

적대적 검토와 수정: `docs/superpowers/plans/2026-10-27-broadening-reviews.md`

원본 16개 production 파일 해시: `docs/superpowers/plans/2026-10-27-broadening-preservation.json`

```powershell
python -m pytest verification -q
```
