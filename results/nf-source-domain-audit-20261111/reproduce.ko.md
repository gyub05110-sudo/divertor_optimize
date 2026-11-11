# 열원·계산 영역 감사 재현

작업 폴더: `D:/DivertorOptimization/worktrees/d3d-report-correction-20261104`.
원본 데이터: `D:/DivertorOptimization/studies/nf-stage5-firefly-20261110-r2`.
새 진단 데이터: `D:/DivertorOptimization/studies/nf-source-domain-audit-20261111`.

원본 5단계의 실행 봉인과 출력은 수정하지 않는다. 아래 검사는 저장된 배열과 원시 NetCDF를 읽으며 열계산·최적화를 새로 실행하지 않는다. Python의 `-O` 옵션은 사용하지 않는다.

```powershell
$env:PYTHONPATH=(Resolve-Path '.local_dependencies/d3d-fixed-equilibrium').Path
$env:PYTHONIOENCODING='utf-8'
$env:PYTHONDONTWRITEBYTECODE='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
& 'C:/Users/FPL-user/miniforge3/python.exe' -B -m verification.source_domain_audit --study D:/DivertorOptimization/studies/nf-stage5-firefly-20261110-r2 --output D:/DivertorOptimization/studies/nf-source-domain-audit-20261111/source-domain.json
& 'C:/Users/FPL-user/miniforge3/python.exe' -B -m pytest verification/test_independent_sources.py verification/test_regional_transport.py verification/test_stage5_inputs.py verification/test_stage5_native.py verification/test_stage5_transport_pin.py -q
```

이번 실행 결과는 7개 OMP 실험, 46개 형상, 276개 전력 장부 검산 통과와 테스트 20개 통과다. `source-domain.json`의 `driver_sha256`는 감사 코드 자체를 식별한다. 각 `raw_sha256`은 검사한 원시 출력의 해시다. 개별 열원 파일, 물질벽, 원래 축소 모델 형상/열 데이터의 해시도 기존 봉인과 대조한다.

이 검사는 현재 7개 평형의 OMP 실험에 한정된다. 이전 IMP 실행의 존재는 과거 기록으로 확인했으며, 현재 7개 평형의 IMP 계산까지 검증했다는 뜻이 아니다.

## 경로 진단 보존본

`native-evidence.zip`에는 실제 사용한 진단 소스·실행 스크립트, 4개 랭크의 4,096개 이력, 최종 NetCDF, 해시 목록, 실패한 진단 시작의 로그와 설명을 함께 보존했다. 원본 FIREFLY/FLARE 라이선스도 포함했다. `native-summary.json`은 압축본 `summary.json`의 동일 사본이다.

재현 명령과 필요한 기존 WSL 실행 환경은 압축본 `README.md`에 있다. 파일은 새 진단 데이터 폴더의 `native` 하위 구조에 대응한다. 개발 중의 생성기를 다시 실행하는 대신 **최종 `instrumented_all` 소스**와 문서에 명시한 빌드 명령을 사용한다. 실제 진단 실행 시간은 317.406초, 누적 CPU 시간은 1,263.340초였다. 원래 열계산과 동일한 시드·형상·열원·수송 조건을 사용했고, 최종 NetCDF의 전체 바이트가 기존 파일과 일치했다.

전체 이력의 영역 판정은 양성 방문 기록이며 연속 경로의 모든 점을 저장하지 않는다. 경계 조건을 변경한 민감도 실험, 새 IMP 최적화 또는 새 형상 간 비교는 수행하지 않았다.
