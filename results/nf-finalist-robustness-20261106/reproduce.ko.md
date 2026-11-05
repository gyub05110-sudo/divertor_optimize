# 최종 후보 검증을 다시 실행하는 방법

저장소 루트에서 원본 DIII-D 입력 및 고정 평형 의존성을 준비한다. 앞선 캠페인의 경로별 원자료를 사용하는 시험에는 [2026-11-05 아카이브 복원 안내](../research-archive-20261105/README.md)가 적용된다. 모든 명령은 새 출력 경로를 사용한다. 1,400 ms의 수치 교차검사와 두 추가 평형의 열계산은 GPU로 실행했다. 수치 교차검사는 네 설정 각각에서 OMP20, χ⊥=0.3의 동일 후보를 CPU로 다시 계산해 피크와 적격성을 확인한다.

```powershell
$trial = 'D:/DivertorOptimization/studies/nf-robustness-replay'
$env:PYTHONPATH = "$PWD/.local_dependencies/d3d-fixed-equilibrium"
$env:OPENBLAS_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:MKL_NUM_THREADS = '1'
foreach ($method in 'cubic','quintic') {
    foreach ($count in 640,1280) {
        python -m verification.nf_finalist_resolution --output "$trial/resolution" --count $count --method $method --backend gpu
        if ($LASTEXITCODE -ne 0) { throw 'profile generation failed' }
    }
}
python -m verification.nf_finalist_decisions --profiles "$trial/resolution" --output "$trial/resolution/decisions.json"
python -m verification.nf_equilibrium_times --output "$trial/times" --phase geometry
python -m verification.nf_equilibrium_times --output "$trial/times" --phase heat
python -m verification.nf_equilibrium_time_audit --study "$trial/times" --output "$trial/time-audit.json"
python -m pytest -q verification
```

[수치 대조 규칙](resolution-protocol.ko.md)과 [고정 후보 벡터](../../verification/fixtures/nf_robust_pairs.json)에는 비교 상대 및 판정법이, [추가 평형 규칙](time-protocol.ko.md)에는 시각·재료 벽·탐색 시드와 실패 처리법이 있다. `resolution-index/` 네 파일에는 이번 실행의 40개 원프로파일 파일명, SHA-256, 실제 사용한 전력·수송 인자와 수치 설정이 있다. 결과 파일의 피크는 저장된 연속 프로파일에서, A₂는 각 셀의 실제 면적과 전력에서 검산했다. [네 수치 설정의 판정](resolution-decisions.json) · [평형과 형상의 독립 기하 검산](time-audit.json)

이번 실행의 대용량 원프로파일은 `D:/DivertorOptimization/studies/nf-finalist-robustness-20261106/resolution/`에 있다. `times/1200/` 및 `times/1600/`의 작은 결과는 이 보고서와 함께 보존한다. 각 단면의 기하 전처리 출력, 정밀 계산 전에 선택한 세 후보 번호, 12개 전수 점수와 전력수지가 포함된다. 최초 도달면, 저온 꼬리, 국소 오차 분해, 무작위 선별 대조에 대한 별도 재현 명령은 각각의 세부 보고서에 있다.
