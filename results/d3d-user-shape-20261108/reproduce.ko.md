# 재현과 증거 검사

비교 설계는 `97ffd16427578839bd58b2b182d71122521a1a8a`에서 먼저 고정했다. 이전 연구 방향으로 돌아가는 체크포인트는 `8e63e34971d379d2e6aea64008e2ce3e7890d185`다. [사전 비교 규칙](protocol.ko.md) · [사용자 실행 안내](../../docs/user-optimizer.ko.md)

## 입력과 환경

`1a3155370c`는 저장소에 포함된 공개 DIII-D 파일이다. `e806f47778`의 실제 바이트는 기존 증거 ZIP 안의 `studies/nf-prospective-inputs-20261107/cases/e806f47778/input-root/input/D3D/d3d/d3d_shot_1a3155370c.parquet`에도 있다. 예전 로더가 요구한 이름으로 보관된 파일이다. 이 멤버를 새 위치에 복사해 임의 이름으로 저장한 뒤 아래 `$other`에 지정한다. SHA-256은 실행 스크립트가 검사한다. 두 선반 실패 입력도 같은 ZIP의 각 사례 디렉터리에 있다. [기존 증거 ZIP](../nf-prospective-search-20261107/evidence.zip) · [공개 원 데이터 고정 버전](https://huggingface.co/datasets/Sophelio/fusion-equilibrium-challenge/tree/1e280905b85f2a6fdde7e06fca8cf3a1edf447cb)

본 계산은 Windows·단일 NVIDIA RTX 3090과 연구용 Python을 사용했다. 각 `result.json`에는 실제 버전과 해시가 있다. 다음은 저장소 루트의 PowerShell에서 쓰는 예시다. `$other`와 `$study`는 자신의 경로로 지정한다. 비교 재실행에는 비어 있는 새 결과 경로를 사용한다.

```powershell
$py = 'C:/Users/FPL-user/miniforge3/python.exe'
$other = 'D:/DivertorOptimization/studies/nf-prospective-inputs-20261107/raw/d3d_shot_e806f47778.parquet'
$study = 'D:/DivertorOptimization/studies/reproduce-input-shape'
$env:PYTHONPATH = (Resolve-Path '.local_dependencies/d3d-fixed-equilibrium').Path
$env:OPENBLAS_NUM_THREADS = $env:MKL_NUM_THREADS = $env:OMP_NUM_THREADS = '1'
$env:DIVERTOR_GPU_RUNTIME = 'D:/DivertorOptimization/runtime/cupy'
$env:CUPY_CACHE_DIR = 'D:/DivertorOptimization/cache/cupy'
$env:CUDA_CACHE_PATH = 'D:/DivertorOptimization/cache/cuda'
$env:MPLCONFIGDIR = 'D:/DivertorOptimization/cache/matplotlib'
$driver = './results/d3d-user-shape-20261108/run-study.ps1'
& $driver -Python $py -OtherEquilibrium $other -Study $study -Phase Search
& $driver -Python $py -OtherEquilibrium $other -Study $study -Phase Replay
& $driver -Python $py -OtherEquilibrium $other -Study $study -Phase Analyze
```

검색은 실제 벽시계 예산을 사용하므로 반환 후보가 실행마다 달라질 수 있다. 원래 형상을 검산하려면 사례별 ZIP에서 `full/`·`floor-only/`를 복원한 후 해당 위치를 다음 명령에 넘긴다. 옮긴 평형 파일이나 기하/정책은 해시가 일치할 때만 허용한다.

```powershell
& $py -B -m verification.user_shape_replay --full "$study/e806f47778/full" --floor "$study/e806f47778/floor-only" --equilibrium $other --output "$study/e806f47778/rebuilt-replay" --backend gpu --levels cubic1280 quintic1280 cubic640
& $py -B -m verification.user_shape_results --full "$study/e806f47778/full" --floor "$study/e806f47778/floor-only" --replay "$study/e806f47778/rebuilt-replay" --output "$study/e806f47778/rebuilt-analysis"
```

## 보존 범위

`entry-checks.zip`은 사용자 진입점 확인 4회의 출력이고, 각 평형 ZIP은 두 탐색의 명령·전 이벤트·원벽 및 선택 형상 R,Z·전 벽 열유속·열원별 자료·재평가 프로토콜과 전 형상 수치 요약이다. 수치 설정마다 수십 MB에 이르는 추가 고정 형상의 원시 배열은 D 드라이브 연구 폴더의 `replay/`에 두었다. 각 파일의 해시는 고정 형상 요약에 있다. 사례별 `analysis/`를 복사한 공개 표·그림과 원시 배열에서 다시 계산하는 검사기를 함께 제공한다.

```powershell
& $py -B results/d3d-user-shape-20261108/package_evidence.py
```

이 명령은 ZIP을 풀지 않고 멤버 목록·크기·SHA-256과 ZIP 자체의 SHA-256을 다시 계산한다. 네 진입점 폴더를 풀어 둔 곳에서 아래 명령으로 원벽·채택 형상과 마감 증거를 검사할 수 있다. 두 사례의 검색 폴더에는 바닥 누적전력 폭을 다시 계산하는 명령도 적용할 수 있다.

```powershell
& $py -B results/d3d-user-shape-20261108/verify_entry.py $study
& $py -B results/d3d-user-shape-20261108/measure_floor.py $study "$study/new-floor-footprint.csv"
& $py -B results/d3d-user-shape-20261108/compare_families.py $study "$study/family-comparison.png"
```
