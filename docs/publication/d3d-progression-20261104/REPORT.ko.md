# D3D 뷰어·OMP 재검수·외부 공개 배포 작업 보고서

2026-11-04 · 대상 브랜치 `codex/selective-candidate-validation`

## 결과

공개 주소: https://d3d-divertor-shape-explorer.sorainpublic.chatgpt.site

기존 로컬 뷰어의 277개 정적 파일(107,192,223 bytes)을 동일한 바이트로 Sites에 배포했다. 접근 정책은 `public`이며 링크를 가진 외부 사용자가 접속할 수 있다. 로컬 PC나 8780 서버의 지속 실행에 의존하지 않는다. 공개 사이트는 버전 1의 저장본이며 이후 로컬 수정은 별도 재배포해야 반영된다.

이번 GitHub 반영은 저장소 밖 D 드라이브에서 진행했던 뷰어·추가 정밀검증·v3 판정·재검수 작업을 함께 보존한다. 수송 solver와 기존 정책 모듈을 교체하거나 새 최적화를 실행한 변경이 아니다.

## 포함한 작업

- `viewers/d3d-progression-20261104/dist/`: 3D/R–Z 형상, 단계별 열분포, q/좌표 조회, leg, 기존·추가·v3 최종 결과와 검증 근거. vendor 파일도 포함해 외부 CDN 의존 없이 제공한다.
- 같은 폴더의 Python/JavaScript: 데이터 추출, 보관 후보 정밀검증, 분포 정책, 적대적·회귀 테스트와 기존 검수 기록. 과거 실행의 재현성을 위해 소스를 바꾸지 않고 반입했다.
- `results/omp-policy-audit-20261104/`: 감사 보고서, 20개 후보 재계산 지표, 원본 보호 해시, 직접 Gaussian 검산, 화면, 변경 전 뷰어 사본과 완료 기록.
- 본 폴더: Sites 배포 ID·버전·공개 주소, 반입 파일 해시 목록, 이번 브랜치의 검수 기록.

정적 뷰어 데이터는 재현 가능한 공개본을 보존하기 위해 포함했다. 원시 NPZ, 경로 캐시, 런타임, 인증 토큰, 배포 저장소의 `.git` 및 업로드 아카이브는 포함하지 않았다.

## 재검수 결론

목표는 고정 D3D·고정 입력·1 MW에서 peak를 낮추고 열을 퍼뜨리는 것이다. 최대점 위치 2 mm만으로 비슷한 두 봉우리의 순위 교체를 탈락시키는 것은 적절하지 않다. 이를 진단값으로 분리한 방향은 유지했다.

v3의 위치별 q 1%, L1/A2 .5%는 연구용 잠정 수치 기준이다. 물리적 허용조건에서 유도되거나 적대적 테스트만으로 정당화된 숫자가 아니다. 20 mm의 19개 후보는 해당 Linf 기준을 초과하지만 목적함수의 계산상 개선이 없는 것은 아니다. 이를 뷰어에서 명확히 구분했다.

| 후보 | cubic1280 peak MW/m² | peak 감소 | A2 증가 | v3 |
|---|---:|---:|---:|---|
| OMP20 4단계, 이전 최종 | 1.444160060 | 47.1150% | 52.4851% | 미통과 |
| OMP20 18단계 | 1.440295123 | 47.2565% | 53.2262% | 미통과 |
| OMP50 34단계 | 1.786940489 | 34.5624% | 11.5371% | 통과 |

20 mm 18단계의 국소 차이는 저장 donor로 Gaussian 합을 직접 계산해 2.693%를 재현했다. 따라서 기존 판정을 삭제하거나 후보가 통과하도록 허용오차를 완화하지 않았다. 18·19단계의 매우 작은 저장 peak 차이도 확정적 우열로 해석하지 않는다. [상세 수치·한계](../../../results/omp-policy-audit-20261104/REPORT.md)를 참고한다.

## 실행과 검수

저장소 루트에서 웹 뷰어만 실행:

```text
python -m http.server 8781 --bind 127.0.0.1 --directory viewers/d3d-progression-20261104/dist
```

브라우저에서 `http://127.0.0.1:8781/`을 연다. 정적 데이터만 사용하므로 Python 수치 라이브러리나 GPU는 필요 없다.

뷰어 폴더에서 회귀 검사:

```text
python -m unittest discover -p "test_*.py" -q
node --test test_math.mjs test_heat_probe.mjs test_stage_profiles.mjs test_qualified_finals.mjs test_distribution_checks.mjs
node verify_distribution_v3.mjs
```

Python 테스트와 재평가 드라이버는 보존된 D 드라이브 연구·원래 checkout 및 수치 의존성을 참조한다. 다른 PC에서 원시 수치 검증까지 재실행하려면 해당 보관 자료와 실행 환경을 복원해야 한다. 경로가 포함된 원본 드라이버를 범용 재실행 패키지라고 주장하지 않는다. JavaScript의 정적 뷰어 테스트는 이 저장소에 포함된 데이터로 실행한다.

원래 감사에서 Python 35개·JavaScript 22개가 통과했고 보호 파일 1,397개가 일치했다. 이번 반입 후 다시 실행한 결과는 `branch-verification.json`과 테스트 로그에 기록한다. `.gitattributes`의 `* -text`를 유지해 JSON 및 소스의 해시를 바꾸는 줄바꿈 변환을 방지한다.

공개 사이트는 native 배포 도구의 성공 및 public 접근 설정을 확인했다. 실제 브라우저에서도 OMP50 peak 1.7869와 1280 경로 분포 로딩을 확인했다. 별도 Python urllib의 비로그인 요청은 Cloudflare 1010으로 차단돼 그 방법의 독립 접속 검사는 완료하지 못했다. 이를 브라우저 장애나 비공개 설정으로 해석하지 않았다.

## 보존 경로와 알려진 한계

- 원본 뷰어: `D:/DivertorOptimization/viewers/d3d-progression-20261104`
- 최초 정밀검증: `D:/DivertorOptimization/studies/omp-stage-qualification-20261104`
- 추가 3개 OMP20 검증: `D:/DivertorOptimization/studies/omp-distribution-extension-20261104`
- 보존된 v3 최종 판정: `D:/DivertorOptimization/studies/omp-distribution-final-20261104`
- 이번 감사 원본: `D:/DivertorOptimization/studies/omp-policy-audit-20261104`
- Sites 소스 checkout: `D:/DivertorOptimization/exports/d3d-progression-public-20261104`

감사 사본의 절대 경로와 체크섬은 당시 원본을 가리키며, 보관된 historical manifest를 새 경로로 재작성하지 않았다. 100 mm·IMP의 새 수치 재평가, 연속 함수 오차 인증, 전역 최적성 증명은 이번 범위에 없다. v3 로더의 비1MW 스케일 인자 문제도 현재 고정1MW 결과와 구분해 상세 감사에 기록했다. 이번 push로 그 frozen 모듈을 수정하지 않는다.

이번 반입 검수: Python 35/35, JavaScript 22/22 통과 및 표시 셀 867개 조회·최종 후보 연결 검사를 통과했다. 반입 345개 파일과 감사 아티팩트의 바이트 해시가 일치했다. CRLF를 인정한 whitespace 검사에서 과거 계획 문서의 말미 빈 줄 2곳과 외부 vendor 들여쓰기 1곳만 남았으며, 원본 해시를 보존하기 위해 수정하지 않았다.
