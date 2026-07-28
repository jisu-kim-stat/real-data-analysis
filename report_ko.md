# CPS·Superconductivity·MEPS 1차 실제 데이터 분석

## 1. 한 줄 결론

세 데이터는 함께 쓰는 것이 좋다. 서로 같은 결과를 반복하는 것이 아니라 논문의 서로 다른 주장을 보여준다.

- **CPS**: 모든 방법이 marginal target을 만족하지만, CQR-TI가 가장 짧고 scale strata별 coverage도 가장 안정적이다.
- **Superconductivity**: SR/ASR score의 full-CDF pivotality가 CQR보다 훨씬 좋고, SR/ASR 구간도 CQR보다 짧다.
- **MEPS 의료비**: Parametric-TI가 명확하게 실패하고, ASR의 비대칭 보정과 CQR의 tail adaptation이 큰 효율 이득을 준다.

가장 중요한 공통 결과는 **SR-TI·ASR-TI·CQR-TI가 세 데이터의 20개 split 모두에서 held-out coverage 0.90 이상을 기록했다**는 것이다. 반면 Parametric-TI의 empirical split-success는 Superconductivity에서 0.80, MEPS에서 0.05였다.

## 2. 데이터와 실험 설정

| 데이터 | \(n\) | \(p\) | 반응·변환 | 주요 특징 |
|---|---:|---:|---|---|
| CPS 2012 | 29,217 | 100 | hourly wage, log scale 학습 후 원척도 복원 | 연속형, 오른쪽 왜도, 교육·경력·성별 이질성 |
| Superconductivity | 21,263 | 81 | critical temperature, 원척도 | 연속형, 비선형, 강한 scale heterogeneity |
| MEPS 2016 Panel 21 | 15,656 | 138 | total medical expenditure, `log1p` 학습 후 원척도 복원 | 18.2%가 0, 극단적 오른쪽 꼬리, 복합 survey |
| MEPS utilization sensitivity | 15,656 | 138 | 의료이용 횟수 합계, `log1p` | CQR 공개 코드가 실제 사용한 response, 27.3%가 0 |

공통 설정:

- \(C=0.90\), \(\alpha=0.05\)
- 각 반복에서 train/calibration/evaluation을 약 1/3씩 분할
- 20 repeated random splits
- mean, scale, lower/upper quantile 모두 gradient boosting 사용
- 300 trees, depth 3, shrinkage 0.04, bag fraction 0.75
- scale은 training split 안에서 별도 20% honest residual로 추정
- squared residual의 상위 0.5%를 winsorize하여 scale learner의 극단값 불안정 완화
- ASR tail scale은 \(\tau=0.10\)
- empirical score quantile은 R `quantile(..., type=1)` 사용

대략적인 calibration level은 다음과 같다.

| 데이터 | \(n_{\rm cal}\) | \(\lambda_\alpha\) | \(C+\lambda_\alpha\) |
|---|---:|---:|---:|
| CPS | 9,739 | 0.0138 | 0.9138 |
| Superconductivity | 7,087 | 0.0161 | 0.9161 |
| MEPS | 5,218 | 0.0188 | 0.9188 |

## 3. 주 결과

아래 coverage와 width는 20개 split의 평균이다. `Empirical success`는 20개 split 중 held-out coverage가 0.90 이상인 비율이다. 이는 유한 evaluation sample에 기초한 진단이며, 참 population content에 대한 theoretical outer probability의 직접 추정치는 아니다.

| 데이터 | 방법 | Coverage | Empirical success | Mean width | Q90 width | Lower miss | Upper miss |
|---|---|---:|---:|---:|---:|---:|---:|
| CPS | Parametric-TI | 0.936 | 1.00 | 39.77 | 59.66 | 0.034 | 0.030 |
| CPS | SR-TI | 0.915 | 1.00 | 37.18 | 63.17 | 0.044 | 0.042 |
| CPS | ASR-TI | 0.915 | 1.00 | 37.68 | 64.13 | 0.045 | 0.040 |
| CPS | CQR-TI | 0.914 | 1.00 | **36.32** | **60.06** | 0.044 | 0.042 |
| Superconductivity | Parametric-TI | 0.907 | 0.80 | **39.41** | **46.58** | 0.054 | 0.039 |
| Superconductivity | SR-TI | 0.917 | 1.00 | 42.61 | 78.61 | 0.044 | 0.039 |
| Superconductivity | ASR-TI | 0.917 | 1.00 | **42.57** | **78.58** | 0.044 | 0.039 |
| Superconductivity | CQR-TI | 0.917 | 1.00 | 45.98 | 79.98 | 0.037 | 0.046 |
| MEPS expenditure | Parametric-TI | 0.884 | 0.05 | 113,044 | 313,210 | 0.096 | 0.020 |
| MEPS expenditure | SR-TI | 0.920 | 1.00 | 80,804 | 179,317 | 0.065 | 0.015 |
| MEPS expenditure | ASR-TI | **0.921** | 1.00 | 56,142 | 125,895 | 0.052 | 0.027 |
| MEPS expenditure | CQR-TI | 0.920 | 1.00 | **13,664** | **29,986** | 0.023 | 0.057 |

굵은 width는 target coverage를 안정적으로 만족한 방법 안에서 해석해야 한다. 예를 들어 MEPS의 Parametric-TI는 넓으면서도 coverage가 부족하므로 효율적인 방법이 아니다.

## 4. 데이터별 해석

### CPS

- 세 PAC 방법의 평균 coverage는 0.914–0.915로 calibration level과 잘 맞는다.
- CQR-TI의 평균 폭 36.32가 가장 짧다.
- SR/ASR의 predicted-scale 1분위 coverage는 약 0.878이고 5분위는 약 0.947이다.
- CQR-TI는 scale 5분위에서 약 0.909–0.916으로 가장 안정적이다.
- 그런데 full-score KS는 SR/ASR가 0.047, CQR가 0.063이다. 즉 CQR의 전체 score CDF가 더 pivotal한 것은 아니지만 target coverage 주변의 relevant quantile alignment는 더 좋을 수 있다.

이 데이터는 논문의 “full pivotality는 편리한 충분조건이지 고정된 content level의 필요조건은 아니다”라는 설명에 잘 맞는다.

### Superconductivity

- Parametric-TI는 평균 coverage 0.907이지만 20개 split 중 4개가 0.90 아래였다.
- 모든 PAC 방법은 20개 split 모두 target을 넘었다.
- SR/ASR 평균 폭은 약 42.6으로 CQR의 46.0보다 약 7.4% 짧다.
- score pivotality 진단은 SR 0.066, ASR 0.066, CQR 0.235로 차이가 매우 크다.
- 이는 **standardized residual score가 CQR score보다 covariate-invariant에 훨씬 가까운 실제 데이터 사례**다.
- 반면 predicted-scale 1분위에서 SR/ASR coverage가 약 0.804로 낮다. 평균 KS가 작더라도 특정 score quantile 또는 특정 covariate region에서 mismatch가 남을 수 있음을 보여준다.

이 데이터는 score-CDF pivotality figure의 주 사례로 가장 좋다.

### MEPS total expenditure

- 반응의 18.2%가 0이고, median은 \$676.5, 99th percentile은 \$54,476.8, maximum은 \$483,812다.
- Parametric-TI coverage는 0.884이고 empirical success는 0.05다.
- Parametric-TI는 lower miss가 0.096으로, 0 또는 작은 지출을 지나치게 많이 제외한다.
- SR-TI는 validity를 회복하지만 평균 폭이 \$80,804로 매우 크다.
- ASR-TI는 coverage를 유지하면서 SR 대비 평균 폭을 약 **30.5% 감소**시킨다.
- CQR-TI는 평균 폭을 \$13,664까지 줄이지만 upper miss가 0.057로 lower miss 0.023보다 크다.
- predicted-scale strata coverage는 CQR가 약 0.911–0.936으로 가장 안정적이다.

MEPS는 parametric failure와 asymmetric adaptation의 실용적 의미를 가장 강하게 보여준다. 다만 zero atom 때문에 SR/ASR의 연속-error conditional theorem을 직접 예증하는 데이터로 쓰면 안 된다.

## 5. Score pivotality 결과

예측 scale 5분위별 score CDF와 전체 score CDF 사이의 empirical KS distance를 요약했다. 각 셀은 `평균 (90th percentile)`이다. 원고의 simulation pivotality figure와 직접 대응되는 값은 괄호 안 90th percentile이다.

| 데이터 | SR-TI | ASR-TI | CQR-TI |
|---|---:|---:|---:|
| CPS | **0.047 (0.090)** | **0.047 (0.090)** | 0.063 (0.115) |
| Superconductivity | **0.066 (0.121)** | **0.066 (0.118)** | 0.235 (0.378) |
| MEPS expenditure | **0.065 (0.101)** | **0.065 (0.096)** | 0.098 (0.170) |

이 결과는 논문의 핵심 구분을 지지한다.

1. 세 score 모두 동일한 PAC calibration을 받아 marginal coverage를 얻는다.
2. score construction에 따라 full-distribution pivotality는 크게 다르다.
3. CQR처럼 full-CDF pivotality가 약해도 target-level conditional coverage는 상대적으로 안정적일 수 있다.

따라서 coverage 표만 보여주는 것보다 pivotality figure와 scale-stratified coverage figure를 같이 넣는 편이 논문의 기여를 훨씬 잘 드러낸다.

## 6. MEPS utilization sensitivity

CQR 공개 코드가 실제로 만든 response는 총 의료비가 아니라 외래·응급·입원 등 의료이용 횟수 합계다.

| 방법 | Coverage | Empirical success | Mean width |
|---|---:|---:|---:|
| Parametric-TI | 0.922 | 1.00 | 21.30 |
| SR-TI | 0.920 | 1.00 | 31.20 |
| ASR-TI | 0.920 | 1.00 | 34.46 |
| CQR-TI | 0.933 | 1.00 | 25.98 |

의료이용 횟수에서는 Parametric-TI도 marginal target을 넘으므로, 본문에서 방법 차이를 보여주는 데는 총 의료비가 더 낫다. 이용횟수는 기존 CQR 문헌과의 연결을 보여주는 appendix sensitivity가 적절하다.

## 7. 원고 구성 권고

본문 Real Data Applications를 다음처럼 구성하는 것이 좋다.

1. **Data and protocol**
   - CPS, Superconductivity, MEPS 설명
   - 공통 split, learner, \(C\), \(\alpha\), transform 명시
2. **Marginal PAC and efficiency**
   - `tables/real_data_table.tex`
   - `figures/coverage_by_method.png`
   - `figures/width_by_method.png`
3. **Conditional diagnostic and pivotality**
   - `figures/coverage_by_predicted_scale.png`
   - `figures/score_pivotality.png`
4. **Appendix**
   - MEPS utilization sensitivity
   - 기존 Happy photo-z sanity check

TSA는 현재 i.i.d. 이론과 random split이 충돌하므로 삭제하거나, chronological/block split을 사용하는 별도 time-series appendix로 보내는 것이 안전하다.

## 8. 반드시 써야 할 한계

- 실제 데이터에서는 \(Y\mid X=x\)의 참 분포를 모르므로 conditional PAC를 직접 검증한 것이 아니다.
- `Empirical success`는 20 repeated splits와 유한 evaluation sample에 근거한 진단이다.
- MEPS 분석은 survey weights를 쓰지 않았으므로 미국 인구가 아니라 관측 sample distribution을 target으로 한다.
- MEPS의 zero atom은 SR/ASR conditional theorem의 연속성 조건과 다르다.
- Superconductivity에는 같은 material family의 dependence가 남아 있을 수 있다.
- 현재 GBM hyperparameter는 공통 고정값이다. 최종 투고 전에는 learner sensitivity 또는 nested tuning을 한 번 확인하는 것이 좋다.
- Parametric-TI는 PDF에 있는 정의를 재현 가능한 형태로 구현한 prototype benchmark다. 원 원고 코드가 있다면 정확히 맞춰 다시 계산해야 한다.

## 9. 현재 결과로 가능한 핵심 문장

> Across all three datasets, the PAC-calibrated procedures attained held-out content above the target in every one of 20 random splits, whereas the normal-theory benchmark failed in 20% of the superconductivity splits and 95% of the MEPS expenditure splits. The score diagnostics reveal a separate phenomenon: SR-TI and ASR-TI produced substantially more covariate-stable score distributions than CQR-TI, especially for superconductivity, even when CQR-TI exhibited stable coverage near the target level. This distinction empirically supports the paper's separation between marginal PAC calibration, full score pivotality, and target-level conditional behavior.
