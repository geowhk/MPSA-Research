# Advanced MPSA Features

방법론 제안 논문의 범위를 벗어나는 고급 MPSA 확장 기능들입니다.

## 📁 파일 구조

```
advanced/
├── MPSA_advanced.R           # 🚀 고급 확장 (베이지안, 다중스케일, 시공간)
└── MPSA_alternatives.R       # 🔬 대안 구현 (적응형, 행렬분해)
```

## 🎯 **심화 연구 포함 이유**

### ✅ **방법론 제안 논문의 명확성 유지**
- 핵심 기여(MPSA 방법론 제안)에 집중
- 논문 분량 및 복잡도 관리
- 독자의 이해도 향상

### ✅ **후속 연구 기반 제공**
- 각 확장은 별도 논문 주제로 활용 가능
- 연구 연속성 확보
- 방법론 생태계 구축

### ✅ **완전성 및 투명성**
- 가능한 모든 확장 방향 제시
- 연구 과정의 완전한 문서화
- 재현 가능한 연구 환경

## 🚀 **MPSA_advanced.R** - 고급 확장 기능

### 📊 **베이지안 MPSA**
```r
compute_bayesian_MPSA(data, W, prior_params = list())
```

#### 핵심 아이디어
- **목적**: MPSA의 불확실성 정량화
- **방법**: 베이지안 추론을 통한 사후분포 추정
- **장점**: 신뢰구간, 모델 비교, 의사결정 지원

#### 구현 특징
- MCMC 기반 사후분포 샘플링
- 다양한 사전분포 선택 옵션
- 베이지안 모델 선택 기준

#### �� 이몽현(2012) 방법과의 비교

```r
# 베이지안 vs 이몽현(2012) 접근 - 정확한 방법론 구현
source("R/mpsa_methods/MPSA.R")
lee_2012_results <- compute_Lee2012_Mahalanobis(data, W)

# 방법론 비교:
# - 베이지안 MPSA: 사후 분포 기반 불확실성 정량화
# - 이몽현(2012): 각 지역과 이웃들의 평균 간 마할라노비스 거리 + 카이제곱 검정
# - 공통점: 다변량 공간 자기상관 측정
# - 차이점: 베이지안은 불확실성 정량화, 이몽현은 이론적 거리 기반 접근

cat("베이지안: 불확실성 정량화 및 계층적 모델링")
cat("이몽현(2012): 공분산 구조를 명시적으로 고려한 이론적 접근")
cat("논문: '마할라노비스 거리를 이용한 다변량 공간 클러스터 분석' (Lee, 2012)")
```

### 🌐 **다중스케일 MPSA**
```r
compute_multiscale_MPSA(data, W_list, scales = c(1, 2, 3))
```

#### 핵심 아이디어
- **목적**: 다양한 공간 스케일에서의 패턴 탐지
- **방법**: 계층적 공간 가중치 행렬 활용
- **장점**: 스케일별 패턴 분해, 다층 분석

#### 구현 특징
- 다중 공간 가중치 행렬 지원
- 스케일별 기여도 분해
- 계층적 클러스터링과의 연계

### ⏰ **시공간 MPSA**
```r
compute_spatiotemporal_MPSA(data_panel, W, temporal_weights)
```

#### 핵심 아이디어
- **목적**: 시공간 데이터에서의 동적 패턴 탐지
- **방법**: 공간-시간 가중치 행렬 확장
- **장점**: 동적 핫스팟 추적, 확산 패턴 분석

#### 구현 특징
- 시공간 가중치 행렬 자동 생성
- 동적 패턴 시각화
- 확산 속도 및 방향 분석

## 🔬 **MPSA_alternatives.R** - 대안 구현 방법

### 🎯 **적응형 MPSA**
```r
compute_adaptive_MPSA(data, W, adaptation_method = "local")
```

#### 핵심 아이디어
- **목적**: 지역별 특성에 맞는 적응적 분석
- **방법**: 국소적 패라미터 조정
- **장점**: 이질성 반영, 세밀한 패턴 탐지

#### 🆕 이몽현(2012) 방법과의 차별점

```r
# 지역별 적응적 모델링 vs 전역 공분산 구조
bayesian_local_effects <- extract_local_effects(bayesian_results)
lee_2012_results <- compute_Lee2012_Mahalanobis(data, W)

# 핵심 차이점:
# - 베이지안: 지역별 적응적 공분산 구조 학습
# - 이몽현(2012): 각 지역과 이웃들의 평균 간 전역 공분산 행렬 사용
# - 베이지안: 계층적 모델로 공간 의존성 명시적 모델링  
# - 이몽현(2012): 4단계 계산 과정 (대상지역 → 이웃선택 → 이웃평균 → 거리계산)

comparison_table <- data.frame(
  Aspect = c("공분산 추정", "공간 의존성", "불확실성", "계산 복잡도"),
  Bayesian = c("지역별 적응적", "계층적 모델링", "사후 분포", "높음 (MCMC)"),
  Lee2012 = c("전역 고정", "이웃 평균 기반", "카이제곱 검정", "중간 (O(n²p²))")
)
```

### 📐 **행렬분해 기반 MPSA**
```r
compute_matrix_decomposition_MPSA(data, W, method = "SVD")
```

#### 핵심 아이디어
- **목적**: 고차원 데이터의 효율적 처리
- **방법**: SVD, NMF 등 행렬분해 기법 활용
- **장점**: 차원축소, 계산 효율성, 잠재요인 해석

#### 구현 특징
- 다양한 분해 방법 지원 (SVD, NMF, ICA)
- 차원축소와 패턴 탐지의 동시 수행
- 잠재요인 해석 기능

### 🧠 **커널 기반 MPSA**
```r
compute_kernel_MPSA(data, W, kernel_type = "RBF")
```

#### 핵심 아이디어
- **목적**: 비선형 패턴의 효과적 탐지
- **방법**: 커널 트릭을 이용한 특징공간 변환
- **장점**: 비선형성 처리, 복잡한 패턴 탐지

## 🆕 **종합 방법론 비교 확장**

### **8개 방법론 확장 비교** (기본 6개 + 고급 2개)
```r
# 기본 6개 방법론
basic_methods <- c("MPSA", "PCA_Moran", "Individual_Moran", 
                   "Euclidean", "LIMSA", "Lee2012_Mahalanobis")

# 🆕 고급 2개 방법론 추가
advanced_methods <- c("Bayesian_MPSA", "Adaptive_MPSA")

# 종합 벤치마크 (8개 방법론)
comprehensive_benchmark(data, W, include_advanced = TRUE)
```

### **특성별 비교 매트릭스**
| 방법 | 이론적 기반 | 불확실성 정량화 | 적응성 | 계산 복잡도 |
|------|-------------|----------------|--------|-------------|
| **MPSA** | RF 이론 | 순열검정 | 중간 | O(n log n) |
| **이몽현(2012)** | 마할라노비스 | 카이제곱 분포 | 낮음 | O(n²p²) |
| **🆕 베이지안 MPSA** | RF + 베이지안 | 사후분포 | 높음 | O(n³) |
| **🆕 적응형 MPSA** | RF + 국소학습 | 지역별 검정 | 매우 높음 | O(n²p) |

## 🎯 **활용 가이드**

### **기본 방법론 제안 논문** (권장)
```r
# 6개 방법론 비교만 사용
source("R/analysis/test_alternatives.R")
benchmark_results <- benchmark_against_traditional_methods(data, W)
```

### **심화 연구 또는 후속 논문**
```r
# 고급 기능 포함 종합 분석
source("R/advanced/MPSA_advanced.R")
source("R/advanced/MPSA_alternatives.R")

# 8개 방법론 종합 비교
comprehensive_results <- comprehensive_benchmark(data, W, include_advanced = TRUE)
```

### **특정 확장 기능 집중 연구**
```r
# 베이지안 MPSA 집중 연구
bayesian_results <- compute_bayesian_MPSA(data, W)
bayesian_analysis <- analyze_bayesian_results(bayesian_results)

# 시공간 MPSA 집중 연구  
panel_data <- load_panel_data()
spatiotemporal_results <- compute_spatiotemporal_MPSA(panel_data, W, temporal_weights)
```

## 🔮 **후속 연구 아이디어**

### **1. 베이지안 MPSA 논문**
- **제목**: "Bayesian Inference for Multivariate Spatial Autocorrelation"
- **기여**: 불확실성 정량화, 모델 선택, 의사결정 지원
- **타겟**: Bayesian Analysis, Spatial Statistics

### **2. 시공간 MPSA 논문**
- **제목**: "Spatiotemporal Extension of MPSA for Dynamic Pattern Detection"
- **기여**: 동적 패턴 탐지, 확산 분석, 예측 모델링
- **타겟**: Journal of Regional Science, Computers & Geosciences

### **3. 🆕 마할라노비스 vs MPSA 심화 비교 논문**
- **제목**: "Comparative Study of Distance-based vs Proximity-based Multivariate Spatial Autocorrelation"
- **기여**: 거리 기반 vs 유사성 기반 방법론 심화 비교
- **타겟**: Geographical Analysis, Spatial Statistics

### **4. 고차원 MPSA 논문**
- **제목**: "High-dimensional Multivariate Spatial Analysis with Matrix Decomposition"
- **기여**: 빅데이터 공간분석, 차원축소와 패턴탐지 융합
- **타겟**: Computational Statistics & Data Analysis

## 📦 **추가 의존성**

### **베이지안 분석용**
```r
install.packages(c("rstan", "brms", "rstanarm"))
```

### **행렬분해용**
```r
install.packages(c("RSpectra", "fastICA", "NMF"))
```

### **시공간 분석용**
```r
install.packages(c("spacetime", "gstat", "stars"))
```

## ⚠️ **주의사항**

### **논문 범위 관리**
- 기본 논문에서는 고급 기능 언급만 (Future Work 섹션)
- 각 확장은 별도 연구 주제로 개발
- 복잡도 증가에 따른 독자 부담 고려

### **계산 자원 요구사항**
- **베이지안 분석**: 32GB RAM 권장 (MCMC)
- **시공간 분석**: 대용량 패널 데이터 처리
- **행렬분해**: GPU 가속 권장 (대규모 행렬)

## 🎯 **기대 효과**

### **학술적 기여**
1. **방법론 생태계 구축**: 기본 → 확장 → 응용의 연구 파이프라인
2. **다학제 융합**: 통계학, 지리학, 컴퓨터과학, 경제학
3. **후속 연구 촉진**: 각 확장별 독립적 연구 주제 제공

### **실용적 가치**
1. **맞춤형 분석**: 데이터 특성에 맞는 방법 선택
2. **확장성**: 다양한 응용 분야 대응
3. **미래 지향성**: 새로운 데이터 유형 및 분석 요구 대응

---

**버전**: 3.2 (🆕 이몽현 2012 + 고급 확장 통합)  
**최종 업데이트**: 2024년  
**상태**: 심화 연구 기반 완성 (8개 방법론 지원) 