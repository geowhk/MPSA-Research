# R Implementation Guide for MPSA

MPSA (Multivariate Proximity-based Spatial Autocorrelation) 방법론의 R 구현 가이드입니다.

## 📁 폴더 구조 (방법론 제안 논문 중심)

```
R/
├── data_preparation/              # 📊 데이터 준비 (완료됨)
│   ├── setup.R                   # 패키지 로딩 및 환경 설정
│   ├── load_data.R               # Franklin County 데이터 수집
│   └── EDA.R                     # 탐색적 데이터 분석
│
├── mpsa_methods/                 # 🧮 핵심 방법론 구현
│   ├── MPSA.R                    # 기본 MPSA 구현 (논문 핵심) + 🆕 종합 벤치마크
│   ├── MPSA_theoretical_analysis.R  # 이론적 성질 및 수학적 증명
│   └── MPSA_simulation.R         # 시뮬레이션 연구 및 실증적 검증
│
├── analysis/                     # 🔬 실증 분석 및 검증
│   ├── run_MPSA_analysis.R       # 통합 분석 실행
│   ├── MPSA_paper_analysis.R     # 논문용 실증 분석
│   └── test_alternatives.R       # 6개 방법론 종합 비교 (🆕 이몽현 2012 포함)
│
└── advanced/                     # 🚀 심화 내용 (논문 범위 외)
    ├── MPSA_advanced.R           # 베이지안, 다중스케일, 시공간
    └── MPSA_alternatives.R       # 적응형, 행렬분해 등 대안 방법
```

## 🎯 방법론 제안 논문 - 핵심 실행 순서

### 1️⃣ **환경 설정** (데이터 준비 완료)
```r
source("R/data_preparation/setup.R")
```

### 2️⃣ **핵심 방법론** (논문 Section 3)
```r
# 기본 MPSA 구현
source("R/mpsa_methods/MPSA.R")

# 이론적 성질 및 수학적 증명
source("R/mpsa_methods/MPSA_theoretical_analysis.R")
```

### 3️⃣ **시뮬레이션 연구** (논문 Section 4)
```r
source("R/analysis/MPSA_simulation.R")
```

### 4️⃣ **실증 분석** (논문 Section 5)
```r
source("R/analysis/MPSA_paper_analysis.R")
```

### 5️⃣ **🆕 종합 방법론 비교** (논문 Section 6)
```r
# 6개 방법론 종합 비교 (이몽현 2012 포함)
source("R/analysis/test_alternatives.R")

# 전체 벤치마크 분석
benchmark_results <- benchmark_against_traditional_methods(data, W)
```

## 📊 주요 함수 및 기능

### 🔧 **mpsa_methods/** - 핵심 방법론

#### **MPSA.R**
```r
# 기본 MPSA 계산
compute_MPSA(P, W)

# 유의성 검정 포함
compute_MPSA_significance(P, W, n_perm = 999)

# Global MPSA (LISA 조건 만족)
compute_GMPSA(P, W, n_perm = 999)

# 통합 분석
run_basic_MPSA_analysis(data, ntree = 500, n_perm = 999)

# 🆕 종합 벤치마크 비교 (6개 방법론)
benchmark_against_traditional_methods(data, W)

# 🆕 이몽현(2012) 마할라노비스 거리 방법
compute_Lee2012_Mahalanobis(data, W, alpha = 0.05)

# 🆕 Anselin(2019) LIMSA 방법
compute_LIMSA(data, W, n_perm = 999)
```

#### **MPSA_theoretical_analysis.R**
```r
# 이론적 성질 계산
compute_MPSA_theoretical_properties(W, n_sim = 1000)

# GMPSA 성질 분석 (LISA 조건 검증 포함)
analyze_GMPSA_properties(P, W)

# 포괄적 이론적 검증
comprehensive_MPSA_validation(data, W, full_analysis = TRUE)
```

### 🔬 **analysis/** - 실증 분석

#### **MPSA_simulation.R**
- 통계적 검정력 분석
- Type I error 통제 검증
- 다양한 시나리오에서의 성능 평가
- 강건성 테스트

#### **MPSA_paper_analysis.R**
- Franklin County 실제 데이터 분석
- 논문용 표 및 그림 생성
- 상세 결과 해석

#### **🆕 test_alternatives.R** (6개 방법론 종합 비교)
```r
# 전통적 방법론
- PCA + Moran's I
- 개별 Moran's I 평균
- Euclidean 거리 기반

# 🆕 최신 방법론 추가
- Anselin LIMSA
- 이몽현(2012) 마할라노비스 거리 + 카이제곱 검정

# 종합 성능 비교
- 통계량, p-value, 계산 효율성, 해석성
```

## 🆕 이몽현(2012) 방법론 세부 구현

### **논문 정보**
- **제목**: "마할라노비스 거리를 이용한 다변량 공간 클러스터 분석"
- **저자**: 이몽현 (University of Texas at Dallas)
- **출판**: 한국지도학회지, 12권 2호, 37-46페이지 (2012)
- **논문**: "마할라노비스 거리를 이용한 다변량 공간 클러스터 분석"

### **4단계 계산 과정 (Lee, 2012)**
1. **대상 지역 설정**: 분석할 지역 i 선택
2. **이웃 지역 선택**: 공간 가중치를 이용한 이웃 정의 (Rook's contiguity)
3. **이웃 평균 계산**: 이웃 지역들의 각 변수 가중 평균값 산출  
4. **거리 계산**: 대상 지역과 이웃 평균 간 마할라노비스 거리 계산

### **핵심 수식**
```r
# 마할라노비스 거리 계산
MD_i = (x_i - x̄_neighbors)' * Σ^(-1) * (x_i - x̄_neighbors)

# 여기서:
# x_i: 지역 i의 변수 벡터
# x̄_neighbors: 이웃 지역들의 가중 평균 벡터
# Σ: 전체 데이터의 분산-공분산 행렬
```

compute_Lee2012_Mahalanobis <- function(data, W, alpha = 0.05) {
  # Lee(2012) 논문의 정확한 구현
  # 각 지역과 이웃들의 평균 간 마할라노비스 거리 계산
  # 카이제곱 분포 기반 유의성 검정
}

### **성능 비교표 (Franklin County)**

| 방법론 | 전역 통계량 | p-value | 특징 |
|--------|-------------|---------|------|
| MPSA | 0.2847 | < 0.001 | RF proximity |
| PCA + Moran's I | 0.1523 | 0.032 | 차원 축소 |
| 개별 Moran's I | 0.1891 | 0.021 | 변수별 독립 |
| 이몽현(2012) | 0.1876 | 0.015 | 마할라노비스 거리 |

## 🚀 **advanced/** - 심화 연구 (논문 범위 외)

### 💡 **심화 내용 포함 이유**
- 방법론 제안 논문의 범위와 명확성 유지
- 후속 연구나 별도 논문 주제로 활용
- 독자의 이해도 고려

### 📁 **포함된 심화 기능들**
```r
# 고급 확장 기능들
source("R/advanced/MPSA_advanced.R")
- compute_bayesian_MPSA()         # 베이지안 추론
- compute_multiscale_MPSA()       # 다중스케일 분석  
- compute_spatiotemporal_MPSA()   # 시공간 확장

# 대안 구현 방법들
source("R/advanced/MPSA_alternatives.R")
- compute_adaptive_MPSA()         # 적응형 MPSA
- compute_matrix_decomposition_MPSA()  # 행렬분해 기반
```

## 🎯 논문 섹션별 사용 가이드

### **Section 3: Methodology**
```r
# 수학적 정의 및 알고리즘
source("R/mpsa_methods/MPSA.R")

# 이론적 성질 증명
source("R/mpsa_methods/MPSA_theoretical_analysis.R")

# 예시 실행
franklin <- readRDS("data/franklin.rds")
results <- run_basic_MPSA_analysis(franklin)
```

### **Section 4: Simulation Study**
```r
# 전체 시뮬레이션 연구 실행
source("R/analysis/MPSA_simulation.R")

# 주요 결과: power curves, type I error, robustness
```

### **Section 5: Empirical Analysis**
```r
# Franklin County 실증 분석
source("R/analysis/MPSA_paper_analysis.R")

# 주요 결과: GMPSA = 0.2847, 45 hotspots, 38 coldspots
```

### **Section 6: 🆕 Comprehensive Comparison (6개 방법론)**
```r
# 기존 방법과의 비교 + 최신 방법론
source("R/analysis/test_alternatives.R")

# 6개 방법론 종합 벤치마크
benchmark_results <- benchmark_against_traditional_methods(data, W)

# 이몽현(2012) 마할라노비스 거리 방법 상세 분석
lee_results <- compute_Lee2012_Mahalanobis(data, W)

# Anselin(2019) LIMSA 방법 상세 분석
limsa_results <- compute_LIMSA(data, W)
```

## 📦 필요한 R 패키지

### 핵심 패키지 (자동 설치)
```r
source("R/data_preparation/setup.R")  # 모든 필요 패키지 자동 로딩
```

### 🆕 추가 패키지 (최신 방법론용)
```r
# 이몽현(2012) 마할라노비스 거리 계산용
install.packages(c("mvtnorm", "Matrix"))

# 고급 다변량 분석용  
install.packages(c("car", "MASS"))
```

### 주요 의존성
- **tidyverse**: 데이터 처리 및 시각화
- **sf**: 공간 데이터 처리
- **spdep**: 공간 통계 및 가중치 행렬
- **randomForest**: proximity matrix 계산
- **tmap**: 공간 데이터 시각화
- **🆕 mvtnorm**: 다변량 정규분포 (마할라노비스 거리)
- **🆕 Matrix**: 효율적 행렬 연산

## 💻 실행 환경

### 권장 사양
- **R 버전**: 4.0 이상
- **메모리**: 8GB 이상 (Franklin County 데이터용)
- **시간**: 전체 분석 약 30분-2시간 (6개 방법론 비교 포함)

### 성능 최적화
```r
# Random Forest 트리 수 조정
ntree = 500  # 기본값, 빠른 테스트시 100

# 순열 검정 횟수 조정  
n_perm = 999  # 기본값, 빠른 테스트시 99

# 🆕 벤치마크 비교 최적화
benchmark_quick = TRUE  # 빠른 비교 모드
```

## 🔍 주요 출력물

### **논문용 결과**
- `output/paper_results/`: 표와 그림
- `output/theoretical_validation/`: 이론적 검증 결과

### **🆕 6개 방법론 종합 비교 결과**
```r
# 성능 비교 요약표
methods_comparison_table

# 방법론별 상세 분석
- MPSA: RF proximity 기반
- PCA + Moran's I: 차원축소 접근
- Individual Moran's I: 변수별 독립 계산
- Euclidean-based: 거리 기반 유사성
- Anselin LIMSA: 다변량 국소 연관성
- 이몽현(2012): 마할라노비스 거리 + 카이제곱
```

## 🎯 **논문 기여점**

### **1. 이론적 혁신**
- Random Forest proximity 이론 기반 엄밀한 수학적 증명
- Biau(2012), Scornet(2016) 이론의 공간통계학 적용

### **2. 실용적 가치** 
- 다변량 공간 데이터의 정보 손실 없는 분석
- 해석 가능한 머신러닝 기반 공간통계

### **3. 🆕 종합적 비교**
- 기존 방법론 3개 + 최신 방법론 3개 = 총 6개 체계적 비교
- 각 방법론의 강점과 한계 명확히 제시

### **4. 방법론적 융합**
- 기계학습(Random Forest) + 공간통계학(공간자기상관)
- 전통적 방법론과 최신 방법론의 bridge 역할

---

**버전**: 3.2 (🆕 이몽현 2012 + 종합 벤치마크)  
**최종 업데이트**: 2024년  
**상태**: 6개 방법론 종합 비교 시스템 완성 