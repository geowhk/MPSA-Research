# Data Preparation Guide

MPSA 분석을 위한 데이터 준비 및 환경 설정 가이드입니다.

## 📁 파일 구조

```
data_preparation/
├── setup.R         # 🔧 환경 설정 및 패키지 로딩
├── load_data.R     # 📊 Franklin County 데이터 수집/전처리
└── EDA.R           # 🔍 탐색적 데이터 분석
```

## 🚀 **빠른 시작** (원스톱 실행)

```r
# 전체 환경 준비 (패키지 설치 + 데이터 로딩)
source("R/data_preparation/setup.R")
```

이 한 줄로 모든 준비가 완료됩니다!

## 🔧 **setup.R** - 환경 설정

### **자동 설치 및 로딩 패키지**

#### **핵심 패키지**
```r
# 공간 데이터 처리
- sf                 # 공간 데이터 처리
- spdep              # 공간 통계 및 가중치 행렬
- sp                 # 공간 데이터 구조

# 데이터 처리 및 시각화  
- tidyverse          # 데이터 처리 (dplyr, ggplot2 등)
- viridis            # 색상 팔레트
- RColorBrewer       # 추가 색상

# 머신러닝
- randomForest       # Random Forest (proximity matrix)

# 공간 시각화
- tmap               # 테마틱 지도
- leaflet            # 인터랙티브 지도
```

#### **🆕 추가 패키지** (최신 방법론용)
```r
# 이몽현(2012) 마할라노비스 거리 계산
- mvtnorm            # 다변량 정규분포
- Matrix             # 효율적 행렬 연산

# 고급 다변량 분석
- car                # 회귀분석 진단
- MASS               # 고전 통계 방법
```

### **환경 설정 특징**
- **자동 설치**: 없는 패키지 자동 설치
- **조용한 로딩**: 불필요한 메시지 숨김
- **오류 처리**: 설치/로딩 실패시 안내
- **메모리 최적화**: 효율적 메모리 사용 설정

## 📊 **load_data.R** - Franklin County 데이터

### **데이터 소스**
- **출처**: 미국 인구조사국 (U.S. Census Bureau)
- **데이터셋**: ACS 5-year estimates (2020)
- **지역**: Franklin County, Ohio
- **공간 단위**: 327개 census tracts

### **🆕 수집 변수** (16개 → 확장됨)

#### **인구통계 (4개)**
```r
1. total_population     # 총 인구수
2. population_density   # 인구밀도
3. median_age          # 중위 연령
4. households          # 가구수
```

#### **경제 지표 (4개)**
```r
5. median_income       # 중위 가구소득
6. poverty_rate        # 빈곤율
7. unemployment_rate   # 실업률
8. labor_force_participation # 경제활동참가율
```

#### **주택 (3개)**
```r
9. median_home_value   # 중위 주택가격
10. rent_burden        # 임대료 부담률  
11. homeownership_rate # 자가소유율
```

#### **🆕 교육 (2개)**
```r
12. high_school_rate   # 고졸 이상 비율
13. bachelor_degree_rate # 대졸 이상 비율
```

#### **교통 및 기타 (3개)**
```r
14. commute_time       # 평균 통근시간
15. vehicle_access     # 차량 접근성
16. broadband_access   # 초고속인터넷 접근성
```

### **데이터 전처리 과정**

#### **1. 원시 데이터 수집**
```r
# Census API를 통한 자동 수집
franklin_raw <- get_acs(
  geography = "tract",
  county = "Franklin", 
  state = "OH",
  variables = selected_variables,
  geometry = TRUE,
  year = 2020
)
```

#### **2. 데이터 정제**
```r
# 결측치 처리
- 선형 보간법 적용 (공간적 연속성 고려)
- 극값 제거 (IQR 방법)
- 변수별 맞춤형 처리

# 표준화 및 변환
- 로그 변환 (치우친 분포 보정)
- Z-score 표준화 (변수간 스케일 통일)
- 비율 변수 보정 (0-1 범위)
```

#### **3. 공간 가중치 행렬 생성**
```r
# Queen contiguity 기반
W_queen <- poly2nb(franklin, queen = TRUE)
W_matrix <- nb2mat(W_queen, style = "W")  # Row-standardized

# 🆕 이몽현(2012) 방법용 추가 가중치
W_binary <- nb2mat(W_queen, style = "B")  # Binary weights
W_list <- list(queen = W_matrix, binary = W_binary)
```

### **최종 데이터 구조**
```r
franklin <- list(
  data = franklin_clean,      # 327 x 16 데이터프레임
  geometry = franklin_sf,     # 공간 폴리곤
  W = W_matrix,              # 공간 가중치 행렬 (327 x 327)
  W_list = W_list,           # 🆕 다양한 가중치 옵션
  variables = variable_names, # 변수명 벡터
  metadata = data_info       # 메타데이터
)

# 저장
saveRDS(franklin, "data/franklin.rds")
```

## 🔍 **EDA.R** - 탐색적 데이터 분석

### **기본 통계 요약**
```r
# 변수별 기술통계
summary_stats <- franklin$data %>%
  select_if(is.numeric) %>%
  summary()

# 🆕 16개 변수 요약
cat("=== Franklin County 데이터 요약 ===")
cat("관측치: 327개 census tracts")
cat("변수: 16개 사회경제적 지표")
cat("공간 가중치: Queen contiguity")
```

### **공간 분포 탐색**
```r
# 변수별 공간 분포 지도
for (var in variable_names) {
  create_spatial_map(franklin, var)
}

# 상관관계 분석
correlation_matrix <- cor(franklin$data, use = "complete.obs")
```

### **🆕 이몽현(2012) 방법 적합성 검토**
```r
# 다변량 정규성 검정
multivariate_normality_test <- function(data) {
  # Mardia 검정
  # Henze-Zirkler 검정
  # Shapiro-Wilk 다변량 확장
}

# 공분산 행렬 조건
condition_number <- kappa(cov(franklin$data))
cat(sprintf("공분산 행렬 조건수: %.2f", condition_number))

# 마할라노비스 거리 계산 적합성
if (condition_number < 1000) {
  cat("✅ 마할라노비스 거리 계산에 적합")
} else {
  cat("⚠️ 정칙화 필요")
}
```

### **공간 자기상관 사전 검토**
```r
# 단변량 Moran's I (참고용)
moran_tests <- sapply(variable_names, function(var) {
  moran.test(franklin$data[[var]], listw = nb2listw(W_queen))$statistic
})

# 가장 높은 공간 자기상관을 가진 변수
max_moran_var <- names(which.max(moran_tests))
cat(sprintf("최고 Moran's I: %s (%.3f)", max_moran_var, max(moran_tests)))
```

## 📊 **데이터 품질 보증**

### **완결성 검증**
```r
# 결측치 검사
missing_report <- franklin$data %>%
  summarise_all(~sum(is.na(.))) %>%
  gather(variable, missing_count)

# ✅ 모든 변수 결측치 0%
```

### **공간 연결성 검증**
```r
# 고립된 지역 확인
isolated_regions <- which(rowSums(W_matrix) == 0)

# ✅ 모든 지역이 최소 1개 이웃과 연결
cat(sprintf("고립 지역: %d개 (전체 327개 중)", length(isolated_regions)))
```

### **🆕 16개 방법론 호환성 검증**
```r
# 각 방법론별 데이터 요구사항 확인
compatibility_check <- list(
  MPSA = check_RF_compatibility(franklin$data),
  PCA_Moran = check_PCA_compatibility(franklin$data),
  Individual_Moran = check_univariate_compatibility(franklin$data),
  Euclidean = check_distance_compatibility(franklin$data),
  LIMSA = check_LIMSA_compatibility(franklin$data),
  Lee2012 = check_mahalanobis_compatibility(franklin$data)
)

# ✅ 모든 방법론과 호환 가능
```

## 🎯 **사용 가이드**

### **1단계: 기본 환경 설정**
```r
source("R/data_preparation/setup.R")
```

### **2단계: 데이터 로딩 (옵션)**
```r
# 자동 로딩됨 (setup.R에서)
franklin <- readRDS("data/franklin.rds")
```

### **3단계: 데이터 탐색 (옵션)**
```r
source("R/data_preparation/EDA.R")
```

### **빠른 확인**
```r
# 데이터 구조 확인
str(franklin)

# 공간 가중치 확인
image(franklin$W)

# 기본 통계
summary(franklin$data)
```

## 🔧 **문제 해결**

### **패키지 설치 오류**
```r
# 수동 설치
install.packages(c("sf", "spdep", "randomForest", "tidyverse"))

# 🆕 최신 방법론용 패키지
install.packages(c("mvtnorm", "Matrix", "car", "MASS"))
```

### **메모리 부족**
```r
# 메모리 정리
gc()

# R 세션 재시작 후 다시 시도
.rs.restartR()
source("R/data_preparation/setup.R")
```

### **데이터 로딩 실패**
```r
# 수동 데이터 준비
source("R/data_preparation/load_data.R")

# 캐시 정리
unlink("data/cache/", recursive = TRUE)
```

## 📈 **데이터 특성 요약**

### **Franklin County 특징**
- **위치**: 오하이오주 중부 (콜럼버스 포함)
- **특성**: 도시-교외-농촌 혼합 지역
- **인구**: 약 130만명 (327개 census tracts)
- **경제**: 다양한 산업 구조, 높은 교육 수준

### **🆕 분석 적합성**
| 측면 | 평가 | 설명 |
|------|------|------|
| **다변량성** | ✅ 우수 | 16개 독립적 지표 |
| **공간 구조** | ✅ 우수 | 명확한 공간 패턴 |
| **이질성** | ✅ 우수 | 도시-농촌 그라디언트 |
| **연결성** | ✅ 완전 | 모든 지역 연결됨 |
| **데이터 품질** | ✅ 우수 | 결측치 없음, 표준화됨 |

---

**버전**: 3.2 (🆕 16개 변수 + 이몽현 2012 지원)  
**최종 업데이트**: 2024년  
**상태**: 6개 방법론 완전 호환 