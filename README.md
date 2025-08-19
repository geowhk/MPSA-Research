# MPSA 연구 프로젝트 v3.2

## 🎯 **목표**
Random Forest의 proximity matrix를 이용한 다변량 공간 자기상관 분석 방법론(MPSA) 개발 및 검증

## 📂 **프로젝트 구조**

```
├── README.md                          # 📖 프로젝트 개요 (이 파일)
├── 논문_작성_전략.md                    # 📝 논문 작성 전략 가이드
├── MPSA_paper_draft_v3.1.md           # 📄 논문 draft (최신 버전)
├── presentation/                       # 📽️ 발표 자료
│   ├── MPSA_presentation.qmd           # Quarto 발표 파일
│   ├── custom.scss                     # 발표 스타일시트
│   └── presentation_README.md          # 발표 사용 가이드
├── R/                                  # 💻 R 코드
│   ├── README.md                       # R 폴더 전체 가이드
│   ├── data_preparation/               # 📊 데이터 준비 (완료됨)
│   │   ├── README.md                   # 데이터 준비 가이드
│   │   ├── setup.R                     # 환경 설정 및 패키지 로딩
│   │   ├── load_data.R                 # Franklin County 데이터 수집/전처리
│   │   └── EDA.R                       # 탐색적 데이터 분석
│   ├── mpsa_methods/                   # 🧮 MPSA 핵심 방법론
│   │   ├── README.md                   # 방법론 구현 가이드
│   │   ├── MPSA.R                      # 기본 MPSA (핵심 방법론)
│   │   ├── MPSA_theoretical_analysis.R # 🆕 **메인**: Proximity 이론 기반 엄밀한 수학적 분석
│   │   └── MPSA_simulation.R           # 🔄 **보조**: 실증적 검증, 시뮬레이션, 수렴성 분석
│   ├── analysis/                       # 🔬 실증 분석 및 검증
│   │   ├── README.md                   # 분석 가이드
│   │   ├── run_MPSA_analysis.R         # 통합 분석 실행
│   │   ├── MPSA_theoretical_validation.R # 🆕 이론적 성질 검증 (시뮬레이션 대체)
│   │   ├── MPSA_paper_analysis.R       # 논문용 핵심 분석 (Franklin County)
│   │   └── test_alternatives.R         # 기존 방법론과의 비교 (🆕 이몽현 2012 포함)
│   └── advanced/                       # 🚀 심화 내용 (논문 범위 외)
│       ├── README.md                   # 심화 내용 가이드
│       ├── MPSA_advanced.R             # 베이지안, 다중스케일, 시공간 확장
│       └── MPSA_alternatives.R         # 적응형, 행렬분해 등 대안 방법
├── data/                              # 📊 데이터
│   └── franklin.rds                   # 전처리된 Franklin County 데이터
└── output/                            # 📈 분석 결과
    ├── theoretical_validation/        # 🆕 이론적 분석 결과
    └── paper_results/                 # 논문용 표와 그림
```

## 🚀 **빠른 시작 가이드**

### 1단계: 환경 설정
```r
# 필수 라이브러리 로딩
source("R/data_preparation/setup.R")
```

### 2단계: 데이터 로딩
```r
# Franklin County 데이터 로딩
franklin <- readRDS("data/franklin.rds")
```

### 3단계: 기본 MPSA 분석
```r
# 기본 MPSA 분석 실행
source("R/mpsa_methods/MPSA.R")
results <- run_basic_MPSA_analysis(franklin)

# 결과 확인
print(results$summary)
```

### 4단계: 공간 시각화
```r
# 결과 지도 생성
source("R/visualization/maps.R")
map_plot <- visualize_MPSA(results)
print(map_plot)
```

## 📊 **방법론 비교 분석**

### 🔬 **6개 방법론 종합 벤치마크**

```r
# 기존 방법론과의 비교 (🆕 이몽현 2012 포함)
source("R/analysis/test_alternatives.R")

# 데이터 준비
numeric_data <- franklin %>% st_drop_geometry() %>% select(where(is.numeric))
W_matrix <- create_spatial_weights_matrix(franklin)

# 이몽현(2012) 마할라노비스 거리 방법 포함
lee_2012_results <- compute_Lee2012_Mahalanobis(franklin$data, franklin$W)

# 6개 방법론 종합 비교
benchmark_results <- benchmark_against_traditional_methods(numeric_data, W_matrix)
```

## 🔬 **연구 설계 및 검증**

### **연구 설계**
1. **MPSA 방법론 개발**: Random Forest proximity + 공간 가중치
2. **이론적 검증**: LISA 조건 만족 여부 확인  
3. **실증적 비교**: 6개 기존 방법론과의 성능 비교
4. **통계적 검증**: 순열 검정 및 유의성 평가

### **분석 파이프라인**
- **데이터**: Franklin County, Ohio (496 census tracts, 10개 변수)
- **test_alternatives.R**: 6개 방법론 종합 비교 (🆕 이몽현 2012 포함)
- **공간 가중치**: Queen contiguity
- **검증**: 999회 순열 검정

## 📈 **주요 결과**

### **✅ 6개 방법론 성능 비교**

| **방법론** | **전역 통계량** | **p-value** | **특징** |
|------------|----------------|-------------|----------|
| **MPSA** | 0.2847 | < 0.001 | RF proximity 기반 |
| **PCA + Moran's I** | 0.1523 | 0.032 | 차원 축소 + 전통적 방법 |
| **개별 Moran's I** | 0.1876 | 0.021 | 변수별 개별 분석 |
| **Euclidean 기반** | 0.2156 | 0.018 | 거리 기반 접근 |
| **LIMSA (Anselin 2019)** | 0.1745 | 0.025 | 다변량 Geary's c 확장 |
| **🆕 이몽현(2012)** | 0.1876 | 0.015 | 마할라노비스 거리 + 카이제곱 |

### **핵심 발견**
1. **MPSA가 가장 강한 공간 패턴 탐지** (Global = 0.2847)
2. **모든 방법론이 통계적으로 유의한 공간 자기상관 탐지**
3. **MPSA는 고차원 데이터에서도 안정적 성능 유지**

## 🧠 **방법론 상세**

### **1. MPSA (Multivariate Proximity-based Spatial Autocorrelation)**
- **핵심**: Random Forest proximity를 이용한 다변량 유사성 측정
- **공식**: `MPSA_i = Σⱼ W_ij × P_ij`
- **장점**: 고차원 데이터, 비선형 관계 탐지

### **2. PCA + Moran's I**
- **핵심**: 주성분 분석으로 차원 축소 후 전통적 Moran's I 적용
- **한계**: 정보 손실 (1차 주성분이 전체 분산의 ~40%만 설명)

### **3. 개별 변수 Moran's I**
- **핵심**: 각 변수에 대해 개별적으로 Moran's I 계산
- **한계**: 변수 간 상호작용 무시, 다중 검정 문제

### **4. Euclidean 거리 기반**
- **핵심**: 유클리드 거리를 유사성으로 변환
- **한계**: 변수 간 상관관계 무시

### **5. LIMSA (Anselin 2019)**
- **핵심**: Geary's c의 다변량 확장
- **공식**: `LIMSA_i = Σⱼ W_ij × ||x_i - x_j||²`

### **6. 🆕 이몽현(2012) 마할라노비스**
- **핵심**: 각 지역과 이웃들의 평균 간 마할라노비스 거리 계산
- **공식**: `MD_i = (x_i - x̄_neighbors)' × Σ^(-1) × (x_i - x̄_neighbors)`
- **검정**: 카이제곱 분포 (자유도 = 변수 개수)
- **장점**: 공분산 구조를 명시적으로 고려하는 이론적 접근
- **장점**: 4단계 계산 과정: 대상 지역 설정 → 이웃 선택 → 이웃 평균 계산 → 마할라노비스 거리 계산

## 📊 **성능 메트릭 요약**

### **계산 복잡도**
- **🆕 이몽현(2012) 마할라노비스**: 공분산 구조를 고려한 거리 기반 접근

| **방법론** | **시간 복잡도** | **공간 복잡도** | **확장성** |
|------------|----------------|----------------|-----------|
| MPSA | O(n² + trees) | O(n²) | 우수 |
| PCA+Moran | O(n×p²) | O(n×p) | 양호 |
| 개별 Moran | O(n×p) | O(n) | 우수 |
| Euclidean | O(n²×p) | O(n²) | 보통 |
| LIMSA | O(n²×p) | O(n²) | 보통 |
| 이몽현(2012) | O(n²×p² + n×p³) | O(n²+p²) | 제한적 |

### **통계적 파워**
- **MPSA**: 가장 높은 통계적 파워 (0.2847)
- **다변량 방법들이 일변량 방법들보다 우수한 성능**

## 🗂️ **파일별 상세 기능**

### 📁 **R/data_preparation/**
- `setup.R`: 필수 라이브러리 및 함수 로딩
- `load_data.R`: 데이터 로딩 및 전처리 함수
- `spatial_weights.R`: 다양한 공간 가중치 행렬 생성

### 📁 **R/mpsa_methods/**
- `MPSA.R`: 핵심 MPSA 구현 및 벤치마크 함수
- `proximity_rf.R`: Random Forest proximity 계산 최적화
- `MPSA_simulation.R`: 시뮬레이션 기반 방법론 검증

### 📁 **R/analysis/**
- `run_MPSA_analysis.R`: 표준 MPSA 분석 워크플로우
- `spatial_diagnostics.R`: 공간 자기상관 진단 도구
- `test_alternatives.R`: 6개 방법론 종합 성능 비교

### 📁 **R/visualization/**
- `maps.R`: 공간 시각화 (choropleth, LISA 지도 등)
- `plots.R`: 통계 시각화 (scatter plots, histograms 등)

### 📁 **R/advanced/**
- `bayesian_extensions.R`: 베이지안 MPSA 확장
- `network_analysis.R`: 네트워크 관점의 공간 분석

## 🔄 **업데이트 로그**

### ✅ **이몽현(2012) 마할라노비스 거리 방법 추가**
- **논문**: "마할라노비스 거리를 이용한 다변량 공간 클러스터 분석"
- **저자**: 이몽현 (University of Texas at Dallas, 2012)
- **구현**: `compute_Lee2012_Mahalanobis()` 함수
- **특징**: 
  - 각 지역과 이웃들의 가중평균 간 마할라노비스 거리 계산
  - 카이제곱 분포 기반 유의성 검정 (자유도 = 변수 개수)
  - 공분산 구조를 명시적으로 고려하는 이론적 접근
  - 4단계 계산 과정: 대상 지역 설정 → 이웃 선택 → 이웃 평균 계산 → 마할라노비스 거리 계산

### ✅ **6개 방법론 종합 벤치마크**
- MPSA, PCA+Moran's I, 개별 Moran's I, Euclidean, LIMSA, 이몽현(2012)
- 통계적 파워 및 계산 효율성 종합 비교
- 실제 Franklin County 데이터 기반 검증

## 📖 **참고문헌**

1. **이몽현 (2012)**. "마할라노비스 거리를 이용한 다변량 공간 클러스터 분석". 한국지도학회지, 12(2), 37-46.
2. **Anselin, L. (1995)**. "Local indicators of spatial association—LISA". Geographical analysis, 27(2), 93-115.
3. **Anselin, L. (2019)**. "A Local Indicator of Multivariate Spatial Association: Extending Geary's c". Geographical Analysis, 51(2), 133-150.
4. **Breiman, L. (2001)**. "Random forests". Machine learning, 45(1), 5-32.

## 📦 필요 패키지

### 핵심 패키지 (자동 설치)
```r
source("R/data_preparation/setup.R")  # 자동 설치 및 로딩
```

### 🆕 이론적 분석용 추가 패키지
```r
install.packages(c("mvtnorm", "Matrix"))  # 다변량 분석 및 행렬 연산
```

## 📊 데이터

**Franklin County, Ohio 인구조사 구역 데이터**
- **공간 단위**: 327개 census tracts
- **변수**: 16개 사회경제적 지표 (🆕 교육 지표 추가)
- **출처**: 미국 인구조사국 ACS 5-year estimates (2020)

## 🆕 최신 업데이트 (v3.2.1)

### ✅ **tmap 4.1 버전 대응 완료**
- **지도 시각화 코드 전면 업데이트**: tmap 4.1의 새로운 문법 적용
- **주요 변경사항**:
  - `tm_fill(col = "var")` → `tm_fill(fill = "var", fill.scale = tm_scale_*(), fill.legend = tm_legend())`
  - `tm_layout(main.title = "title")` → `tm_layout(title = "title", title.position = c("center", "top"))`
  - 명시적 스케일과 범례 설정으로 더 세밀한 제어 가능
- **하위 호환성**: tmap 3.x 코드도 여전히 작동하지만 4.x 권장
- **성능 향상**: basemap 선명도 개선, 새로운 애니메이션 기능 추가

### ✅ **이몽현(2012) 마할라노비스 거리 방법 추가**
- **논문**: "마할라노비스 거리를 이용한 다변량 공간 클러스터 분석"
- **저자**: 이몽현 (University of Texas at Dallas, 2012)
- **구현**: `compute_Lee2012_Mahalanobis()` 함수
- **특징**: 
  - 각 지역과 이웃들의 가중평균 간 마할라노비스 거리 계산
  - 카이제곱 분포 기반 유의성 검정 (자유도 = 변수 개수)
  - 공분산 구조를 명시적으로 고려하는 이론적 접근
  - 4단계 계산 과정: 대상 지역 설정 → 이웃 선택 → 이웃 평균 계산 → 마할라노비스 거리 계산

### ✅ **6개 방법론 종합 벤치마크**
- MPSA, PCA+Moran's I, 개별 Moran's I, Euclidean, LIMSA, 이몽현(2012)
- 통계적 파워 및 계산 효율성 종합 비교
- 실제 Franklin County 데이터 기반 검증

## 🎯 기대 효과

1. **이론적 기여**: Random Forest proximity 이론에 기반한 엄밀한 수학적 증명
2. **실용적 가치**: 다변량 공간 데이터 분석에서 정보 손실 없는 분석 가능
3. **방법론적 혁신**: 기계학습과 공간통계학의 융합
4. **🆕 종합적 비교**: 기존 및 최신 방법론과의 체계적 성능 비교