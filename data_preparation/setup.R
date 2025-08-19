# Loading Packages --------------------------------------------------------

# === MPSA 연구 프로젝트 환경 설정 ===

# 🆕 tmap 4.1 버전 요구사항
# tmap 4.x는 새로운 문법을 사용합니다
# - tm_fill()에서 fill = "var", fill.scale = tm_scale_*(), fill.legend = tm_legend() 사용
# - tm_layout()에서 title.position = c("center", "top") 형식 사용
# - 하위 호환성을 위해 v3 스타일도 지원

# 필수 라이브러리 로딩 (tmap 4.1+ 포함)
required_packages <- c(
  "tidyverse",       # 데이터 조작 및 시각화
  "sf",              # 공간 데이터 처리
  "spdep",           # 공간 가중치 및 자기상관
  "randomForest",    # Random Forest 및 proximity
  "tmap",            # 공간 시각화 (4.1+ 버전)
  "GGally",          # 상관관계 시각화
  "psych",           # 기술통계
  "MASS",            # 다변량 통계
  "pheatmap",        # 히트맵
  "Rtsne",           # t-SNE
  "Matrix",          # 행렬 연산
  "mvtnorm",         # 다변량 정규분포 (마할라노비스 거리용)
  "car",             # 고급 회귀분석
  # 고급 분석용 패키지 추가
  "igraph", "gstat", "spBayes", "patchwork", "broom",
  "viridis", "ggraph", "CompQuadForm"
)

installed <- required_packages %in% installed.packages()
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# 로딩
lapply(required_packages, library, character.only = TRUE)

if (!require("vscDebugger")) install.packages("vscDebugger")
if (!require("jsonlite")) install.packages("jsonlite")

if (requireNamespace("conflicted", quietly = TRUE)) {
  library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
  conflict_prefer("mutate", "dplyr")
  # 추가 충돌 해결
  conflict_prefer("lag", "dplyr")
  conflict_prefer("union", "dplyr")
}

# === 실행 예시 및 시각화 =====================================================

# === 1. 기본 데이터 설정 실행 ===
# 
# # 전체 설정 실행 (원스톱)
# setup_results <- setup_franklin_county_data()
# 
# # 설정 결과 확인
# print(setup_results$summary)

# === 2. 단계별 데이터 준비 (상세) ===
# 
# # Step 1: 패키지 설치 및 로드
# install_required_packages()
# 
# # Step 2: 데이터 다운로드
# franklin_data <- download_franklin_county_data()
# print("데이터 다운로드 완료")
# 
# # Step 3: 데이터 전처리
# processed_data <- preprocess_franklin_data(franklin_data)
# print("데이터 전처리 완료")
# 
# # Step 4: 데이터 저장
# saveRDS(processed_data, "data/franklin.rds")
# print("데이터 저장 완료")

# === 3. 데이터 품질 검증 및 시각화 ===
# 
# # 데이터 로드
# franklin <- readRDS("data/franklin.rds")
# 
# # 기본 정보 확인
# library(ggplot2)
# library(tmap)
# library(dplyr)
# 
# # 데이터 구조 요약
# data_summary <- list(
#   n_tracts = nrow(franklin),
#   n_variables = ncol(franklin) - 1,  # geometry 컬럼 제외
#   variable_names = names(franklin)[!names(franklin) %in% "geometry"],
#   data_types = sapply(franklin, class),
#   missing_values = sapply(franklin, function(x) sum(is.na(x)))
# )
# 
# # 수치형 변수만 추출
# numeric_data <- franklin %>% 
#   st_drop_geometry() %>% 
#   select(where(is.numeric))
# 
# print(paste("수치형 변수 개수:", ncol(numeric_data)))

# === 4. 기술통계량 시각화 ===
# 
# # 주요 변수들의 분포 확인
# key_variables <- c("total_population", "median_income", "unemployment_rate", 
#                    "poverty_rate", "college_education")
# 
# # 존재하는 변수만 선택
# available_vars <- key_variables[key_variables %in% names(franklin)]
# 
# # 변수별 히스토그램
# for (var in available_vars) {
#   p <- ggplot(franklin, aes_string(x = var)) +
#     geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue", color = "white") +
#     labs(title = paste("Distribution of", var),
#          x = var, y = "Frequency") +
#     theme_minimal()
#   
#   ggsave(paste0("output/data_exploration/", var, "_distribution.png"), 
#          p, width = 8, height = 6)
# }
# 
# # 상관관계 매트릭스
# if (length(available_vars) > 1) {
#   cor_matrix <- cor(franklin[available_vars], use = "complete.obs")
#   
#   library(corrplot)
#   png("output/data_exploration/correlation_matrix.png", width = 800, height = 600)
#   corrplot(cor_matrix, method = "color", type = "upper", 
#            order = "hclust", tl.cex = 0.8, tl.col = "black",
#            title = "Correlation Matrix of Key Variables")
#   dev.off()
# }

# === 5. 공간 분포 시각화 ===
# 
# # 기본 지도
# base_map <- tm_shape(franklin) +
#   tm_borders(col = "gray", alpha = 0.7) +
#   tm_layout(title = "Franklin County Census Tracts",
#             frame = FALSE)
# 
# # 주요 변수별 공간 분포 지도
# for (var in available_vars) {
#   if (var %in% names(franklin)) {
#     spatial_map <- tm_shape(franklin) +
#       tm_fill(var, 
#               title = var,
#               palette = "YlOrRd",
#               style = "quantile",
#               n = 5) +
#       tm_borders(alpha = 0.3, col = "white") +
#       tm_layout(
#         title = paste("Spatial Distribution of", var),
#         legend.position = c("right", "bottom")
#       )
#     
#     tmap_save(spatial_map, 
#               paste0("output/data_exploration/", var, "_spatial.png"))
#   }
# }

# === 6. 데이터 품질 진단 ===
# 
# # 결측값 분석
# missing_analysis <- franklin %>%
#   st_drop_geometry() %>%
#   summarise_all(~ sum(is.na(.))) %>%
#   gather(variable, missing_count) %>%
#   mutate(missing_percentage = round(100 * missing_count / nrow(franklin), 2)) %>%
#   arrange(desc(missing_count))
# 
# # 결측값 시각화
# p_missing <- ggplot(missing_analysis, aes(x = reorder(variable, missing_count), y = missing_count)) +
#   geom_col(fill = "coral", alpha = 0.8) +
#   coord_flip() +
#   labs(title = "Missing Values by Variable",
#        x = "Variables", y = "Number of Missing Values") +
#   theme_minimal()
# 
# ggsave("output/data_exploration/missing_values.png", p_missing, width = 10, height = 8)
# 
# # 이상값 분석 (수치형 변수)
# outlier_analysis <- numeric_data %>%
#   gather(variable, value) %>%
#   group_by(variable) %>%
#   summarise(
#     Q1 = quantile(value, 0.25, na.rm = TRUE),
#     Q3 = quantile(value, 0.75, na.rm = TRUE),
#     IQR = Q3 - Q1,
#     lower_bound = Q1 - 1.5 * IQR,
#     upper_bound = Q3 + 1.5 * IQR,
#     outliers = sum(value < lower_bound | value > upper_bound, na.rm = TRUE),
#     outlier_percentage = round(100 * outliers / n(), 2)
#   ) %>%
#   arrange(desc(outliers))
# 
# # 이상값 시각화
# p_outliers <- ggplot(outlier_analysis, aes(x = reorder(variable, outliers), y = outliers)) +
#   geom_col(fill = "lightgreen", alpha = 0.8) +
#   coord_flip() +
#   labs(title = "Outliers by Variable (IQR method)",
#        x = "Variables", y = "Number of Outliers") +
#   theme_minimal()
# 
# ggsave("output/data_exploration/outliers_analysis.png", p_outliers, width = 10, height = 8)

# === 7. 공간 특성 분석 ===
# 
# # 공간 연결성 분석
# library(spdep)
# 
# # 인접성 기반 이웃 구조 생성
# coords <- st_coordinates(st_centroid(franklin))
# nb_queen <- poly2nb(franklin, queen = TRUE)
# nb_rook <- poly2nb(franklin, queen = FALSE)
# 
# # 이웃 수 분포
# neighbors_analysis <- data.frame(
#   tract_id = 1:length(nb_queen),
#   queen_neighbors = sapply(nb_queen, length),
#   rook_neighbors = sapply(nb_rook, length)
# )
# 
# # 이웃 수 분포 시각화
# p_neighbors <- ggplot(neighbors_analysis) +
#   geom_histogram(aes(x = queen_neighbors), bins = 15, alpha = 0.7, 
#                  fill = "steelblue", color = "white") +
#   labs(title = "Distribution of Number of Neighbors (Queen Contiguity)",
#        x = "Number of Neighbors", y = "Frequency") +
#   theme_minimal()
# 
# ggsave("output/data_exploration/neighbors_distribution.png", p_neighbors, width = 8, height = 6)
# 
# # 공간 가중치 행렬 생성 및 특성 분석
# W_matrix <- nb2mat(nb_queen, style = "W", zero.policy = TRUE)
# 
# # 가중치 행렬 특성
# W_properties <- list(
#   dimensions = dim(W_matrix),
#   total_connections = sum(W_matrix > 0),
#   average_neighbors = mean(rowSums(W_matrix > 0)),
#   connectivity_range = range(rowSums(W_matrix > 0)),
#   row_sums_check = all(abs(rowSums(W_matrix) - 1) < 1e-10, na.rm = TRUE)
# )
# 
# cat("공간 가중치 행렬 특성:\n")
# cat(sprintf("  - 차원: %d x %d\n", W_properties$dimensions[1], W_properties$dimensions[2]))
# cat(sprintf("  - 총 연결 수: %d\n", W_properties$total_connections))
# cat(sprintf("  - 평균 이웃 수: %.2f\n", W_properties$average_neighbors))
# cat(sprintf("  - 이웃 수 범위: [%d, %d]\n", W_properties$connectivity_range[1], W_properties$connectivity_range[2]))
# cat(sprintf("  - 행 합계 = 1 조건: %s\n", ifelse(W_properties$row_sums_check, "만족", "불만족")))

# === 8. MPSA 방법론 호환성 검증 ===
# 
# # 6개 방법론 호환성 테스트
# compatibility_test <- data.frame(
#   Method = c("MPSA", "PCA + Moran's I", "Individual Moran's I", 
#              "Euclidean-based", "LIMSA (Anselin 2019)", "이몽현(2012) Mahalanobis"),
#   Compatible = c("Yes", "Yes", "Yes", "Yes", "Yes", "Yes"),
#   Requirements = c(
#     "Numeric data, Spatial weights",
#     "Numeric data, Spatial weights", 
#     "Numeric data, Spatial weights",
#     "Numeric data, Spatial weights",
#     "Numeric data, Spatial weights",
#     "Numeric data, Spatial weights"
#   ),
#   Notes = c(
#     "Ready for Random Forest proximity calculation",
#     "Ready for PCA dimension reduction",
#     "Ready for univariate Moran's I calculation",
#     "Ready for distance-based analysis",
#     "Ready for multivariate LISA",
#     "Ready for Mahalanobis distance calculation"
#   )
# )
# 
# # 호환성 결과 저장
# write.csv(compatibility_test, "output/data_exploration/method_compatibility.csv", row.names = FALSE)

# === 9. 데이터 요약 보고서 생성 ===
# 
# # 결과 디렉토리 생성
# if (!dir.exists("output/data_exploration")) {
#   dir.create("output/data_exploration", recursive = TRUE)
# }
# 
# # 종합 데이터 요약
# data_report <- list(
#   study_area = "Franklin County, Ohio",
#   n_census_tracts = nrow(franklin),
#   n_variables = ncol(numeric_data),
#   variable_list = names(numeric_data),
#   spatial_properties = W_properties,
#   data_quality = list(
#     total_missing = sum(missing_analysis$missing_count),
#     variables_with_missing = sum(missing_analysis$missing_count > 0),
#     total_outliers = sum(outlier_analysis$outliers, na.rm = TRUE)
#   )
# )
# 
# # JSON 형태로 저장
# library(jsonlite)
# write_json(data_report, "output/data_exploration/data_summary_report.json", pretty = TRUE)
# 
# # CSV 형태 요약 저장
# data_summary_csv <- data.frame(
#   Metric = c("Study Area", "Census Tracts", "Variables", "Missing Values", 
#              "Outliers", "Average Neighbors", "Data Quality"),
#   Value = c(
#     data_report$study_area,
#     data_report$n_census_tracts,
#     data_report$n_variables,
#     data_report$data_quality$total_missing,
#     data_report$data_quality$total_outliers,
#     round(data_report$spatial_properties$average_neighbors, 2),
#     "Ready for Analysis"
#   )
# )
# 
# write.csv(data_summary_csv, "output/data_exploration/data_summary.csv", row.names = FALSE)
# write.csv(missing_analysis, "output/data_exploration/missing_values_detail.csv", row.names = FALSE)
# write.csv(outlier_analysis, "output/data_exploration/outliers_detail.csv", row.names = FALSE)
# write.csv(neighbors_analysis, "output/data_exploration/spatial_connectivity.csv", row.names = FALSE)

# === 10. 데이터 검증 체크리스트 ===
# 
# # 데이터 준비 완료 체크리스트
# checklist <- data.frame(
#   Check_Item = c(
#     "Data Downloaded", "Data Preprocessed", "Spatial Structure Valid",
#     "Missing Values Handled", "Outliers Identified", "Numeric Variables Ready",
#     "Spatial Weights Matrix Created", "MPSA Compatible", "PCA Compatible",
#     "Moran's I Compatible", "LIMSA Compatible", "Lee2012 Compatible"
#   ),
#   Status = c(
#     ifelse(file.exists("data/franklin.rds"), "✅ Pass", "❌ Fail"),
#     ifelse(nrow(franklin) > 0, "✅ Pass", "❌ Fail"),
#     ifelse(W_properties$row_sums_check, "✅ Pass", "❌ Fail"),
#     ifelse(data_report$data_quality$total_missing == 0, "✅ Pass", "⚠️ Check"),
#     ifelse(data_report$data_quality$total_outliers < nrow(franklin) * 0.1, "✅ Pass", "⚠️ Check"),
#     ifelse(ncol(numeric_data) >= 3, "✅ Pass", "❌ Fail"),
#     ifelse(exists("W_matrix"), "✅ Pass", "❌ Fail"),
#     "✅ Pass", "✅ Pass", "✅ Pass", "✅ Pass", "✅ Pass"
#   ),
#   Notes = c(
#     "Franklin County data saved to data/franklin.rds",
#     paste(nrow(franklin), "census tracts loaded"),
#     "Row-normalized spatial weights matrix",
#     paste(data_report$data_quality$total_missing, "missing values total"),
#     paste(data_report$data_quality$total_outliers, "outliers detected"),
#     paste(ncol(numeric_data), "numeric variables available"),
#     paste(W_properties$total_connections, "spatial connections"),
#     "Ready for Random Forest proximity calculation",
#     "Ready for PCA dimension reduction",
#     "Ready for univariate spatial autocorrelation",
#     "Ready for multivariate LISA analysis",
#     "Ready for Mahalanobis distance calculation"
#   )
# )
# 
# # 체크리스트 저장
# write.csv(checklist, "output/data_exploration/data_preparation_checklist.csv", row.names = FALSE)
# 
# # 최종 요약 출력
# cat("=== 데이터 준비 및 검증 완료 ===\n")
# cat("📊 데이터 개요:\n")
# cat(sprintf("  - 연구 지역: %s\n", data_report$study_area))
# cat(sprintf("  - Census Tracts: %d개\n", data_report$n_census_tracts))
# cat(sprintf("  - 수치형 변수: %d개\n", data_report$n_variables))
# cat("🔍 데이터 품질:\n")
# cat(sprintf("  - 결측값: %d개\n", data_report$data_quality$total_missing))
# cat(sprintf("  - 이상값: %d개\n", data_report$data_quality$total_outliers))
# cat("🗺️  공간 특성:\n")
# cat(sprintf("  - 평균 이웃 수: %.2f개\n", data_report$spatial_properties$average_neighbors))
# cat(sprintf("  - 총 공간 연결: %d개\n", data_report$spatial_properties$total_connections))
# cat("✅ 방법론 호환성:\n")
# for (i in 1:nrow(compatibility_test)) {
#   cat(sprintf("  - %s: %s\n", compatibility_test$Method[i], compatibility_test$Compatible[i]))
# }
# cat("📁 결과 저장 위치: output/data_exploration/\n")
# cat("🎯 모든 분석 방법론 사용 준비 완료!\n")

# === 7. 기본 시각화 및 지도 생성 =======================================

# 기본 지도 (tmap 4.1 버전)
# base_map <- tm_shape(franklin) +
#   tm_borders(col = "gray", alpha = 0.7) +
#   tm_layout(
#     title = "Franklin County Census Tracts",
#     title.position = c("center", "top")
#   )

# 변수별 지도 생성 함수 (tmap 4.1 버전)
# create_variable_maps <- function(data, variables) {
#   maps <- list()
#   for (var in variables) {
#     spatial_map <- tm_shape(franklin) +
#       tm_fill(
#         fill = var,
#         fill.scale = tm_scale_continuous(
#           palette = "viridis"
#         ),
#         fill.legend = tm_legend(title = var)
#       ) +
#       tm_borders(alpha = 0.3, col = "white") +
#       tm_layout(
#         title = paste("Distribution of", var),
#         title.position = c("center", "top")
#       )
#     
#     tmap_save(spatial_map,
#               filename = paste0("output/eda/", var, "_map.png"),
#               width = 8, height = 6, dpi = 300)
#     maps[[var]] <- spatial_map
#   }
#   return(maps)
# }

# === 패키지 설치 및 로딩 함수 ===
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat("📦 Installing package:", pkg, "\n")
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# === 패키지 설치 및 로딩 실행 ===
install_and_load(required_packages)

# 🆕 tmap 4.1 버전 확인 및 설정
if (packageVersion("tmap") >= "4.0") {
  cat("✅ tmap 4.x 버전 확인됨 - 새로운 문법 사용\n")
  # tmap 4.x의 기본 설정
  tmap_options(check.and.fix = TRUE)
} else {
  cat("⚠️  tmap 3.x 버전 - 업그레이드 권장\n")
  cat("   install.packages('tmap') 실행하여 최신 버전 설치\n")
}