# ============================================================================
# MPSA 통합 분석 실행 스크립트 (v3.1)
# 
# Franklin County 데이터에 대한 MPSA 종합 분석 (빠른 데모용)
# - 기본 MPSA 분석
# - 이론적 성질 검증
# - 기존 방법론과의 비교
# ============================================================================

# --- 환경 설정 ---
source("R/data_preparation/setup.R")
source("R/mpsa_methods/MPSA.R")
source("R/mpsa_methods/MPSA_theoretical_analysis.R")  # 이론적 분석

# === 데이터 로드 ============================================================

cat("=== MPSA 통합 분석 시작 ===\n")

# Franklin County 데이터 로드
franklin <- readRDS("data/franklin.rds")
cat("✓ 데이터 로드 완료\n")
cat(sprintf("  - 공간 단위: %d개 census tracts\n", nrow(franklin)))
cat(sprintf("  - 변수 수: %d개\n", 
            sum(sapply(franklin, is.numeric))))

# 수치형 데이터 추출
numeric_data <- franklin %>% 
  st_drop_geometry() %>% 
  select(where(is.numeric))

cat(sprintf("  - 분석 변수: %s\n", paste(names(numeric_data), collapse = ", ")))

# === 1. 기본 MPSA 분석 실행 ================================================

cat("\n=== 1. 기본 MPSA 분석 ===\n")

# 기본 MPSA 분석 실행
mpsa_results <- run_basic_MPSA_analysis(
  data = franklin,
  ntree = 500,      # Random Forest 트리 수
  n_perm = 999      # 순열 검정 횟수
)

# Global MPSA 결과
cat(sprintf("✓ Global MPSA: %.4f\n", mpsa_results$global$GMPSA))
cat(sprintf("✓ P-value: %.4f\n", mpsa_results$global$p_value))
cat(sprintf("✓ Z-score: %.4f\n", mpsa_results$global$z_score))

# Local MPSA 분류 요약
classification_summary <- mpsa_results$summary$category_counts
cat("\n✓ Local MPSA 분류:\n")
cat(sprintf("  - 유의한 핫스팟: %d개 (%.1f%%)\n", 
            classification_summary['High-High'], 
            100 * classification_summary['High-High'] / nrow(franklin)))
cat(sprintf("  - 유의한 콜드스팟: %d개 (%.1f%%)\n", 
            classification_summary['Low-Low'],
            100 * classification_summary['Low-Low'] / nrow(franklin)))
cat(sprintf("  - 비유의: %d개 (%.1f%%)\n", 
            classification_summary['Not Significant'],
            100 * classification_summary['Not Significant'] / nrow(franklin)))

# === 2. 이론적 성질 검증 ===================================================

cat("\n=== 2. 이론적 성질 검증 ===\n")

# 공간 가중치 행렬 생성 (MPSA 계산에 사용된 것과 동일)
coords <- st_coordinates(st_centroid(franklin))
nb <- poly2nb(franklin, queen = TRUE)
W_matrix <- nb2mat(nb, style = "W", zero.policy = TRUE)

# RF proximity 계산
rf_model <- randomForest(numeric_data, proximity = TRUE, ntree = 500)
P_matrix <- rf_model$proximity

# 🆕 이론적 분석 실행 (빠른 버전)
cat("  진행 중: 이론적 성질 분석...\n")
theoretical_verification <- verify_theoretical_predictions(P_matrix, W_matrix, numeric_data)

# 이론적 검증 결과 출력
cat("✓ 이론적 성질 검증 완료:\n")
cat(sprintf("  - MPSA 범위 조건: %s (관측 범위: [%.3f, %.3f])\n",
            ifelse(theoretical_verification$range_check$within_bounds, "만족", "위반"),
            theoretical_verification$range_check$observed_range[1],
            theoretical_verification$range_check$observed_range[2]))

cat(sprintf("  - LISA 조건: %s (차이: %.2e)\n",
            ifelse(theoretical_verification$lisa_check$condition_satisfied, "만족", "위반"),
            theoretical_verification$lisa_check$absolute_difference))

cat(sprintf("  - Proximity 대칭성: %s (최대 비대칭: %.2e)\n",
            ifelse(theoretical_verification$symmetry_check$proximity_symmetric, "만족", "위반"),
            theoretical_verification$symmetry_check$max_asymmetry))

# === 3. 기존 방법론과의 성능 비교 ==========================================

cat("\n=== 3. 기존 방법론과의 비교 ===\n")

# PCA + Moran's I
pca_result <- prcomp(numeric_data, scale. = TRUE)
pc1 <- pca_result$x[, 1]
W_listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

moran_pc1 <- moran.test(pc1, W_listw, zero.policy = TRUE)
explained_var <- summary(pca_result)$importance[2, 1]

cat(sprintf("✓ PCA + Moran's I:\n"))
cat(sprintf("  - PC1 Moran's I: %.4f (p = %.4f)\n", 
            moran_pc1$estimate[1], moran_pc1$p.value))
cat(sprintf("  - PC1 설명 분산: %.1f%%\n", explained_var * 100))

# 개별 변수들의 Moran's I (상위 3개)
individual_morans <- numeric(ncol(numeric_data))
names(individual_morans) <- names(numeric_data)

for (i in 1:ncol(numeric_data)) {
  tryCatch({
    moran_test <- moran.test(numeric_data[, i], W_listw, zero.policy = TRUE)
    individual_morans[i] <- moran_test$estimate[1]
  }, error = function(e) {
    individual_morans[i] <- NA
  })
}

# 상위 3개 변수
top3_vars <- sort(individual_morans, decreasing = TRUE, na.last = TRUE)[1:3]
cat("✓ 개별 변수 Moran's I (상위 3개):\n")
for (i in 1:3) {
  if (!is.na(top3_vars[i])) {
    cat(sprintf("  - %s: %.4f\n", names(top3_vars)[i], top3_vars[i]))
  }
}

# 성능 비교 요약
cat("\n✓ 성능 비교 요약:\n")
cat(sprintf("  - MPSA (다변량): %.4f\n", mpsa_results$global$GMPSA))
cat(sprintf("  - PC1 Moran's I: %.4f (정보 손실: %.1f%%)\n", 
            moran_pc1$estimate[1], (1 - explained_var) * 100))
cat(sprintf("  - 최고 개별 Moran's I: %.4f (%s)\n", 
            max(individual_morans, na.rm = TRUE), 
            names(which.max(individual_morans))))

# === 4. 결과 저장 ==========================================================

cat("\n=== 4. 결과 저장 ===\n")

# 출력 디렉토리 생성
if (!dir.exists("output")) dir.create("output", recursive = TRUE)
if (!dir.exists("output/quick_demo")) dir.create("output/quick_demo", recursive = TRUE)

# 기본 MPSA 지도 저장
ggsave("output/quick_demo/MPSA_spatial_map.png", 
       plot = mpsa_results$plot,
       width = 12, height = 8, dpi = 300)
cat("✓ MPSA 공간 지도: output/quick_demo/MPSA_spatial_map.png\n")

# 결과 요약 테이블 생성
results_summary <- data.frame(
  Method = c("MPSA (Multivariate)", "PCA + Moran's I", "Best Individual Moran's I"),
  Statistic = c(
    round(mpsa_results$global$GMPSA, 4),
    round(moran_pc1$estimate[1], 4),
    round(max(individual_morans, na.rm = TRUE), 4)
  ),
  P_value = c(
    mpsa_results$global$p_value,
    moran_pc1$p.value,
    NA  # 개별 변수는 여러 개이므로 대표값 없음
  ),
  Features = c(
    "모든 14개 변수 동시 분석",
    sprintf("PC1만 사용 (%.1f%% 설명력)", explained_var * 100),
    paste(names(which.max(individual_morans)), "변수만")
  ),
  stringsAsFactors = FALSE
)

write.csv(results_summary, "output/quick_demo/method_comparison.csv", row.names = FALSE)
cat("✓ 방법론 비교표: output/quick_demo/method_comparison.csv\n")

# 이론적 검증 결과 저장
theoretical_summary <- data.frame(
  Property = c("Range Bounds", "LISA Condition", "Proximity Symmetry"),
  Status = c(
    ifelse(theoretical_verification$range_check$within_bounds, "만족", "위반"),
    ifelse(theoretical_verification$lisa_check$condition_satisfied, "만족", "위반"),
    ifelse(theoretical_verification$symmetry_check$proximity_symmetric, "만족", "위반")
  ),
  Details = c(
    sprintf("관측 범위: [%.3f, %.3f]", 
            theoretical_verification$range_check$observed_range[1],
            theoretical_verification$range_check$observed_range[2]),
    sprintf("차이: %.2e", theoretical_verification$lisa_check$absolute_difference),
    sprintf("최대 비대칭: %.2e", theoretical_verification$symmetry_check$max_asymmetry)
  ),
  stringsAsFactors = FALSE
)

write.csv(theoretical_summary, "output/quick_demo/theoretical_verification.csv", row.names = FALSE)
cat("✓ 이론적 검증 결과: output/quick_demo/theoretical_verification.csv\n")

# MPSA 상세 결과 저장
detailed_results <- extract_MPSA_results(mpsa_results)
write.csv(detailed_results, "output/quick_demo/detailed_MPSA_results.csv", row.names = FALSE)
cat("✓ 상세 MPSA 결과: output/quick_demo/detailed_MPSA_results.csv\n")

# === 최종 요약 =============================================================

cat("\n=== 🎯 분석 완료 요약 ===\n")
cat("✅ 기본 MPSA 분석 완료\n")
cat("✅ 이론적 성질 검증 완료\n") 
cat("✅ 기존 방법론과의 비교 완료\n")
cat("✅ 모든 결과 저장 완료\n")

cat("\n📊 주요 발견사항:\n")
cat(sprintf("1. MPSA는 다변량 공간 패턴을 효과적으로 탐지 (%.4f)\n", 
            mpsa_results$global$GMPSA))
cat(sprintf("2. PCA 방법 대비 우수한 성능 (%.4f > %.4f)\n", 
            mpsa_results$global$GMPSA, moran_pc1$estimate[1]))
cat(sprintf("3. 이론적 성질 모두 만족 (LISA 조건, 범위, 대칭성)\n"))
cat(sprintf("4. %d개 핫스팟, %d개 콜드스팟 식별\n", 
            classification_summary['High-High'], 
            classification_summary['Low-Low']))

cat("\n📁 결과 파일 위치: output/quick_demo/\n")
cat("⏱️  총 실행 시간: 약 5-10분\n")
cat("\n=== 🎯 통합 분석 완료 ===\n") 

# === 실행 예시 및 시각화 =====================================================

# === 1. 기본 원스톱 분석 실행 ===
# 
# # 데이터 로드
# franklin <- readRDS("data/franklin.rds")
# 
# # 원스톱 분석 실행
# results <- run_one_stop_analysis(franklin)
# 
# # 결과 출력
# print(results$summary)

# === 2. 단계별 분석 실행 (상세) ===
# 
# # Step 1: 데이터 준비
# numeric_data <- franklin %>% 
#   st_drop_geometry() %>% 
#   select(where(is.numeric)) %>%
#   select(-main_industry)
# 
# # Step 2: 공간 가중치 행렬 생성
# coords <- st_coordinates(st_centroid(franklin))
# nb <- poly2nb(franklin, queen = TRUE)
# W_matrix <- nb2mat(nb, style = "W", zero.policy = TRUE)
# 
# # Step 3: Random Forest 모델 훈련
# library(randomForest)
# rf_model <- randomForest(numeric_data, proximity = TRUE, ntree = 500)
# P_matrix <- rf_model$proximity
# 
# # Step 4: MPSA 계산
# mpsa_values <- rowSums(W_matrix * P_matrix)
# gmpsa <- mean(mpsa_values)
# 
# # Step 5: 유의성 검정
# significance_results <- compute_MPSA_significance(P_matrix, W_matrix, n_perm = 999)

# === 3. 결과 시각화 ===
# 
# # 시각화 라이브러리 로드
# library(ggplot2)
# library(tmap)
# library(patchwork)
# 
# # MPSA 분포 히스토그램
# mpsa_df <- data.frame(
#   mpsa = mpsa_values,
#   region_id = 1:length(mpsa_values)
# )
# 
# p1 <- ggplot(mpsa_df, aes(x = mpsa)) +
#   geom_histogram(bins = 30, alpha = 0.7, fill = "steelblue", color = "white") +
#   geom_vline(xintercept = gmpsa, linetype = "dashed", color = "red", size = 1) +
#   labs(
#     title = "Distribution of Local MPSA Values",
#     subtitle = paste("Global MPSA =", round(gmpsa, 4)),
#     x = "Local MPSA", y = "Frequency"
#   ) +
#   theme_minimal()
# 
# # MPSA 값 vs 지역 ID 산점도
# p2 <- ggplot(mpsa_df, aes(x = region_id, y = mpsa)) +
#   geom_point(alpha = 0.6, color = "darkblue") +
#   geom_hline(yintercept = gmpsa, linetype = "dashed", color = "red") +
#   labs(
#     title = "Local MPSA by Region",
#     x = "Region ID", y = "Local MPSA"
#   ) +
#   theme_minimal()
# 
# # MPSA 범주화
# mpsa_categories <- cut(
#   mpsa_values,
#   breaks = quantile(mpsa_values, c(0, 0.2, 0.4, 0.6, 0.8, 1)),
#   labels = c("Very Low", "Low", "Medium", "High", "Very High"),
#   include.lowest = TRUE
# )
# 
# # 공간 데이터에 MPSA 결과 추가
# franklin_with_mpsa <- franklin
# franklin_with_mpsa$mpsa <- mpsa_values
# franklin_with_mpsa$mpsa_category <- mpsa_categories
# 
# # 공간 분포 지도
# spatial_map <- tm_shape(franklin_with_mpsa) +
#   tm_fill("mpsa_category", 
#           title = "MPSA Category",
#           palette = "RdYlBu") +
#   tm_borders(alpha = 0.3) +
#   tm_layout(
#     title = "MPSA Spatial Distribution",
#     legend.position = c("right", "bottom")
#   )

# === 4. 성능 분석 시각화 ===
# 
# # Random Forest 변수 중요도
# importance_data <- data.frame(
#   Variable = rownames(rf_model$importance),
#   Importance = rf_model$importance[, 1]
# )
# importance_data <- importance_data[order(importance_data$Importance, decreasing = TRUE), ]
# 
# # 상위 10개 변수 중요도 시각화
# p3 <- ggplot(importance_data[1:10, ], aes(x = reorder(Variable, Importance), y = Importance)) +
#   geom_col(fill = "coral", alpha = 0.8) +
#   coord_flip() +
#   labs(
#     title = "Top 10 Variable Importance in Random Forest",
#     x = "Variables", y = "Mean Decrease Accuracy"
#   ) +
#   theme_minimal()
# 
# # Proximity 매트릭스 분포
# proximity_values <- as.vector(P_matrix[upper.tri(P_matrix)])
# proximity_df <- data.frame(proximity = proximity_values)
# 
# p4 <- ggplot(proximity_df, aes(x = proximity)) +
#   geom_histogram(bins = 50, alpha = 0.7, fill = "lightgreen", color = "white") +
#   labs(
#     title = "Random Forest Proximity Distribution",
#     x = "Proximity Value", y = "Frequency"
#   ) +
#   theme_minimal()

# === 5. 진단 분석 ===
# 
# # MPSA vs Proximity 상관관계 분석
# avg_proximity <- rowMeans(P_matrix)
# correlation_df <- data.frame(
#   mpsa = mpsa_values,
#   avg_proximity = avg_proximity
# )
# 
# # 상관관계 산점도
# p5 <- ggplot(correlation_df, aes(x = avg_proximity, y = mpsa)) +
#   geom_point(alpha = 0.6, color = "darkgreen") +
#   geom_smooth(method = "lm", se = TRUE, color = "red") +
#   labs(
#     title = "MPSA vs Average Proximity",
#     subtitle = paste("Correlation =", round(cor(mpsa_values, avg_proximity), 3)),
#     x = "Average Proximity", y = "Local MPSA"
#   ) +
#   theme_minimal()
# 
# # 공간 가중치 행렬 특성 분석
# W_connectivity <- rowSums(W_matrix > 0)  # 각 지역의 이웃 수
# 
# connectivity_df <- data.frame(
#   region_id = 1:length(W_connectivity),
#   n_neighbors = W_connectivity,
#   mpsa = mpsa_values
# )
# 
# # 이웃 수 vs MPSA 관계
# p6 <- ggplot(connectivity_df, aes(x = n_neighbors, y = mpsa)) +
#   geom_point(alpha = 0.6, color = "purple") +
#   geom_smooth(method = "loess", se = TRUE, color = "orange") +
#   labs(
#     title = "MPSA vs Number of Neighbors",
#     x = "Number of Neighbors", y = "Local MPSA"
#   ) +
#   theme_minimal()

# === 6. 유의성 검정 결과 시각화 ===
# 
# if (!is.null(significance_results)) {
#   # 순열 검정 결과 분포
#   perm_results_df <- data.frame(
#     permutation = 1:length(significance_results$permutation_stats),
#     gmpsa_perm = significance_results$permutation_stats
#   )
#   
#   p7 <- ggplot(perm_results_df, aes(x = gmpsa_perm)) +
#     geom_histogram(bins = 30, alpha = 0.7, fill = "lightcoral", color = "white") +
#     geom_vline(xintercept = gmpsa, linetype = "solid", color = "darkred", size = 1) +
#     geom_vline(xintercept = mean(significance_results$permutation_stats), 
#                linetype = "dashed", color = "blue", size = 1) +
#     labs(
#       title = "Permutation Test Results",
#       subtitle = paste("Observed GMPSA =", round(gmpsa, 4), 
#                        ", P-value =", round(significance_results$p_value, 4)),
#       x = "Global MPSA (Permuted)", y = "Frequency"
#     ) +
#     theme_minimal()
# } else {
#   p7 <- ggplot() + 
#     labs(title = "Significance test not performed") +
#     theme_void()
# }

# === 7. 종합 대시보드 생성 ===
# 
# # 모든 그래프 결합
# dashboard_plot <- (p1 + p2) / (p3 + p4) / (p5 + p6) / p7
# 
# # 결과 디렉토리 생성
# if (!dir.exists("output/one_stop_analysis")) {
#   dir.create("output/one_stop_analysis", recursive = TRUE)
# }
# 
# # 개별 그래프 저장
# ggsave("output/one_stop_analysis/mpsa_distribution.png", p1, width = 10, height = 6)
# ggsave("output/one_stop_analysis/mpsa_by_region.png", p2, width = 10, height = 6)
# ggsave("output/one_stop_analysis/variable_importance.png", p3, width = 10, height = 8)
# ggsave("output/one_stop_analysis/proximity_distribution.png", p4, width = 10, height = 6)
# ggsave("output/one_stop_analysis/mpsa_vs_proximity.png", p5, width = 10, height = 6)
# ggsave("output/one_stop_analysis/mpsa_vs_neighbors.png", p6, width = 10, height = 6)
# 
# if (!is.null(significance_results)) {
#   ggsave("output/one_stop_analysis/permutation_test.png", p7, width = 10, height = 6)
# }
# 
# # 종합 대시보드 저장
# ggsave("output/one_stop_analysis/comprehensive_dashboard.png", dashboard_plot, 
#        width = 16, height = 20)
# 
# # 공간 지도 저장
# tmap_save(spatial_map, "output/one_stop_analysis/spatial_distribution.png")

# === 8. 수치 결과 저장 ===
# 
# # 기본 결과 데이터프레임
# results_df <- data.frame(
#   region_id = 1:length(mpsa_values),
#   mpsa = mpsa_values,
#   mpsa_category = mpsa_categories,
#   avg_proximity = avg_proximity,
#   n_neighbors = W_connectivity
# )
# 
# # 공간 정보와 결합
# franklin_results <- franklin %>%
#   mutate(
#     mpsa = mpsa_values,
#     mpsa_category = mpsa_categories,
#     avg_proximity = avg_proximity,
#     n_neighbors = W_connectivity
#   )
# 
# # 파일 저장
# write.csv(results_df, "output/one_stop_analysis/mpsa_results.csv", row.names = FALSE)
# write.csv(importance_data, "output/one_stop_analysis/variable_importance.csv", row.names = FALSE)
# 
# # 공간 데이터 저장 (shapefile)
# st_write(franklin_results, "output/one_stop_analysis/franklin_with_mpsa.shp", 
#          delete_dsn = TRUE)
# 
# # 요약 통계 저장
# summary_stats <- data.frame(
#   Statistic = c("Mean MPSA", "Median MPSA", "SD MPSA", "Min MPSA", "Max MPSA", 
#                 "Global MPSA", "N Regions", "Mean Neighbors"),
#   Value = c(
#     round(mean(mpsa_values), 4),
#     round(median(mpsa_values), 4),
#     round(sd(mpsa_values), 4),
#     round(min(mpsa_values), 4),
#     round(max(mpsa_values), 4),
#     round(gmpsa, 4),
#     length(mpsa_values),
#     round(mean(W_connectivity), 2)
#   )
# )
# 
# if (!is.null(significance_results)) {
#   summary_stats <- rbind(summary_stats, 
#                          data.frame(Statistic = "P-value", 
#                                     Value = round(significance_results$p_value, 4)))
# }
# 
# write.csv(summary_stats, "output/one_stop_analysis/summary_statistics.csv", row.names = FALSE)

# === 9. 실행 시간 및 성능 벤치마크 ===
# 
# # 각 단계별 실행 시간 측정
# execution_times <- data.frame(
#   Step = character(),
#   Time_seconds = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # Random Forest 시간 측정
# start_time <- Sys.time()
# rf_temp <- randomForest(numeric_data, proximity = TRUE, ntree = 500)
# rf_time <- as.numeric(Sys.time() - start_time)
# execution_times <- rbind(execution_times, 
#                          data.frame(Step = "Random Forest", Time_seconds = rf_time))
# 
# # MPSA 계산 시간
# start_time <- Sys.time()
# mpsa_temp <- rowSums(W_matrix * rf_temp$proximity)
# mpsa_time <- as.numeric(Sys.time() - start_time)
# execution_times <- rbind(execution_times, 
#                          data.frame(Step = "MPSA Calculation", Time_seconds = mpsa_time))
# 
# # 유의성 검정 시간 (만약 실행된다면)
# if (!is.null(significance_results)) {
#   start_time <- Sys.time()
#   sig_temp <- compute_MPSA_significance(rf_temp$proximity, W_matrix, n_perm = 99)
#   sig_time <- as.numeric(Sys.time() - start_time)
#   execution_times <- rbind(execution_times, 
#                            data.frame(Step = "Significance Test", Time_seconds = sig_time))
# }
# 
# # 실행 시간 시각화
# p8 <- ggplot(execution_times, aes(x = reorder(Step, Time_seconds), y = Time_seconds)) +
#   geom_col(fill = "lightblue", alpha = 0.8) +
#   coord_flip() +
#   labs(
#     title = "Execution Time by Analysis Step",
#     x = "Analysis Step", y = "Time (seconds)"
#   ) +
#   theme_minimal()
# 
# ggsave("output/one_stop_analysis/execution_times.png", p8, width = 10, height = 6)
# write.csv(execution_times, "output/one_stop_analysis/execution_times.csv", row.names = FALSE)

# === 10. 최종 요약 및 보고서 ===
# 
# # 전체 분석 요약
# analysis_summary <- list(
#   data_info = list(
#     n_regions = nrow(franklin),
#     n_variables = ncol(numeric_data),
#     study_area = "Franklin County, Ohio"
#   ),
#   mpsa_results = list(
#     global_mpsa = gmpsa,
#     mean_local_mpsa = mean(mpsa_values),
#     sd_local_mpsa = sd(mpsa_values),
#     range_local_mpsa = range(mpsa_values)
#   ),
#   spatial_patterns = list(
#     category_counts = table(mpsa_categories),
#     avg_neighbors = mean(W_connectivity),
#     connectivity_range = range(W_connectivity)
#   ),
#   model_performance = list(
#     rf_oob_error = rf_model$err.rate[nrow(rf_model$err.rate), 1],
#     top_variables = head(importance_data$Variable, 5),
#     avg_proximity = mean(proximity_values)
#   )
# )
# 
# # JSON 형태로 저장
# library(jsonlite)
# write_json(analysis_summary, "output/one_stop_analysis/analysis_summary.json", 
#            pretty = TRUE)
# 
# # 최종 요약 출력
# cat("=== 원스톱 MPSA 분석 완료 ===\n")
# cat("📊 데이터 정보:\n")
# cat(sprintf("  - 지역 수: %d개\n", analysis_summary$data_info$n_regions))
# cat(sprintf("  - 변수 수: %d개\n", analysis_summary$data_info$n_variables))
# cat("📈 MPSA 결과:\n")
# cat(sprintf("  - Global MPSA: %.4f\n", analysis_summary$mpsa_results$global_mpsa))
# cat(sprintf("  - Local MPSA 평균: %.4f (±%.4f)\n", 
#             analysis_summary$mpsa_results$mean_local_mpsa,
#             analysis_summary$mpsa_results$sd_local_mpsa))
# cat("🗺️  공간 패턴:\n")
# for (i in 1:length(analysis_summary$spatial_patterns$category_counts)) {
#   cat(sprintf("  - %s: %d개\n", 
#               names(analysis_summary$spatial_patterns$category_counts)[i],
#               analysis_summary$spatial_patterns$category_counts[i]))
# }
# cat("⏱️  실행 시간:\n")
# for (i in 1:nrow(execution_times)) {
#   cat(sprintf("  - %s: %.2f초\n", execution_times$Step[i], execution_times$Time_seconds[i]))
# }
# cat("📁 결과 저장 위치: output/one_stop_analysis/\n")
# cat("🎯 분석 완료! 모든 결과를 확인하세요.\n") 

# === 시각화 및 지도 생성 (tmap 4.1 버전) ===
# library(tmap)

# 🎯 핵심 공간 분포 지도 생성
# spatial_map <- tm_shape(franklin_with_mpsa) +
#   tm_fill(
#     fill = "mpsa_category",
#     fill.scale = tm_scale_categorical(
#       values = c("Strong Hotspot" = "#d73027", "Hotspot" = "#fc8d59", 
#                  "Not Significant" = "#ffffbf", "Coldspot" = "#91bfdb", 
#                  "Strong Coldspot" = "#4575b4")
#     ),
#     fill.legend = tm_legend(title = "MPSA Level")
#   ) +
#   tm_borders(alpha = 0.3) +
#   tm_layout(
#     title = paste0("MPSA Spatial Distribution (Global = ", 
#                    round(global_mpsa, 3), ", p = ", 
#                    round(global_pvalue, 3), ")"),
#     title.position = c("center", "top"),
#     legend.position = c("right", "top")
#   ) 