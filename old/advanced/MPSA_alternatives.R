# ============================================================================
# Key Alternative LIRF Structures
# 
# 논문에서 robustness 검증용으로 사용할 핵심 대안들
# ============================================================================

# --- 환경 설정 ---
source("R/data_preparation/setup.R")

# === 1. Adaptive LIRF (Moran-optimized) ===================================

#' Adaptive LIRF with Moran's I optimization
#' 
#' @description 공간 자기상관을 최대화하는 학습된 가중치를 사용
#' @param P proximity matrix
#' @param W spatial weight matrix
#' @param X original data (for optimization)
#' @return list with LIRF values and learned weights
compute_LIRF_adaptive <- function(P, W, X) {
  n <- nrow(P)
  
  # Moran's I 최대화 목적함수
  objective <- function(theta) {
    theta_matrix <- matrix(rep(theta, each = n), n, n)
    lirf_weighted <- rowSums(W * theta_matrix * P)
    
    # Moran's I 계산
    lirf_centered <- lirf_weighted - mean(lirf_weighted)
    moran_num <- sum(W * outer(lirf_centered, lirf_centered))
    moran_denom <- sum(lirf_centered^2)
    
    -moran_num / moran_denom  # 최소화 문제로 변환
  }
  
  # 최적화 실행
  init_theta <- rep(1, n)
  opt_result <- optim(init_theta, objective, method = "L-BFGS-B",
                     lower = 0.1, upper = 10)
  
  theta_opt <- opt_result$par
  
  # 학습된 가중치로 LIRF 계산
  theta_matrix <- matrix(rep(theta_opt, each = n), n, n)
  lirf_adaptive <- rowSums(W * theta_matrix * P)
  
  return(list(
    LIRF = lirf_adaptive,
    weights = theta_opt,
    optimization = opt_result
  ))
}

# === 2. Matrix Factorization LIRF ==========================================

#' Matrix Factorization LIRF
#' 
#' @description SVD로 노이즈를 제거하고 주요 패턴만 추출
#' @param P proximity matrix
#' @param W spatial weight matrix
#' @param k number of components to keep
#' @return list with LIRF values and reconstructed matrix
compute_LIRF_factorization <- function(P, W, k = 10) {
  n <- nrow(P)
  
  # SVD 분해
  svd_result <- svd(P)
  
  # 상위 k개 성분으로 재구성
  P_reconstructed <- svd_result$u[, 1:k] %*% 
                     diag(svd_result$d[1:k]) %*% 
                     t(svd_result$v[, 1:k])
  
  # 0-1 범위로 제한
  P_reconstructed <- pmax(0, pmin(1, P_reconstructed))
  
  # LIRF 계산
  lirf_mf <- rowSums(W * P_reconstructed)
  
  return(list(
    LIRF = lirf_mf,
    P_reconstructed = P_reconstructed,
    explained_variance = cumsum(svd_result$d^2) / sum(svd_result$d^2)
  ))
}

# === 3. 대안들과 기본 LIRF 비교 ============================================

#' Compare LIRF alternatives with basic LIRF
#' 
#' @param data sf object with spatial data
#' @return comparison results
compare_LIRF_methods <- function(data) {
  
  # 데이터 준비
  rf_data <- data |> st_drop_geometry() |> select(where(is.numeric))
  
  # 공간 가중치
  nb <- poly2nb(data, queen = TRUE)
  W <- nb2mat(nb, style = "W", zero.policy = TRUE)
  
  # Random Forest
  rf <- randomForest(rf_data, proximity = TRUE, ntree = 500)
  P <- rf$proximity
  
  cat("Computing LIRF alternatives...\n")
  
  # 1. 기본 LIRF (from LIRF.R)
  source("R/lirf_methods/LIRF.R")
  basic_lirf <- compute_LIRF(P, W)
  
  # 2. Adaptive LIRF
  adaptive_result <- compute_LIRF_adaptive(P, W, rf_data)
  
  # 3. Matrix Factorization LIRF
  mf_result <- compute_LIRF_factorization(P, W, k = 10)
  
  # 성능 비교
  results <- list(
    basic = basic_lirf,
    adaptive = adaptive_result$LIRF,
    factorization = mf_result$LIRF
  )
  
  # Moran's I 계산
  moran_stats <- sapply(results, function(lirf) {
    moran.test(lirf, mat2listw(W))$statistic
  })
  
  # 기본 통계
  performance <- data.frame(
    method = names(results),
    mean = sapply(results, mean),
    sd = sapply(results, sd),
    moran_I = moran_stats
  )
  
  cat("\nPerformance Comparison:\n")
  print(performance)
  
  return(list(
    results = results,
    performance = performance,
    adaptive_details = adaptive_result,
    factorization_details = mf_result
  ))
}

# === 4. 시각화 비교 ========================================================

#' Visualize LIRF method comparison
#' 
#' @param comparison_results results from compare_LIRF_methods
#' @param data original spatial data
visualize_LIRF_comparison <- function(comparison_results, data) {
  library(tmap)
  library(dplyr)
  
  # 데이터에 결과 추가
  data$LIRF_basic <- comparison_results$results$basic
  data$LIRF_adaptive <- comparison_results$results$adaptive
  data$LIRF_factorization <- comparison_results$results$factorization
  
  # 지도 생성
  tm1 <- tm_shape(data) +
    tm_fill(
      fill = "LIRF_basic",
      fill.scale = tm_scale_continuous(
        palette = "viridis"
      ),
      fill.legend = tm_legend(title = "LIRF (Basic)")
    ) +
    tm_borders(alpha = 0.3) +
    tm_layout(title = "Basic LIRF", title.position = c("center", "top"))
  
  tm2 <- tm_shape(data) +
    tm_fill(
      fill = "LIRF_adaptive",
      fill.scale = tm_scale_continuous(
        palette = "plasma"
      ),
      fill.legend = tm_legend(title = "LIRF (Adaptive)")
    ) +
    tm_borders(alpha = 0.3) +
    tm_layout(title = "Adaptive LIRF (Moran-optimized)", title.position = c("center", "top"))
  
  tm3 <- tm_shape(data) +
    tm_fill(
      fill = "LIRF_factorization",
      fill.scale = tm_scale_continuous(
        palette = "cividis"
      ),
      fill.legend = tm_legend(title = "LIRF (Matrix)")
    ) +
    tm_borders(alpha = 0.3) +
    tm_layout(title = "Matrix Factorization LIRF", title.position = c("center", "top"))
  
  return(tmap_arrange(tm1, tm2, tm3, ncol = 3))
}

# === 실행 예시 및 시각화 =====================================================

# === 1. 대안적 MPSA 방법론 실행 ===
# 
# # 데이터 로드
# franklin <- readRDS("data/franklin.rds")
# numeric_data <- franklin %>% 
#   st_drop_geometry() %>% 
#   select(where(is.numeric)) %>%
#   select(-main_industry)
# 
# # 공간 가중치 행렬 생성
# coords <- st_coordinates(st_centroid(franklin))
# nb <- poly2nb(franklin, queen = TRUE)
# W_matrix <- nb2mat(nb, style = "W", zero.policy = TRUE)
# 
# # 다중스케일 MPSA 실행
# multiscale_results <- compute_multiscale_MPSA(numeric_data, W_matrix)
# print("다중스케일 MPSA 완료")
# 
# # 시간적 MPSA 실행 (시뮬레이션 데이터 사용)
# temporal_data <- generate_temporal_data(numeric_data, n_time_periods = 5)
# temporal_results <- compute_temporal_MPSA(temporal_data, W_matrix)
# print("시간적 MPSA 완료")

# === 2. 방법론별 결과 비교 시각화 ===
# 
# # 기본 MPSA 계산
# library(randomForest)
# rf_basic <- randomForest(numeric_data, proximity = TRUE, ntree = 500)
# mpsa_basic <- rowSums(W_matrix * rf_basic$proximity)
# 
# # 비교 데이터 준비
# comparison_df <- data.frame(
#   region_id = 1:length(mpsa_basic),
#   basic_mpsa = mpsa_basic,
#   local_scale = multiscale_results$local_scale,
#   regional_scale = multiscale_results$regional_scale,
#   global_scale = multiscale_results$global_scale,
#   temporal_trend = temporal_results$temporal_trend,
#   temporal_variability = temporal_results$temporal_variability
# )
# 
# # 분포 비교 시각화
# library(ggplot2)
# library(tidyr)
# 
# methods_long <- comparison_df %>%
#   gather(method, value, -region_id) %>%
#   mutate(method = factor(method, 
#                         levels = c("basic_mpsa", "local_scale", "regional_scale", 
#                                   "global_scale", "temporal_trend", "temporal_variability"),
#                         labels = c("Basic MPSA", "Local Scale", "Regional Scale", 
#                                   "Global Scale", "Temporal Trend", "Temporal Variability")))
# 
# p1 <- ggplot(methods_long, aes(x = value, fill = method)) +
#   geom_histogram(bins = 25, alpha = 0.7) +
#   facet_wrap(~method, scales = "free", ncol = 2) +
#   labs(title = "Distribution Comparison of Alternative MPSA Methods",
#        x = "Value", y = "Frequency") +
#   theme_minimal() +
#   theme(legend.position = "none")

# === 3. 다중스케일 분석 시각화 ===
# 
# # 스케일별 공간 패턴
# franklin_multiscale <- franklin
# franklin_multiscale$local_scale <- multiscale_results$local_scale
# franklin_multiscale$regional_scale <- multiscale_results$regional_scale
# franklin_multiscale$global_scale <- multiscale_results$global_scale
# 
# # 스케일별 범주화
# franklin_multiscale$local_category <- cut(multiscale_results$local_scale, 
#                                           breaks = 5, labels = c("VL", "L", "M", "H", "VH"))
# franklin_multiscale$regional_category <- cut(multiscale_results$regional_scale, 
#                                              breaks = 5, labels = c("VL", "L", "M", "H", "VH"))
# franklin_multiscale$global_category <- cut(multiscale_results$global_scale, 
#                                            breaks = 5, labels = c("VL", "L", "M", "H", "VH"))
# 
# # === 다중스케일 MPSA 시각화 (tmap 4.1 버전) ===
# # library(tmap)
# 
# # 지역 스케일 지도
# # map_local <- tm_shape(franklin_multiscale) +
# #   tm_fill(
# #     fill = "local_category", 
# #     fill.scale = tm_scale_categorical(
# #       values = "Reds"
# #     ),
# #     fill.legend = tm_legend(title = "Local Scale")
# #   ) +
# #   tm_borders(alpha = 0.3) +
# #   tm_layout(title = "Local Scale MPSA", title.position = c("center", "top"))
# 
# # 지역 스케일 지도  
# # map_regional <- tm_shape(franklin_multiscale) +
# #   tm_fill(
# #     fill = "regional_category", 
# #     fill.scale = tm_scale_categorical(
# #       values = "Blues"
# #     ),
# #     fill.legend = tm_legend(title = "Regional Scale")
# #   ) +
# #   tm_borders(alpha = 0.3) +
# #   tm_layout(title = "Regional Scale MPSA", title.position = c("center", "top"))
# 
# # 전역 스케일 지도
# # map_global <- tm_shape(franklin_multiscale) +
# #   tm_fill(
# #     fill = "global_category", 
# #     fill.scale = tm_scale_categorical(
# #       values = "Greens"
# #     ),
# #     fill.legend = tm_legend(title = "Global Scale")
# #   ) +
# #   tm_borders(alpha = 0.3) +
# #   tm_layout(title = "Global Scale MPSA", title.position = c("center", "top"))
# 
# # 다중스케일 지도 결합
# # multiscale_maps <- tmap_arrange(map_local, map_regional, map_global, ncol = 3)

# === 4. 시간적 변화 패턴 분석 ===
# 
# # 시간적 트렌드 시각화
# p2 <- ggplot(comparison_df, aes(x = region_id, y = temporal_trend)) +
#   geom_line(color = "darkblue", alpha = 0.7) +
#   geom_smooth(method = "loess", se = TRUE, color = "red") +
#   labs(title = "Temporal Trend in MPSA",
#        subtitle = "Positive values indicate increasing spatial autocorrelation",
#        x = "Region ID", y = "Temporal Trend") +
#   theme_minimal()
# 
# # 시간적 변동성 시각화
# p3 <- ggplot(comparison_df, aes(x = region_id, y = temporal_variability)) +
#   geom_point(alpha = 0.6, color = "darkgreen") +
#   geom_smooth(method = "loess", se = TRUE, color = "orange") +
#   labs(title = "Temporal Variability in MPSA",
#        subtitle = "Higher values indicate more temporal instability",
#        x = "Region ID", y = "Temporal Variability") +
#   theme_minimal()
# 
# # 트렌드 vs 변동성 산점도
# p4 <- ggplot(comparison_df, aes(x = temporal_trend, y = temporal_variability)) +
#   geom_point(alpha = 0.6, color = "purple") +
#   geom_smooth(method = "lm", se = TRUE, color = "red") +
#   labs(title = "Temporal Trend vs Variability",
#        subtitle = paste("Correlation:", 
#                         round(cor(comparison_df$temporal_trend, 
#                                  comparison_df$temporal_variability), 3)),
#        x = "Temporal Trend", y = "Temporal Variability") +
#   theme_minimal()

# === 5. 스케일 간 상관관계 분석 ===
# 
# # 스케일별 상관관계 매트릭스
# scale_data <- comparison_df[, c("basic_mpsa", "local_scale", "regional_scale", "global_scale")]
# scale_cor <- cor(scale_data, use = "complete.obs")
# 
# # 상관관계 히트맵
# library(corrplot)
# png("output/alternative_mpsa/scale_correlation.png", width = 600, height = 600)
# corrplot(scale_cor, method = "color", type = "upper", 
#          order = "hclust", tl.cex = 1, tl.col = "black",
#          title = "Correlation Between Different Scales")
# dev.off()
# 
# # 스케일별 산점도 매트릭스
# library(GGally)
# p5 <- ggpairs(scale_data, 
#               columnLabels = c("Basic", "Local", "Regional", "Global"),
#               title = "Pairwise Comparisons of Different Scales") +
#   theme_minimal()

# === 6. 공간 계층구조 분석 ===
# 
# # 계층별 분산 분해
# hierarchical_analysis <- analyze_spatial_hierarchy(multiscale_results)
# 
# # 분산 분해 시각화
# variance_df <- data.frame(
#   Scale = c("Local", "Regional", "Global", "Residual"),
#   Variance_Explained = hierarchical_analysis$variance_explained,
#   Percentage = round(100 * hierarchical_analysis$variance_explained / 
#                     sum(hierarchical_analysis$variance_explained), 1)
# )
# 
# p6 <- ggplot(variance_df, aes(x = reorder(Scale, Variance_Explained), y = Variance_Explained)) +
#   geom_col(fill = "steelblue", alpha = 0.8) +
#   geom_text(aes(label = paste0(Percentage, "%")), 
#             hjust = -0.1, size = 3) +
#   coord_flip() +
#   labs(title = "Variance Decomposition by Spatial Scale",
#        x = "Spatial Scale", y = "Variance Explained") +
#   theme_minimal()

# === 7. 시간적 패턴의 공간 분포 ===
# 
# # 시간적 결과를 공간 데이터에 추가
# franklin_temporal <- franklin
# franklin_temporal$temporal_trend <- temporal_results$temporal_trend
# franklin_temporal$temporal_variability <- temporal_results$temporal_variability
# 
# # 시간적 패턴 범주화
# franklin_temporal$trend_category <- cut(temporal_results$temporal_trend,
#                                        breaks = c(-Inf, -0.01, 0.01, Inf),
#                                        labels = c("Decreasing", "Stable", "Increasing"))
# 
# franklin_temporal$variability_category <- cut(temporal_results$temporal_variability,
#                                              breaks = quantile(temporal_results$temporal_variability, c(0, 0.33, 0.67, 1)),
#                                              labels = c("Low", "Medium", "High"),
#                                              include.lowest = TRUE)
# 
# === 시공간 트렌드 시각화 (tmap 4.1 버전) ===
# trend_map <- tm_shape(franklin_temporal) +
#   tm_fill(
#     fill = "trend_category",
#     fill.scale = tm_scale_categorical(
#       values = c("Increasing" = "#d73027", "Stable" = "#ffffbf", 
#                  "Decreasing" = "#4575b4")
#     ),
#     fill.legend = tm_legend(title = "MPSA Trend")
#   ) +
#   tm_borders(alpha = 0.3) +
#   tm_layout(
#     title = "Spatiotemporal MPSA Trends",
#     title.position = c("center", "top")
#   )

# === 8. 방법론별 성능 비교 ===
# 
# # 각 방법론의 탐지 능력 평가
# detection_performance <- data.frame(
#   Method = c("Basic MPSA", "Local Scale", "Regional Scale", "Global Scale", 
#              "Temporal Trend", "Temporal Variability"),
#   Mean_Value = c(
#     mean(mpsa_basic),
#     mean(multiscale_results$local_scale),
#     mean(multiscale_results$regional_scale),
#     mean(multiscale_results$global_scale),
#     mean(temporal_results$temporal_trend),
#     mean(temporal_results$temporal_variability)
#   ),
#   SD_Value = c(
#     sd(mpsa_basic),
#     sd(multiscale_results$local_scale),
#     sd(multiscale_results$regional_scale),
#     sd(multiscale_results$global_scale),
#     sd(temporal_results$temporal_trend),
#     sd(temporal_results$temporal_variability)
#   ),
#   Detection_Range = c(
#     diff(range(mpsa_basic)),
#     diff(range(multiscale_results$local_scale)),
#     diff(range(multiscale_results$regional_scale)),
#     diff(range(multiscale_results$global_scale)),
#     diff(range(temporal_results$temporal_trend)),
#     diff(range(temporal_results$temporal_variability))
#   )
# )
# 
# # 성능 비교 시각화
# p7 <- ggplot(detection_performance, aes(x = reorder(Method, Detection_Range), y = Detection_Range)) +
#   geom_col(fill = "lightcoral", alpha = 0.8) +
#   coord_flip() +
#   labs(title = "Detection Range by Method",
#        subtitle = "Larger range indicates better discriminatory power",
#        x = "Method", y = "Detection Range") +
#   theme_minimal()

# === 9. 결과 저장 및 종합 분석 ===
# 
# # 결과 디렉토리 생성
# if (!dir.exists("output/alternative_mpsa")) {
#   dir.create("output/alternative_mpsa", recursive = TRUE)
# }
# 
# # 모든 시각화 저장
# library(patchwork)
# combined_temporal <- (p2 + p3) / p4
# combined_performance <- (p6 + p7)
# 
# ggsave("output/alternative_mpsa/method_distributions.png", p1, width = 12, height = 10)
# ggsave("output/alternative_mpsa/temporal_patterns.png", combined_temporal, width = 16, height = 10)
# ggsave("output/alternative_mpsa/pairwise_comparisons.png", p5, width = 12, height = 10)
# ggsave("output/alternative_mpsa/performance_analysis.png", combined_performance, width = 16, height = 8)
# 
# # 지도 저장
# tmap_save(multiscale_maps, "output/alternative_mpsa/multiscale_spatial.png")
# tmap_save(temporal_maps, "output/alternative_mpsa/temporal_spatial.png")
# 
# # 수치 결과 저장
# write.csv(comparison_df, "output/alternative_mpsa/comprehensive_comparison.csv", row.names = FALSE)
# write.csv(variance_df, "output/alternative_mpsa/variance_decomposition.csv", row.names = FALSE)
# write.csv(detection_performance, "output/alternative_mpsa/performance_metrics.csv", row.names = FALSE)
# 
# # 공간 데이터 저장
# st_write(franklin_multiscale, "output/alternative_mpsa/multiscale_results.shp", delete_dsn = TRUE)
# st_write(franklin_temporal, "output/alternative_mpsa/temporal_results.shp", delete_dsn = TRUE)

# === 10. 방법론 특성 요약 및 권고사항 ===
# 
# # 각 방법론의 특성 요약
# method_characteristics <- data.frame(
#   Method = c("Basic MPSA", "Multiscale Local", "Multiscale Regional", 
#              "Multiscale Global", "Temporal Trend", "Temporal Variability"),
#   Spatial_Focus = c("Single Scale", "Neighborhood", "Community", "Overall", "Change", "Stability"),
#   Temporal_Aspect = c("Static", "Static", "Static", "Static", "Dynamic", "Dynamic"),
#   Computational_Cost = c("Medium", "High", "Very High", "Very High", "High", "High"),
#   Best_Use_Case = c(
#     "General spatial autocorrelation",
#     "Fine-scale hotspot detection",
#     "Community-level pattern analysis", 
#     "Overall spatial structure",
#     "Monitoring spatial change",
#     "Identifying unstable regions"
#   ),
#   Interpretation = c(
#     "Standard proximity-based autocorrelation",
#     "Local neighborhood similarities",
#     "Broader regional clustering",
#     "Global spatial organization",
#     "Direction of spatial change",
#     "Consistency of spatial patterns"
#   )
# )
# 
# # 상관관계 요약
# all_correlations <- cor(comparison_df[, -1], use = "complete.obs")
# correlation_summary <- data.frame(
#   Method_Pair = c("Basic-Local", "Basic-Regional", "Basic-Global", 
#                   "Local-Regional", "Regional-Global", "Trend-Variability"),
#   Correlation = c(
#     all_correlations["basic_mpsa", "local_scale"],
#     all_correlations["basic_mpsa", "regional_scale"],
#     all_correlations["basic_mpsa", "global_scale"],
#     all_correlations["local_scale", "regional_scale"],
#     all_correlations["regional_scale", "global_scale"],
#     all_correlations["temporal_trend", "temporal_variability"]
#   )
# )
# 
# # 종합 평가 저장
# library(jsonlite)
# comprehensive_evaluation <- list(
#   method_characteristics = method_characteristics,
#   correlation_summary = correlation_summary,
#   performance_metrics = detection_performance,
#   variance_decomposition = variance_df
# )
# 
# write_json(comprehensive_evaluation, "output/alternative_mpsa/comprehensive_evaluation.json", pretty = TRUE)
# write.csv(method_characteristics, "output/alternative_mpsa/method_characteristics.csv", row.names = FALSE)
# write.csv(correlation_summary, "output/alternative_mpsa/correlation_summary.csv", row.names = FALSE)
# 
# # 최종 권고사항 생성
# recommendations <- data.frame(
#   Research_Question = c(
#     "General spatial pattern detection",
#     "Fine-scale hotspot identification", 
#     "Regional clustering analysis",
#     "Global spatial organization",
#     "Monitoring spatial change over time",
#     "Identifying spatially unstable regions"
#   ),
#   Recommended_Method = c(
#     "Basic MPSA",
#     "Multiscale Local MPSA",
#     "Multiscale Regional MPSA", 
#     "Multiscale Global MPSA",
#     "Temporal Trend MPSA",
#     "Temporal Variability MPSA"
#   ),
#   Expected_Outcome = c(
#     "Overall spatial autocorrelation pattern",
#     "Precise local hotspot boundaries",
#     "Community-level spatial clusters",
#     "System-wide spatial structure",
#     "Trend direction and magnitude",
#     "Temporal stability assessment"
#   )
# )
# 
# write.csv(recommendations, "output/alternative_mpsa/method_recommendations.csv", row.names = FALSE)
# 
# # 최종 요약 출력
# cat("=== 대안적 MPSA 방법론 분석 완료 ===\n")
# cat("📊 방법론별 특성:\n")
# for (i in 1:nrow(method_characteristics)) {
#   cat(sprintf("  - %s: %s\n", 
#               method_characteristics$Method[i], 
#               method_characteristics$Best_Use_Case[i]))
# }
# cat("🔗 주요 상관관계:\n")
# for (i in 1:nrow(correlation_summary)) {
#   cat(sprintf("  - %s: %.3f\n", 
#               correlation_summary$Method_Pair[i], 
#               correlation_summary$Correlation[i]))
# }
# cat("📈 성능 지표:\n")
# cat(sprintf("  - 최고 탐지 범위: %s (%.3f)\n", 
#             detection_performance$Method[which.max(detection_performance$Detection_Range)],
#             max(detection_performance$Detection_Range)))
# cat("📁 결과 저장 위치: output/alternative_mpsa/\n")
# cat("💡 권고사항: method_recommendations.csv 참조\n")
# cat("🎉 대안적 방법론 분석 완료!\n") 