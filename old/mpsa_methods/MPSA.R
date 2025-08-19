# ============================================================================
# Multivariate Proximity-based Spatial Autocorrelation (MPSA)
# 
# 기본 MPSA 구현 - 논문의 주요 방법론
# ============================================================================

# --- 환경 설정 ---
source("R/data_preparation/setup.R")

franklin <- readRDS("data/franklin.rds")
x <- franklin |> 
  st_drop_geometry() |> 
  dplyr::select(where(is.numeric))

synthetic <- as.data.frame(lapply(x, function(col) sample(col, replace = TRUE)))
data_all <- rbind(x, synthetic)
y_all <- factor(c(rep("real", nrow(x)), rep("synthetic", nrow(synthetic))))

rf <- randomForest(x = data_all,
                   y = y_all,
                   proximity = TRUE,
                   ntree = 500)
P <- rf$proximity
P <- P[1:327, 1:327]

nb <- poly2nb(franklin, queen = TRUE)
W <- nb2mat(nb, style = "W", zero.policy = TRUE)

# === 1. 기본 MPSA 계산 ======================================================

#' Compute MPSA (Multivariate Proximity-based Spatial Autocorrelation)
#' 
#' @description MPSA_i = Σ_j W_ij * P_ij
#' @param P proximity matrix from random forest
#' @param W spatial weight matrix (row-standardized)
#' @return MPSA values
compute_MPSA <- function(P, W) {
  rowSums(W * P)
}

#' Compute MPSA with significance testing
#' 
#' @param P proximity matrix
#' @param W spatial weight matrix
#' @param n_perm number of permutations for significance testing
#' @param alpha significance level
#' @param method permutation method ("random", "spatial", "conditional")
#' @return list with MPSA values, p-values, and categories
compute_MPSA_significance <- function(P, W, n_perm = 999, alpha = 0.05, 
                                    method = "random") {
  n <- nrow(P)
  observed_mpsa <- compute_MPSA(P, W)
  
  # 개선된 순열 검정
  perm_matrix <- matrix(NA, nrow = n, ncol = n_perm)
  
  for (b in 1:n_perm) {
    if (method == "random") {
      # 완전 무작위 순열
      perm_idx <- sample(1:n)
      P_perm <- P[perm_idx, perm_idx]
      
    } else if (method == "spatial") {
      # 공간 구조 보존 순열 (제한된 순열)
      # 이웃들끼리만 위치 교환
      perm_idx <- 1:n
      for (i in 1:n) {
        neighbors <- which(W[i,] > 0)
        if (length(neighbors) > 1) {
          # 확률적으로 이웃과 위치 교환
          if (runif(1) < 0.3) {
            swap_with <- sample(neighbors, 1)
            temp <- perm_idx[i]
            perm_idx[i] <- perm_idx[swap_with]
            perm_idx[swap_with] <- temp
          }
        }
      }
      P_perm <- P[perm_idx, perm_idx]
      
    } else if (method == "conditional") {
      # 조건부 순열 (RF proximity의 분포 보존)
      # Bootstrap에 기반한 proximity 재샘플링
      P_perm <- P[sample(1:n, replace = TRUE), sample(1:n, replace = TRUE)]
    }
    
    perm_matrix[, b] <- compute_MPSA(P_perm, W)
  }
  
  # 개선된 p값 계산
  p_values <- sapply(1:n, function(i) {
    obs <- observed_mpsa[i]
    null_dist <- perm_matrix[i, ]
    
    # 정확한 p-value (양측검정)
    p_greater <- mean(null_dist >= obs)
    p_less <- mean(null_dist <= obs)
    
    # 더 conservative한 양측검정
    p_two_sided <- 2 * min(p_greater, p_less)
    return(p_two_sided)
  })
  
  # 다중 검정 보정 (FDR)
  p_adjusted <- p.adjust(p_values, method = "fdr")
  
  # 효과 크기 기반 분류
  expected_mpsa <- apply(perm_matrix, 1, mean)
  sd_mpsa <- apply(perm_matrix, 1, sd)
  effect_size <- (observed_mpsa - expected_mpsa) / sd_mpsa
  
  # 5단계 분류 체계
  category <- cut(effect_size,
                 breaks = c(-Inf, -2, -1, 1, 2, Inf),
                 labels = c("Strong Coldspot", "Coldspot", "Not Significant", 
                           "Hotspot", "Strong Hotspot"),
                 include.lowest = TRUE)
  
  # 유의성과 효과 크기 모두 고려
  final_category <- as.character(category)
  final_category[p_adjusted >= alpha] <- "Not Significant"
  
  return(list(
    MPSA = observed_mpsa,
    p_raw = p_values,
    p_adjusted = p_adjusted,
    effect_size = effect_size,
    category = factor(final_category, 
                     levels = c("Strong Hotspot", "Hotspot", "Not Significant", 
                               "Coldspot", "Strong Coldspot")),
    expected = expected_mpsa,
    sd_null = sd_mpsa,
    method = method
  ))
}


# === 2. GMPSA (Global Multivariate Proximity-based Spatial Autocorrelation) ===

#' Compute GMPSA - Global spatial autocorrelation measure
#' 
#' @description GMPSA following Anselin's LISA conditions: mean of all Local MPSA values
#' @param P proximity matrix
#' @param W spatial weight matrix (row-standardized)
#' @param n_perm number of permutations
#' @return GMPSA statistic and p-value
#' @details GMPSA = (1/n) * Σᵢ MPSA_i = (1/n) * Σᵢ Σⱼ W_ij * P_ij
#'          This satisfies Anselin's LISA condition: Σᵢ MPSA_i = n * GMPSA
compute_GMPSA <- function(P, W, n_perm = 999) {
  n <- nrow(P)
  
  # Calculate Local MPSA values
  local_mpsa <- rowSums(W * P)
  
  # GMPSA = mean of all Local MPSA values (satisfies LISA condition)
  GMPSA <- mean(local_mpsa)
  
  # Alternative: standardized version similar to Moran's I
  W_sum <- sum(W)
  P_mean <- mean(P)
  P_centered <- P - P_mean
  
  GMPSA_standardized <- (n / W_sum) * sum(W * P_centered) / sum(P_centered^2)
  
  # Permutation testing for significance
  GMPSA_perm <- numeric(n_perm)
  GMPSA_std_perm <- numeric(n_perm)
  
  for (i in 1:n_perm) {
    perm_idx <- sample(1:n)
    P_perm <- P[perm_idx, perm_idx]
    local_mpsa_perm <- rowSums(W * P_perm)
    
    # Basic GMPSA
    GMPSA_perm[i] <- mean(local_mpsa_perm)
    
    # Standardized GMPSA
    P_perm_mean <- mean(P_perm)
    P_perm_centered <- P_perm - P_perm_mean
    GMPSA_std_perm[i] <- (n / W_sum) * sum(W * P_perm_centered) / sum(P_perm_centered^2)
  }
  
  # p-values
  p_value_basic <- (sum(abs(GMPSA_perm) >= abs(GMPSA)) + 1) / (n_perm + 1)
  p_value_std <- (sum(abs(GMPSA_std_perm) >= abs(GMPSA_standardized)) + 1) / (n_perm + 1)
  
  # Verification of LISA condition
  lisa_sum <- sum(local_mpsa)
  lisa_expected <- n * GMPSA
  lisa_condition_satisfied <- abs(lisa_sum - lisa_expected) < 1e-10
  
  return(list(
    GMPSA = GMPSA,
    GMPSA_standardized = GMPSA_standardized,
    GMPSA_perm = GMPSA_perm,
    GMPSA_std_perm = GMPSA_std_perm,
    p_value = p_value_basic,
    p_value_standardized = p_value_std,
    expected = mean(GMPSA_perm),
    variance = var(GMPSA_perm),
    z_score = (GMPSA - mean(GMPSA_perm)) / sd(GMPSA_perm),
    local_mpsa = local_mpsa,
    lisa_condition_satisfied = lisa_condition_satisfied,
    lisa_verification = list(
      sum_local = lisa_sum,
      n_times_global = lisa_expected,
      difference = abs(lisa_sum - lisa_expected)
    )
  ))
}

# === 3. 기본 MPSA 분석 함수 ===============================================

#' Run basic MPSA analysis
#' 
#' @param data sf object with spatial data
#' @param ntree number of trees for random forest
#' @param n_perm number of permutations
#' @return basic MPSA analysis results
run_basic_MPSA_analysis <- function(data, ntree = 500, n_perm = 999) {
  
  cat("=== MPSA Analysis ===\n")
  cat(sprintf("Data: %d spatial units\n", nrow(data)))
  
  # 1. 데이터 준비
  rf_data <- data |> 
    st_drop_geometry() |> 
    dplyr::select(where(is.numeric), "main_industry")
  
  cat(sprintf("Variables: %d numeric variables\n", ncol(rf_data)))
  
  # 2. Random Forest 및 Proximity 계산
  cat("\nComputing Random Forest proximity matrix...\n")
  
  set.seed(42)
  synthetic <- as.data.frame(lapply(rf_data, function(col) sample(col, replace = TRUE)))
  data_all <- rbind(rf_data, synthetic)
  y_all <- factor(c(rep("real", nrow(rf_data)), rep("synthetic", nrow(synthetic))))
  
  rf_model <- randomForest(
    x = data_all,
    y = y_all,
    proximity = TRUE,
    ntree = ntree,
    importance = TRUE
  )
  
  P <- rf_model$proximity[1:nrow(rf_data), 1:nrow(rf_data)]
  
  # 3. 공간 가중치 행렬 구성
  cat("Creating spatial weight matrix...\n")
  nb <- poly2nb(data, queen = TRUE)
  W <- nb2listw(nb, style = "W", zero.policy = TRUE)
  W_matrix <- listw2mat(W)
  
  # 4. Global MPSA 계산
  cat("Computing Global MPSA...\n")
  global_result <- compute_GMPSA(P, W_matrix, n_perm = n_perm)
  
  # 별도로 GMPSA permutation 결과 추출
  global_perm_result <- list(
    observed = global_result$GMPSA,
    permuted = global_result$GMPSA_perm,
    p_value = global_result$p_value,
    expected = global_result$expected,
    variance = global_result$variance,
    z_score = global_result$z_score
  )
  
  cat(sprintf("Global MPSA: %.4f (p = %.3f)\n", 
              global_result$GMPSA, global_result$p_value))
  
  # 5. Local MPSA 계산
  cat("\nComputing Local MPSA with significance testing...\n")
  local_result <- compute_MPSA_significance(P, W_matrix, n_perm = n_perm)
  
  # 6. 결과 요약
  category_summary <- table(local_result$category)
  cat("\nMPSA Classification Summary:\n")
  print(category_summary)
  
  # 7. 결과를 원 데이터에 결합
  results_data <- data |> 
    mutate(
      MPSA = local_result$MPSA,
      MPSA_pvalue = local_result$p_adjusted,
      MPSA_category = factor(local_result$category, levels = c("Strong Hotspot", 
                                                               "Hotspot",
                                                               "Not Significant",
                                                               "Coldspot",
                                                               "Strong Coldspot")),
      MPSA_effect_size = local_result$effect_size,
      MPSA_expected = local_result$expected
    )
  

  
  # 8. 기본 시각화 생성 (tmap 4.1 버전)
  map_plot <- tm_shape(results_data) +
    tm_polygons(
      fill = "MPSA_category",
      fill.scale = tm_scale_categorical(values = c(
        "Strong Hotspot" = "#d73027",
        "Hotspot" = "#fc8d59",
        "Not Significant" = "#ffffbf", 
        "Coldspot" = "#91bfdb",
        "Strong Coldspot" = "#4575b4"
        )
      ),
      col = "black",
      fill.legend = tm_legend(title = "MPSA Category",
                              position = c(0.85, 0.25),
                              bg.color = "white"
                              )
    ) +
    tm_title("Local MPSA Results") +
    tm_compass(type = "arrow",
               position = c("left", "bottom")) +
    tm_scalebar(position = c("left", "bottom"),
                bg.color = "white")

  
  # 9. 반환 결과
  return(list(
    data = results_data,
    global = global_result,
    local = local_result,
    global_perm_result = global_perm_result,
    proximity = P,
    weights = W_matrix,
    rf_model = rf_model,
    plot = map_plot,
    summary = list(
      n_units = nrow(data),
      n_variables = ncol(rf_data),
      global_mpsa = global_result$GMPSA,
      global_pvalue = global_result$p_value,
      category_counts = category_summary
    )
  ))
}

# === 4. 전통적 방법론과의 벤치마크 비교 (공통 함수) ======================

#' Benchmark Against Traditional Methods
#' 
#' @description 전통적 방법들과의 성능 비교를 위한 공통 함수
#' @param data 분석 데이터 (numeric variables only)
#' @param W 공간 가중치 행렬
#' @return 벤치마크 비교 결과
benchmark_against_traditional_methods <- function(data, W, style = "W") {
  
  n <- nrow(data)
  p <- ncol(data)
  
  benchmark_results <- list()
  
  # 1. MPSA
  rf <- randomForest(data, proximity = TRUE, ntree = 500)
  P <- rf$proximity
  mpsa_values <- rowSums(W * P)
  gmpsa <- mean(mpsa_values)
  
  # 2. PCA + Moran's I
  pca_result <- prcomp(data, scale. = TRUE)
  pc1_moran <- tryCatch({
    W_listw <- mat2listw(W)
    moran.test(pca_result$x[, 1], W_listw, zero.policy = TRUE)
  }, error = function(e) list(statistic = NA, p.value = NA))
  
  # 3. 개별 변수 Moran's I
  individual_morans <- data.frame(
    variable = colnames(data),
    moran_statistic = NA,
    moran_pvalue = NA
  )
  
  for (j in 1:p) {
    moran_result <- tryCatch({
      W_listw <- mat2listw(W)
      moran.test(data[, j], W_listw, zero.policy = TRUE)
    }, error = function(e) list(statistic = NA, p.value = NA))
    
    individual_morans$moran_statistic[j] <- moran_result$statistic
    individual_morans$moran_pvalue[j] <- moran_result$p.value
  }
  
  # 4. 다변량 거리 기반 방법
  # Euclidean distance-based spatial autocorrelation
  euclidean_dist <- as.matrix(dist(data))
  similarity_euclidean <- 1 / (1 + euclidean_dist)  # Convert to similarity
  diag(similarity_euclidean) <- 1
  mpsa_euclidean <- rowSums(W * similarity_euclidean)
  gmpsa_euclidean <- mean(mpsa_euclidean)
  
  # 5. 🆕 Anselin(2019) LIMSA - Local Indicator of Multivariate Spatial Association
  # Extending Geary's c to multivariate case
  limsa_results <- compute_LIMSA(data, W)
  
  # 6. 🆕 이몽현(2012) 마할라노비스 거리 기반 방법
  # Mahalanobis distance-based spatial autocorrelation with chi-square test
  lee_2012_results <- compute_Lee2012_Mahalanobis(data, W)
  
  benchmark_results <- list(
    MPSA = list(
      global = gmpsa,
      local = mpsa_values,
      method = "Random Forest Proximity"
    ),
    PCA_Moran = list(
      global = pc1_moran$statistic,
      p_value = pc1_moran$p.value,
      variance_explained = summary(pca_result)$importance[2, 1],
      method = "PCA + Moran's I"
    ),
    Individual_Moran = list(
      results = individual_morans,
      mean_statistic = mean(individual_morans$moran_statistic, na.rm = TRUE),
      method = "Individual Variable Moran's I"
    ),
    Euclidean_based = list(
      global = gmpsa_euclidean,
      local = mpsa_euclidean,
      method = "Euclidean Distance-based"
    ),
    LIMSA_Anselin2019 = list(
      global = limsa_results$global_limsa,
      local = limsa_results$local_limsa,
      p_value = limsa_results$p_value,
      method = "Anselin(2019) LIMSA - Multivariate Geary's c Extension"
    ),
    Lee2012_Mahalanobis = list(
      global = lee_2012_results$global_stat,
      local = lee_2012_results$local_stats,
      p_value = lee_2012_results$global_p_value,
      local_p_values = lee_2012_results$local_p_values,
      significant_regions = lee_2012_results$significant_regions,
      method = "이몽현(2012) Mahalanobis Distance with Chi-square Test"
    )
  )
  
  return(benchmark_results)
}

#' Compute LIMSA (Local Indicator of Multivariate Spatial Association)
#' 
#' @description Anselin(2019)의 LIMSA 구현 - Geary's c의 다변량 확장
#' @param data 다변량 데이터
#' @param W 공간 가중치 행렬
#' @param n_perm 순열 검정 횟수
#' @return LIMSA 결과
#' @details 
#' LIMSA_i = Σ_j W_ij * ||x_i - x_j||²
#' 여기서 ||x_i - x_j||²는 i와 j 간의 다변량 거리제곱
compute_LIMSA <- function(data, W, n_perm = 999) {
  n <- nrow(data)
  p <- ncol(data)
  
  # 데이터 표준화
  data_scaled <- scale(data)
  
  # Local LIMSA 계산
  local_limsa <- numeric(n)
  
  for (i in 1:n) {
    # i번째 관측값과 이웃들 간의 거리제곱 계산
    neighbors <- which(W[i, ] > 0)
    if (length(neighbors) > 0) {
      # 다변량 거리제곱 계산
      distances_sq <- apply(data_scaled[neighbors, , drop = FALSE], 1, function(x_j) {
        sum((data_scaled[i, ] - x_j)^2)
      })
      
      # 가중 평균으로 LIMSA 계산
      local_limsa[i] <- sum(W[i, neighbors] * distances_sq)
    } else {
      local_limsa[i] <- 0
    }
  }
  
  # Global LIMSA (평균)
  global_limsa <- mean(local_limsa)
  
  # 순열 검정을 통한 유의성 평가
  if (n_perm > 0) {
    perm_global <- numeric(n_perm)
    
    for (perm in 1:n_perm) {
      # 데이터 순열
      perm_idx <- sample(1:n)
      data_perm <- data_scaled[perm_idx, ]
      
      # 순열된 데이터로 LIMSA 계산
      perm_local <- numeric(n)
      for (i in 1:n) {
        neighbors <- which(W[i, ] > 0)
        if (length(neighbors) > 0) {
          distances_sq <- apply(data_perm[neighbors, , drop = FALSE], 1, function(x_j) {
            sum((data_perm[i, ] - x_j)^2)
          })
          perm_local[i] <- sum(W[i, neighbors] * distances_sq)
        }
      }
      perm_global[perm] <- mean(perm_local)
    }
    
    # p-value 계산 (낮은 LIMSA가 높은 유사성을 의미하므로 단측검정)
    p_value <- (sum(perm_global <= global_limsa) + 1) / (n_perm + 1)
  } else {
    p_value <- NA
  }
  
  return(list(
    local_limsa = local_limsa,
    global_limsa = global_limsa,
    p_value = p_value,
    expected_null = ifelse(n_perm > 0, mean(perm_global), NA),
    method = "Anselin(2019) LIMSA"
  ))
}

#' Compute 이몽현(2012) Mahalanobis Distance-based Spatial Autocorrelation
#' 
#' @description 이몽현(2012)이 제안한 한 지역과 주변 지역의 마할라노비스 거리를 이용하고 
#'              카이제곱 분포로 검정하는 방법
#' @param data 다변량 데이터
#' @param W 공간 가중치 행렬
#' @param alpha 유의수준 (기본값: 0.05)
#' @return 이몽현(2012) 방법 결과
#' @details 
#' 이 방법은 각 지역과 그 이웃 지역들의 평균 간의 마할라노비스 거리를 계산하고,
#' 이 거리들이 카이제곱 분포를 따른다고 가정하여 통계적 유의성을 검정합니다.
#' 
#' 단계 (Lee, 2012):
#' 1. 대상 지역 (i) 설정
#' 2. 이웃 지역들 (j) 선택
#' 3. 이웃 지역들의 각 변수 평균값 계산  
#' 4. 대상 지역과 이웃 지역들의 평균 간 마할라노비스 거리 계산
#' 
#' MD_i = sqrt((x_i - x̄_neighbors)' * C^(-1) * (x_i - x̄_neighbors))
compute_Lee2012_Mahalanobis <- function(data, W, alpha = 0.05) {
  n <- nrow(data)
  p <- ncol(data)
  
  # 데이터 표준화
  data_scaled <- scale(data)
  
  # 전체 데이터의 공분산 행렬 계산
  cov_matrix <- cov(data_scaled)
  
  # 공분산 행렬의 역행렬 계산 (특이값 분해를 이용한 안정적 계산)
  cov_inv <- tryCatch({
    solve(cov_matrix)
  }, error = function(e) {
    # 특이값 분해를 이용한 pseudo-inverse
    svd_result <- svd(cov_matrix)
    # 작은 특이값 제거 (수치적 안정성)
    tolerance <- 1e-10
    d_inv <- ifelse(svd_result$d > tolerance, 1/svd_result$d, 0)
    svd_result$v %*% diag(d_inv) %*% t(svd_result$u)
  })
  
  # 각 지역의 로컬 마할라노비스 거리 기반 통계량 계산
  local_stats <- numeric(n)
  local_p_values <- numeric(n)
  
  for (i in 1:n) {
    # i번째 지역의 이웃들 찾기
    neighbors <- which(W[i, ] > 0)
    
    if (length(neighbors) > 0) {
      # Lee(2012) 방법: 이웃들의 가중 평균 계산
      weights <- W[i, neighbors] / sum(W[i, neighbors])  # 정규화된 가중치
      
      # 이웃들의 가중 평균 벡터 계산
      neighbor_mean <- apply(data_scaled[neighbors, , drop = FALSE], 2, function(x) {
        sum(weights * x)
      })
      
      # 해당 지역과 이웃들의 평균 간의 차이 벡터
      diff_vector <- data_scaled[i, ] - neighbor_mean
      
      # 마할라노비스 거리의 제곱 계산
      local_stats[i] <- as.numeric(t(diff_vector) %*% cov_inv %*% diff_vector)
      
      # 카이제곱 검정 (자유도 = 변수 개수)
      # 높은 거리 = 낮은 유사성 = 음의 공간자기상관을 의미
      # 낮은 거리 = 높은 유사성 = 양의 공간자기상관을 의미
      local_p_values[i] <- 1 - pchisq(local_stats[i], df = p)
      
    } else {
      # 이웃이 없는 경우
      local_stats[i] <- 0
      local_p_values[i] <- 1
    }
  }
  
  # 전역 통계량 계산
  global_stat <- mean(local_stats)
  
  # 전역 유의성 검정
  # 각 지역의 통계량들이 독립적이지 않으므로, 보수적인 접근 사용
  # 유의한 지역의 비율을 기반으로 전체적 유의성 평가
  n_significant <- sum(local_p_values < alpha)
  expected_significant <- n * alpha  # 귀무가설 하에서 기대되는 유의한 지역 수
  
  # 이항검정을 이용한 전역 p-value
  global_p_value <- 1 - pbinom(n_significant - 1, n, alpha)
  
  # 유의한 지역들 분류
  significant_regions <- data.frame(
    region_id = 1:n,
    mahalanobis_stat = local_stats,
    p_value = local_p_values,
    is_significant = local_p_values < alpha,
    spatial_pattern = ifelse(local_p_values < alpha,
                           ifelse(local_stats < qchisq(0.5, df = p), "High Similarity (Hotspot)", 
                                  "Low Similarity (Coldspot)"),
                           "Not Significant")
  )
  
  # 다중검정 보정 (FDR)
  local_p_adjusted <- p.adjust(local_p_values, method = "fdr")
  significant_regions$p_adjusted <- local_p_adjusted
  significant_regions$is_significant_adjusted <- local_p_adjusted < alpha
  
  return(list(
    local_stats = local_stats,
    local_p_values = local_p_values,
    local_p_adjusted = local_p_adjusted,
    global_stat = global_stat,
    global_p_value = global_p_value,
    significant_regions = significant_regions,
    n_significant = n_significant,
    n_significant_adjusted = sum(local_p_adjusted < alpha),
    cov_matrix = cov_matrix,
    method = "이몽현(2012) Mahalanobis Distance with Chi-square Test",
    description = paste0(
      "각 지역과 이웃 지역들의 평균 간 마할라노비스 거리를 이용한 공간자기상관 분석. ",
      "카이제곱 분포(df=", p, ")로 유의성 검정. ",
      "논문 원문: '마할라노비스 거리를 이용한 다변량 공간 클러스터 분석' (Lee, 2012). ",
      "낮은 거리 = 높은 유사성(Hotspot), 높은 거리 = 낮은 유사성(Coldspot)"
    )
  ))
}

# === 5. 보조 함수들 =======================================================

#' Extract MPSA values and categories from results
#' 
#' @param mpsa_results results from run_basic_MPSA_analysis()
#' @return data frame with MPSA statistics
extract_MPSA_results <- function(mpsa_results) {
  data.frame(
    MPSA = mpsa_results$local$MPSA,
    p_value = mpsa_results$local$p_adjusted,
    category = mpsa_results$local$category,
    effect_size = mpsa_results$local$effect_size,
    stringsAsFactors = FALSE
  )
}

#' Create MPSA summary table
#' 
#' @param mpsa_results results from run_basic_MPSA_analysis()
#' @return summary statistics table
create_MPSA_summary <- function(mpsa_results) {
  stats <- mpsa_results$local
  
  summary_table <- data.frame(
    Statistic = c("Mean MPSA", "Median MPSA", "SD MPSA", "Min MPSA", "Max MPSA",
                  "Global MPSA", "P-value (Global)", "Significant Units (α=0.05)"),
    Value = c(
      round(mean(stats$MPSA), 4),
      round(median(stats$MPSA), 4), 
      round(sd(stats$MPSA), 4),
      round(min(stats$MPSA), 4),
      round(max(stats$MPSA), 4),
      round(mpsa_results$global$GMPSA, 4),
      round(mpsa_results$global$p_value, 4),
      paste0(sum(stats$p_adjusted < 0.05), " / ", length(stats$p_adjusted))
    ),
    stringsAsFactors = FALSE
  )
  
  return(summary_table)
}

# === 실행 예시 =============================================================
results <- run_basic_MPSA_analysis(franklin)
map_plot <- results$plot
map_plot
tmap_save(map_plot, "output/Local MPSA Results.png", dpi = 900)

# 히스토그램 시각화
df_local <- data.frame(MPSA = results$local$MPSA)

# MPSA 히스토그램
MPSA_plot <- ggplot(df_local, aes(x = MPSA)) +
  geom_histogram(
    bins = 10,
    fill = "steelblue",
    color = "white",
    alpha = 0.85
  ) +
  labs(
    title = "Distribution of Local MPSA Values",
    x = "Local MPSA",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)

ggsave("output/MPSA_hist.png", MPSA_plot, dpi = 900)

df_local$p_adjusted <- results$local$p_adjusted

p_plot <- ggplot(df_local, aes(x = p_adjusted)) +
  geom_histogram(
    bins = 10,
    fill = "tomato",
    color = "white",
    alpha = 0.85
  ) +
  labs(
    title = "Distribution of Adjusted p-values (FDR)",
    x = "Adjusted p-value",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14)
ggsave("output/Pvalue_hist.png", p_plot, dpi = 900)

# 핫스팟 지역 파악
pmsa_results <- extract_MPSA_results(results)
franklin_results <- franklin |> 
  cbind(pmsa_results)
tmap_mode("view")
tm_shape(franklin_results) +
  tm_polygons(fill = "p_value")

df_box <- franklin_results |> 
  dplyr::select(MPSA, category, p_value)

boxplot <- ggplot(df_box, aes(x = category, y = MPSA, fill = category)) +
  geom_boxplot(alpha = 0.8, color = "black") +
  scale_fill_manual(values = c(
    "Strong Hotspot" = "#d73027",
    "Hotspot" = "#fc8d59",
    "Not Significant" = "#ffffbf", 
    "Coldspot" = "#91bfdb",
    "Strong Coldspot" = "#4575b4"
  )) +
  labs(
    title = "Local MPSA Distribution by Category",
    x = "MPSA Category",
    y = "Local MPSA"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1)

ggsave("output/boxplot.png", boxplot, dpi = 900)

boxplot2 <- ggplot(df_box, aes(x = category, y = p_value, fill = category)) +
  geom_boxplot(alpha = 0.8, color = "black") +
  scale_fill_manual(values = c(
    "Strong Hotspot" = "#d73027",
    "Hotspot" = "#fc8d59",
    "Not Significant" = "#ffffbf", 
    "Coldspot" = "#91bfdb",
    "Strong Coldspot" = "#4575b4"
  )) +
  labs(
    title = "Adjusted P-value Distribution by Category",
    x = "MPSA Category",
    y = "Adjusted P-value"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
  coord_cartesian(ylim = c(-0.25, 1))

ggsave("output/p_boxplot.png", boxplot2, dpi = 900)

# 변수 중요도 시각화

imp <- as.data.frame(results$rf_model$importance)
imp$Variable <- rownames(imp)
imp_sd <- as.data.frame(results$rf_model$importanceSD)
imp_sd$Variable <- rownames(imp_sd)

imp_plot <- merge(imp, imp_sd, by = "Variable", suffixes = c("", "_sd"))

ggplot(imp_plot, aes(x = reorder(Variable, MeanDecreaseAccuracy), 
                     y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = MeanDecreaseAccuracy - MeanDecreaseAccuracy_sd,
                    ymax = MeanDecreaseAccuracy + MeanDecreaseAccuracy_sd),
                width = 0.2) +
  coord_flip() +
  labs(x = "Variable", y = "Mean Decrease in Accuracy",
       title = "Variable Importance with Standard Deviation")


# 벤치마크 비교 예시
numeric_data <- franklin |>
  st_drop_geometry() |>
  dplyr::select(where(is.numeric))
benchmark_results <- benchmark_against_traditional_methods(numeric_data, W)

randomForest:::randomForest.default
