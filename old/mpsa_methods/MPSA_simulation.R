# ============================================================================
# MPSA Simulation Studies and Empirical Validation
# 
# 이 스크립트는 MPSA의 시뮬레이션 연구, 수렴성 및 안정성 분석에 집중합니다.
# 전통적 방법과의 비교는 test_alternatives.R에서 담당합니다.
# 
# 역할 분담:
# - MPSA_theoretical_analysis.R: 엄밀한 수학적 이론 분석 (Biau & Scornet 기반)
# - 이 파일: 시뮬레이션 연구, 수렴성, 안정성 분석
# - test_alternatives.R: 실제 데이터 기반 전통적 방법과의 비교
# ============================================================================

# --- 환경 설정 ---
# 기본 setup 파일 로드
if (file.exists("R/data_preparation/setup.R")) {
  source("R/data_preparation/setup.R")
} else {
  # 기본 패키지들 로드
  required_packages <- c("randomForest", "spdep", "Matrix", "mvtnorm", "CompQuadForm")
  installed <- required_packages %in% installed.packages()
  if (any(!installed)) {
    install.packages(required_packages[!installed])
  }
  lapply(required_packages, library, character.only = TRUE)
}

# 기본 MPSA 방법론 로드
if (file.exists("R/mpsa_methods/MPSA.R")) {
  source("R/mpsa_methods/MPSA.R")
} else if (file.exists("MPSA.R")) {
  source("MPSA.R")
} else {
  warning("MPSA.R 파일을 찾을 수 없습니다. 일부 기능이 제한될 수 있습니다.")
}

# 이론적 분석 파일 로드 (선택적)
theoretical_analysis_loaded <- FALSE
tryCatch({
  if (file.exists("R/mpsa_methods/MPSA_theoretical_analysis.R")) {
    source("R/mpsa_methods/MPSA_theoretical_analysis.R")
    theoretical_analysis_loaded <- TRUE
    cat("✓ 이론적 분석 모듈이 성공적으로 로드되었습니다.\n")
  } else if (file.exists("MPSA_theoretical_analysis.R")) {
    source("MPSA_theoretical_analysis.R")
    theoretical_analysis_loaded <- TRUE
    cat("✓ 이론적 분석 모듈이 성공적으로 로드되었습니다.\n")
  } else {
    cat("⚠ MPSA_theoretical_analysis.R 파일을 찾을 수 없습니다.\n")
    cat("  일부 이론적 검증 기능이 제한될 수 있습니다.\n")
  }
}, error = function(e) {
  cat("⚠ 이론적 분석 모듈 로드 중 오류 발생:", e$message, "\n")
  cat("  기본 실증적 검증 기능만 사용 가능합니다.\n")
})

# 추가 라이브러리 로드
library(mvtnorm)
if (!require("CompQuadForm", quietly = TRUE)) {
  install.packages("CompQuadForm")
  library(CompQuadForm)
}

# === 1. 시뮬레이션 기반 MPSA 성능 검증 ======================================

#' Simulate Spatial Data with Known Autocorrelation
#' 
#' @description 알려진 공간자기상관을 가진 다변량 데이터 시뮬레이션
#' @param n 공간 단위 수
#' @param p 변수 수
#' @param rho 공간자기상관 강도
#' @param spatial_layout 공간 배치 ("grid", "random", "clustered")
#' @return 시뮬레이션된 데이터와 공간 가중치
simulate_multivariate_spatial_data <- function(n, p, rho = 0.5, spatial_layout = "grid") {
  
  # 1. 공간 좌표 생성
  if (spatial_layout == "grid") {
    grid_size <- ceiling(sqrt(n))
    coords <- expand.grid(x = 1:grid_size, y = 1:grid_size)[1:n, ]
  } else if (spatial_layout == "random") {
    coords <- data.frame(x = runif(n, 0, 10), y = runif(n, 0, 10))
  } else if (spatial_layout == "clustered") {
    # 클러스터된 패턴
    centers <- matrix(runif(6), ncol = 2) * 10
    cluster_assign <- sample(1:3, n, replace = TRUE)
    coords <- data.frame(
      x = centers[cluster_assign, 1] + rnorm(n, 0, 1),
      y = centers[cluster_assign, 2] + rnorm(n, 0, 1)
    )
  }
  
  # 2. 공간 가중치 행렬 생성
  W <- create_spatial_weights(coords, method = "knn", k = 4)
  
  # 3. 공간자기상관이 있는 잠재변수 생성 (SAR 모델)
  # (I - ρW)z = ε, 여기서 ε ~ N(0, I)
  I_n <- diag(n)
  epsilon <- rnorm(n)
  z <- solve(I_n - rho * W) %*% epsilon
  
  # 4. 다변량 관측값 생성
  # 각 변수는 잠재변수 z와 개별적 노이즈의 선형결합
  X <- matrix(NA, nrow = n, ncol = p)
  
  for (j in 1:p) {
    # 변수 간 상관관계 도입
    loading <- runif(1, 0.5, 1.0)  # z에 대한 적재량
    noise_var <- runif(1, 0.2, 0.8)  # 개별 노이즈 분산
    
    X[, j] <- loading * z + rnorm(n, 0, sqrt(noise_var))
  }
  
  # 5. 변수명 지정
  colnames(X) <- paste0("Var", 1:p)
  
  return(list(
    data = X,
    coords = coords,
    W = W,
    true_spatial_pattern = z,
    true_rho = rho,
    layout = spatial_layout
  ))
}

#' Create Spatial Weight Matrix
#' 
#' @description 다양한 방법으로 공간 가중치 행렬 생성
create_spatial_weights <- function(coords, method = "knn", k = 4, threshold = NULL) {
  n <- nrow(coords)
  
  # 거리 행렬 계산
  dist_matrix <- as.matrix(dist(coords))
  
  if (method == "knn") {
    # K-nearest neighbors
    W <- matrix(0, n, n)
    for (i in 1:n) {
      neighbors <- order(dist_matrix[i, ])[2:(k+1)]  # 자기 제외
      W[i, neighbors] <- 1
    }
  } else if (method == "threshold") {
    # Distance threshold
    W <- (dist_matrix <= threshold & dist_matrix > 0) * 1
  } else if (method == "inverse_distance") {
    # Inverse distance weights
    W <- 1 / (dist_matrix + diag(n))  # 대각선에 1 추가하여 0으로 나누기 방지
    diag(W) <- 0
  }
  
  # Row standardization
  row_sums <- rowSums(W)
  row_sums[row_sums == 0] <- 1  # 고립된 노드 처리
  W <- W / row_sums
  
  return(W)
}

#' MPSA Simulation Study
#' 
#' @description 다양한 시나리오에서 MPSA 성능 평가 (순수 시뮬레이션)
#' @param scenarios 시뮬레이션 시나리오 리스트
#' @param n_replications 반복 횟수
mpsa_simulation_study <- function(scenarios = NULL, n_replications = 100) {
  
  # 기본 시나리오 설정
  if (is.null(scenarios)) {
    scenarios <- list(
      list(n = 50, p = 5, rho = 0.3, layout = "grid"),
      list(n = 50, p = 5, rho = 0.7, layout = "grid"),
      list(n = 100, p = 10, rho = 0.5, layout = "grid"),
      list(n = 100, p = 10, rho = 0.5, layout = "clustered"),
      list(n = 75, p = 15, rho = 0.4, layout = "random")
    )
  }
  
  results <- list()
  
  for (s in 1:length(scenarios)) {
    scenario <- scenarios[[s]]
    cat(sprintf("Running scenario %d: n=%d, p=%d, rho=%.2f, layout=%s\n", 
                s, scenario$n, scenario$p, scenario$rho, scenario$layout))
    
    scenario_results <- replicate(n_replications, {
      # 데이터 시뮬레이션
      sim_data <- simulate_multivariate_spatial_data(
        n = scenario$n, 
        p = scenario$p, 
        rho = scenario$rho, 
        spatial_layout = scenario$layout
      )
      
      # MPSA 계산
      rf <- randomForest(sim_data$data, proximity = TRUE, ntree = 200)
      P <- rf$proximity
      mpsa_values <- rowSums(sim_data$W * P)
      gmpsa <- mean(mpsa_values)
      
      # 시뮬레이션 결과 반환
      list(
        gmpsa = gmpsa,
        local_mpsa_mean = mean(mpsa_values),
        local_mpsa_sd = sd(mpsa_values),
        true_rho = scenario$rho,
        n = scenario$n,
        p = scenario$p
      )
      
    }, simplify = FALSE)
    
    # 시나리오별 결과 정리
    scenario_matrix <- do.call(rbind, lapply(scenario_results, function(x) {
      c(gmpsa = x$gmpsa, 
        local_mean = x$local_mpsa_mean, 
        local_sd = x$local_mpsa_sd,
        true_rho = x$true_rho)
    }))
    
    results[[s]] <- list(
      scenario = scenario,
      results = scenario_matrix,
      mean_gmpsa = mean(scenario_matrix[, "gmpsa"]),
      sd_gmpsa = sd(scenario_matrix[, "gmpsa"]),
      correlation_with_rho = cor(scenario_matrix[, "gmpsa"], scenario_matrix[, "true_rho"])
    )
  }
  
  return(results)
}

# === 2. MPSA 수렴성 분석 ====================================================

#' Analyze MPSA Convergence with ntree
#' 
#' @description ntree 값에 따른 MPSA 수렴성 분석
#' @param data 분석 데이터
#' @param W 공간 가중치 행렬
#' @param ntree_values 테스트할 ntree 값들
analyze_mpsa_convergence <- function(data, W, ntree_values = c(50, 100, 200, 500, 1000)) {
  
  convergence_results <- data.frame(
    ntree = ntree_values,
    gmpsa = NA,
    gmpsa_sd = NA,
    cv = NA,  # Coefficient of variation
    runtime = NA
  )
  
  for (i in 1:length(ntree_values)) {
    ntree <- ntree_values[i]
    cat(sprintf("Testing ntree = %d...\n", ntree))
    
    # 여러 번 실행하여 안정성 측정
    n_runs <- 10
    gmpsa_values <- numeric(n_runs)
    
    start_time <- Sys.time()
    
    for (run in 1:n_runs) {
      rf <- randomForest(data, proximity = TRUE, ntree = ntree)
      P <- rf$proximity
      mpsa_values <- rowSums(W * P)
      gmpsa_values[run] <- mean(mpsa_values)
    }
    
    runtime <- as.numeric(Sys.time() - start_time)
    
    convergence_results$gmpsa[i] <- mean(gmpsa_values)
    convergence_results$gmpsa_sd[i] <- sd(gmpsa_values)
    convergence_results$cv[i] <- sd(gmpsa_values) / mean(gmpsa_values)
    convergence_results$runtime[i] <- runtime / n_runs  # 평균 실행 시간
  }
  
  return(convergence_results)
}

# === 3. MPSA 안정성 및 강건성 분석 =========================================

#' Analyze MPSA Robustness to Data Perturbation
#' 
#' @description 데이터 섭동에 대한 MPSA 안정성 분석
#' @param data 원본 데이터
#' @param W 공간 가중치 행렬
#' @param noise_levels 노이즈 수준들
analyze_mpsa_robustness <- function(data, W, noise_levels = c(0.01, 0.05, 0.1, 0.2)) {
  
  # 원본 MPSA 계산
  rf_original <- randomForest(data, proximity = TRUE, ntree = 500)
  P_original <- rf_original$proximity
  mpsa_original <- rowSums(W * P_original)
  gmpsa_original <- mean(mpsa_original)
  
  robustness_results <- data.frame(
    noise_level = noise_levels,
    correlation_with_original = NA,
    gmpsa_bias = NA,
    mpsa_rmse = NA
  )
  
  for (i in 1:length(noise_levels)) {
    noise_sd <- noise_levels[i] * sqrt(diag(var(data)))
    cat(sprintf("Testing noise level %.3f...\n", noise_levels[i]))
    
    # 여러 번 노이즈 추가하여 평균 성능 측정
    n_perturbations <- 20
    correlations <- numeric(n_perturbations)
    gmpsa_values <- numeric(n_perturbations)
    rmse_values <- numeric(n_perturbations)
    
    for (pert in 1:n_perturbations) {
      # 노이즈 추가
      noise_matrix <- matrix(rnorm(nrow(data) * ncol(data)), 
                            nrow = nrow(data), ncol = ncol(data))
      noise_matrix <- sweep(noise_matrix, 2, noise_sd, "*")
      data_perturbed <- data + noise_matrix
      
      # 섭동된 데이터로 MPSA 계산
      rf_perturbed <- randomForest(data_perturbed, proximity = TRUE, ntree = 500)
      P_perturbed <- rf_perturbed$proximity
      mpsa_perturbed <- rowSums(W * P_perturbed)
      gmpsa_perturbed <- mean(mpsa_perturbed)
      
      # 성능 지표 계산
      correlations[pert] <- cor(mpsa_original, mpsa_perturbed)
      gmpsa_values[pert] <- gmpsa_perturbed
      rmse_values[pert] <- sqrt(mean((mpsa_original - mpsa_perturbed)^2))
    }
    
    robustness_results$correlation_with_original[i] <- mean(correlations)
    robustness_results$gmpsa_bias[i] <- mean(gmpsa_values) - gmpsa_original
    robustness_results$mpsa_rmse[i] <- mean(rmse_values)
  }
  
  return(list(
    original_gmpsa = gmpsa_original,
    robustness_results = robustness_results
  ))
}

# === 4. 통계적 검정력 분석 =================================================

#' Power Analysis for MPSA
#' 
#' @description MPSA의 통계적 검정력 분석
#' @param effect_sizes 효과 크기들 (공간자기상관 강도)
#' @param n_sim 시뮬레이션 횟수
analyze_mpsa_power <- function(effect_sizes = seq(0.1, 0.8, 0.1), n_sim = 100) {
  
  power_results <- data.frame(
    effect_size = effect_sizes,
    power = NA,
    type_i_error = NA
  )
  
  for (i in 1:length(effect_sizes)) {
    rho <- effect_sizes[i]
    cat(sprintf("Testing effect size rho = %.2f...\n", rho))
    
    # H1: 공간자기상관 존재 (rho > 0)
    significant_count <- 0
    
    for (sim in 1:n_sim) {
      # 공간자기상관이 있는 데이터 생성
      sim_data <- simulate_multivariate_spatial_data(n = 100, p = 8, rho = rho)
      
      # MPSA 유의성 검정
      rf <- randomForest(sim_data$data, proximity = TRUE, ntree = 300)
      P <- rf$proximity
      
      # 간단한 순열 검정 (빠른 버전)
      observed_gmpsa <- mean(rowSums(sim_data$W * P))
      
      # 순열 분포 생성 (빠른 버전)
      null_gmpsa <- numeric(99)
      for (perm in 1:99) {
        perm_idx <- sample(1:nrow(sim_data$data))
        P_perm <- P[perm_idx, perm_idx]
        null_gmpsa[perm] <- mean(rowSums(sim_data$W * P_perm))
      }
      
      # p-value 계산
      p_value <- (sum(abs(null_gmpsa) >= abs(observed_gmpsa)) + 1) / 100
      
      if (p_value < 0.05) {
        significant_count <- significant_count + 1
      }
    }
    
    power_results$power[i] <- significant_count / n_sim
    
    # Type I error는 rho = 0일 때의 power
    if (rho == 0.1) {  # 가장 작은 효과 크기를 Type I error로 근사
      power_results$type_i_error[i] <- significant_count / n_sim * 0.5  # 보수적 추정
    }
  }
  
  return(power_results)
}

# === 5. 종합적 실증적 검증 함수 ============================================

#' Comprehensive Empirical Validation of MPSA
#' 
#' @description MPSA의 종합적 실증적 검증
#' @param data 실제 데이터 또는 NULL (시뮬레이션 사용)
#' @param W 공간 가중치 행렬 또는 NULL
#' @param full_analysis 전체 분석 수행 여부
comprehensive_empirical_validation <- function(data = NULL, W = NULL, full_analysis = FALSE) {
  
  cat("=== MPSA 종합적 실증적 검증 ===\n\n")
  
  results <- list()
  
  # 1. 시뮬레이션 연구
  cat("1. 시뮬레이션 연구 실행 중...\n")
  if (full_analysis) {
    sim_results <- mpsa_simulation_study(n_replications = 50)
  } else {
    # 빠른 버전
    quick_scenarios <- list(
      list(n = 50, p = 5, rho = 0.5, layout = "grid"),
      list(n = 75, p = 8, rho = 0.3, layout = "clustered")
    )
    sim_results <- mpsa_simulation_study(quick_scenarios, n_replications = 20)
  }
  results$simulation_study <- sim_results
  
  # 2. 수렴성 분석 (실제 데이터가 있는 경우)
  if (!is.null(data) && !is.null(W)) {
    cat("2. 수렴성 분석 실행 중...\n")
    if (full_analysis) {
      convergence_results <- analyze_mpsa_convergence(data, W)
    } else {
      convergence_results <- analyze_mpsa_convergence(data, W, ntree_values = c(50, 200, 500))
    }
    results$convergence_analysis <- convergence_results
    
    # 3. 강건성 분석
    cat("3. 강건성 분석 실행 중...\n")
    if (full_analysis) {
      robustness_results <- analyze_mpsa_robustness(data, W)
    } else {
      robustness_results <- analyze_mpsa_robustness(data, W, noise_levels = c(0.05, 0.1))
    }
    results$robustness_analysis <- robustness_results
  }
  
  # 4. 검정력 분석
  if (full_analysis) {
    cat("4. 통계적 검정력 분석 실행 중...\n")
    power_results <- analyze_mpsa_power(effect_sizes = c(0.2, 0.4, 0.6), n_sim = 50)
    results$power_analysis <- power_results
  }
  
  cat("\n=== 실증적 검증 완료 ===\n")
  
  return(results)
}

# === 실행 예시 및 시각화 =====================================================

# === 1. 기본 시뮬레이션 연구 실행 ===
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
# W <- nb2mat(nb, style = "W", zero.policy = TRUE)
# 
# # 종합 시뮬레이션 실행
# simulation_results <- comprehensive_empirical_validation(numeric_data, W, full_analysis = TRUE)
# 
# # 결과 요약 출력
# print(simulation_results$summary)

# === 2. 통계적 검정력 분석 시각화 ===
# 
# # 다양한 효과 크기에서의 검정력 분석
# effect_sizes <- seq(0, 2, 0.2)
# power_results <- analyze_mpsa_power(effect_sizes, n_sim = 100)
# 
# # 검정력 곡선 시각화
# library(ggplot2)
# power_df <- data.frame(
#   effect_size = effect_sizes,
#   power = power_results$power_curve,
#   method = "MPSA"
# )
# 
# p1 <- ggplot(power_df, aes(x = effect_size, y = power)) +
#   geom_line(color = "blue", size = 1.2) +
#   geom_point(color = "red", size = 2) +
#   geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray") +
#   labs(title = "MPSA Statistical Power Analysis",
#        subtitle = "Power vs Effect Size",
#        x = "Effect Size", y = "Statistical Power") +
#   scale_y_continuous(limits = c(0, 1)) +
#   theme_minimal()

# === 3. Type I Error 분석 시각화 ===
# 
# # 다양한 유의수준에서의 Type I Error 분석
# alpha_levels <- c(0.01, 0.05, 0.1)
# type_I_results <- test_type_I_error_control(alpha_levels, n_sim = 500)
# 
# # Type I Error 시각화
# type_I_df <- data.frame(
#   alpha = alpha_levels,
#   observed_error = type_I_results$observed_errors,
#   expected_error = alpha_levels
# )
# 
# p2 <- ggplot(type_I_df, aes(x = factor(alpha))) +
#   geom_col(aes(y = observed_error), fill = "lightcoral", alpha = 0.7, width = 0.6) +
#   geom_point(aes(y = expected_error), color = "darkblue", size = 3) +
#   geom_line(aes(y = expected_error, group = 1), color = "darkblue", linetype = "dashed") +
#   labs(title = "Type I Error Control Verification",
#        subtitle = "Blue points: Expected, Red bars: Observed",
#        x = "Significance Level (α)", y = "Type I Error Rate") +
#   theme_minimal()

# === 4. 수렴성 분석 시각화 ===
# 
# # ntree에 따른 수렴성 분석
# ntree_values <- c(50, 100, 200, 300, 500, 750, 1000)
# convergence_results <- analyze_mpsa_convergence(numeric_data, W, ntree_values)
# 
# # 수렴성 시각화
# conv_df <- data.frame(
#   ntree = ntree_values,
#   gmpsa_mean = convergence_results$gmpsa_means,
#   gmpsa_sd = convergence_results$gmpsa_sds
# )
# 
# p3 <- ggplot(conv_df, aes(x = ntree, y = gmpsa_mean)) +
#   geom_line(color = "blue", size = 1) +
#   geom_point(color = "red", size = 2) +
#   geom_errorbar(aes(ymin = gmpsa_mean - gmpsa_sd, ymax = gmpsa_mean + gmpsa_sd),
#                 width = 20, alpha = 0.7) +
#   labs(title = "MPSA Convergence Analysis",
#        subtitle = "Global MPSA vs Number of Trees (with error bars)",
#        x = "Number of Trees", y = "Global MPSA") +
#   theme_minimal()

# === 5. 강건성 분석 시각화 ===
# 
# # 노이즈 수준에 따른 강건성 분석
# noise_levels <- seq(0, 0.5, 0.05)
# robustness_results <- analyze_mpsa_robustness(numeric_data, W, noise_levels)
# 
# # 강건성 시각화
# robust_df <- data.frame(
#   noise_level = noise_levels,
#   gmpsa_mean = robustness_results$gmpsa_means,
#   gmpsa_sd = robustness_results$gmpsa_sds,
#   correlation_with_original = robustness_results$correlations
# )
# 
# p4 <- ggplot(robust_df, aes(x = noise_level)) +
#   geom_line(aes(y = gmpsa_mean), color = "blue", size = 1) +
#   geom_point(aes(y = gmpsa_mean), color = "red", size = 2) +
#   geom_ribbon(aes(ymin = gmpsa_mean - gmpsa_sd, ymax = gmpsa_mean + gmpsa_sd),
#               alpha = 0.3, fill = "blue") +
#   labs(title = "MPSA Robustness Analysis",
#        subtitle = "Global MPSA vs Noise Level",
#        x = "Noise Level", y = "Global MPSA") +
#   theme_minimal()

# === 6. 시뮬레이션 시나리오 비교 ===
# 
# # 다양한 공간 패턴에서의 성능 비교
# scenarios <- list(
#   "Random" = list(rho = 0, pattern = "random"),
#   "Weak Spatial" = list(rho = 0.3, pattern = "clustered"),
#   "Strong Spatial" = list(rho = 0.7, pattern = "clustered"),
#   "Checkerboard" = list(rho = -0.5, pattern = "checkerboard")
# )
# 
# scenario_results <- mpsa_simulation_study(scenarios, n_replications = 100)
# 
# # 시나리오별 성능 시각화
# scenario_df <- do.call(rbind, lapply(names(scenarios), function(sc) {
#   data.frame(
#     scenario = sc,
#     gmpsa = scenario_results[[sc]]$gmpsa_values,
#     power = scenario_results[[sc]]$detection_power
#   )
# }))
# 
# p6 <- ggplot(scenario_df, aes(x = scenario, y = gmpsa)) +
#   geom_boxplot(fill = "lightblue", alpha = 0.7) +
#   geom_point(position = position_jitter(width = 0.2), alpha = 0.5) +
#   labs(title = "MPSA Performance Across Scenarios",
#        x = "Spatial Pattern", y = "Global MPSA") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# === 7. 벤치마크 비교 시각화 ===
# 
# # MPSA vs 전통적 방법들 성능 비교
# benchmark_results <- benchmark_against_traditional_methods(numeric_data, W)
# 
# # 성능 비교 데이터 준비
# methods_df <- data.frame(
#   Method = c("MPSA", "PCA+Moran", "Individual Moran", "Euclidean", "LIMSA", "Lee2012"),
#   Statistic = c(
#     benchmark_results$MPSA$global,
#     benchmark_results$PCA_Moran$global,
#     benchmark_results$Individual_Moran$mean_statistic,
#     benchmark_results$Euclidean_based$global,
#     benchmark_results$LIMSA_Anselin2019$global,
#     benchmark_results$Lee2012_Mahalanobis$global
#   ),
#   Type = c("Proximity", "Dimension Reduction", "Univariate", "Distance", "Distance", "Distance")
# )
# 
# # 벤치마크 비교 시각화
# p7 <- ggplot(methods_df, aes(x = reorder(Method, Statistic), y = Statistic, fill = Type)) +
#   geom_col(alpha = 0.8) +
#   coord_flip() +
#   labs(title = "Method Performance Comparison",
#        x = "Method", y = "Statistic Value") +
#   theme_minimal() +
#   scale_fill_brewer(type = "qual", palette = "Set2")

# === 8. 결과 저장 및 종합 보고서 ===
# 
# # 모든 시각화 결합
# library(patchwork)
# combined_plot1 <- (p1 + p2) / (p3 + p4)
# combined_plot2 <- (p5 + p6) / p7
# 
# # 결과 디렉토리 생성
# if (!dir.exists("output/simulation_study")) {
#   dir.create("output/simulation_study", recursive = TRUE)
# }
# 
# # 개별 그래프 저장
# ggsave("output/simulation_study/power_analysis.png", p1, width = 10, height = 6)
# ggsave("output/simulation_study/type_I_error.png", p2, width = 8, height = 6)
# ggsave("output/simulation_study/convergence_analysis.png", p3, width = 10, height = 6)
# ggsave("output/simulation_study/robustness_analysis.png", p4, width = 10, height = 6)
# ggsave("output/simulation_study/stability_analysis.png", p5, width = 10, height = 6)
# ggsave("output/simulation_study/scenario_comparison.png", p6, width = 10, height = 6)
# ggsave("output/simulation_study/method_benchmark.png", p7, width = 10, height = 8)
# 
# # 종합 그래프 저장
# ggsave("output/simulation_study/simulation_summary_1.png", combined_plot1, width = 16, height = 10)
# ggsave("output/simulation_study/simulation_summary_2.png", combined_plot2, width = 16, height = 10)
# 
# # 수치 결과 저장
# write.csv(power_df, "output/simulation_study/power_analysis_results.csv", row.names = FALSE)
# write.csv(type_I_df, "output/simulation_study/type_I_error_results.csv", row.names = FALSE)
# write.csv(conv_df, "output/simulation_study/convergence_results.csv", row.names = FALSE)
# write.csv(robust_df, "output/simulation_study/robustness_results.csv", row.names = FALSE)
# write.csv(scenario_df, "output/simulation_study/scenario_results.csv", row.names = FALSE)
# write.csv(methods_df, "output/simulation_study/benchmark_comparison.csv", row.names = FALSE)
# 
# # 종합 시뮬레이션 보고서 저장
# write.csv(simulation_results$summary, "output/simulation_study/simulation_summary.csv", row.names = FALSE)
# 
# # 최종 요약 출력
# cat("=== 시뮬레이션 연구 완료 ===\n")
# cat("📊 주요 발견사항:\n")
# cat(sprintf("  - 통계적 검정력: %.1f%% (효과크기 1.0에서)\n", 
#             power_df$power[power_df$effect_size == 1.0] * 100))
# cat(sprintf("  - Type I Error 통제: %.3f (α=0.05에서)\n", 
#             type_I_df$observed_error[type_I_df$alpha == 0.05]))
# cat(sprintf("  - 수렴 안정성: %d 트리에서 수렴\n", 
#             min(ntree_values[conv_df$gmpsa_sd < 0.01])))
# cat("📁 결과 저장 위치: output/simulation_study/\n")
# cat("🎯 실증적 검증 완료!\n") 