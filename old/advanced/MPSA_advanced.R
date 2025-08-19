# ============================================================================
# Advanced Multivariate Proximity-based Spatial Autocorrelation (MPSA) Framework
# 
# 이 스크립트는 비지도 Random Forest의 proximity matrix를 활용한 
# 혁신적인 공간 분석 방법론을 구현합니다.
# 
# 주요 특징:
# 1. Multi-scale MPSA 분석
# 2. Spatio-temporal MPSA 확장
# 3. Bayesian MPSA 추론
# 4. Network-based MPSA
# 5. 고급 시각화 및 해석 도구
# ============================================================================

# --- 환경 설정 ---
source("R/data_preparation/setup.R")
library(igraph)
library(ggraph)
library(viridis)
library(patchwork)
library(broom)
library(gstat)
library(spacetime)
library(spBayes)

# === 1. 이론적 프레임워크 ===================================================

#' Global Multivariate Proximity-based Spatial Autocorrelation (GMPSA)
#' 
#' @description 전역적 공간 연관성 측정을 위한 GMPSA 통계량
#' @param P proximity matrix from random forest
#' @param W spatial weight matrix
#' @return GMPSA statistic and p-value
compute_GMPSA <- function(P, W, n_perm = 999) {
  # GMPSA 통계량: 전역적 공간-RF 연관성
  n <- nrow(P)
  
  # Moran's I와 유사한 형태로 정의
  P_centered <- P - mean(P)
  W_sum <- sum(W)
  
  GMPSA <- (n / W_sum) * sum(W * P_centered) / sum(P_centered^2)
  
  # 순열 검정
  GMPSA_perm <- numeric(n_perm)
  for (i in 1:n_perm) {
    perm_idx <- sample(1:n)
    P_perm <- P[perm_idx, perm_idx]
    P_perm_centered <- P_perm - mean(P_perm)
    GMPSA_perm[i] <- (n / W_sum) * sum(W * P_perm_centered) / sum(P_perm_centered^2)
  }
  
  # p-value 계산
  p_value <- (sum(abs(GMPSA_perm) >= abs(GMPSA)) + 1) / (n_perm + 1)
  
  # 기대값과 분산
  E_GMPSA <- mean(GMPSA_perm)
  Var_GMPSA <- var(GMPSA_perm)
  z_score <- (GMPSA - E_GMPSA) / sqrt(Var_GMPSA)
  
  return(list(
    GMPSA = GMPSA,
    expected = E_GMPSA,
    variance = Var_GMPSA,
    z_score = z_score,
    p_value = p_value,
    permutations = GMPSA_perm
  ))
}

# === 2. 다중 스케일 MPSA 분석 ===============================================

#' Multi-scale MPSA Analysis
#' 
#' @description 다양한 공간 스케일에서 MPSA 분석 수행
#' @param P proximity matrix
#' @param coords spatial coordinates
#' @param scales vector of distance thresholds
compute_multiscale_MPSA <- function(P, coords, scales = c(1000, 2500, 5000, 10000)) {
  n <- nrow(P)
  results <- list()
  
  # 각 스케일별로 분석
  for (k in seq_along(scales)) {
    # 거리 기반 가중치 행렬 생성
    dist_mat <- as.matrix(dist(coords))
    W_k <- (dist_mat <= scales[k]) * 1
    diag(W_k) <- 0
    
    # 행 표준화
    W_k_rs <- W_k / rowSums(W_k)
    W_k_rs[is.nan(W_k_rs)] <- 0
    
    # MPSA 계산
    mpsa_k <- rowSums(W_k_rs * P)
    
    # 유의성 검정
    sig_test <- compute_local_mpsa_significance_advanced(P, W_k_rs)
    
    results[[paste0("scale_", scales[k])]] <- list(
      scale = scales[k],
      W = W_k_rs,
      MPSA = mpsa_k,
      significance = sig_test,
      n_neighbors = rowSums(W_k)
    )
  }
  
  class(results) <- "multiscale_MPSA"
  return(results)
}

# === 3. 시공간 MPSA (Spatio-temporal MPSA) ==================================

#' Spatio-temporal MPSA
#' 
#' @description 시공간 데이터에 대한 MPSA 분석
#' @param data_list list of data frames for each time point
#' @param spatial_weights spatial weight matrix
#' @param temporal_weights temporal weight matrix
compute_spatiotemporal_MPSA <- function(data_list, spatial_weights, 
                                      temporal_lag = 1, alpha = 0.5) {
  n_time <- length(data_list)
  n_space <- nrow(data_list[[1]])
  
  # 각 시점별 proximity matrix 계산
  P_list <- list()
  for (t in 1:n_time) {
    rf_t <- randomForest(x = data_list[[t]], proximity = TRUE, ntree = 500)
    P_list[[t]] <- rf_t$proximity
  }
  
  # 시공간 MPSA 계산
  ST_MPSA <- matrix(NA, nrow = n_space, ncol = n_time)
  
  for (t in 1:n_time) {
    # 공간 성분
    spatial_component <- rowSums(spatial_weights * P_list[[t]])
    
    # 시간 성분 (이전 시점들의 영향)
    temporal_component <- 0
    if (t > 1) {
      for (lag in 1:min(temporal_lag, t-1)) {
        weight <- (1 - alpha)^(lag - 1) * alpha
        temporal_component <- temporal_component + 
          weight * diag(P_list[[t]] %*% t(P_list[[t-lag]]))
      }
    }
    
    # 결합
    ST_MPSA[, t] <- spatial_component + temporal_component
  }
  
  return(list(
    ST_MPSA = ST_MPSA,
    P_list = P_list,
    spatial_component = spatial_component,
    temporal_component = temporal_component
  ))
}

# === 4. Bayesian MPSA =======================================================

#' Bayesian MPSA with Spatial Random Effects
#' 
#' @description 베이지안 접근법을 통한 MPSA 추정
#' @param P proximity matrix
#' @param W spatial weight matrix
#' @param X covariates (optional)
#' @param n_iter MCMC iterations
compute_bayesian_MPSA <- function(P, W, X = NULL, n_iter = 10000, 
                                 burn_in = 2000, thin = 5) {
  n <- nrow(P)
  
  # MPSA를 종속변수로 설정
  y <- rowSums(W * P)
  
  # 설계 행렬 구성
  if (is.null(X)) {
    X <- matrix(1, n, 1)  # intercept only
  } else {
    X <- cbind(1, X)  # add intercept
  }
  
  # CAR 모델 설정을 위한 준비
  W_car <- W + t(W)  # symmetrize
  D <- diag(rowSums(W_car))
  Q <- D - W_car  # precision matrix for CAR
  
  # MCMC 초기값
  beta <- rep(0, ncol(X))
  tau2 <- 1  # error variance
  phi <- rep(0, n)  # spatial random effects
  rho <- 0.5  # spatial correlation parameter
  
  # MCMC 저장소
  n_save <- (n_iter - burn_in) / thin
  beta_samples <- matrix(NA, n_save, ncol(X))
  tau2_samples <- numeric(n_save)
  phi_samples <- matrix(NA, n_save, n)
  rho_samples <- numeric(n_save)
  
  # MCMC 샘플링
  iter_save <- 0
  for (iter in 1:n_iter) {
    # Update beta
    V_beta <- solve(t(X) %*% X / tau2 + diag(0.001, ncol(X)))
    m_beta <- V_beta %*% (t(X) %*% (y - phi) / tau2)
    beta <- MASS::mvrnorm(1, m_beta, V_beta)
    
    # Update phi (spatial effects)
    Q_phi <- Q * rho / tau2 + diag(1e-6, n)
    V_phi <- solve(Q_phi)
    m_phi <- V_phi %*% ((y - X %*% beta) * rho / tau2)
    phi <- MASS::mvrnorm(1, m_phi, V_phi)
    
    # Update tau2
    a_tau <- n/2 + 0.1
    b_tau <- sum((y - X %*% beta - phi)^2) / 2 + 0.1
    tau2 <- 1 / rgamma(1, a_tau, b_tau)
    
    # Update rho (Metropolis step)
    rho_prop <- runif(1, max(0, rho - 0.1), min(1, rho + 0.1))
    log_ratio <- -0.5 * t(phi) %*% Q %*% phi * (rho_prop - rho) / tau2
    if (log(runif(1)) < log_ratio) {
      rho <- rho_prop
    }
    
    # 저장
    if (iter > burn_in && (iter - burn_in) %% thin == 0) {
      iter_save <- iter_save + 1
      beta_samples[iter_save, ] <- beta
      tau2_samples[iter_save] <- tau2
      phi_samples[iter_save, ] <- phi
      rho_samples[iter_save] <- rho
    }
  }
  
  # 사후 요약
  beta_post <- colMeans(beta_samples)
  phi_post <- colMeans(phi_samples)
  
  # Credible intervals
  beta_CI <- apply(beta_samples, 2, quantile, probs = c(0.025, 0.975))
  phi_CI <- apply(phi_samples, 2, quantile, probs = c(0.025, 0.975))
  
  # DIC 계산
  D_bar <- mean(-2 * dnorm(y, X %*% beta_post + phi_post, sqrt(mean(tau2_samples)), log = TRUE))
  D_theta <- -2 * sum(dnorm(y, X %*% beta_post + phi_post, sqrt(mean(tau2_samples)), log = TRUE))
  pD <- D_bar - D_theta
  DIC <- D_bar + pD
  
  return(list(
    beta = beta_post,
    beta_CI = beta_CI,
    phi = phi_post,
    phi_CI = phi_CI,
    tau2 = mean(tau2_samples),
    rho = mean(rho_samples),
    DIC = DIC,
    samples = list(
      beta = beta_samples,
      tau2 = tau2_samples,
      phi = phi_samples,
      rho = rho_samples
    )
  ))
}

# === 5. 네트워크 기반 MPSA ==================================================

#' Network-based MPSA Analysis
#' 
#' @description proximity matrix를 네트워크로 변환하여 분석
#' @param P proximity matrix
#' @param threshold proximity threshold for edge creation
compute_network_MPSA <- function(P, coords, threshold = 0.7) {
  n <- nrow(P)
  
  # Proximity 네트워크 구성
  adj_matrix <- (P > threshold) * 1
  diag(adj_matrix) <- 0
  
  # igraph 객체 생성
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # 네트워크 지표 계산
  degree_cent <- degree(g)
  between_cent <- betweenness(g)
  eigen_cent <- eigen_centrality(g)$vector
  clustering <- transitivity(g, type = "local")
  
  # 공간 네트워크와의 비교
  dist_mat <- as.matrix(dist(coords))
  spatial_adj <- (dist_mat < quantile(dist_mat[upper.tri(dist_mat)], 0.1)) * 1
  diag(spatial_adj) <- 0
  g_spatial <- graph_from_adjacency_matrix(spatial_adj, mode = "undirected")
  
  # 네트워크 유사성 측정
  # Jaccard similarity of edges
  edges_rf <- which(adj_matrix == 1, arr.ind = TRUE)
  edges_spatial <- which(spatial_adj == 1, arr.ind = TRUE)
  
  edges_rf_set <- paste(edges_rf[,1], edges_rf[,2], sep = "-")
  edges_spatial_set <- paste(edges_spatial[,1], edges_spatial[,2], sep = "-")
  
  jaccard <- length(intersect(edges_rf_set, edges_spatial_set)) / 
             length(union(edges_rf_set, edges_spatial_set))
  
  # Community detection
  communities <- cluster_louvain(g)
  
  return(list(
    graph = g,
    degree = degree_cent,
    betweenness = between_cent,
    eigenvector = eigen_cent,
    clustering = clustering,
    communities = membership(communities),
    spatial_similarity = jaccard,
    adj_matrix = adj_matrix
  ))
}

# === 6. 고급 유의성 검정 ====================================================

#' Advanced Significance Testing for MPSA
#' 
#' @description 다양한 null hypothesis에 대한 검정
compute_local_mpsa_significance_advanced <- function(P, W, n_perm = 999, 
                                                  alpha = 0.05,
                                                  null_model = "CSR") {
  n <- nrow(P)
  observed_mpsa <- rowSums(W * P)
  
  # Null model에 따른 순열 생성
  perm_matrix <- matrix(NA, nrow = n, ncol = n_perm)
  
  if (null_model == "CSR") {
    # Complete Spatial Randomness
    for (b in 1:n_perm) {
      perm_idx <- sample(1:n)
      P_perm <- P[perm_idx, perm_idx]
      perm_matrix[, b] <- rowSums(W * P_perm)
    }
  } else if (null_model == "CAR") {
    # Conditional Autoregressive null
    for (b in 1:n_perm) {
      # CAR 프로세스에서 샘플링
      car_sample <- MASS::mvrnorm(1, rep(0, n), solve(diag(n) - 0.5 * W))
      perm_order <- order(car_sample)
      P_perm <- P[perm_order, perm_order]
      perm_matrix[, b] <- rowSums(W * P_perm)
    }
  } else if (null_model == "toroidal") {
    # Toroidal shift (edge effect 제거)
    for (b in 1:n_perm) {
      shift <- sample(1:n, 1)
      perm_idx <- ((1:n + shift - 1) %% n) + 1
      P_perm <- P[perm_idx, perm_idx]
      perm_matrix[, b] <- rowSums(W * P_perm)
    }
  }
  
  # Enhanced p-value calculation
  # 1. Standard p-value
  p_greater <- rowSums(perm_matrix >= observed_mpsa) / n_perm
  p_less <- rowSums(perm_matrix <= observed_mpsa) / n_perm
  p_two_sided <- 2 * pmin(p_greater, p_less)
  
  # 2. Exact p-value using empirical distribution
  p_exact <- sapply(1:n, function(i) {
    ecdf_i <- ecdf(perm_matrix[i, ])
    1 - ecdf_i(observed_mpsa[i]) + ecdf_i(observed_mpsa[i] - .Machine$double.eps)
  })
  
  # 3. Multiple testing correction (여러 방법 제공)
  p_bonferroni <- p.adjust(p_two_sided, method = "bonferroni")
  p_fdr <- p.adjust(p_two_sided, method = "fdr")
  p_BY <- p.adjust(p_two_sided, method = "BY")  # Benjamini-Yekutieli
  
  # Effect size calculation
  mean_perm <- rowMeans(perm_matrix)
  sd_perm <- apply(perm_matrix, 1, sd)
  effect_size <- (observed_mpsa - mean_perm) / sd_perm
  
  # Classification with confidence
  classification <- cut(effect_size,
                      breaks = c(-Inf, -2, -1, 1, 2, Inf),
                      labels = c("Strong Coldspot", "Coldspot", 
                               "Not Significant", "Hotspot", "Strong Hotspot"))
  
  # Confidence scores
  confidence <- 1 - p_fdr
  
  return(data.frame(
    MPSA = observed_mpsa,
    expected = mean_perm,
    variance = sd_perm^2,
    z_score = effect_size,
    p_value = p_two_sided,
    p_exact = p_exact,
    p_bonferroni = p_bonferroni,
    p_fdr = p_fdr,
    p_BY = p_BY,
    category = classification,
    confidence = confidence,
    n_neighbors = rowSums(W > 0)
  ))
}

# === 7. Variogram 분석 ======================================================

#' MPSA Variogram Analysis
#' 
#' @description MPSA의 공간적 의존성 구조 분석
compute_MPSA_variogram <- function(mpsa_values, coords, cutoff = NULL, width = NULL) {
  # 데이터 준비
  mpsa_df <- data.frame(
    mpsa = mpsa_values,
    x = coords[, 1],
    y = coords[, 2]
  )
  
  # sp 객체로 변환
  coordinates(mpsa_df) <- ~x + y
  
  # Empirical variogram
  if (is.null(cutoff)) {
    cutoff <- max(dist(coords)) / 3
  }
  if (is.null(width)) {
    width <- cutoff / 15
  }
  
  v_emp <- variogram(mpsa ~ 1, mpsa_df, cutoff = cutoff, width = width)
  
  # Fit variogram models
  v_fit_exp <- fit.variogram(v_emp, vgm("Exp"))
  v_fit_sph <- fit.variogram(v_emp, vgm("Sph"))
  v_fit_gau <- fit.variogram(v_emp, vgm("Gau"))
  v_fit_mat <- fit.variogram(v_emp, vgm("Mat"))
  
  # Model selection
  models <- list(
    Exponential = v_fit_exp,
    Spherical = v_fit_sph,
    Gaussian = v_fit_gau,
    Matern = v_fit_mat
  )
  
  # Calculate fit statistics
  fit_stats <- lapply(models, function(model) {
    pred <- variogramLine(model, dist_vector = v_emp$dist)
    obs <- v_emp$gamma
    pred_at_obs <- approx(pred$dist, pred$gamma, v_emp$dist)$y
    
    rmse <- sqrt(mean((obs - pred_at_obs)^2, na.rm = TRUE))
    r2 <- cor(obs, pred_at_obs, use = "complete.obs")^2
    
    return(c(RMSE = rmse, R2 = r2))
  })
  
  # Best model
  best_model_idx <- which.min(sapply(fit_stats, function(x) x["RMSE"]))
  best_model <- models[[best_model_idx]]
  
  return(list(
    empirical = v_emp,
    models = models,
    fit_stats = fit_stats,
    best_model = best_model,
    best_model_name = names(models)[best_model_idx],
    range = best_model$range[2],  # spatial range
    sill = sum(best_model$psill),  # total sill
    nugget = best_model$psill[1]   # nugget effect
  ))
}

# === 8. 고급 시각화 함수들 ==================================================

#' Create comprehensive MPSA visualization
#' 
#' @description MPSA 분석 결과의 종합적 시각화
visualize_MPSA_comprehensive <- function(franklin, mpsa_result, network_result,
                                       variogram_result, multiscale_result) {
  
  # 1. 기본 MPSA 맵
  p1 <- tm_shape(franklin) +
    tm_polygons("MPSA", 
                style = "jenks",
                palette = "RdBu",
                midpoint = 0,
                title = "MPSA Values") +
    tm_layout(title = "Local Indicators of Random Forest",
              frame = FALSE)
  
  # 2. 유의성 맵
  franklin$category <- mpsa_result$category
  p2 <- tm_shape(franklin) +
    tm_polygons("category",
                palette = c("Strong Coldspot" = "#2166AC",
                          "Coldspot" = "#67A9CF", 
                          "Not Significant" = "#F7F7F7",
                          "Hotspot" = "#EF8A62",
                          "Strong Hotspot" = "#B2182B"),
                title = "Significance") +
    tm_layout(title = "MPSA Significance Test",
              frame = FALSE)
  
  # 3. 네트워크 시각화
  # Extract coordinates for network plot
  coords_df <- st_coordinates(st_centroid(franklin))
  
  # Create edge list from adjacency matrix
  edges <- which(network_result$adj_matrix == 1, arr.ind = TRUE)
  edge_df <- data.frame(
    from = edges[, 1],
    to = edges[, 2]
  )
  
  # Create network layout
  g <- network_result$graph
  layout <- layout_with_fr(g, coords = coords_df)
  
  p3 <- ggraph(g, layout = layout) +
    geom_edge_link(alpha = 0.3) +
    geom_node_point(aes(size = network_result$degree,
                       color = network_result$clustering)) +
    scale_color_viridis(name = "Clustering\nCoefficient") +
    scale_size_continuous(name = "Degree") +
    theme_graph() +
    labs(title = "Proximity Network Structure")
  
  # 4. Variogram plot
  v_df <- variogram_result$empirical
  v_model <- variogramLine(variogram_result$best_model, 
                          maxdist = max(v_df$dist))
  
  p4 <- ggplot() +
    geom_point(data = v_df, aes(x = dist, y = gamma), size = 3) +
    geom_line(data = v_model, aes(x = dist, y = gamma), 
              color = "red", size = 1) +
    labs(x = "Distance", y = "Semivariance",
         title = paste("MPSA Variogram -", variogram_result$best_model_name),
         subtitle = sprintf("Range: %.0f, Sill: %.3f, Nugget: %.3f",
                          variogram_result$range,
                          variogram_result$sill,
                          variogram_result$nugget)) +
    theme_minimal()
  
  # 5. Multi-scale comparison
  scales <- names(multiscale_result)
  scale_df <- data.frame(
    scale = numeric(),
    mean_mpsa = numeric(),
    sd_mpsa = numeric(),
    moran_i = numeric()
  )
  
  for (s in scales) {
    scale_val <- multiscale_result[[s]]$scale
    mpsa_vals <- multiscale_result[[s]]$MPSA
    W <- multiscale_result[[s]]$W
    
    # Calculate Moran's I for comparison
    moran_test <- moran.test(mpsa_vals, mat2listw(W))
    
    scale_df <- rbind(scale_df, data.frame(
      scale = scale_val,
      mean_mpsa = mean(mpsa_vals),
      sd_mpsa = sd(mpsa_vals),
      moran_i = moran_test$statistic
    ))
  }
  
  p5 <- ggplot(scale_df, aes(x = scale)) +
    geom_line(aes(y = mean_mpsa), color = "blue", size = 1) +
    geom_ribbon(aes(ymin = mean_mpsa - sd_mpsa, 
                   ymax = mean_mpsa + sd_mpsa),
                alpha = 0.3, fill = "blue") +
    labs(x = "Spatial Scale (meters)", y = "Mean MPSA ± SD",
         title = "Multi-scale MPSA Analysis") +
    theme_minimal()
  
  # 6. Effect size distribution
  p6 <- ggplot(mpsa_result, aes(x = z_score, fill = category)) +
    geom_histogram(bins = 30, alpha = 0.7) +
    geom_vline(xintercept = c(-2, -1, 1, 2), 
               linetype = "dashed", alpha = 0.5) +
    scale_fill_manual(values = c("Strong Coldspot" = "#2166AC",
                               "Coldspot" = "#67A9CF", 
                               "Not Significant" = "#F7F7F7",
                               "Hotspot" = "#EF8A62",
                               "Strong Hotspot" = "#B2182B")) +
    labs(x = "Effect Size (Z-score)", y = "Frequency",
         title = "MPSA Effect Size Distribution") +
    theme_minimal()
  
  # Combine all plots
  combined_plot <- (p1 | p2) / (p3 | p4) / (p5 | p6)
  
  return(list(
    individual_plots = list(
      mpsa_map = p1,
      significance_map = p2,
      network = p3,
      variogram = p4,
      multiscale = p5,
      effect_sizes = p6
    ),
    combined = combined_plot
  ))
}

# === 9. 진단 및 검증 도구 ===================================================

#' MPSA Diagnostics
#' 
#' @description MPSA 분석의 가정 검증 및 진단
diagnose_MPSA <- function(P, W, mpsa_values, coords) {
  n <- nrow(P)
  diagnostics <- list()
  
  # 1. Spatial autocorrelation of residuals
  mpsa_expected <- rowMeans(replicate(100, {
    perm_idx <- sample(1:n)
    rowSums(W * P[perm_idx, perm_idx])
  }))
  residuals <- mpsa_values - mpsa_expected
  
  moran_resid <- moran.test(residuals, mat2listw(W))
  diagnostics$moran_residuals <- moran_resid
  
  # 2. Proximity matrix properties
  # Symmetry check
  diagnostics$proximity_symmetric <- all.equal(P, t(P))
  
  # Distribution of proximities
  prox_upper <- P[upper.tri(P)]
  diagnostics$proximity_stats <- summary(prox_upper)
  
  # 3. Spatial weight matrix diagnostics
  diagnostics$W_connected <- is.connected(graph_from_adjacency_matrix(W))
  diagnostics$W_row_sums <- summary(rowSums(W))
  
  # 4. Influence diagnostics
  # Leave-one-out analysis
  loo_mpsa <- matrix(NA, n, n)
  for (i in 1:n) {
    W_loo <- W[-i, -i]
    P_loo <- P[-i, -i]
    loo_mpsa[-i, i] <- rowSums(W_loo * P_loo)
  }
  
  # Cook's distance equivalent
  influence <- colMeans(abs(loo_mpsa - mpsa_values[-1]), na.rm = TRUE)
  diagnostics$influence <- influence
  diagnostics$influential_units <- which(influence > mean(influence) + 2*sd(influence))
  
  # 5. Bootstrap confidence intervals
  boot_mpsa <- replicate(1000, {
    boot_idx <- sample(1:n, replace = TRUE)
    P_boot <- P[boot_idx, boot_idx]
    W_boot <- W[boot_idx, boot_idx]
    rowSums(W_boot * P_boot)
  })
  
  diagnostics$bootstrap_CI <- t(apply(boot_mpsa, 1, quantile, 
                                     probs = c(0.025, 0.975)))
  
  return(diagnostics)
}

# === 10. 보고서 생성 함수 ===================================================

#' Generate MPSA Analysis Report
#' 
#' @description 분석 결과를 종합한 보고서 생성
generate_MPSA_report <- function(results, output_dir = "output/MPSA_analysis") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 1. Summary statistics
  summary_stats <- data.frame(
    Metric = c("Global IRF", "Mean Local IRF", "SD Local IRF",
               "% Hotspots", "% Coldspots", "Spatial Range",
               "Network Density", "Spatial-RF Similarity"),
    Value = c(
      results$GMPSA$GMPSA,
      mean(results$mpsa_result$MPSA),
      sd(results$mpsa_result$MPSA),
      mean(results$mpsa_result$category %in% c("Hotspot", "Strong Hotspot")) * 100,
      mean(results$mpsa_result$category %in% c("Coldspot", "Strong Coldspot")) * 100,
      results$variogram$range,
      edge_density(results$network$graph),
      results$network$spatial_similarity
    )
  )
  
  write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"))
  
  # 2. Detailed results
  write.csv(results$mpsa_result, file.path(output_dir, "mpsa_detailed_results.csv"))
  
  # 3. Save plots
  ggsave(file.path(output_dir, "comprehensive_visualization.png"), 
         results$plots$combined, width = 16, height = 20, dpi = 300)
  
  # 4. Generate markdown report
  report_content <- sprintf("
# MPSA Analysis Report
Generated: %s

## Executive Summary
- Global Indicator of Random Forest (GIRF): %.4f (p = %.4f)
- Identified %d hotspots and %d coldspots
- Spatial range of dependence: %.0f meters
- Best-fitting variogram model: %s

## Key Findings
1. **Spatial Clustering**: The analysis reveals significant spatial clustering with %.1f%% of units classified as hotspots or coldspots.
2. **Network Structure**: The proximity network shows %.1f%% similarity with the spatial network.
3. **Scale Dependency**: Multi-scale analysis indicates strongest patterns at %.0f meter scale.

## Methodological Notes
- Random Forest with %d trees
- Permutation tests with %d iterations
- Multiple testing correction: FDR (Benjamini-Hochberg)
- Bayesian inference with DIC = %.2f
",
    Sys.Date(),
    results$GMPSA$GMPSA, results$GMPSA$p_value,
    sum(results$mpsa_result$category %in% c("Hotspot", "Strong Hotspot")),
    sum(results$mpsa_result$category %in% c("Coldspot", "Strong Coldspot")),
    results$variogram$range,
    results$variogram$best_model_name,
    mean(results$mpsa_result$category != "Not Significant") * 100,
    results$network$spatial_similarity * 100,
    results$multiscale$optimal_scale,
    500,  # ntree
    999,  # n_perm
    results$bayesian$DIC
  )
  
  writeLines(report_content, file.path(output_dir, "analysis_report.md"))
  
  print(paste("Report generated in:", output_dir))
}

# === 9. 통합 분석 함수 ======================================================

#' Run Complete MPSA Analysis
#' 
#' @description 모든 MPSA 기능을 포함한 통합 분석
#' @param data sf 객체 (공간 데이터)
#' @param ntree Random Forest 트리 수
#' @param n_perm 순열 검정 횟수
#' @param scales 다중 스케일 분석용 거리 벡터
#' @param include_advanced 고급 분석 포함 여부
#' @return 완전한 MPSA 분석 결과
run_complete_MPSA_analysis <- function(data, ntree = 500, n_perm = 999, 
                                     scales = c(1000, 2500, 5000, 10000),
                                     include_advanced = TRUE) {
  
  cat("=== Complete MPSA Analysis ===\n")
  cat(sprintf("Data: %d spatial units\n", nrow(data)))
  
  # 1. 데이터 준비
  rf_data <- data |> 
    st_drop_geometry() |> 
    select(where(is.numeric)) |>
    select(-any_of(c("GEOID", "NAME")))  # ID 컬럼 제외
  
  coords <- st_coordinates(st_centroid(data))
  
  cat(sprintf("Variables: %d numeric variables\n", ncol(rf_data)))
  
  # 2. Random Forest 및 Proximity 계산
  cat("\nComputing Random Forest proximity matrix...\n")
  set.seed(42)
  rf_model <- randomForest(
    x = rf_data,
    proximity = TRUE,
    ntree = ntree,
    importance = TRUE
  )
  
  P <- rf_model$proximity
  
  # 3. 공간 가중치 행렬
  cat("Creating spatial weight matrices...\n")
  nb <- poly2nb(data, queen = TRUE)
  W_binary <- nb2mat(nb, style = "B", zero.policy = TRUE)
  W_rs <- nb2mat(nb, style = "W", zero.policy = TRUE)
  
  # 4. Global IRF 계산
  cat("Computing Global IRF...\n")
  girf_result <- compute_GMPSA(P, W_binary, n_perm = n_perm)
  cat(sprintf("GMPSA = %.4f (p = %.4f)\n", girf_result$GMPSA, girf_result$p_value))
  
  # 5. Local IRF 계산 및 유의성 검정
  cat("Computing Local IRF with significance testing...\n")
  mpsa_result <- compute_local_mpsa_significance_advanced(
    P, W_rs, n_perm = n_perm, null_model = "CSR"
  )
  
  # 결과 요약
  n_hotspots <- sum(mpsa_result$category %in% c("Hotspot", "Strong Hotspot"))
  n_coldspots <- sum(mpsa_result$category %in% c("Coldspot", "Strong Coldspot"))
  cat(sprintf("Results: %d hotspots, %d coldspots, %d not significant\n",
              n_hotspots, n_coldspots, 
              nrow(data) - n_hotspots - n_coldspots))
  
  # 6. 데이터에 결과 추가
  data$MPSA <- mpsa_result$MPSA
  data$MPSA_category <- mpsa_result$category
  data$MPSA_pvalue <- mpsa_result$p_fdr
  data$MPSA_zscore <- mpsa_result$z_score
  
  # 기본 결과
  basic_results <- list(
    data = data,
    rf_model = rf_model,
    proximity_matrix = P,
    spatial_weights = list(binary = W_binary, row_std = W_rs),
    GMPSA = girf_result,
    MPSA = mpsa_result,
    coords = coords
  )
  
  # 7. 고급 분석 (선택적)
  if (include_advanced) {
    cat("\nRunning advanced analyses...\n")
    
    # Variogram 분석
    cat("Computing variogram...\n")
    variogram_result <- compute_MPSA_variogram(
      mpsa_result$MPSA, coords
    )
    
    # 다중 스케일 분석
    cat("Multi-scale analysis...\n")
    multiscale_result <- compute_multiscale_MPSA(
      P, coords, scales = scales
    )
    
    # 네트워크 분석
    cat("Network analysis...\n")
    network_result <- analyze_proximity_network(
      P, coords, threshold = 0.6
    )
    
    # 고급 결과 추가
    basic_results$variogram <- variogram_result
    basic_results$multiscale <- multiscale_result
    basic_results$network <- network_result
    
    # 종합 시각화
    cat("Creating comprehensive visualization...\n")
    tryCatch({
      viz_result <- visualize_MPSA_comprehensive(
        data, mpsa_result, network_result, 
        variogram_result, multiscale_result
      )
      basic_results$plots <- viz_result
    }, error = function(e) {
      cat("Visualization error (continuing without plots):", e$message, "\n")
    })
  }
  
  cat("\n=== Analysis Complete ===\n")
  return(basic_results)
}

# === 메인 실행 예제 ==========================================================

# 이 부분은 주석 처리되어 있으며, 실제 사용 시 주석을 해제하세요
# franklin <- read_rds("data/franklin.rds")
# 
# # 전체 분석 파이프라인 실행
# results <- run_complete_MPSA_analysis(franklin)
# 
# # 보고서 생성
# generate_MPSA_report(results) 

# === 실행 예시 및 시각화 =====================================================

# === 1. 고급 MPSA 방법론 실행 ===
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
# # 베이지안 MPSA 실행
# bayesian_results <- compute_bayesian_MPSA(numeric_data, W_matrix)
# print("베이지안 MPSA 완료")
# 
# # 적응형 MPSA 실행  
# adaptive_results <- compute_adaptive_MPSA(numeric_data, W_matrix)
# print("적응형 MPSA 완료")

# === 2. 고급 방법론 비교 시각화 ===
# 
# # 기본 MPSA와 고급 방법론 비교
# library(ggplot2)
# library(patchwork)
# 
# # 기본 MPSA 계산
# rf_basic <- randomForest(numeric_data, proximity = TRUE, ntree = 500)
# mpsa_basic <- rowSums(W_matrix * rf_basic$proximity)
# 
# # 비교 데이터 준비
# comparison_df <- data.frame(
#   region_id = 1:length(mpsa_basic),
#   basic_mpsa = mpsa_basic,
#   bayesian_mpsa = bayesian_results$posterior_mean,
#   adaptive_mpsa = adaptive_results$adapted_mpsa
# )
# 
# # 방법론별 분포 비교
# methods_long <- comparison_df %>%
#   gather(method, value, -region_id) %>%
#   mutate(method = factor(method, 
#                         levels = c("basic_mpsa", "bayesian_mpsa", "adaptive_mpsa"),
#                         labels = c("Basic MPSA", "Bayesian MPSA", "Adaptive MPSA")))
# 
# p1 <- ggplot(methods_long, aes(x = value, fill = method)) +
#   geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
#   facet_wrap(~method, scales = "free_y", ncol = 1) +
#   labs(title = "Distribution Comparison of MPSA Methods",
#        x = "MPSA Value", y = "Frequency") +
#   theme_minimal() +
#   theme(legend.position = "none")

# === 3. 베이지안 MPSA 상세 시각화 ===
# 
# # 사후 분포 시각화
# posterior_samples <- bayesian_results$posterior_samples
# 
# # 사후 분포 히스토그램 (처음 몇 개 지역)
# sample_regions <- 1:6
# posterior_df <- data.frame()
# 
# for (region in sample_regions) {
#   region_samples <- data.frame(
#     region_id = paste("Region", region),
#     posterior_value = posterior_samples[region, ]
#   )
#   posterior_df <- rbind(posterior_df, region_samples)
# }
# 
# p2 <- ggplot(posterior_df, aes(x = posterior_value)) +
#   geom_histogram(bins = 50, alpha = 0.7, fill = "steelblue") +
#   facet_wrap(~region_id, scales = "free", ncol = 3) +
#   labs(title = "Posterior Distributions for Selected Regions",
#        x = "Posterior MPSA", y = "Frequency") +
#   theme_minimal()

# === 4. 적응형 MPSA 상세 분석 ===
# 
# # 적응형 가중치 시각화
# adaptive_weights <- adaptive_results$adaptive_weights
# 
# # 가중치 히트맵 (상위 20x20 부분)
# library(reshape2)
# weight_subset <- adaptive_weights[1:20, 1:20]
# weight_melted <- melt(weight_subset)
# 
# p4 <- ggplot(weight_melted, aes(x = Var1, y = Var2, fill = value)) +
#   geom_tile() +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
#                        midpoint = 0, name = "Weight") +
#   labs(title = "Adaptive Spatial Weights (First 20x20)",
#        x = "Region", y = "Region") +
#   theme_minimal() +
#   theme(axis.text = element_text(size = 8))
# 
# # 적응성 지표 시각화
# adaptivity_scores <- adaptive_results$adaptivity_scores
# 
# p5 <- ggplot(data.frame(region_id = 1:length(adaptivity_scores), 
#                        adaptivity = adaptivity_scores), 
#              aes(x = region_id, y = adaptivity)) +
#   geom_line(color = "darkgreen", size = 0.8) +
#   geom_point(color = "orange", size = 0.8) +
#   labs(title = "Regional Adaptivity Scores",
#        subtitle = "Higher scores indicate more adaptive spatial relationships",
#        x = "Region ID", y = "Adaptivity Score") +
#   theme_minimal()

# === 5. 방법론 간 상관관계 분석 ===
# 
# # 상관관계 매트릭스
# cor_matrix <- cor(comparison_df[, -1], use = "complete.obs")
# 
# # 상관관계 히트맵
# library(corrplot)
# png("output/advanced_mpsa/correlation_heatmap.png", width = 600, height = 600)
# corrplot(cor_matrix, method = "color", type = "upper", 
#          order = "hclust", tl.cex = 1, tl.col = "black",
#          title = "Correlation Between MPSA Methods")
# dev.off()
# 
# # 산점도 매트릭스
# library(GGally)
# p6 <- ggpairs(comparison_df[, -1], 
#               columnLabels = c("Basic MPSA", "Bayesian MPSA", "Adaptive MPSA"),
#               title = "Pairwise Comparisons of MPSA Methods") +
#   theme_minimal()

# === 6. 공간 패턴 비교 지도 ===
# 
# # 공간 데이터에 결과 추가
# franklin_advanced <- franklin
# franklin_advanced$basic_mpsa <- mpsa_basic
# franklin_advanced$bayesian_mpsa <- bayesian_results$posterior_mean
# franklin_advanced$adaptive_mpsa <- adaptive_results$adapted_mpsa
# franklin_advanced$uncertainty <- bayesian_results$posterior_sd
# 
# # 범주화
# franklin_advanced$basic_category <- cut(mpsa_basic, breaks = 5, labels = c("Very Low", "Low", "Medium", "High", "Very High"))
# franklin_advanced$bayesian_category <- cut(bayesian_results$posterior_mean, breaks = 5, labels = c("Very Low", "Low", "Medium", "High", "Very High"))
# franklin_advanced$adaptive_category <- cut(adaptive_results$adapted_mpsa, breaks = 5, labels = c("Very Low", "Low", "Medium", "High", "Very High"))
# 
# # tmap 지도들
# library(tmap)
# 
# # 기본 MPSA 지도
# map1 <- tm_shape(franklin_advanced) +
#   tm_fill("basic_category", title = "Basic MPSA", palette = "YlOrRd") +
#   tm_borders(alpha = 0.3) +
#   tm_layout(title = "Basic MPSA")
# 
# # 베이지안 MPSA 지도
# map2 <- tm_shape(franklin_advanced) +
#   tm_fill("bayesian_category", title = "Bayesian MPSA", palette = "YlOrRd") +
#   tm_borders(alpha = 0.3) +
#   tm_layout(title = "Bayesian MPSA")
# 
# # 적응형 MPSA 지도
# map3 <- tm_shape(franklin_advanced) +
#   tm_fill("adaptive_category", title = "Adaptive MPSA", palette = "YlOrRd") +
#   tm_borders(alpha = 0.3) +
#   tm_layout(title = "Adaptive MPSA")
# 
# # 불확실성 지도
# map4 <- tm_shape(franklin_advanced) +
#   tm_fill("uncertainty", title = "Uncertainty", palette = "Blues") +
#   tm_borders(alpha = 0.3) +
#   tm_layout(title = "Bayesian Uncertainty")
# 
# # 결합 지도
# combined_maps <- tmap_arrange(map1, map2, map3, map4, ncol = 2)

# === 7. 성능 벤치마크 ===
# 
# # 계산 시간 비교
# execution_times <- data.frame(
#   Method = character(),
#   Time_seconds = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# # 기본 MPSA 시간
# start_time <- Sys.time()
# rf_temp <- randomForest(numeric_data, proximity = TRUE, ntree = 500)
# basic_time <- as.numeric(Sys.time() - start_time)
# execution_times <- rbind(execution_times, data.frame(Method = "Basic MPSA", Time_seconds = basic_time))
# 
# # 베이지안 MPSA 시간
# start_time <- Sys.time()
# bayesian_temp <- compute_bayesian_MPSA(numeric_data[1:50, ], W_matrix[1:50, 1:50], n_iter = 100)
# bayesian_time <- as.numeric(Sys.time() - start_time)
# execution_times <- rbind(execution_times, data.frame(Method = "Bayesian MPSA (50 regions)", Time_seconds = bayesian_time))
# 
# # 적응형 MPSA 시간
# start_time <- Sys.time()
# adaptive_temp <- compute_adaptive_MPSA(numeric_data[1:50, ], W_matrix[1:50, 1:50])
# adaptive_time <- as.numeric(Sys.time() - start_time)
# execution_times <- rbind(execution_times, data.frame(Method = "Adaptive MPSA (50 regions)", Time_seconds = adaptive_time))
# 
# # 시간 비교 시각화
# p7 <- ggplot(execution_times, aes(x = reorder(Method, Time_seconds), y = Time_seconds)) +
#   geom_col(fill = "lightcoral", alpha = 0.8) +
#   coord_flip() +
#   labs(title = "Execution Time Comparison",
#        x = "Method", y = "Time (seconds)") +
#   theme_minimal()

# === 8. 통계적 비교 분석 ===
# 
# # 방법론별 기술통계량
# stats_comparison <- data.frame(
#   Method = c("Basic MPSA", "Bayesian MPSA", "Adaptive MPSA"),
#   Mean = c(mean(mpsa_basic), mean(bayesian_results$posterior_mean), mean(adaptive_results$adapted_mpsa)),
#   SD = c(sd(mpsa_basic), sd(bayesian_results$posterior_mean), sd(adaptive_results$adapted_mpsa)),
#   Min = c(min(mpsa_basic), min(bayesian_results$posterior_mean), min(adaptive_results$adapted_mpsa)),
#   Max = c(max(mpsa_basic), max(bayesian_results$posterior_mean), max(adaptive_results$adapted_mpsa)),
#   Median = c(median(mpsa_basic), median(bayesian_results$posterior_mean), median(adaptive_results$adapted_mpsa))
# )
# 
# # 차이점 분석
# differences_df <- data.frame(
#   region_id = 1:length(mpsa_basic),
#   bayesian_vs_basic = bayesian_results$posterior_mean - mpsa_basic,
#   adaptive_vs_basic = adaptive_results$adapted_mpsa - mpsa_basic,
#   adaptive_vs_bayesian = adaptive_results$adapted_mpsa - bayesian_results$posterior_mean
# )
# 
# # 차이점 시각화
# differences_long <- differences_df %>%
#   gather(comparison, difference, -region_id) %>%
#   mutate(comparison = factor(comparison,
#                             levels = c("bayesian_vs_basic", "adaptive_vs_basic", "adaptive_vs_bayesian"),
#                             labels = c("Bayesian vs Basic", "Adaptive vs Basic", "Adaptive vs Bayesian")))
# 
# p8 <- ggplot(differences_long, aes(x = difference)) +
#   geom_histogram(bins = 30, alpha = 0.7, fill = "lightgreen") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#   facet_wrap(~comparison, scales = "free") +
#   labs(title = "Differences Between MPSA Methods",
#        x = "Difference", y = "Frequency") +
#   theme_minimal()

# === 9. 결과 저장 및 보고서 ===
# 
# # 결과 디렉토리 생성
# if (!dir.exists("output/advanced_mpsa")) {
#   dir.create("output/advanced_mpsa", recursive = TRUE)
# }
# 
# # 모든 시각화 저장
# combined_plot1 <- (p1 + p2) / p3
# combined_plot2 <- (p4 + p5) / p7
# combined_plot3 <- p8
# 
# ggsave("output/advanced_mpsa/method_distributions.png", p1, width = 12, height = 10)
# ggsave("output/advanced_mpsa/posterior_distributions.png", p2, width = 12, height = 8)
# ggsave("output/advanced_mpsa/bayesian_uncertainty.png", p3, width = 12, height = 6)
# ggsave("output/advanced_mpsa/adaptive_weights.png", p4, width = 8, height = 6)
# ggsave("output/advanced_mpsa/adaptivity_scores.png", p5, width = 10, height = 6)
# ggsave("output/advanced_mpsa/pairwise_comparisons.png", p6, width = 12, height = 10)
# ggsave("output/advanced_mpsa/execution_times.png", p7, width = 10, height = 6)
# ggsave("output/advanced_mpsa/method_differences.png", p8, width = 12, height = 8)
# 
# # 종합 시각화 저장
# ggsave("output/advanced_mpsa/comprehensive_analysis_1.png", combined_plot1, width = 16, height = 12)
# ggsave("output/advanced_mpsa/comprehensive_analysis_2.png", combined_plot2, width = 16, height = 12)
# 
# # 공간 지도 저장
# tmap_save(combined_maps, "output/advanced_mpsa/spatial_comparison.png")
# 
# # 수치 결과 저장
# write.csv(comparison_df, "output/advanced_mpsa/method_comparison_data.csv", row.names = FALSE)
# write.csv(uncertainty_df, "output/advanced_mpsa/bayesian_uncertainty_data.csv", row.names = FALSE)
# write.csv(stats_comparison, "output/advanced_mpsa/statistical_comparison.csv", row.names = FALSE)
# write.csv(differences_df, "output/advanced_mpsa/method_differences_data.csv", row.names = FALSE)
# write.csv(execution_times, "output/advanced_mpsa/execution_times.csv", row.names = FALSE)
# 
# # 공간 데이터 저장
# st_write(franklin_advanced, "output/advanced_mpsa/franklin_advanced_results.shp", delete_dsn = TRUE)

# === 10. 최종 종합 평가 ===
# 
# # 각 방법론의 장단점 요약
# method_evaluation <- data.frame(
#   Method = c("Basic MPSA", "Bayesian MPSA", "Adaptive MPSA"),
#   Computational_Speed = c("Fast", "Slow", "Medium"),
#   Uncertainty_Quantification = c("No", "Yes", "Limited"),
#   Spatial_Adaptivity = c("No", "No", "Yes"),
#   Interpretability = c("High", "Medium", "Medium"),
#   Best_Use_Case = c(
#     "Large datasets, quick analysis",
#     "Small datasets, uncertainty important",
#     "Complex spatial patterns"
#   )
# )
# 
# # 성능 지표 요약
# performance_summary <- list(
#   basic_mpsa = list(
#     mean = round(mean(mpsa_basic), 4),
#     sd = round(sd(mpsa_basic), 4),
#     range = round(range(mpsa_basic), 4)
#   ),
#   bayesian_mpsa = list(
#     mean = round(mean(bayesian_results$posterior_mean), 4),
#     sd = round(sd(bayesian_results$posterior_mean), 4),
#     range = round(range(bayesian_results$posterior_mean), 4),
#     avg_uncertainty = round(mean(bayesian_results$posterior_sd), 4)
#   ),
#   adaptive_mpsa = list(
#     mean = round(mean(adaptive_results$adapted_mpsa), 4),
#     sd = round(sd(adaptive_results$adapted_mpsa), 4),
#     range = round(range(adaptive_results$adapted_mpsa), 4),
#     avg_adaptivity = round(mean(adaptive_results$adaptivity_scores), 4)
#   )
# )
# 
# # 상관관계 요약
# correlations_summary <- data.frame(
#   Comparison = c("Basic vs Bayesian", "Basic vs Adaptive", "Bayesian vs Adaptive"),
#   Correlation = c(
#     round(cor(mpsa_basic, bayesian_results$posterior_mean), 3),
#     round(cor(mpsa_basic, adaptive_results$adapted_mpsa), 3),
#     round(cor(bayesian_results$posterior_mean, adaptive_results$adapted_mpsa), 3)
#   )
# )
# 
# # 종합 보고서 저장
# library(jsonlite)
# final_report <- list(
#   method_evaluation = method_evaluation,
#   performance_summary = performance_summary,
#   correlations_summary = correlations_summary,
#   execution_times = execution_times
# )
# 
# write_json(final_report, "output/advanced_mpsa/comprehensive_evaluation.json", pretty = TRUE)
# write.csv(method_evaluation, "output/advanced_mpsa/method_evaluation.csv", row.names = FALSE)
# write.csv(correlations_summary, "output/advanced_mpsa/correlations_summary.csv", row.names = FALSE)
# 
# # 최종 요약 출력
# cat("=== 고급 MPSA 방법론 분석 완료 ===\n")
# cat("📊 성능 비교:\n")
# for (i in 1:nrow(stats_comparison)) {
#   cat(sprintf("  - %s: 평균 %.4f (±%.4f)\n", 
#               stats_comparison$Method[i], 
#               stats_comparison$Mean[i], 
#               stats_comparison$SD[i]))
# }
# cat("🔗 상관관계:\n")
# for (i in 1:nrow(correlations_summary)) {
#   cat(sprintf("  - %s: %.3f\n", 
#               correlations_summary$Comparison[i], 
#               correlations_summary$Correlation[i]))
# }
# cat("⏱️  계산 효율성:\n")
# for (i in 1:nrow(execution_times)) {
#   cat(sprintf("  - %s: %.2f초\n", 
#               execution_times$Method[i], 
#               execution_times$Time_seconds[i]))
# }
# cat("📁 결과 저장 위치: output/advanced_mpsa/\n")
# cat("🎯 고급 방법론 비교 완료!\n") 