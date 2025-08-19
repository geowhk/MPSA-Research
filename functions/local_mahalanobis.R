mahalanobis_local <- function(data, W, alpha = 0.05) {
  n <- nrow(data)
  p <- ncol(data)
  
  # 데이터 정규화
  data_scaled <- scale(data)
  
  # 공분산 행렬 및 역행렬 계산
  cov_matrix <- cov(data_scaled)
  cov_inv <- tryCatch({
    solve(cov_matrix)
  }, error = function(e) {
    svd_result <- svd(cov_matrix)
    tolerance <- 1e-10
    d_inv <- ifelse(svd_result$d > tolerance, 1 / svd_result$d, 0)
    svd_result$v %*% diag(d_inv) %*% t(svd_result$u)
  })
  
  # 통계량 및 유의확률
  local_stats <- numeric(n)
  local_p_values <- numeric(n)
  
  for (i in 1:n) {
    neighbors <- which(W[i, ] > 0)
    
    if (length(neighbors) > 0) {
      weights <- W[i, neighbors] / sum(W[i, neighbors])
      neighbor_mean <- colSums(t(t(data_scaled[neighbors, , drop = FALSE]) * weights))
      
      diff <- data_scaled[i, ] - neighbor_mean
      stat <- as.numeric(t(diff) %*% cov_inv %*% diff)
      local_stats[i] <- stat
      local_p_values[i] <- 1 - pchisq(stat, df = p)
    } else {
      local_stats[i] <- 0
      local_p_values[i] <- 1
    }
  }
  
  # 효과 크기 기반 정규화: Z-score
  expected <- mean(local_stats)
  sd_stat <- sd(local_stats)
  effect_size <- (local_stats - expected) / sd_stat
  
  # 5단계 분류 체계
  category <- cut(effect_size,
                  breaks = c(-Inf, -2, -1, 1, 2, Inf),
                  labels = c("Strong Coldspot", "Coldspot", "Not Significant",
                             "Hotspot", "Strong Hotspot"),
                  include.lowest = TRUE)
  
  # 유의하지 않은 지역은 Not Significant로 통일
  final_category <- as.character(category)
  final_category[local_p_values >= alpha] <- "Not Significant"
  
  result <- list(
    local_stat = local_stats,
    p_value = local_p_values,
    effect_size = effect_size,
    category = factor(final_category, levels = c("Strong Coldspot", "Coldspot", "Not Significant",
                                                 "Hotspot", "Strong Hotspot"))
  )
  
  return(result)
}