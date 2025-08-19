local_mpsa <- function(data, W, ntree, alpha = 0.05, n_perm = 999, seed = 42) {
  library(randomForest)
  
  set.seed(seed)
  n <- nrow(data)
  
  # 1. RF proximity matrix
  synthetic <- as.data.frame(lapply(data, function(col) sample(col, replace = TRUE)))
  data_all <- rbind(data, synthetic)
  y_all <- factor(c(rep("real", nrow(data)), rep("synthetic", nrow(synthetic))))
  
  rf <- randomForest(x = data_all,
                     y = y_all,
                     proximity = TRUE,
                     ntree = 500)
  P <- rf$proximity[1:nrow(x), 1:nrow(x)]
  
  # 2. 정규화된 Local MPSA (row-wise inner product normalized)
  P_mean <- mean(P)
  P_centered <- P - P_mean
  W_sum <- sum(W)
  
  observed <- (n^2 / W_sum) * rowSums(W * P_centered) / sum(P_centered^2)
  
  # 3. 조건부 순열을 통한 귀무분포 생성
  perm_matrix <- matrix(NA, nrow = n, ncol = n_perm)
  upper_idx <- which(upper.tri(P), arr.ind = TRUE)
  
  for (b in 1:n_perm) {
    P_perm <- matrix(0, n, n)
    vals <- sample(P[upper_idx], replace = FALSE)
    P_perm[upper_idx] <- vals
    P_perm <- P_perm + t(P_perm)
    diag(P_perm) <- 1
    
    P_perm_centered <- P_perm - mean(P_perm)
    perm_matrix[, b] <- (n^2 / W_sum) * rowSums(W * P_perm_centered) / sum(P_perm_centered^2)
  }
  
  # 4. p-value 계산
  p_values <- sapply(1:n, function(i) {
    null_dist <- perm_matrix[i, ]
    obs <- observed[i]
    p_greater <- mean(null_dist >= obs)
    p_less <- mean(null_dist <= obs)
    2 * min(p_greater, p_less)
  })
  
  # 5. FDR 보정 및 효과 크기
  p_adjusted <- p.adjust(p_values, method = "fdr")
  expected <- apply(perm_matrix, 1, mean)
  sd_vals <- apply(perm_matrix, 1, sd)
  effect_size <- (observed - expected) / sd_vals
  
  # 6. 5단계 분류
  category <- cut(effect_size,
                  breaks = c(-Inf, -2, -1, 1, 2, Inf),
                  labels = c("Strong Coldspot", "Coldspot", "Not Significant", 
                             "Hotspot", "Strong Hotspot"),
                  include.lowest = TRUE)
  category[p_adjusted >= alpha] <- "Not Significant"
  
  # 결과 반환
  result <- list(
    local_mpsa = observed,
    p_value = p_values,
    p_adjusted = p_adjusted,
    effect_size = effect_size,
    category = factor(category, 
                      levels = c("Strong Coldspot", "Coldspot", "Not Significant",
                                 "Hotspot", "Strong Hotspot"))
  )
  
  return(result)
}
