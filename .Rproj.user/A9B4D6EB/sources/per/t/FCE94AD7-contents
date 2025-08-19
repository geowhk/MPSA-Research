global_mpsa <- function(data, W, ntree = 500, n_perm = 999, seed = NULL) {
  library(randomForest)
  
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(data)
  W_sum <- sum(W)
  
  # 1. Random Forest proximity matrix
  synthetic <- as.data.frame(lapply(data, function(col) sample(col, replace = TRUE)))
  data_all <- rbind(data, synthetic)
  y_all <- factor(c(rep("real", nrow(data)), rep("synthetic", nrow(synthetic))))
  
  rf <- randomForest(x = data_all,
                     y = y_all,
                     proximity = TRUE,
                     ntree = 500)
  P <- rf$proximity[1:nrow(x), 1:nrow(x)]
  
  # 2. 중심화
  P_mean <- mean(P)
  P_centered <- P - P_mean
  denom <- sum(P_centered^2)
  
  # 3. 정규화된 Global MPSA 계산
  observed_mpsa <- (n / W_sum) * (sum(W * P_centered) / denom)
  
  # 4. 조건부 순열: 상삼각 성분만 섞고 대칭화
  permuted_mpsa <- numeric(n_perm)
  upper_idx <- which(upper.tri(P), arr.ind = TRUE)
  
  for (b in 1:n_perm) {
    P_perm <- matrix(0, n, n)
    vals <- sample(P[upper_idx], replace = FALSE)
    P_perm[upper_idx] <- vals
    P_perm <- P_perm + t(P_perm)
    diag(P_perm) <- 1
    
    P_perm_centered <- P_perm - mean(P_perm)
    permuted_mpsa[b] <- (n / W_sum) * (sum(W * P_perm_centered) / sum(P_perm_centered^2))
  }
  
  # 5. 양측검정 p-value 및 z-score
  p_greater <- mean(permuted_mpsa >= observed_mpsa)
  p_less <- mean(permuted_mpsa <= observed_mpsa)
  p_value <- 2 * min(p_greater, p_less)
  
  z_score <- (observed_mpsa - mean(permuted_mpsa)) / sd(permuted_mpsa)
  
  return(list(
    observed_mpsa = observed_mpsa,
    p_value = p_value,
    z_score = z_score,
    null_distribution = permuted_mpsa
  ))
}
