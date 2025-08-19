
#' Compute Anselin(2019) Local Multivariate Geary's C
#'
#' @param data A numeric matrix or data frame with rows as regions and columns as variables
#' @param W Spatial weight matrix (n x n), row-standardized
#' @param alpha Significance level (default: 0.05)
#'
#' @return A list with components: local_stat, p_value, adjusted_p_value, category
multivariate_local_geary <- function(data, W, alpha = 0.05, n_perm = 999) {
  n <- nrow(data)
  data <- scale(data)
  g_values <- numeric(n)

  for (i in 1:n) {
    neighbors <- which(W[i, ] > 0)
    if (length(neighbors) > 0) {
      diffs <- sweep(data[neighbors, , drop = FALSE], 2, data[i, ])
      g_values[i] <- sum(W[i, neighbors] * rowSums(diffs^2))
    } else {
      g_values[i] <- NA
    }
  }

  # Permutation test
  perm_matrix <- matrix(NA, nrow = n, ncol = n_perm)
  for (b in 1:n_perm) {
    data_perm <- data[sample(1:n), ]
    for (i in 1:n) {
      neighbors <- which(W[i, ] > 0)
      if (length(neighbors) > 0) {
        diffs <- sweep(data_perm[neighbors, , drop = FALSE], 2, data_perm[i, ])
        perm_matrix[i, b] <- sum(W[i, neighbors] * rowSums(diffs^2))
      } else {
        perm_matrix[i, b] <- NA
      }
    }
  }

  p_values <- sapply(1:n, function(i) {
    null_dist <- perm_matrix[i, ]
    obs <- g_values[i]
    if (is.na(obs)) return(1)
    p_greater <- mean(null_dist >= obs, na.rm = TRUE)
    p_less <- mean(null_dist <= obs, na.rm = TRUE)
    2 * min(p_greater, p_less)
  })

  # Adjust for multiple testing
  p_adjusted <- p.adjust(p_values, method = "fdr")

  # Effect size: standardized to mean and SD of null distribution
  expected <- apply(perm_matrix, 1, mean, na.rm = TRUE)
  sd_val <- apply(perm_matrix, 1, sd, na.rm = TRUE)
  effect_size <- (g_values - expected) / sd_val

  category <- cut(effect_size,
                  breaks = c(-Inf, -2, -1, 1, 2, Inf),
                  labels = c("Strong Coldspot", "Coldspot", "Not Significant", 
                             "Hotspot", "Strong Hotspot"),
                  include.lowest = TRUE)
  final_category <- as.character(category)
  final_category[p_adjusted >= alpha] <- "Not Significant"

  return(data.frame(
    region_id = 1:n,
    local_stat = g_values,
    p_value = p_values,
    adjusted_p_value = p_adjusted,
    category = final_category
  ))
}
