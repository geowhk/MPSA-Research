# --- 환경 설정 ---
source("R/data_preparation/setup.R")


# 데이터 준비 ------------------------------------------------------------------

franklin <- readRDS("data/franklin.rds")

# P matrix
x <- franklin |> 
  st_drop_geometry() |> 
  dplyr::select(where(is.numeric), main_industry)

synthetic <- as.data.frame(lapply(x, function(col) sample(col, replace = TRUE)))
data_all <- rbind(x, synthetic)
y_all <- factor(c(rep("real", nrow(x)), rep("synthetic", nrow(synthetic))))

rf <- randomForest(x = data_all,
                   y = y_all,
                   proximity = TRUE,
                   ntree = 500)
P <- rf$proximity[1:nrow(x), 1:nrow(x)]

# W matrix
nb <- poly2nb(franklin, queen = TRUE)
nb_self <- include.self(nb)

W <- nb2mat(nb, style = "B", zero.policy = TRUE) # 자기 자신 이웃 포함 X
W_rs <- nb2mat(nb, style = "W", zero.policy = TRUE) # 자기 자신 이웃 포함 X, 행표준화
V <- nb2mat(nb_self, style = "B", zero.policy = TRUE) # 자기 자신 이웃 포함
V_rs <- nb2mat(nb_self, style = "W", zero.policy = TRUE) # 자기 자신 이웃 포함, 행표준화


# Global MPSA(Standardized) -------------------------------------------------------------

gmpsa <- function(P, W) {
  P_mean <- mean(P)
  P_centered <- P - P_mean
  
  W_sum <- sum(W)
  
  gmpsa <- (nrow(P) / W_sum) * sum(W * P_centered) / sum(P_centered^2)
  
  return(gmpsa)
}

gmpsa_W <- gmpsa(P, W)
gmpsa_W_rs <- gmpsa(P, W_rs)
gmpsa_V <- gmpsa(P, V)
gmpsa_V_rs <- gmpsa(P, V_rs)

# significance test
mpsa_permutation_test <- function(P, W, n_perm = 999, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(P)
  
  # 중심화 함수
  center_matrix <- function(P) {
    P - mean(P)
  }
  
  # 전역 MPSA 계산 함수
  compute_global_mpsa <- function(P_centered, W) {
    W_sum <- sum(W)
    numerator <- sum(W * P_centered)
    denominator <- sum(P_centered^2)
    (n / W_sum) * (numerator / denominator)
  }
  
  # 순열 함수: 상삼각만 섞고 대칭 복원, diag = 1 유지
  permute_P_matrix <- function(P) {
    P_perm <- matrix(0, n, n)
    upper_idx <- which(upper.tri(P))
    permuted_vals <- sample(P[upper_idx])
    P_perm[upper_idx] <- permuted_vals
    P_perm <- P_perm + t(P_perm)
    diag(P_perm) <- 1
    return(P_perm)
  }
  
  # 관측값 계산
  P_centered_obs <- center_matrix(P)
  observed_mpsa <- compute_global_mpsa(P_centered_obs, W)
  
  # 순열 분포 계산
  permuted_mpsas <- numeric(n_perm)
  for (i in 1:n_perm) {
    P_perm <- permute_P_matrix(P)
    P_centered_perm <- center_matrix(P_perm)
    permuted_mpsas[i] <- compute_global_mpsa(P_centered_perm, W)
  }
  
  # 단측 (우측) p-value
  p_value <- (sum(permuted_mpsas >= observed_mpsa) + 1) / (n_perm + 1)
  
  # 결과 반환
  list(
    observed_mpsa = observed_mpsa,
    permuted_mpsas = permuted_mpsas,
    p_value = p_value
  )
}

W_result <- mpsa_permutation_test(P, W)
W_rs_result <- mpsa_permutation_test(P, W_rs)
V_result <- mpsa_permutation_test(P, V)
V_rs_result <- mpsa_permutation_test(P, V_rs)

cat("MPSA(W)", W_result$observed_mpsa, "/ p-value(W)", W_result$p_value, "\n",
    "MPSA(W_rs)", W_rs_result$observed_mpsa, "/ p-value(W_rs)", W_rs_result$p_value, "\n",
    "MPSA(V)", V_result$observed_mpsa, "/ p-value(V)", V_result$p_value, "\n",
    "MPSA(V_rs)", V_rs_result$observed_mpsa, "/ p-value(V_rs)", V_rs_result$p_value, "\n")

plot_mpsa_histogram <- function(result, buffer = 0.01) {
  
  # 추출
  perm <- result$permuted_mpsas
  obs <- result$observed_mpsa
  
  # 데이터 프레임 생성
  df <- data.frame(mpsa = perm)
  
  # x축 범위 설정 (관측값 포함해서 여유 공간)
  x_max <- max(c(perm, obs)) + buffer
  x_min <- min(perm)
  
  # 시각화
  ggplot(df, aes(x = mpsa)) +
    geom_histogram(
      bins = 50,
      fill = "lightgray",
      color = "white"
    ) +
    geom_density(color = "gray", linewidth = 1) +
    geom_vline(
      xintercept = obs,
      color = "red",
      linewidth = 1
    ) +
    annotate(
      "text",
      x = obs,
      y = max(table(cut(perm, 30))) * 0.95,
      label = paste0("Observed = ", round(obs, 4)),
      hjust = -0.1,
      color = "red"
    ) +
    labs(
      title = "Global MPSA Permutation Distribution",
      x = "MPSA",
      y = "Frequency"
    ) +
    coord_cartesian(xlim = c(x_min, x_max)) +
    theme_minimal(base_size = 14)
}
global_W_plot <- plot_mpsa_histogram(W_result)
global_W_rs_plot <- plot_mpsa_histogram(W_rs_result)
global_V_plot <- plot_mpsa_histogram(V_result)
global_V_rs_plot <- plot_mpsa_histogram(V_rs_result)

ggsave("output/global(standardized)/global_W_plot.png", plot = global_W_plot, dpi = 900)
ggsave("output/global(standardized)/global_W_rs_plot.png", plot = global_W_rs_plot, dpi = 900)
ggsave("output/global(standardized)/global_V_plot.png", plot = global_V_plot, dpi = 900)
ggsave("output/global(standardized)/global_V_rs_plot.png", plot = global_V_rs_plot, dpi = 900)

# Local MPSA(Standardized) ------------------------------------------------

local_mpsa <- function(P, W) {
  n <- nrow(P)
  
  # 중심화
  P_mean <- mean(P)
  P_centered <- P - P_mean
  
  # 분모: 전체 분산
  denom <- sum(P_centered^2)
  
  # 공간 가중치 총합
  W_sum <- sum(W)
  
  # 지역별 LISA 계산
  local_vals <- numeric(n)
  for (i in 1:n) {
    local_vals[i] <- sum(W[i, ] * P_centered[i, ])
  }
  
  # 정규화 계수: n^2 / sum(W)
  local_mpsa_vals <- (n^2 / W_sum) * local_vals / denom
  return(local_mpsa_vals)
}

local_mpsa(P, W)

# significance test
local_mpsa_significance_test <- function(P, W, n_perm = 999, alpha = 0.05) {
  n <- nrow(P)
  
  # 1. Local MPSA 계산 (정규화된 버전)
  P_mean <- mean(P)
  P_centered <- P - P_mean
  local_mpsa <- (n^2) * rowSums(W * P_centered) / sum(P_centered^2)
  
  # 2. 조건부 순열을 통한 perm_matrix 생성
  perm_matrix <- matrix(NA, nrow = n, ncol = n_perm)
  upper_idx <- which(upper.tri(P), arr.ind = TRUE)
  
  for (b in 1:n_perm) {
    vals <- sample(P[upper_idx], replace = FALSE)
    P_perm <- matrix(0, n, n)
    P_perm[upper_idx] <- vals
    P_perm <- P_perm + t(P_perm)
    diag(P_perm) <- 1
    
    P_perm_centered <- P_perm - mean(P_perm)
    perm_matrix[, b] <- (n^2) * rowSums(W * P_perm_centered) / sum(P_perm_centered^2)
  }
  
  # 3. p값 계산 (양측검정)
  p_values <- sapply(1:n, function(i) {
    obs <- local_mpsa[i]
    null_dist <- perm_matrix[i, ]
    p_greater <- mean(null_dist >= obs)
    p_less <- mean(null_dist <= obs)
    p_two_sided <- 2 * min(p_greater, p_less)
    return(p_two_sided)
  })
  
  # 4. 다중검정 보정
  p_adjusted <- p.adjust(p_values, method = "fdr")
  
  # 5. 효과 크기 기반 분류
  expected_mpsa <- apply(perm_matrix, 1, mean)
  sd_mpsa <- apply(perm_matrix, 1, sd)
  effect_size <- (local_mpsa - expected_mpsa) / sd_mpsa
  
  # 6. 5단계 분류 체계
  category <- cut(effect_size,
                  breaks = c(-Inf, -2, -1, 1, 2, Inf),
                  labels = c("Strong Coldspot", "Coldspot", "Not Significant", 
                             "Hotspot", "Strong Hotspot"),
                  include.lowest = TRUE)
  
  # 7. 유의성 필터링
  final_category <- as.character(category)
  final_category[p_adjusted >= alpha] <- "Not Significant"
  
  # 결과 반환
  return(list(
    local_mpsa = local_mpsa,
    p_value = p_values,
    p_adjusted = p_adjusted,
    effect_size = effect_size,
    category = final_category
  ))
}

local_W_result <- local_mpsa_significance_test(P, W)
local_W_rs_result <- local_mpsa_significance_test(P, W_rs)
local_V_result <- local_mpsa_significance_test(P, V)
local_V_rs_result <- local_mpsa_significance_test(P, V_rs)

franklin$W_category <- local_W_result$category
franklin$W_rs_category <- local_W_rs_result$category
franklin$V_category <- local_V_result$category
franklin$V_rs_category <- local_V_rs_result$category


local_W_map <- tm_shape(franklin) +
  tm_polygons(
    fill = "W_category",
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
  tm_title("Local MPSA Results(diag = 0, binary)") +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white")

tmap_save(tm = local_W_map, "output/local(standardized)/local_W_map.png", dpi = 900)

local_W_rs_map <- tm_shape(franklin) +
  tm_polygons(
    fill = "W_rs_category",
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
  tm_title("Local MPSA Results(diag = 0, row-standardized)") +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white")

tmap_save(tm = local_W_rs_map, "output/local(standardized)/local_W_rs_map.png", dpi = 900)

local_V_map <- tm_shape(franklin) +
  tm_polygons(
    fill = "V_category",
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
  tm_title("Local MPSA Results(diag = 1, binary)") +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white")

tmap_save(tm = local_V_map, "output/local(standardized)/local_V_map.png", dpi = 900)

local_V_rs_map <- tm_shape(franklin) +
  tm_polygons(
    fill = "V_rs_category",
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
  tm_title("Local MPSA Results(diag = 1, row-standardized)") +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white")

tmap_save(tm = local_V_rs_map, "output/local(standardized)/local_V_rs_map.png", dpi = 900)
