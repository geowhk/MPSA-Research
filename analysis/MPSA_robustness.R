# MPSA의 robustness 탐구

source("data_preparation/setup.R")

franklin <- readRDS("data/franklin.rds")
nb <- franklin |> 
  poly2nb(queen = TRUE)
nb_self <- include.self(nb)
W <- nb2mat(nb, style = "W", zero.policy = TRUE)

data <- franklin |>
  st_drop_geometry() |> 
  dplyr::select(where(is.numeric), main_industry)


# 0. 논문용 지도 제작 ------------------------------------------------------------
# 0-1. 3개 변수의 Moran Significance Map

# 분석할 변수 목록
vars <- c("med_income", "pct_renter", "unemployment_rate")

# 각 변수에 대해 Local Moran's I 및 hotspot/coldspot 분류
for (var in vars) {
  x <- scale(data[[var]])[, 1]  # 정규화
  
  # Local Moran's I 계산
  local_moran <- localmoran(x, listw = mat2listw(W, style = "W"), zero.policy = TRUE)
  Ii <- local_moran[, "Ii"]
  pval <- local_moran[, "Pr(z != E(Ii))"]
  
  # 유의 수준 0.05 기준으로 hotspot/coldspot 분류
  category <- rep("Not Significant", length(pval))
  category[x >= 0 & Ii > 0 & pval < 0.05] <- "Hotspot"
  category[x <  0 & Ii > 0 & pval < 0.05] <- "Coldspot"
  category[x >= 0 & Ii < 0 & pval < 0.05] <- "Spatial Outlier High-Low"
  category[x <  0 & Ii < 0 & pval < 0.05] <- "Spatial Outlier Low-High"
  
  # 변수명_category 형식으로 결과 저장
  franklin[[paste0(var, "_category")]] <- factor(category,
                                             levels = c("Hotspot", "Coldspot",
                                                        "Spatial Outlier High-Low",
                                                        "Spatial Outlier Low-High",
                                                        "Not Significant"))
}

# Local MPSA 계산
franklin_result <- local_mpsa(data, W, ntree = 500)
franklin$category <- franklin_result$category


# 시각화
tm1 <- tm_shape(franklin) +
  tm_polygons(
    fill = "med_income_category",
    fill.scale = tm_scale_categorical(values = c(
      "Hotspot" = "#d73027",
      "Spatial Outlier High-Low" = "#fc8d59",
      "Spatial Outlier Low-High" = "lightgreen",
      "Coldspot" = "#4575b4",
      "Not Significant" = "lightgray"
    )
    ),
    col = "black",
    fill.legend = tm_legend(title = "Types",
                            position = c(0.825, 0.25),
                            bg.color = "white"
    )
  ) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white") +
  tm_credits("(A)", position = c("left", "top"), size = 3)

tm2 <- tm_shape(franklin) +
  tm_polygons(
    fill = "pct_renter_category",
    fill.scale = tm_scale_categorical(values = c(
      "Hotspot" = "#d73027",
      "Spatial Outlier High-Low" = "#fc8d59",
      "Spatial Outlier Low-High" = "lightgreen",
      "Coldspot" = "#4575b4",
      "Not Significant" = "lightgray"
    )
    ),
    col = "black",
    fill.legend = tm_legend(title = "Types",
                            position = c(0.825, 0.25),
                            bg.color = "white"
    )
  ) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white") +
  tm_credits("(B)", position = c("left", "top"), size = 3)

tm3 <- tm_shape(franklin) +
  tm_polygons(
    fill = "unemployment_rate_category",
    fill.scale = tm_scale_categorical(values = c(
      "Hotspot" = "#d73027",
      "Spatial Outlier High-Low" = "#fc8d59",
      "Spatial Outlier Low-High" = "lightgreen",
      "Coldspot" = "#4575b4",
      "Not Significant" = "lightgray"
    )
    ),
    col = "black",
    fill.legend = tm_legend(title = "Types",
                            position = c(0.825, 0.25),
                            bg.color = "white"
    )
  ) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white") +
  tm_credits("(C)", position = c("left", "top"), size = 3)

tm4 <- tm_shape(franklin) +
  tm_polygons(
    fill = "category",
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
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white") +
  tm_credits("(D)", position = c("left", "top"), size = 3)

tmap_arrange(tm1, tm2, tm3, tm4, ncol = 2, nrow = 2)

tmap_save(tm1,
          filename = "output/robustness/med_income.png",
          width = 8, height = 7, units = "in", dpi = 300)
tmap_save(tm2,
          filename = "output/robustness/pct_renter.png",
          width = 8, height = 7, units = "in", dpi = 300)
tmap_save(tm3,
          filename = "output/robustness/unemp.png",
          width = 8, height = 7, units = "in", dpi = 300)
tmap_save(tm4,
          filename = "output/robustness/local_mpsa.png",
          width = 8, height = 7, units = "in", dpi = 300)
tmap_save(tmap_arrange(tm1, tm2, tm3, tm4, ncol = 2, nrow = 2),
          filename = "output/robustness/data_description.png",
          width = 17, height = 15, units = "in", dpi = 300)


# 1. Seed 미고정 반복 실험 -------------------------------------------------------
# 1-1. Global MPSA
# 각 ntree에 따라 시드 고정 없이 100번 MPSA 계산
source("functions/global_MPSA.R")

# 설정
ntree_values <- c(50, 100, 500, 1000)
n_iter <- 100

# 결과 저장용 데이터프레임 초기화
mpsa_results <- data.frame()

# 반복 실행
for (ntree in ntree_values) {
  
  for (i in 1:n_iter) {
    cat("ntree =", ntree, "/", i, "\n")
    mpsa_value <- global_mpsa(data = data, W = W, ntree = ntree)
    mpsa_results <- rbind(mpsa_results, data.frame(
      ntree = ntree,
      iteration = i,
      MPSA = mpsa_value$observed_mpsa,
      p_value = mpsa_value$p_value
    ))
  }
}


# 요약 확인
summary_stats <- mpsa_results |> 
  group_by(ntree) |> 
  summarise(mean_MPSA = mean(MPSA),
            sd_MPSA = sd(MPSA),
            cv_MPSA = sd(MPSA) / mean(MPSA))

print(summary_stats)

mpsa_results$ntree <- factor(mpsa_results$ntree)

write_rds(summary_stats, "output/robustness/global_mpsa_summary.rds")
write_rds(mpsa_results, "output/robustness/global_mpsa_robustnee.rds")

# 시각화: x축에 ntree 값을 둔 boxplot
global_mpsa_robustness <- ggplot(mpsa_results, aes(x = ntree, y = MPSA)) +
  geom_boxplot(fill = "lightblue", outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.3, color = "black", size = 1) +
  labs(x = "ntree",
       y = "MPSA") +
  theme_bw(base_size = 13)
global_mpsa_robustness
ggsave("output/robustness/global_mpsa_boxplot.png", global_mpsa_robustness, dpi = 900)

# 1-2. Local MPSA
source("functions/local_MPSA.R")

# 설정
ntree_values <- c(50, 100, 500, 1000)
n_iter <- 100
n_region <- nrow(data)
local_mpsa_array <- array(NA, dim = c(n_region, n_iter, length(ntree_values)),
                          dimnames = list(paste0("region_", 1:n_region), NULL, paste0("ntree_", ntree_values)))

# 반복 실행 (시드 미고정)
for (k in seq_along(ntree_values)) {
  ntree <- ntree_values[k]
  
  for (i in 1:n_iter) {
    cat("ntree =", ntree, "/", i, "\n")  # 실행 진행 상황 출력
    result <- local_mpsa(data, W, ntree = ntree, seed = NULL)  # 시드 고정 X
    local_mpsa_array[, i, k] <- result$local_mpsa
  }
}

# CV 계산
cv_df <- data.frame(region_id = 1:n_region)

for (k in seq_along(ntree_values)) {
  ntree <- ntree_values[k]
  mpsa_mat <- local_mpsa_array[, , k]
  mean_vals <- rowMeans(mpsa_mat)
  sd_vals <- apply(mpsa_mat, 1, sd)
  cv_vals <- sd_vals / abs(mean_vals)  # CV 정의
  cv_df[[paste0("CV_ntree_", ntree)]] <- cv_vals
}

write_rds(cv_df, "output/local_mpsa_cv.rds")

franklin$cv_50 <- cv_df$CV_ntree_50
franklin$cv_100 <- cv_df$CV_ntree_100
franklin$cv_500 <- cv_df$CV_ntree_500
franklin$cv_1000 <- cv_df$CV_ntree_1000

local_tm1 <- tm_shape(franklin) +
  tm_polygons(
    fill = "cv_50",
    fill.scale = tm_scale_continuous(
      n = 5,
      limits = c(0.03, 0.6)
    ),
    col = "black",
    fill.legend = tm_legend(title = "Local MPSA CV",
                            orientation = "portrait",
                            title.size = 1.2,
                            text.size = 1,
                            position = c(0.83, 0.3),
                            item.height = 1,
                            frame = FALSE,
                            bg = FALSE
    )
  ) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white") +
  tm_credits("(A)", position = c("left", "top"), size = 3)

local_tm2 <- tm_shape(franklin) +
  tm_polygons(
    fill = "cv_100",
    fill.scale = tm_scale_continuous(
      n = 5,
      limits = c(0.03, 0.6)
    ),
    col = "black",
    fill.legend = tm_legend(title = "Local MPSA CV",
                            orientation = "portrait",
                            title.size = 1.2,
                            text.size = 1,
                            position = c(0.83, 0.3),
                            item.height = 1,
                            frame = FALSE,
                            bg = FALSE
    )
  ) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white") +
  tm_credits("(B)", position = c("left", "top"), size = 3)

local_tm3 <- tm_shape(franklin) +
  tm_polygons(
    fill = "cv_500",
    fill.scale = tm_scale_continuous(
      n = 5,
      limits = c(0.03, 0.6)
    ),
    col = "black",
    fill.legend = tm_legend(title = "Local MPSA CV",
                            orientation = "portrait",
                            title.size = 1.2,
                            text.size = 1,
                            position = c(0.83, 0.3),
                            item.height = 1,
                            frame = FALSE,
                            bg = FALSE
    )
  ) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white") +
  tm_credits("(C)", position = c("left", "top"), size = 3)

local_tm4 <- tm_shape(franklin) +
  tm_polygons(
    fill = "cv_1000",
    fill.scale = tm_scale_continuous(
      n = 5,
      limits = c(0.03, 0.6)
    ),
    col = "black",
    fill.legend = tm_legend(title = "Local MPSA CV",
                            orientation = "portrait",
                            title.size = 1.2,
                            text.size = 1,
                            position = c(0.83, 0.3),
                            item.height = 1,
                            frame = FALSE,
                            bg = FALSE
    )
  ) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white") +
  tm_credits("(D)", position = c("left", "top"), size = 3)
local_tm4

tmap_save(tmap_arrange(local_tm1, local_tm2, local_tm3, local_tm4,
                       ncol = 2, 
                       nrow = 2),
          filename = "output/robustness/local_mpsa_cv.png",
          width = 23, height = 15, units = "in", dpi = 300)

summary(franklin$cv_50)
summary(franklin$cv_100)
summary(franklin$cv_500)
summary(franklin$cv_1000)


# 2. Franklin 변수별 Moran's I 계산 --------------------------------------------
source("functions/global_MPSA.R")
compare_variablewise_moransI_and_mpsa <- function(data, W, ntree = 500, n_perm = 999, seed = NULL) {
  library(spdep)
  library(randomForest)
  
  if (!is.null(seed)) set.seed(seed)
  
  listw <- mat2listw(W, style = "W")
  
  numeric_vars <- data[, sapply(data, is.numeric)]
  var_names <- colnames(numeric_vars)
  
  morans_I_values <- numeric(length(var_names))
  morans_pvals <- numeric(length(var_names))
  mpsa_values <- numeric(length(var_names))
  mpsa_pvals <- numeric(length(var_names))
  
  for (i in seq_along(var_names)) {
    var <- numeric_vars[[i]]
    
    if (any(is.na(var))) {
      morans_I_values[i] <- NA
      morans_pvals[i] <- NA
      mpsa_values[i] <- NA
      mpsa_pvals[i] <- NA
      next
    }
    
    # Moran's I (순열 기반)
    moran_result <- moran.mc(x = var, listw = listw, nsim = n_perm)
    morans_I_values[i] <- moran_result$statistic
    morans_pvals[i] <- moran_result$p.value
    
    # Global MPSA
    rf_model <- randomForest(x = data.frame(var), proximity = TRUE, ntree = ntree, importance = FALSE)
    P <- rf_model$proximity
    
    n <- length(var)
    P_centered <- P - mean(P)
    W_sum <- sum(W)
    denom <- sum(P_centered^2)
    observed_mpsa <- (n / W_sum) * (sum(W * P_centered) / denom)
    
    # Global MPSA 순열 검정
    upper_idx <- which(upper.tri(P), arr.ind = TRUE)
    permuted_mpsa <- numeric(n_perm)
    
    for (b in 1:n_perm) {
      P_perm <- matrix(0, n, n)
      vals <- sample(P[upper_idx], replace = FALSE)
      P_perm[upper_idx] <- vals
      P_perm <- P_perm + t(P_perm)
      diag(P_perm) <- 1
      P_perm_centered <- P_perm - mean(P_perm)
      permuted_mpsa[b] <- (n / W_sum) * (sum(W * P_perm_centered) / sum(P_perm_centered^2))
    }
    
    p_greater <- mean(permuted_mpsa >= observed_mpsa)
    p_less <- mean(permuted_mpsa <= observed_mpsa)
    p_value <- 2 * min(p_greater, p_less)
    
    mpsa_values[i] <- observed_mpsa
    mpsa_pvals[i] <- p_value
  }
  
  result <- data.frame(
    variable = var_names,
    morans_I = morans_I_values,
    morans_p = morans_pvals,
    global_mpsa = mpsa_values,
    mpsa_p = mpsa_pvals
  )
  
  return(result)
}


result <- compare_variablewise_moransI_and_mpsa(data, W)
result <- result |> 
  dplyr::filter(variable != "main_industry")
ggplot(result, aes(x = morans_I, y = global_mpsa)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Moran's I vs. Global MPSA",
       x = "Moran's I",
       y = "Global MPSA") +
  theme_bw()
cor(result$morans_I, result$global_mpsa, use = "complete.obs", method = "spearman")
result
