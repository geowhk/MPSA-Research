# 격자 데이터를 활용한 시뮬레이션 연구
# Anselin(1995)의 SAR 기반 데이터 생성 방식 사용

library(spdep)
library(sf)
library(tidyverse)
library(tmap)

# 함수 정의(SAR 및 격자 생성) ------------------------------------------------------

simulate_sar <- function(nrow = 10, ncol = 10, rho = 0.7, epsilon = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow * ncol
  nb <- cell2nb(nrow, ncol, type = "queen")
  W <- nb2mat(nb, style = "W")
  
  # 외부에서 ε를 주지 않았다면 생성
  if (is.null(epsilon)) {
    epsilon <- rnorm(n)
  } else {
    if (length(epsilon) != n) stop("Length of epsilon must be equal to nrow * ncol")
  }
  
  y <- solve(diag(n) - rho * W) %*% epsilon
  return(as.numeric(y))
}


create_grid_sf_multi <- function(nrow = 10, ncol = 10, values_list = list()) {
  require(sf)
  
  bbox <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = ncol, ymax = nrow)))
  grid <- st_make_grid(x = bbox, n = c(ncol, nrow), what = "polygons", square = TRUE)
  grid_sf <- st_sf(geometry = grid)
  
  if (length(values_list) > 0) {
    for (i in seq_along(values_list)) {
      grid_sf[[paste0("var", i)]] <- values_list[[i]]
    }
  }
  
  return(grid_sf)
}


# 시뮬레이션용 데이터 생성 -----------------------------------------------------------

# Simulation 1: 모든 변수에서 동일한 공간 자기상관 존재
s1_var1 <- simulate_sar(rho = 0.7, seed = 1)
s1_var2 <- simulate_sar(rho = 0.7, seed = 2)
s1_var3 <- simulate_sar(rho = 0.7, seed = 3)

grid_s1 <- create_grid_sf_multi(values_list = list(s1_var1, s1_var2, s1_var3))
qtm(grid_s1, fill = "var1")
qtm(grid_s1, fill = "var2")
qtm(grid_s1, fill = "var3")


# Simulation 2: 모든 변수가 서로 다른 공간적 자기상관 정도
shared_epsilon <- rnorm(100)
s2_var1 <- simulate_sar(rho = 0.3, epsilon = shared_epsilon)
s2_var2 <- simulate_sar(rho = 0.7, epsilon = shared_epsilon)
s2_var3 <- simulate_sar(rho = 0.9, epsilon = shared_epsilon)

grid_s2 <- create_grid_sf_multi(values_list = list(s2_var1, s2_var2, s2_var3))
qtm(grid_s2, fill = "var1")
qtm(grid_s2, fill = "var2")
qtm(grid_s2, fill = "var3")

# Simulation 3: 하나의 변수만 공간적 자기상관이 나타나고, 나머지 변수는 랜덤
s3_var1 <- simulate_sar(rho = 0.7, seed = 123)
s3_var2 <- rnorm(100)
s3_var3 <- rnorm(100)

grid_s3 <- create_grid_sf_multi(values_list = list(s3_var1, s3_var2, s3_var3))
qtm(grid_s3, fill = "var1")
qtm(grid_s3, fill = "var2")
qtm(grid_s3, fill = "var3")


# Simulation 4: 변수들 간 선형 조합
# 원 SAR 변수 생성
set.seed(123)
base1 <- simulate_sar(rho = 0.7, seed = 1)
base2 <- simulate_sar(rho = 0.7, seed = 2)
base3 <- simulate_sar(rho = 0.7, seed = 3)

# 선형 조합된 3개 변수
s4_var1 <- 0.6 * base1 + 0.4 * base2
s4_var2 <- 0.5 * base2 + 0.5 * base3
s4_var3 <- 0.3 * base1 + 0.3 * base2 + 0.4 * base3

# 다변량 grid 생성
grid_s4 <- create_grid_sf_multi(
  nrow = 10, ncol = 10,
  values_list = list(s4_var1, s4_var2, s4_var3)
)

# Simulation 5: 변수들 간 비선형 조합
set.seed(456)
base1 <- simulate_sar(rho = 0.7, seed = 11)
base2 <- simulate_sar(rho = 0.7, seed = 22)
base3 <- simulate_sar(rho = 0.7, seed = 33)

# 비선형 조합된 3개 변수
s5_var1 <- sin(base1) + log(abs(base2) + 1)
s5_var2 <- base2^2 - sqrt(abs(base3))
s5_var3 <- exp(-abs(base1 - base3))

# 다변량 grid 생성
grid_s5 <- create_grid_sf_multi(
  nrow = 10, ncol = 10,
  values_list = list(s5_var1, s5_var2, s5_var3)
)


# 시뮬레이션 수행 ----------------------------------------------------------------

source("functions/local_MPSA.R")
source("functions/pca_local_moran.R")
source("functions/multivariate_local_geary.R")
source("functions/local_mahalanobis.R")

nb = poly2nb(grid_s1, queen = TRUE)
W = nb2mat(nb, style = "W")

plot_comparison_facets <- function(sf_data, W, alpha = 0.05, ntree = 100, n_perm = 999) {
  data <- sf_data %>% st_drop_geometry() %>% select(starts_with("var"))
  
  # 1. 각 통계량 계산
  res_mpsa <- local_mpsa(data, W, ntree = ntree, alpha = alpha, n_perm = n_perm)
  res_pca <- pca_local_moran(data, W, alpha = alpha)
  res_geary <- multivariate_local_geary(data, W, alpha = alpha, n_perm = n_perm)
  res_maha <- mahalanobis_local(data, W, alpha = alpha)
  
  # 2. 결과 병합
  category_df <- tibble(
    Local_MPSA = res_mpsa$category,
    PCA_Moran = res_pca$category,
    Local_Geary = res_geary$category,
    Mahalanobis = res_maha$category
  ) %>%
    mutate(id = row_number()) %>%
    pivot_longer(-id, names_to = "Method", values_to = "Category")
  
  sf_data$id <- 1:nrow(sf_data)
  plot_data <- left_join(sf_data, category_df, by = "id")
  
  # 3. 시각화
  p <- ggplot(plot_data) +
    geom_sf(aes(fill = Category), color = "white", size = 0.1) +
    facet_wrap(~Method, ncol = 2) +
    scale_fill_manual(values = c(
      "Strong Coldspot" = "#0571b0",
      "Coldspot" = "#92c5de",
      "Not Significant" = "grey80",
      "Hotspot" = "#f4a582",
      "Strong Hotspot" = "#ca0020"
    )) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(fill = "Category", title = "Local Spatial Association by Method")
  
  return(p)
}

plot_comparison_facets(grid_s1, W)
plot_comparison_facets(grid_s2, W)
plot_comparison_facets(grid_s3, W)
plot_comparison_facets(grid_s4, W)
plot_comparison_facets(grid_s5, W)
