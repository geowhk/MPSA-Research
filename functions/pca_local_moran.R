# 다변량 데이터에 대해 PCA를 적용한 뒤 첫 번째 주성분에 대해 Local Moran's I 계산
pca_local_moran <- function(data, W, alpha = 0.05){
  library(spdep)
  
  # 데이터 표준화 및 PCA
  data_scaled <- scale(data)
  pc1 <- prcomp(data_scaled)$x[, 1]
  
  # Local Moran's I 계산
  lm_result <- localmoran(pc1, listw = mat2listw(W, style = "W"))
  Ii <- lm_result[, 1]
  pval <- lm_result[, 5]
  
  # 효과크기 기반 분류
  category <- cut(Ii,
                  breaks = c(-Inf, -2, -1, 1, 2, Inf),
                  labels = c("Strong Coldspot", "Coldspot", "Not Significant", 
                             "Hotspot", "Strong Hotspot"),
                  include.lowest = TRUE)
  
  category[pval >= alpha] <- "Not Significant"
  
  return(list(category = category))
}