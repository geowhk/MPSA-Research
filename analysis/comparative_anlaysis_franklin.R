source("data_preparation/setup.R")

# 데이터 불러오기 ----------------------------------------------------------------

franklin <- readRDS("data/franklin.rds") |> 
  select(where(is.numeric), -main_industry) # categorical data 제외
  
# 공간 가중치 행렬 생성(자기 자신 포함X, 행표준화 행렬로 통일)
nb <- franklin |> 
  poly2nb(queen = TRUE)
nb_self <- include.self(nb)
W <- nb2mat(nb, style = "W", zero.policy = TRUE) # 자기 자신 이웃 포함, 행표준화

data <- franklin |> 
  st_drop_geometry()

# functions for local SAS calculation -------------------------------------
source("functions/pca_local_moran.R")
source("functions/multivariate_local_geary.R")
source("functions/local_mahalanobis.R")
source("functions/local_MPSA.R")

pca_result <- pca_local_moran(data, W)
geary_result <- multivariate_local_geary(data, W)
mahalanobis_result <- mahalanobis_local(data, W)
MPSA_result <- local_mpsa(data, W, ntree = 500)

# hotspot visualization ---------------------------------------------------

franklin$pca_category <- pca_result$category
franklin$geary_category <- geary_result$category
franklin$mahalanobis_category <- mahalanobis_result$category
franklin$mpsa_category <- MPSA_result$category

hotspot_pca <- tm_shape(franklin) +
  tm_polygons(
    fill = "pca_category",
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

hotspot_geary <- tm_shape(franklin) +
  tm_polygons(
    fill = "geary_category",
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

hotspot_mahalanobis <- tm_shape(franklin) +
  tm_polygons(
    fill = "mahalanobis_category",
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
  tm_credits("(C)", position = c("left", "top"), size = 3)

hotspot_mpsa <- tm_shape(franklin) +
  tm_polygons(
    fill = "mpsa_category",
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

tmap_save(tm = hotspot_pca,
          "output/comparative_analysis(franklin)/hotspot_pca.png",
          dpi = 900)
tmap_save(tm = hotspot_geary,
          "output/comparative_analysis(franklin)/hotspot_geary.png",
          dpi = 900)
tmap_save(tm = hotspot_mahalanobis,
          "output/comparative_analysis(franklin)/hotspot_mahalanobis.png",
          dpi = 900)
tmap_save(tm = hotspot_mpsa,
          "output/comparative_analysis(franklin)/hotspot_mpsa.png",
          dpi = 900)
tmap_save(tmap_arrange(hotspot_pca, hotspot_geary, hotspot_mahalanobis, hotspot_mpsa,
                       ncol = 2, 
                       nrow = 2),
          filename = "output/comparative_analysis(franklin)/comparative_analysis.png",
          width = 23, height = 15, units = "in", dpi = 300)
