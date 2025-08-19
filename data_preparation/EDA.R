# --- 환경 설정 ---
source("R/data_preparation/setup.R")

# --- 데이터 불러오기 ---
franklin <- read_rds("data/franklin.rds")

# 기초 통계량 확인
franklin <- franklin |> 
  st_drop_geometry() |> 
  select(-GEOID, -NAME)

describe(franklin)

desc_df <- describe(franklin)

write_xlsx(desc_df, "output/기초통계량.xlsx")

# 변수 간 산점도 행렬
ggpairs <- ggpairs(franklin, 
        upper = list(continuous = wrap("smooth", method = "loess", se = FALSE)),
        lower = list(continuous = wrap("points", alpha = 0.5)))

ggsave("output/ggpairs_plot.png", plot = ggpairs, dpi = 900)

# 변수 간 상관행렬
png("corrplot_plot.png", width = 800, height = 800, res = 150)
corrplot(cor(franklin), method = "color", tl.cex = 0.8)
dev.off()

# PCA
pca <- prcomp(franklin, scale. = TRUE)
autoplot(pca, data = franklin, loadings = TRUE, loadings.label = TRUE)
