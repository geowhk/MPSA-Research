# 하이퍼파라미터에 따른 P 행렬의 분포 탐구

library(randomForest)
library(sf)
library(tidyverse)

make_p_matrix <- function(x, maxnode) {
  synthetic <- as.data.frame(lapply(x, function(col) sample(col, replace = TRUE)))
  data_all <- rbind(x, synthetic)
  y_all <- factor(c(rep("real", nrow(x)), rep("synthetic", nrow(synthetic))))
  
  rf <- randomForest(x = data_all,
                     y = y_all,
                     proximity = TRUE,
                     maxnode = maxnode)
  P <- rf$proximity[1:nrow(x), 1:nrow(x)]
  
  return(P)
}

franklin <- readRDS("data/franklin.rds")

x <- franklin |> 
  st_drop_geometry() |> 
  dplyr::select(where(is.numeric), main_industry)

P <- make_p_matrix(x, 3)
hist(P)



# maxnode에 따른 P행렬 분포 시각화 --------------------------------------------------

maxnode <- 50
n_iter <- 100
zero_counts <- c()
bin_breaks <- seq(0, 1, by = 0.1)
bin_labels <- paste0("[", head(bin_breaks, -1), ",", tail(bin_breaks, -1), ")")

# 반복 실행 결과 저장
hist_data <- data.frame()

for (i in 1:n_iter) {
  P <- make_p_matrix(x, maxnode = maxnode)
  pvec <- P[lower.tri(P)]
  
  # 히스토그램 bin 집계
  bin_count <- as.data.frame(table(cut(pvec, breaks = bin_breaks, right = FALSE)))
  names(bin_count) <- c("Bin", "Count")
  bin_count$Iteration <- i
  
  hist_data <- rbind(hist_data, bin_count)
  
  zero_counts[i] <- sum(pvec == 0)
}

mean_zero <- round(mean(zero_counts))
mean_zero_percent <- round(mean(zero_counts) / length(pvec) * 100, 2)

# Bin 레벨 정렬
hist_data$Bin <- factor(hist_data$Bin, levels = bin_labels)

# boxplot 시각화
p <- ggplot(hist_data, aes(x = Bin, y = Count)) +
  geom_boxplot(fill = "skyblue", outlier.alpha = 0.3) +
  theme_bw() +
  labs(title = paste0("Distribution of Proximity Histogram Bin Frequencies across Iterations (maxnode =", maxnode, ")" ),
       subtitle = paste0("Average number of zero elements in P: ", mean_zero, " (", mean_zero_percent, "%)"),
       x = "Proximity Bin", y = "Frequency (Count)")
ggsave("output/Pmatrix/maxnode50.png", p, dpi = 900)
