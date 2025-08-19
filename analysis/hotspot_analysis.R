# 핫스팟과 주변 지역 비교 -----------------------------------------------------------

franklin_scaled <- franklin |> 
  mutate(across(where(is.numeric), scale))

# 1. 관심 지역 선택
target <- franklin_scaled |> filter(GEOID == "39049000730")  # 예시 GEOID

# 2. 주변 지역 선택 (인접하는 1차 이웃)
neighbors <- franklin_scaled |> 
  filter(st_touches(geometry, target$geometry, sparse = FALSE)[, 1])

# 3. 변수 선택
all_vars <- franklin_scaled |> 
  st_drop_geometry() |> 
  dplyr::select(where(is.numeric)) |> 
  colnames()

# 4. group 부여 및 결합
df_compare <- bind_rows(
  target |> mutate(group = "Target"),
  neighbors |> mutate(group = "Neighbor")
) |> 
  st_drop_geometry() |> 
  dplyr::select(group, all_of(all_vars))

# 5. 평균값 요약
summary_stats <- df_compare |> 
  group_by(group) |> 
  summarise(across(everything(), mean, na.rm = TRUE)) |> 
  pivot_longer(-group, names_to = "variable", values_to = "value") |> 
  pivot_wider(names_from = group, values_from = value)

# ✅ 6. Z-score 표준화
summary_scaled <- summary_stats

# 7. fmsb 패키지용 구조: max, min, rows
radar_data <- summary_scaled[, c("Target", "Neighbor")] |> 
  as.data.frame()
rownames(radar_data) <- summary_stats$variable

max_val <- max(radar_data, na.rm = TRUE)
min_val <- min(radar_data, na.rm = TRUE)

radar_input <- rbind(
  max = rep(max_val, nrow(radar_data)),
  min = rep(min_val, nrow(radar_data)),
  t(radar_data)
)
colnames(radar_input) <- summary_stats$variable

# 8. 레이더 차트 시각화
fmsb::radarchart(
  as.data.frame(radar_input),
  axistype = 1,
  axislabcol = NA,
  pcol = c("red", "blue"),
  plwd = 2,
  plty = 1,
  title = "Target vs Neighbor: Variable Profile (Z-score standardized)",
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  vlcex = 0.8
)
legend("bottomright", legend = c("Target", "Neighbor"), col = c("red", "blue"), lty = 1, lwd = 2)


tmap_mode("plot")

hotspot_analysis <- tm_shape(franklin_results) +
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
    fill_alpha = 0.3,
    fill.legend = tm_legend(title = "MPSA Category",
                            position = c(0.85, 0.25),
                            bg.color = "white"
    )
  ) +
  tm_title("Local MPSA Results") +
  tm_shape(neighbors) +
  tm_polygons(fill = "#d73027",
              col = "black",
              lwd = 2) +
  tm_shape(target) +
  tm_polygons(fill = "#d73027",
              col = "black",
              lwd = 2) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white")

tmap_save(hotspot_analysis, "output/hotspot_analysis.png", dpi = 900)


# 유의x 지역 주변 ---------------------------------------------------------------

# 1. 관심 지역 선택
target <- franklin_scaled |> filter(GEOID == "39049007114")  # 예시 GEOID

# 2. 주변 지역 선택 (인접하는 1차 이웃)
neighbors <- franklin_scaled |> 
  filter(st_touches(geometry, target$geometry, sparse = FALSE)[, 1])

# 3. 변수 선택
all_vars <- franklin_scaled |> 
  st_drop_geometry() |> 
  dplyr::select(where(is.numeric)) |> 
  colnames()

# 4. group 부여 및 결합
df_compare <- bind_rows(
  target |> mutate(group = "Target"),
  neighbors |> mutate(group = "Neighbor")
) |> 
  st_drop_geometry() |> 
  dplyr::select(group, all_of(all_vars))

# 5. 평균값 요약
summary_stats <- df_compare |> 
  group_by(group) |> 
  summarise(across(everything(), mean, na.rm = TRUE)) |> 
  pivot_longer(-group, names_to = "variable", values_to = "value") |> 
  pivot_wider(names_from = group, values_from = value)

# ✅ 6. Z-score 표준화
summary_scaled <- summary_stats

# 7. fmsb 패키지용 구조: max, min, rows
radar_data <- summary_scaled[, c("Target", "Neighbor")] |> 
  as.data.frame()
rownames(radar_data) <- summary_stats$variable

max_val <- max(radar_data, na.rm = TRUE)
min_val <- min(radar_data, na.rm = TRUE)

radar_input <- rbind(
  max = rep(max_val, nrow(radar_data)),
  min = rep(min_val, nrow(radar_data)),
  t(radar_data)
)
colnames(radar_input) <- summary_stats$variable

# 8. 레이더 차트 시각화
fmsb::radarchart(
  as.data.frame(radar_input),
  axistype = 1,
  axislabcol = NA,
  pcol = c("red", "blue"),
  plwd = 2,
  plty = 1,
  title = "Target vs Neighbor: Variable Profile (Z-score standardized)",
  cglcol = "grey",
  cglty = 1,
  axislabcol = "grey",
  vlcex = 0.8
)
legend("bottomright", legend = c("Target", "Neighbor"), col = c("red", "blue"), lty = 1, lwd = 2)


tmap_mode("plot")

nsig_analysis <- tm_shape(franklin_results) +
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
    fill_alpha = 0.3,
    fill.legend = tm_legend(title = "MPSA Category",
                            position = c(0.85, 0.25),
                            bg.color = "white"
    )
  ) +
  tm_title("Local MPSA Results") +
  tm_shape(neighbors) +
  tm_polygons(fill = "#ffffbf",
              col = "black",
              lwd = 2) +
  tm_shape(target) +
  tm_polygons(fill = "#ffffbf",
              col = "black",
              lwd = 2) +
  tm_compass(type = "arrow",
             position = c("left", "bottom")) +
  tm_scalebar(position = c("left", "bottom"),
              bg.color = "white")

tmap_save(nsig_analysis, "output/nsig_analysis.png", dpi = 900)

qtm(franklin, fill = "unemployment_rate")
