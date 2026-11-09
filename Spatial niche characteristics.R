
# 加载必要的库
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(reshape2)

# 设置桌面路径
desktop_path <- "C:/Users/Administrator/Desktop/生态位宽度及重叠度计算/"

# 读取 CSV 文件
file_path <- "C:/Users/Administrator/Desktop/生态位宽度及重叠度计算/生态位计算.csv"

# 尝试读取文件，指定编码
data <- read_csv(file_path, locale = locale(encoding = "GBK"))

# 检查数据是否正确读取
print("读取的数据：")
print(data)

# 检查列名并清理
colnames(data) <- make.names(colnames(data), unique = TRUE)
print("清理后的列名：")
print(colnames(data))

# 删除无效的列名（如 ...38、...39 等）
invalid_columns <- colnames(data)[grepl("\\.\\.", colnames(data))]
if(length(invalid_columns) > 0) {
  data <- data[, !colnames(data) %in% invalid_columns]
}

# 将数据转换为长格式，并处理重复值（取平均）
data_long <- data %>%
  pivot_longer(cols = -空间间隔, names_to = "Species", values_to = "Abundance") %>%
  group_by(空间间隔, Species) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

print("转换为长格式后的数据（前20行）：")
print(head(data_long, 20))

# 计算每个物种的总个体数
total_abundance <- data_long %>%
  group_by(Species) %>%
  summarise(Total = sum(Abundance, na.rm = TRUE)) %>%
  filter(Total > 0)

# 过滤掉总丰度为0的物种
data_long <- data_long %>%
  filter(Species %in% total_abundance$Species)

# 计算每个物种在每个空间间隔的相对丰度
data_long <- data_long %>%
  left_join(total_abundance, by = "Species") %>%
  mutate(Relative_Abundance = Abundance / Total) %>%
  select(-Total, -Abundance)

print("计算相对丰度后的数据（前20行）：")
print(head(data_long, 20))

# 计算 Levin 指数 (Bi)
s <- length(unique(data_long$空间间隔))

levin_indices <- data_long %>%
  group_by(Species) %>%
  summarise(
    Levin_Index_Bi = 1 / (s * sum(Relative_Abundance^2, na.rm = TRUE)),
    .groups = "drop"
  )

print("Levin 指数 (Bi)：")
print(levin_indices)

# ========== 计算 Pianka 指数 (O_ij) ==========
# 将数据转换为宽格式（物种 x 空间间隔矩阵）
data_wide <- data_long %>%
  pivot_wider(
    names_from = 空间间隔,
    values_from = Relative_Abundance,
    values_fill = 0
  )

# 转换为矩阵（去掉 Species 列）
species_names <- data_wide$Species
data_matrix <- as.matrix(data_wide[, -1])
rownames(data_matrix) <- species_names

print("宽格式数据矩阵维度：")
print(dim(data_matrix))

# 计算 Pianka 指数的函数
calculate_pianka <- function(matrix_data) {
  n_species <- nrow(matrix_data)
  pianka_matrix <- matrix(0, nrow = n_species, ncol = n_species)
  rownames(pianka_matrix) <- rownames(matrix_data)
  colnames(pianka_matrix) <- rownames(matrix_data)
  
  for(i in 1:n_species) {
    for(j in 1:n_species) {
      if(i == j) {
        pianka_matrix[i, j] <- 1
      } else {
        numerator <- sum(matrix_data[i, ] * matrix_data[j, ], na.rm = TRUE)
        denominator <- sqrt(
          sum(matrix_data[i, ]^2, na.rm = TRUE) *
            sum(matrix_data[j, ]^2, na.rm = TRUE)
        )
        
        if(denominator > 0) {
          pianka_matrix[i, j] <- numerator / denominator
        } else {
          pianka_matrix[i, j] <- 0
        }
      }
    }
  }
  return(pianka_matrix)
}

# 计算 Pianka 指数矩阵
pianka_matrix <- calculate_pianka(data_matrix)

print("Pianka指数矩阵（前5x5）：")
print(pianka_matrix[1:min(5, nrow(pianka_matrix)), 1:min(5, ncol(pianka_matrix))])

# ========== 物种排序 ==========
# 按 Pianka 相似性做层次聚类排序，更容易形成块状结构
dist_mat <- as.dist(1 - pianka_matrix)
hc <- hclust(dist_mat, method = "average")
species_order <- rownames(pianka_matrix)[hc$order]

# 按新顺序重排矩阵
pianka_matrix <- pianka_matrix[species_order, species_order]

# ========== 导出完整 Pianka 指数矩阵 ==========
pianka_output <- file.path(desktop_path, "Pianka指数矩阵.csv")
write.csv(pianka_matrix, pianka_output, fileEncoding = "GBK")
cat("Pianka指数矩阵已保存至：", pianka_output, "\n")

# ========== 转长格式 ==========
pianka_long <- melt(
  pianka_matrix,
  varnames = c("Species1", "Species2"),
  value.name = "Pianka_Index"
)

# 添加行列编号
pianka_long <- pianka_long %>%
  mutate(
    Species1 = factor(Species1, levels = species_order),
    Species2 = factor(Species2, levels = species_order),
    row_id = as.numeric(Species2),   # 行
    col_id = as.numeric(Species1)    # 列
  )

# ========== 导出各物种间 Pianka 指数表 ==========
# 只保留一次配对，不包含对角线
pianka_pairwise_table <- pianka_long %>%
  filter(row_id > col_id) %>%
  select(Species1, Species2, Pianka_Index) %>%
  arrange(desc(Pianka_Index))

pairwise_output <- file.path(desktop_path, "各物种间Pianka指数表.csv")
write.csv(
  pianka_pairwise_table,
  pairwise_output,
  row.names = FALSE,
  fileEncoding = "GBK"
)
cat("各物种间Pianka指数表已保存至：", pairwise_output, "\n")

# ========== 只保留三角区用于作图 ==========
# 保留含对角线的下三角
pianka_triangle_plot <- pianka_long %>%
  filter(row_id >= col_id)

# ========== 创建生态位重叠等级分类变量 ==========
pianka_triangle_plot <- pianka_triangle_plot %>%
  mutate(
    Overlap_Level = case_when(
      Pianka_Index >= 0.60 ~ "High (≥0.60)",
      Pianka_Index >= 0.40 ~ "Moderate (0.40-0.60)",
      TRUE ~ "Low (<0.40)"
    )
  ) %>%
  mutate(
    Overlap_Level = factor(
      Overlap_Level,
      levels = c("Low (<0.40)", "Moderate (0.40-0.60)", "High (≥0.60)")
    )
  )

print("生态位重叠等级分布：")
print(table(pianka_triangle_plot$Overlap_Level))

# ========== 动态设置尺寸 ==========
n_species <- nrow(pianka_matrix)

if(n_species <= 10) {
  tile_text_size <- 4
  fig_width <- 8
  fig_height <- 8
} else if(n_species <= 20) {
  tile_text_size <- 4
  fig_width <- 10
  fig_height <- 10
} else {
  tile_text_size <- 3.5
  fig_width <- 12
  fig_height <- 12
}

# ========== 构建左侧标签和底部标签 ==========
left_label_df <- data.frame(
  x = 0.35,
  y = 1:n_species,
  label = species_order
)

bottom_label_df <- data.frame(
  x = 1:n_species,
  y = n_species + 0.6,
  label = species_order
)

# ========== 绘图 ==========
p <- ggplot() +
  # 三角热图主体
  geom_tile(
    data = pianka_triangle_plot,
    aes(x = col_id, y = row_id, fill = Overlap_Level),
    color = "white",
    linewidth = 0.4
  ) +
  
  # 左侧物种名称
  geom_text(
    data = left_label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    size = tile_text_size,
    colour = "blue",
    fontface = "plain"
  ) +
  
  # 底部物种名称
  geom_text(
    data = bottom_label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    angle = 90,
    hjust = 1,
    vjust = 0.5,
    size = tile_text_size,
    colour = "blue",
    fontface = "plain"
  ) +
  
  scale_fill_manual(
    values = c(
      "Low (<0.40)" = "#3b5bdb",
      "Moderate (0.40-0.60)" = "#f1c40f",
      "High (≥0.60)" = "#ff2d2d"
    ),
    name = "Pianka index (Oij)"
  ) +
  
  scale_x_continuous(
    limits = c(0.2, n_species + 0.9),
    expand = c(0, 0)
  ) +
  scale_y_reverse(
    limits = c(n_species + 1.2, 0.5),
    expand = c(0, 0)
  ) +
  
  coord_fixed(clip = "off") +
  
  labs(
    title = "Spring fish community",
    x = NULL,
    y = NULL
  ) +
  
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.margin = margin(20, 30, 60, 60)
  )

# ========== 保存图形 ==========
output_png <- file.path(desktop_path, "生态位重叠三角热图_白色背景.png")
png(output_png, width = fig_width, height = fig_height, units = "in", res = 300)
print(p)
dev.off()
cat("PNG热图已保存至：", output_png, "\n")

output_pdf <- file.path(desktop_path, "生态位重叠三角热图_白色背景_可编辑.pdf")
ggsave(
  filename = output_pdf,
  plot = p,
  device = cairo_pdf,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 300
)
cat("PDF热图（文字可编辑）已保存至：", output_pdf, "\n")

# 保存 Levin 指数结果
levin_output <- file.path(desktop_path, "Levin指数.csv")
write.csv(levin_indices, levin_output, row.names = FALSE, fileEncoding = "GBK")
cat("Levin指数已保存至：", levin_output, "\n")