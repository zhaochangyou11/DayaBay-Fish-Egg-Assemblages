# ============================================
# Mantel检验与相关性热图分析（可调节距离版本）
# ============================================

# 清除环境
rm(list = ls())

# 加载必要的包
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(readxl)
library(reshape2)

# ============================================
# 1. 读取数据
# ============================================

# 设置文件路径（请根据实际情况修改）
file_path <- "C:/Users/Administrator/Desktop/大亚湾鱼卵及鱼类群落春秋两季动态变化研究/结果表示图/Mantel检验/"

# 读取物种数据（A组和B组）
A <- read_excel(paste0(file_path, "A组物种集.xlsx"), sheet = 1)
B <- read_excel(paste0(file_path, "B组物种集.xlsx"), sheet = 1)
env_data <- read_excel(paste0(file_path, "环境因子.xlsx"), sheet = 1)

# ============================================
# 2. 数据准备和样本匹配
# ============================================

# 处理A组数据
A <- as.data.frame(A)
rownames(A) <- A$Species
A <- A[, -1]

# 处理B组数据（删除空白行）
B <- as.data.frame(B)
B <- B[!is.na(B$Species) & B$Species != "", ]
rownames(B) <- B$Species
B <- B[, -1]

# 处理环境因子数据
env_data <- as.data.frame(env_data)
rownames(env_data) <- env_data$Stations
env_data <- env_data[, -1]

# 提取样本名称
A_samples <- colnames(A)
B_samples <- colnames(B)
env_samples <- rownames(env_data)

# 找出共同样本
common_samples <- intersect(intersect(A_samples, B_samples), env_samples)
cat("共同样本数:", length(common_samples), ":", paste(common_samples, collapse = ", "), "\n")

# ============================================
# 3. 准备物种数据
# ============================================

# 转置数据：行为样本，列为物种
A_t <- as.data.frame(t(A[, common_samples]))
B_t <- as.data.frame(t(B[, common_samples]))

# 给列名添加前缀
colnames(A_t) <- paste0("A_", colnames(A_t))
colnames(B_t) <- paste0("B_", colnames(B_t))

# 合并数据
species_combined <- cbind(A_t, B_t)

cat("A群落物种数:", ncol(A_t), "\n")
cat("B群落物种数:", ncol(B_t), "\n")

# ============================================
# 4. 准备环境因子数据（保持原始名称）
# ============================================

env_factors <- env_data[common_samples, ]

# 查看实际列名
cat("环境因子实际列名:\n")
print(colnames(env_factors))

# 根据实际列名选择（使用模糊匹配）
target_vars <- c("pH", "DO", "SST", "SSS", "Chl", "SRP", "SiO", "WD", "SDD", "SBT", "SBS")

selected_cols <- c()
for(var in target_vars) {
  matched <- grep(var, colnames(env_factors), value = TRUE, ignore.case = TRUE)
  if(length(matched) > 0) {
    selected_cols <- c(selected_cols, matched[1])
  }
}

selected_cols <- unique(selected_cols)
cat("选中的环境因子:", paste(selected_cols, collapse = ", "), "\n")

# 提取选定的环境因子
env_factors_selected <- env_factors[, selected_cols]

# 保持原始列名，不做任何修改
cat("最终列名（保持原始名称）:", paste(colnames(env_factors_selected), collapse = ", "), "\n")

# ============================================
# 5. 执行Mantel检验（逐个环境因子）- 修复版
# ============================================

# 检查缺失值
cat("A_t 中的缺失值数量:", sum(is.na(A_t)), "\n")
cat("B_t 中的缺失值数量:", sum(is.na(B_t)), "\n")

# 方法1：删除含有缺失值的样本（推荐，如果缺失不多）
# A_t <- A_t[complete.cases(A_t), ]
# B_t <- B_t[complete.cases(B_t), ]

# 方法2：将缺失值替换为0（假设缺失表示该物种未被检测到）
A_t[is.na(A_t)] <- 0
B_t[is.na(B_t)] <- 0

# 方法3：如果某物种全部为0或NA，删除该物种列
A_t <- A_t[, colSums(A_t, na.rm = TRUE) > 0]
B_t <- B_t[, colSums(B_t, na.rm = TRUE) > 0]

cat("处理后 A群落物种数:", ncol(A_t), "\n")
cat("处理后 B群落物种数:", ncol(B_t), "\n")

# 计算物种Bray-Curtis距离矩阵
A_dist <- vegdist(A_t, method = "bray")
B_dist <- vegdist(B_t, method = "bray")

# 对每个环境因子执行Mantel检验
env_names <- colnames(env_factors_selected)
n_env <- length(env_names)

mantel_results <- data.frame()

for(i in 1:n_env) {
  # 提取单个环境因子
  env_single <- env_factors_selected[, i]
  
  # 检查环境因子是否有缺失值
  if(any(is.na(env_single))) {
    cat("警告:", env_names[i], "存在缺失值，已删除对应样本\n")
    valid_idx <- !is.na(env_single)
    env_single <- env_single[valid_idx]
    A_dist_sub <- as.dist(as.matrix(A_dist)[valid_idx, valid_idx])
    B_dist_sub <- as.dist(as.matrix(B_dist)[valid_idx, valid_idx])
  } else {
    A_dist_sub <- A_dist
    B_dist_sub <- B_dist
  }
  
  env_dist <- dist(as.matrix(env_single), method = "euclidean")
  
  # A群落Mantel检验
  mantel_A <- mantel(A_dist_sub, env_dist, permutations = 999)
  # B群落Mantel检验
  mantel_B <- mantel(B_dist_sub, env_dist, permutations = 999)
  
  # 记录结果
  mantel_results <- rbind(mantel_results, data.frame(
    env_var = env_names[i],
    env_index = i,
    r_A = mantel_A$statistic,
    p_A = mantel_A$signif,
    r_B = mantel_B$statistic,
    p_B = mantel_B$signif
  ))
}

# 分类r和p值
mantel_results <- mantel_results %>%
  mutate(
    rd_A = cut(r_A, breaks = c(-Inf, 0.2, 0.4, Inf),
               labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    pd_A = cut(p_A, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
    rd_B = cut(r_B, breaks = c(-Inf, 0.2, 0.4, Inf),
               labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    pd_B = cut(p_B, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  )

print(mantel_results)

# ============================================
# 6. 计算环境因子相关性矩阵
# ============================================

env_cor <- cor(env_factors_selected, method = "pearson")

# 转换为长格式（只保留下三角）
env_cor_melt <- melt(env_cor)
env_cor_melt$Var1 <- as.character(env_cor_melt$Var1)
env_cor_melt$Var2 <- as.character(env_cor_melt$Var2)

# 只保留下三角（不包括对角线）
env_cor_lower <- env_cor_melt %>%
  filter(as.numeric(factor(Var1, levels = env_names)) > 
           as.numeric(factor(Var2, levels = env_names)))

# ============================================
# 7. 创建经典Mantel热图（可调节距离版本）
# ============================================

# 定义颜色
color_r <- c("< 0.2" = 0.5, "0.2 - 0.4" = 1.5, ">= 0.4" = 2.5)
color_p <- c("< 0.01" = "#F37600", "0.01 - 0.05" = "#3A9B51", ">= 0.05" = "#A9A9A9")

# ========== 可调节参数区域 ==========
# 这些参数控制节点间距，可根据需要调整
x_env_node  <- n_env + 0.2    # 环境因子节点x位置（左侧，固定）
x_community <- n_env + 1.5    # 群落节点x位置（右侧，调大=水平距离增大）
community_y_A <- n_env + 1.5  # A群落节点y位置（调大=向上移）
community_y_B <- -0.5         # B群落节点y位置（调小=向下移）
# ====================================

# 创建基础热图
p <- ggplot() +
  # 绘制下三角相关性方块
  geom_tile(data = env_cor_lower, 
            aes(x = Var2, y = Var1, fill = value),
            color = "white", size = 0.5) +
  
  # 设置颜色渐变（红蓝配色）
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "RdBu"),
    limits = c(-1, 1),
    name = "Pearson's r"
  ) +
  
  # 设置坐标轴
  scale_x_discrete(limits = env_names, position = "bottom",
                   expand = expansion(add = c(0.5, 2))) +
  scale_y_discrete(limits = rev(env_names),
                   expand = expansion(add = c(0.5, 0.5))) +
  
  # 主题设置
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                               size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(20, 150, 20, 20)  # 增大右边距防止截断
  )

# 添加Mantel检验连线（使用可调节参数）
y_pos <- seq(n_env, 1, by = -1)

for(i in 1:n_env) {
  y <- y_pos[i]
  
  # 获取当前环境因子的颜色和线宽
  line_col_A <- color_p[as.character(mantel_results$pd_A[i])]
  line_size_A <- color_r[as.character(mantel_results$rd_A[i])]
  line_col_B <- color_p[as.character(mantel_results$pd_B[i])]
  line_size_B <- color_r[as.character(mantel_results$rd_B[i])]
  
  # A群落连线（从环境因子到A群落节点）
  p <- p + annotate("segment", 
                    x = x_env_node, xend = x_community,
                    y = y, yend = community_y_A,
                    color = line_col_A, 
                    size = line_size_A,
                    alpha = 0.8)
  
  # B群落连线（从环境因子到B群落节点）
  p <- p + annotate("segment", 
                    x = x_env_node, xend = x_community,
                    y = y, yend = community_y_B,
                    color = line_col_B, 
                    size = line_size_B,
                    alpha = 0.8)
  
  # 环境因子节点（红点）
  p <- p + annotate("point", 
                    x = x_env_node, y = y,
                    color = "black", fill = "#DC143C", 
                    shape = 21, size = 3, stroke = 0.5)
}

# 添加A群落节点（顶部）- 放在循环外面，只添加一次
p <- p + annotate("point", 
                  x = x_community, y = community_y_A,
                  color = "black", fill = "#DC143C", 
                  shape = 21, size = 5, stroke = 0.8) +
  annotate("text", x = x_community + 0.5, y = community_y_A, 
           label = "A community", size = 4, hjust = 0, fontface = "bold")

# 添加B群落节点（底部）- 放在循环外面，只添加一次
p <- p + annotate("point", 
                  x = x_community, y = community_y_B,
                  color = "black", fill = "#DC143C", 
                  shape = 21, size = 5, stroke = 0.8) +
  annotate("text", x = x_community + 0.5, y = community_y_B, 
           label = "B community", size = 4, hjust = 0, fontface = "bold")

# 添加图例（位置根据新坐标自动调整）
legend_x <- x_community + 2.5  # 基于群落节点位置动态调整

# Mantel's r 图例
p <- p + 
  annotate("text", x = legend_x, y = n_env - 1, 
           label = "Mantel's r", fontface = "bold", size = 4, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 2, yend = n_env - 2, size = 0.5, color = "grey35") +
  annotate("text", x = legend_x + 1, y = n_env - 2, 
           label = "< 0.2", size = 3.5, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 3, yend = n_env - 3, size = 1.5, color = "grey35") +
  annotate("text", x = legend_x + 1, y = n_env - 3, 
           label = "0.2 - 0.4", size = 3.5, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 4, yend = n_env - 4, size = 2.5, color = "grey35") +
  annotate("text", x = legend_x + 1, y = n_env - 4, 
           label = ">= 0.4", size = 3.5, hjust = 0)

# Mantel's p 图例
p <- p + 
  annotate("text", x = legend_x, y = n_env - 6, 
           label = "Mantel's p", fontface = "bold", size = 4, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 7, yend = n_env - 7, size = 1.5, color = "#F37600") +
  annotate("text", x = legend_x + 1, y = n_env - 7, 
           label = "< 0.01", size = 3.5, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 8, yend = n_env - 8, size = 1.5, color = "#3A9B51") +
  annotate("text", x = legend_x + 1, y = n_env - 8, 
           label = "0.01 - 0.05", size = 3.5, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 9, yend = n_env - 9, size = 1.5, color = "#A9A9A9") +
  annotate("text", x = legend_x + 1, y = n_env - 9, 
           label = ">= 0.05", size = 3.5, hjust = 0)

# 动态调整绘图范围
xlim_max <- legend_x + 3       # 确保图例有足够空间
ylim_min <- community_y_B - 1  # 底部留白
ylim_max <- community_y_A + 1  # 顶部留白

p <- p + coord_cartesian(xlim = c(0.5, xlim_max), 
                         ylim = c(ylim_min, ylim_max),
                         clip = "off")

print(p)

# ============================================
# 8. 保存结果
# ============================================

# 保存为PDF（矢量格式）
ggsave(paste0(file_path, "Mantel_heatmap.pdf"), 
       plot = p, width = 14, height = 12, device = cairo_pdf)

# 保存为SVG（推荐用于Illustrator编辑）
ggsave(paste0(file_path, "Mantel_heatmap.svg"), 
       plot = p, width = 14, height = 12, device = "svg")

# 保存为EPS
ggsave(paste0(file_path, "Mantel_heatmap.eps"), 
       plot = p, width = 14, height = 12, device = cairo_ps)

# 保存为高分辨率TIFF
ggsave(paste0(file_path, "Mantel_heatmap.tiff"), 
       plot = p, width = 14, height = 12, dpi = 600, compression = "lzw")

# 保存Mantel结果
write.csv(mantel_results, paste0(file_path, "mantel_results.csv"), row.names = FALSE)

cat("\n分析完成！已生成文件：\n")
cat("1. Mantel_heatmap.pdf\n")
cat("2. Mantel_heatmap.svg（推荐用于编辑）\n")
cat("3. Mantel_heatmap.eps\n")
cat("4. Mantel_heatmap.tiff\n")
cat("5. mantel_results.csv\n")




# ============================================
# Mantel检验与相关性热图分析（可调节距离版本）
# ============================================

# 清除环境
rm(list = ls())

# 加载必要的包
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(readxl)
library(reshape2)

# ============================================
# 1. 读取数据
# ============================================

# 设置文件路径（请根据实际情况修改）
file_path <- "C:/Users/Administrator/Desktop/大亚湾鱼卵及鱼类群落春秋两季动态变化研究/结果表示图/Mantel检验/"

# 读取物种数据（A组和B组）
A <- read_excel(paste0(file_path, "A组物种集.xlsx"), sheet = 1)
B <- read_excel(paste0(file_path, "B组物种集.xlsx"), sheet = 1)
env_data <- read_excel(paste0(file_path, "环境因子.xlsx"), sheet = 1)

# ============================================
# 2. 数据准备和样本匹配
# ============================================

# 处理A组数据
A <- as.data.frame(A)
rownames(A) <- A$Species
A <- A[, -1]

# 处理B组数据（删除空白行）
B <- as.data.frame(B)
B <- B[!is.na(B$Species) & B$Species != "", ]
rownames(B) <- B$Species
B <- B[, -1]

# 处理环境因子数据
env_data <- as.data.frame(env_data)
rownames(env_data) <- env_data$Stations
env_data <- env_data[, -1]

# 提取样本名称
A_samples <- colnames(A)
B_samples <- colnames(B)
env_samples <- rownames(env_data)

# 找出共同样本
common_samples <- intersect(intersect(A_samples, B_samples), env_samples)
cat("共同样本数:", length(common_samples), ":", paste(common_samples, collapse = ", "), "\n")

# ============================================
# 3. 准备物种数据
# ============================================

# 转置数据：行为样本，列为物种
A_t <- as.data.frame(t(A[, common_samples]))
B_t <- as.data.frame(t(B[, common_samples]))

# 给列名添加前缀
colnames(A_t) <- paste0("A_", colnames(A_t))
colnames(B_t) <- paste0("B_", colnames(B_t))

# 合并数据
species_combined <- cbind(A_t, B_t)

cat("A群落物种数:", ncol(A_t), "\n")
cat("B群落物种数:", ncol(B_t), "\n")

# ============================================
# 4. 准备环境因子数据（保持原始名称）
# ============================================

env_factors <- env_data[common_samples, ]

# 查看实际列名
cat("环境因子实际列名:\n")
print(colnames(env_factors))

# 根据实际列名选择（使用模糊匹配）
target_vars <- c("pH", "DO", "SST", "SSS", "Chl", "SRP", "SiO", "WD", "SDD", "SBT", "SBS")

selected_cols <- c()
for(var in target_vars) {
  matched <- grep(var, colnames(env_factors), value = TRUE, ignore.case = TRUE)
  if(length(matched) > 0) {
    selected_cols <- c(selected_cols, matched[1])
  }
}

selected_cols <- unique(selected_cols)
cat("选中的环境因子:", paste(selected_cols, collapse = ", "), "\n")

# 提取选定的环境因子
env_factors_selected <- env_factors[, selected_cols]

# 保持原始列名，不做任何修改
cat("最终列名（保持原始名称）:", paste(colnames(env_factors_selected), collapse = ", "), "\n")

# ============================================
# 5. 执行Mantel检验（逐个环境因子）- 修复版
# ============================================

# 检查缺失值
cat("A_t 中的缺失值数量:", sum(is.na(A_t)), "\n")
cat("B_t 中的缺失值数量:", sum(is.na(B_t)), "\n")

# 方法1：删除含有缺失值的样本（推荐，如果缺失不多）
# A_t <- A_t[complete.cases(A_t), ]
# B_t <- B_t[complete.cases(B_t), ]

# 方法2：将缺失值替换为0（假设缺失表示该物种未被检测到）
A_t[is.na(A_t)] <- 0
B_t[is.na(B_t)] <- 0

# 方法3：如果某物种全部为0或NA，删除该物种列
A_t <- A_t[, colSums(A_t, na.rm = TRUE) > 0]
B_t <- B_t[, colSums(B_t, na.rm = TRUE) > 0]

cat("处理后 A群落物种数:", ncol(A_t), "\n")
cat("处理后 B群落物种数:", ncol(B_t), "\n")

# 计算物种Bray-Curtis距离矩阵
A_dist <- vegdist(A_t, method = "bray")
B_dist <- vegdist(B_t, method = "bray")

# 对每个环境因子执行Mantel检验
env_names <- colnames(env_factors_selected)
n_env <- length(env_names)

mantel_results <- data.frame()

for(i in 1:n_env) {
  # 提取单个环境因子
  env_single <- env_factors_selected[, i]
  
  # 检查环境因子是否有缺失值
  if(any(is.na(env_single))) {
    cat("警告:", env_names[i], "存在缺失值，已删除对应样本\n")
    valid_idx <- !is.na(env_single)
    env_single <- env_single[valid_idx]
    A_dist_sub <- as.dist(as.matrix(A_dist)[valid_idx, valid_idx])
    B_dist_sub <- as.dist(as.matrix(B_dist)[valid_idx, valid_idx])
  } else {
    A_dist_sub <- A_dist
    B_dist_sub <- B_dist
  }
  
  env_dist <- dist(as.matrix(env_single), method = "euclidean")
  
  # A群落Mantel检验
  mantel_A <- mantel(A_dist_sub, env_dist, permutations = 999)
  # B群落Mantel检验
  mantel_B <- mantel(B_dist_sub, env_dist, permutations = 999)
  
  # 记录结果
  mantel_results <- rbind(mantel_results, data.frame(
    env_var = env_names[i],
    env_index = i,
    r_A = mantel_A$statistic,
    p_A = mantel_A$signif,
    r_B = mantel_B$statistic,
    p_B = mantel_B$signif
  ))
}

# 分类r和p值
mantel_results <- mantel_results %>%
  mutate(
    rd_A = cut(r_A, breaks = c(-Inf, 0.2, 0.4, Inf),
               labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    pd_A = cut(p_A, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
    rd_B = cut(r_B, breaks = c(-Inf, 0.2, 0.4, Inf),
               labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")),
    pd_B = cut(p_B, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
  )

print(mantel_results)

# ============================================
# 6. 计算环境因子相关性矩阵及显著性
# ============================================

# ========== 可调节阈值区域 ==========
cor_method <- "pearson"     # 相关性方法：pearson / spearman
r_threshold <- 0.00         # 相关系数阈值；只有 |r| >= 该值时才显示星号
p_threshold <- 0.05         # 显著性阈值；只有 p < 该值时才显示星号
# 星号规则：
# *** p < 0.001
# **  p < 0.01
# *   p < 0.05
# ====================================

# 计算环境因子相关系数矩阵
env_cor <- cor(
  env_factors_selected,
  method = cor_method,
  use = "pairwise.complete.obs"
)

# 计算环境因子相关性显著性 p 值矩阵
env_p <- matrix(NA, nrow = n_env, ncol = n_env)
rownames(env_p) <- env_names
colnames(env_p) <- env_names

for(i in 1:n_env) {
  for(j in 1:n_env) {
    if(i == j) {
      env_p[i, j] <- NA
    } else {
      test_result <- cor.test(
        env_factors_selected[[i]],
        env_factors_selected[[j]],
        method = cor_method,
        use = "pairwise.complete.obs"
      )
      env_p[i, j] <- test_result$p.value
    }
  }
}

# 转换为长格式
env_cor_melt <- melt(
  env_cor,
  varnames = c("Var1", "Var2"),
  value.name = "value"
)

env_p_melt <- melt(
  env_p,
  varnames = c("Var1", "Var2"),
  value.name = "p_value"
)

env_cor_melt$Var1 <- as.character(env_cor_melt$Var1)
env_cor_melt$Var2 <- as.character(env_cor_melt$Var2)
env_p_melt$Var1 <- as.character(env_p_melt$Var1)
env_p_melt$Var2 <- as.character(env_p_melt$Var2)

# 合并相关系数和显著性p值，并生成星号
env_cor_melt <- env_cor_melt %>%
  left_join(env_p_melt, by = c("Var1", "Var2")) %>%
  mutate(
    sig_label = case_when(
      !is.na(p_value) & abs(value) >= r_threshold & p_value < 0.001 ~ "***",
      !is.na(p_value) & abs(value) >= r_threshold & p_value < 0.01  ~ "**",
      !is.na(p_value) & abs(value) >= r_threshold & p_value < p_threshold ~ "*",
      TRUE ~ ""
    )
  )

# 只保留下三角（不包括对角线）
env_cor_lower <- env_cor_melt %>%
  filter(as.numeric(factor(Var1, levels = env_names)) > 
           as.numeric(factor(Var2, levels = env_names)))

# 导出环境因子相关性及显著性结果
write.csv(
  env_cor_melt,
  paste0(file_path, "environment_correlation_significance.csv"),
  row.names = FALSE
)

# ============================================
# 7. 创建经典Mantel热图（可调节距离版本）
# ============================================

# 定义颜色
color_r <- c("< 0.2" = 0.5, "0.2 - 0.4" = 1.5, ">= 0.4" = 2.5)
color_p <- c("< 0.01" = "#F37600", "0.01 - 0.05" = "#3A9B51", ">= 0.05" = "#A9A9A9")

# ========== 可调节参数区域 ==========
# 这些参数控制节点间距，可根据需要调整
x_env_node  <- n_env + 0.2    # 环境因子节点x位置（左侧，固定）
x_community <- n_env + 1.5    # 群落节点x位置（右侧，调大=水平距离增大）
community_y_A <- n_env + 1.5  # A群落节点y位置（调大=向上移）
community_y_B <- -0.5         # B群落节点y位置（调小=向下移）
# ====================================

# 创建基础热图
p <- ggplot() +
  # 绘制下三角相关性方块
  geom_tile(data = env_cor_lower, 
            aes(x = Var2, y = Var1, fill = value),
            color = "white", size = 0.5) +
  
  # 添加环境因子相关性显著性星号
  geom_text(
    data = env_cor_lower,
    aes(x = Var2, y = Var1, label = sig_label),
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  
  # 设置颜色渐变（红蓝配色）
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "RdBu"),
    limits = c(-1, 1),
    name = "Pearson's r"
  ) +
  
  # 设置坐标轴
  scale_x_discrete(limits = env_names, position = "bottom",
                   expand = expansion(add = c(0.5, 2))) +
  scale_y_discrete(limits = rev(env_names),
                   expand = expansion(add = c(0.5, 0.5))) +
  
  # 主题设置
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, 
                               size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(20, 150, 20, 20)
  )

# 添加Mantel检验连线（使用可调节参数）
y_pos <- seq(n_env, 1, by = -1)

for(i in 1:n_env) {
  y <- y_pos[i]
  
  # 获取当前环境因子的颜色和线宽
  line_col_A <- color_p[as.character(mantel_results$pd_A[i])]
  line_size_A <- color_r[as.character(mantel_results$rd_A[i])]
  line_col_B <- color_p[as.character(mantel_results$pd_B[i])]
  line_size_B <- color_r[as.character(mantel_results$rd_B[i])]
  
  # A群落连线（从环境因子到A群落节点）
  p <- p + annotate("segment", 
                    x = x_env_node, xend = x_community,
                    y = y, yend = community_y_A,
                    color = line_col_A, 
                    size = line_size_A,
                    alpha = 0.8)
  
  # B群落连线（从环境因子到B群落节点）
  p <- p + annotate("segment", 
                    x = x_env_node, xend = x_community,
                    y = y, yend = community_y_B,
                    color = line_col_B, 
                    size = line_size_B,
                    alpha = 0.8)
  
  # 环境因子节点（红点）
  p <- p + annotate("point", 
                    x = x_env_node, y = y,
                    color = "black", fill = "#DC143C", 
                    shape = 21, size = 3, stroke = 0.5)
}

# 添加A群落节点（顶部）- 放在循环外面，只添加一次
p <- p + annotate("point", 
                  x = x_community, y = community_y_A,
                  color = "black", fill = "#DC143C", 
                  shape = 21, size = 5, stroke = 0.8) +
  annotate("text", x = x_community + 0.5, y = community_y_A, 
           label = "A community", size = 4, hjust = 0, fontface = "bold")

# 添加B群落节点（底部）- 放在循环外面，只添加一次
p <- p + annotate("point", 
                  x = x_community, y = community_y_B,
                  color = "black", fill = "#DC143C", 
                  shape = 21, size = 5, stroke = 0.8) +
  annotate("text", x = x_community + 0.5, y = community_y_B, 
           label = "B community", size = 4, hjust = 0, fontface = "bold")

# 添加图例（位置根据新坐标自动调整）
legend_x <- x_community + 2.5

# Mantel's r 图例
p <- p + 
  annotate("text", x = legend_x, y = n_env - 1, 
           label = "Mantel's r", fontface = "bold", size = 4, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 2, yend = n_env - 2, size = 0.5, color = "grey35") +
  annotate("text", x = legend_x + 1, y = n_env - 2, 
           label = "< 0.2", size = 3.5, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 3, yend = n_env - 3, size = 1.5, color = "grey35") +
  annotate("text", x = legend_x + 1, y = n_env - 3, 
           label = "0.2 - 0.4", size = 3.5, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 4, yend = n_env - 4, size = 2.5, color = "grey35") +
  annotate("text", x = legend_x + 1, y = n_env - 4, 
           label = ">= 0.4", size = 3.5, hjust = 0)

# Mantel's p 图例
p <- p + 
  annotate("text", x = legend_x, y = n_env - 6, 
           label = "Mantel's p", fontface = "bold", size = 4, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 7, yend = n_env - 7, size = 1.5, color = "#F37600") +
  annotate("text", x = legend_x + 1, y = n_env - 7, 
           label = "< 0.01", size = 3.5, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 8, yend = n_env - 8, size = 1.5, color = "#3A9B51") +
  annotate("text", x = legend_x + 1, y = n_env - 8, 
           label = "0.01 - 0.05", size = 3.5, hjust = 0) +
  annotate("segment", x = legend_x, xend = legend_x + 0.8, 
           y = n_env - 9, yend = n_env - 9, size = 1.5, color = "#A9A9A9") +
  annotate("text", x = legend_x + 1, y = n_env - 9, 
           label = ">= 0.05", size = 3.5, hjust = 0)

# 动态调整绘图范围
xlim_max <- legend_x + 3
ylim_min <- community_y_B - 1
ylim_max <- community_y_A + 1

p <- p + coord_cartesian(xlim = c(0.5, xlim_max), 
                         ylim = c(ylim_min, ylim_max),
                         clip = "off")

print(p)

# ============================================
# 8. 保存结果
# ============================================

# 保存为PDF（矢量格式）
ggsave(paste0(file_path, "Mantel_heatmap.pdf"), 
       plot = p, width = 14, height = 12, device = cairo_pdf)

# 保存为SVG（推荐用于Illustrator编辑）
ggsave(paste0(file_path, "Mantel_heatmap.svg"), 
       plot = p, width = 14, height = 12, device = "svg")

# 保存为EPS
ggsave(paste0(file_path, "Mantel_heatmap.eps"), 
       plot = p, width = 14, height = 12, device = cairo_ps)

# 保存为高分辨率TIFF
ggsave(paste0(file_path, "Mantel_heatmap.tiff"), 
       plot = p, width = 14, height = 12, dpi = 600, compression = "lzw")

# 保存Mantel结果
write.csv(mantel_results, paste0(file_path, "mantel_results.csv"), row.names = FALSE)

cat("\n分析完成！已生成文件：\n")
cat("1. Mantel_heatmap.pdf\n")
cat("2. Mantel_heatmap.svg（推荐用于编辑）\n")
cat("3. Mantel_heatmap.eps\n")
cat("4. Mantel_heatmap.tiff\n")
cat("5. mantel_results.csv\n")
cat("6. environment_correlation_significance.csv\n")