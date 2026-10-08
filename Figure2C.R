##### 先取 miRNA family 的 average

# 设置工作目录
setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\1seedregionInmiRNA\\miRNAFamily")

# 加载必要的包
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(tibble)     # 需要用到 recode

# ===== 1. 读取数据 =====
file_path <- "miRNA_seed_similarity_analysis.csv"

# 关键改动 ①：把空字符串、空格也当成 NA 读进来
data <- read.csv(file_path, na.strings = c("NA", "", " ", "  "))

# 保险起见，把所有字符列的空白字符串统一转成 NA
data <- data %>%
  mutate(across(where(is.character), ~ na_if(trimws(.x), "")))

# ===== 2. 数据清洗 + 按 family 聚合 =====
# 关键改动 ②：filter 里加上 family 和 category 的非 NA 判断
data_family_avg <- data %>%
  filter(
    !is.na(mam_avg_consistant),
    !is.na(non_mam_avg_consistant),
    !is.na(family),
    !is.na(category)
  ) %>%
  group_by(family, category) %>%
  summarise(
    mam_avg_sim     = mean(mam_avg_consistant),
    non_mam_avg_sim = mean(non_mam_avg_consistant),
    .groups = "drop"
  )

# ===== 3. 宽格式 → 长格式 =====
data_long <- data_family_avg %>%
  pivot_longer(
    cols = c(mam_avg_sim, non_mam_avg_sim),
    names_to = "Group_Type",
    values_to = "Similarity"
  )

# ===== 4. 优化因子标签和顺序 =====

# 4.1 Group_Type（分面标签）
data_long$Group_Type <- factor(
  data_long$Group_Type,
  levels = c("mam_avg_sim", "non_mam_avg_sim"),
  labels = c("Mammals", "Non-Mammals")
)

# 4.2 category：先重命名值，再设 levels（顺序不能反！）
data_long$category <- dplyr::recode(
  data_long$category,
  "Viruslike" = "virus-like hmiRNAs",
  "Others"    = "other hmiRNAs"
)

data_long$category <- factor(
  data_long$category,
  levels = c("virus-like hmiRNAs", "other hmiRNAs")
)

# ===== 5. 画图 =====
p <- ggplot(data_long, aes(x = category, y = Similarity, fill = category)) +
  geom_violin(
    trim      = FALSE,
    alpha     = 0.7,
    width     = 0.8,
    color     = "black",
    linewidth = 0.15
  ) +
  
  # 显示每个数据点
  geom_jitter(
    aes(color = category),
    width = 0.06,
    size  = 1.3,
    alpha = 0.65,
    shape = 16
  ) +
  
  facet_wrap(~Group_Type) +
  
  stat_compare_means(
    method   = "wilcox.test",
    label    = "p.signif",
    label.x  = 1.5,
    size     = 5,
    hide.ns  = TRUE,
    symnum.args = list(
      cutpoints = c(0, 0.05, 1),
      symbols   = c("*", "ns")
    )
  ) +
  
  # 小提琴颜色
  scale_fill_manual(
    values = c(
      "virus-like hmiRNAs" = "#4E79A7",
      "other hmiRNAs"      = "gray"
    )
  ) +
  
  # 散点颜色
  scale_color_manual(
    values = c(
      "virus-like hmiRNAs" = "black",
      "other hmiRNAs"      = "black"
    )
  ) +
  
  labs(
    title = "",
    x     = "",
    y     = "Average seed consistency\nin each miRNA family",
    fill  = "Category"
  ) +
  
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, 0.2)
  ) +
  
  theme_classic() +
  
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "none",
    axis.text       = element_text(color = "black", size = 10),
    
    axis.line = element_line(
      color = "black",
      linewidth = 0.15
    ),
    
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.15
    ),
    
    strip.background = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.15
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    )
  )
print(p)

data_export <- data_long %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%   # 数值列保留 3 位
  arrange(Group_Type, category)

readr::write_excel_csv(data_export, "data_long.csv")