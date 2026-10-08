
library(ggplot2)
library(reshape2)  # 用于数据转换

library(rrvgo)
library(RColorBrewer)
library(tibble) 

setwd("D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/")

vir  <- read.csv("miRNABP/virus_BP_representative_terms.csv")  %>%
  mutate(group = "Virus")
oth  <- read.csv("miRNABP/others_BP_representative_terms.csv") %>%
  mutate(group = "Others")
dat <- dplyr::bind_rows(vir, oth) %>%
  dplyr::select(parentTerm, score, group) %>%
  dplyr::group_by(parentTerm, group) %>%
  dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

View(dat)
mat <- dat %>%
  pivot_wider(names_from = parentTerm,
              values_from = score,
              values_fill = 0) %>%
  column_to_rownames("group") %>%
  as.matrix()

mat[mat == 0] <- NA


tag <- mat %>%               # 2×n 矩阵
  as.data.frame() %>%
  tibble::rownames_to_column("group") %>%
  tidyr::pivot_longer(-group, names_to = "term", values_to = "score") %>%
  group_by(term) %>%
  summarise(
    both   = all(!is.na(score)),        # 两行都有值
    onlyO  = !is.na(score[group=="Others"]) & is.na(score[group=="Virus"]),
    onlyV  = is.na(score[group=="Others"]) & !is.na(score[group=="Virus"]),
    .groups = "drop"
  ) %>%
  mutate(
    order = case_when(
      both  ~ 1,
      onlyO ~ 2,
      onlyV ~ 3,
      TRUE  ~ 4
    )
  ) %>%
  arrange(order) %>%
  pull(term)                 # 拿到排好序的通路名

mat_sorted <- mat[, tag]

# ===== 1. 转置矩阵，交换 x / y 轴 =====
mat_t <- t(mat_sorted)

# ===== 2. 重命名分组（转置后 group 变成列名）=====
# 保险写法：按原名字映射，避免列顺序问题
colnames(mat_t) <- ifelse(
  colnames(mat_t) == "Virus",    "virus-like hmiRNAs",
  ifelse(colnames(mat_t) == "Others", "other hmiRNAs",
         colnames(mat_t))
)

# 可选：指定列顺序（想先显示哪个就先放哪个）
mat_t <- mat_t[, c( "virus-like hmiRNAs","other hmiRNAs")]

# ===== 3. 0.15 mm 换算成 grid 的 lwd 单位 =====
# grid 里 1 lwd ≈ 1/96 英寸 ≈ 0.2646 mm
lw <- 0.15 * 96 / 25.4        # ≈ 0.567

blue_palette <- colorRampPalette(c("#A3C9D5", "#276C9E"))(30)

# ===== 4. 画图 =====
library(ComplexHeatmap)

Heatmap(
  mat_t,
  name = "average of -log10(qvalue)",
  na_col = "white",
  rect_gp = gpar(col = "white", lwd = lw),        # 单元格白线 0.15mm
  col = blue_palette,
  border = TRUE,
  border_gp = gpar(col = "black", lwd = lw),      # 外框 0.15mm
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "right",
  column_names_side = "bottom",
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8),           # 列名 8pt
  row_names_gp    = gpar(fontsize = 8),           # 行名 8pt
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 8),               # 图例标题 8pt
    labels_gp = gpar(fontsize = 8)                # 图例刻度 8pt
  )
)

mat_t %>%
  as.data.frame() %>%
  tibble::rownames_to_column("parentTerm") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
  readr::write_excel_csv("Figure 2E.csv")
