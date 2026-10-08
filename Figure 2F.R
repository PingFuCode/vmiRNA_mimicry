# 加载必要的包
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\6chr")
##################################################重新画倍富集图############################
library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(circlize)

# 读取病毒模拟hmiRNA数据
file_path <- "D:/0work/0wholetransctiptome/2sRNAminic/6chr/0006mer_target_with_hsaID_chr.csv"
data <- read_csv(file_path)

# 处理染色体列：拆分含分号的行，并过滤掉以H开头的染色体
expanded_data <- data %>%
  separate_rows(Chromosome, sep = ";") %>%
  filter(Chromosome != "", !str_starts(Chromosome, "H"))

# 统计每个病毒在每个染色体上的出现次数 (1)
count_data <- expanded_data %>%
  count(abbreviation, Chromosome, name = "Count")

# 转换为宽格式（病毒为行，染色体为列）
count_matrix <- count_data %>%
  pivot_wider(
    id_cols = abbreviation,
    names_from = Chromosome,
    values_from = Count,
    values_fill = 0
  ) %>%
  column_to_rownames("abbreviation") %>%
  as.matrix()

# 读取hmiRNA分布数据 (2)
hmiRNA_dist <- read_delim("hmiRNADistribution.txt", delim = "\t") %>%
  select(Chr, miRNA) %>%
  distinct() %>%
  filter(Chr != "", !str_starts(Chr, "H"))

# 计算每条染色体的hmiRNA总数
chr_total_miRNAs <- hmiRNA_dist %>%
  count(Chr, name = "total_miRNAs")

# 计算总hmiRNA数用于比例计算
total_all_miRNAs <- sum(chr_total_miRNAs$total_miRNAs)

# 计算每条染色体的hmiRNA比例 (2)
chr_proportions <- chr_total_miRNAs %>%
  mutate(proportion = total_miRNAs / total_all_miRNAs) %>%
  select(Chr, proportion)

# 计算Fold enrichment
# 首先计算每个病毒在每个染色体的比例 (1)
virus_chr_proportions <- count_data %>%
  group_by(abbreviation) %>%
  mutate(total_virus_miRNAs = sum(Count)) %>%
  ungroup() %>%
  mutate(proportion_virus = Count / total_virus_miRNAs) %>%
  select(abbreviation, Chromosome, proportion_virus)

# 转换为宽格式并计算Fold enrichment
fold_enrichment_matrix <- virus_chr_proportions %>%
  left_join(chr_proportions, by = c("Chromosome" = "Chr")) %>%
  mutate(fold_enrichment = proportion_virus / proportion) %>%
  select(abbreviation, Chromosome, fold_enrichment) %>%
  pivot_wider(
    id_cols = abbreviation,
    names_from = Chromosome,
    values_from = fold_enrichment,
    values_fill = 0
  ) %>%
  column_to_rownames("abbreviation") %>%
  as.matrix() %>%
  round(4)

write.csv(fold_enrichment_matrix, "Fold.csv")


chr.order <- c(paste0( 1:22), "X")
chr.order <- intersect(chr.order, colnames(fold_enrichment_matrix))

fold_enrichment_matrix <- fold_enrichment_matrix[, chr.order]

# 绘制热图
col_fun <- colorRamp2(c(0, max(fold_enrichment_matrix, na.rm = TRUE)), 
                      c("white", "#276C9E"))

ht <- Heatmap(
  fold_enrichment_matrix,
  name = "Fold enrichment",
  col = col_fun,
  na_col = "white",
  rect_gp = gpar(col = NA, lwd = 0),
  border = TRUE,
  border_gp = gpar(col = "black", lwd = 1),
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_side = "left",
  column_names_rot = 0,
  row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 11))
 # cell_fun = function(j, i, x, y, width, height, fill) {
  #  if(fold_enrichment_matrix[i, j] > 0.0) {
  ##    grid.text(sprintf("%.2f", fold_enrichment_matrix[i, j]),
    #            x, y, gp = gpar(fontsize = 9))
   # }
 # }

ht
# 保存热图
pdf("Fold_enrichment_heatmap.pdf", width = 7, height = 7)
draw(ht)
dev.off()


mat_export <- fold_enrichment_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("rowname") %>%          # 行名变成第一列
  mutate(across(where(is.numeric), ~ round(.x, 3)))  # 数值保留3位

readr::write_excel_csv(mat_export, "Figure 2F.csv")





# Fisher精确检验
# 准备数据
all_chromosomes <- unique(c(colnames(fold_enrichment_matrix), chr_proportions$Chr))

# 创建包含所有miRNA信息的数据框
all_miRNAs_info <- expanded_data %>%
  select(abbreviation, Chromosome) %>%
  mutate(is_mimic = "Mimic") %>%
  bind_rows(
    hmiRNA_dist %>%
      rename(abbreviation = miRNA, Chromosome = Chr) %>%
      mutate(is_mimic = "Non-mimic")
  )

# 初始化结果数据框
results <- data.frame(
  Virus = character(),
  Chr = character(),
  Mimic_On = integer(),
  Mimic_Not = integer(),
  NonMimic_On = integer(),
  NonMimic_Not = integer(),
  p.value = numeric(),
  stringsAsFactors = FALSE
)

# 对每个病毒和每个染色体进行Fisher检验
for (virus in rownames(fold_enrichment_matrix)) {
  for (chr in all_chromosomes) {
    A <- sum(all_miRNAs_info$is_mimic == "Mimic" & 
               all_miRNAs_info$abbreviation == virus & 
               all_miRNAs_info$Chromosome == chr, na.rm = TRUE)
    
    B <- sum(all_miRNAs_info$is_mimic == "Mimic" & 
               all_miRNAs_info$abbreviation == virus & 
               (all_miRNAs_info$Chromosome != chr | is.na(all_miRNAs_info$Chromosome)), na.rm = TRUE)
    
    C <- sum(all_miRNAs_info$is_mimic == "Non-mimic" & 
               all_miRNAs_info$Chromosome == chr, na.rm = TRUE)
    
    D <- sum(all_miRNAs_info$is_mimic == "Non-mimic" & 
               (all_miRNAs_info$Chromosome != chr | is.na(all_miRNAs_info$Chromosome)), na.rm = TRUE)
    
    contingency_table <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE,
                                dimnames = list(c("On chromosome", "Not on chromosome"),
                                                c("Mimic", "Non-mimic")))
    
    fisher_test <- fisher.test(contingency_table)
    
    results <- rbind(results, data.frame(
      Virus = virus,
      Chr = chr,
      Mimic_On = A,
      Mimic_Not = B,
      NonMimic_On = C,
      NonMimic_Not = D,
      p.value = fisher_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

# ===== 新增：BH 校正 + 显著性标记（基于校正后的 p 值）=====
results <- results %>%
  mutate(
    p.adj = p.adjust(p.value, method = "BH"),   # BH / FDR 校正
    p.value = round(p.value, 3),
    p.adj   = round(p.adj, 3),
    sig = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

# 写出带标记的文件
write.csv(results, "Fold_pvalue.csv", row.names = FALSE)











x <- c(51,14,9,2,29,336,90,9,72,100,6,3,3,4,38,21,4,6,40,8,3,16)
y <- c(22,9,6,2,16,23,20,8,22,19,5,3,3,3,18,14,4,5,19,4,3,11)

cor.test(x, y, method = "pearson")


















# 读取数据
file_path <- "D:/0work/0wholetransctiptome/2sRNAminic/6chr/0006mer_target_with_hsaID_chr.csv"
data <- read_csv(file_path)

# 处理染色体列：拆分含分号的行，并过滤掉以H开头的染色体
expanded_data <- data %>%
  separate_rows(Chromosome, sep = ";") %>%
  filter(Chromosome != "", 
         !str_starts(Chromosome, "H"))  # 过滤掉以H开头的染色体

# 统计每个病毒在每个染色体上的出现次数
count_data <- expanded_data %>%
  count(abbreviation, Chromosome, name = "Count")

# 转换为宽格式（病毒为行，染色体为列）
count_matrix <- count_data %>%
  pivot_wider(
    id_cols = abbreviation,
    names_from = Chromosome,
    values_from = Count,
    values_fill = 0
  ) %>%
  column_to_rownames("abbreviation") %>%
  as.matrix()
View(count_matrix)
# 归一化处理：按行（病毒）计算比例
normalized_matrix <- apply(count_matrix, 1, function(x) {
  if(sum(x) > 0) round(x / sum(x), 3) else x
})

# 转置矩阵，使行仍然代表病毒，列代表染色体
normalized_matrix <- t(normalized_matrix)

# 定义染色体排序顺序
chr_order <- c(as.character(1:22), "X", "Y")

# 提取当前存在的染色体
existing_chr <- colnames(normalized_matrix)

# 按照1-22,X,Y的顺序排序染色体
sorted_chr <- chr_order[chr_order %in% existing_chr]

# 对矩阵列进行重新排序
normalized_matrix_sorted <- normalized_matrix[, sorted_chr, drop = FALSE]

# 导出归一化矩阵
output_file <- "D:/0work/0wholetransctiptome/2sRNAminic/6chr/000virus_chr.csv"
write.csv(normalized_matrix_sorted, output_file)








#####画图
head(normalized_matrix_sorted)
#rowSums(normalized_matrix_sorted)
row_sum_normalize <- function(x) {
  x / sum(x)
}

# 对每一行应用归一化
normalized_matrix_final <- t(apply(normalized_matrix_sorted, 1, row_sum_normalize))
head(rowSums(normalized_matrix_final))  # 应该都接近1

mat <- as.matrix(normalized_matrix_sorted)   # 保证是数值矩阵
rownames(normalized_matrix_final) <- rownames(mat)
#normalized_matrix_final
# 设置热图颜色
col_fun <- colorRamp2(c(0, max(normalized_matrix_final)), c("white", "#276C9E"))


ht <- Heatmap(
  normalized_matrix_final,
      name = "Row-sum normalization\nProportion",
      col = col_fun,
      na_col = "white",
      rect_gp = gpar(col = NA, lwd = 0),  # 👈 去掉网格线
      #rect_gp = gpar(col = "grey90", lwd = 0.5),# 👈 去掉网格线
      border = TRUE,
      border_gp = gpar(col = "black", lwd = 1),
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "left",
      column_names_rot = 0,
      row_names_gp = gpar(fontsize = 11),
      column_names_gp = gpar(fontsize = 11),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(normalized_matrix_final[i, j] > 0.0) {
    grid.text(sprintf("%.2f", normalized_matrix_final[i, j]),
              x, y, gp = gpar(fontsize = 9))}
      })
      
      # 绘制热图
#pdf("D:/0work/0wholetransctiptome/2sRNAminic/6chr/virus_chr_heatmap.pdf", 
         # width = 10, height = 7)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "bottom")
#dev.off()