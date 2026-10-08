######################## GSE184320 HIV-1 皮肤和PBMC样本分析 ########################

# 设置工作目录
setwd("D:/0work/0wholetransctiptome/2sRNAminic/17HIV1/GEO_HIV_bulk")

# 加载必要的包
library(limma)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)

# 1. 读取原始counts数据 ------------------------------------------------------------
counts_data <- read.delim(
  "GSE184320_Gene_Raw_Count_Bulk_RNA-seq_Smart-seq2.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 查看数据结构
cat("原始数据维度：", dim(counts_data), "\n")
cat("列名：", colnames(counts_data), "\n")

# 2. 数据准备和RPM计算 ------------------------------------------------------------
# 提取基因信息
gene_info <- counts_data[, c("HGNC", "Ensembl_ID")]
counts_matrix <- counts_data[, -c(1, 2)]  # 移除HGNC和Ensembl_ID列

# 将counts矩阵转换为数值型
counts_matrix <- as.matrix(counts_matrix)
mode(counts_matrix) <- "numeric"

# RPM计算函数
calculate_rpm <- function(counts_df) {
  # 对每个样本，计算每百万reads的计数
  rpm_df <- sweep(counts_df, 2, colSums(counts_df) / 1e6, FUN = "/")
  return(rpm_df)
}

# 计算RPM
rpm_matrix <- calculate_rpm(counts_matrix)

# 设置行名为基因名
rownames(rpm_matrix) <- gene_info$HGNC

# 添加一个小常数避免log2(0)
rpm_matrix_log <- log2(rpm_matrix + 1)

# 查看RPM统计
cat("\nRPM矩阵维度：", dim(rpm_matrix), "\n")
cat("RPM矩阵范围：", range(rpm_matrix), "\n")

# 3. 构建样本分组信息 ------------------------------------------------------------
sample_names <- colnames(rpm_matrix)
sample_group <- data.frame(
  Sample = sample_names,
  Group = ifelse(grepl("^HC_", sample_names), "Healthy", "HIV"),
  Type = ifelse(grepl("SKIN$", sample_names), "Skin", "PBMC"),
  stringsAsFactors = FALSE
)

# 解析样本信息
sample_group$ID <- gsub("^.*_(SKIN|PBMC)$", "", sample_group$Sample)
sample_group$ID <- gsub("_", "", sample_group$ID)

cat("\n样本分组信息：\n")
print(sample_group)

# 4. 过滤低表达基因 ------------------------------------------------------------
# 至少在2个样本中RPM > 1（对应原始counts）
keep <- rowSums(rpm_matrix > 1, na.rm = TRUE) >= 2
rpm_matrix_filtered <- rpm_matrix[keep, ]
rpm_matrix_log_filtered <- rpm_matrix_log[keep, ]

cat("\n过滤前基因数：", nrow(rpm_matrix), "\n")
cat("过滤后基因数：", nrow(rpm_matrix_filtered), "\n")

# 5. 读取HIV-1共有靶基因信息 ----------------------------------------------------
target_info <- read.csv(
  "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 提取HIV-1共有靶基因
hiv1_target_genes <- target_info %>%
  dplyr::filter(abbreviation == "HIV-1") %>%
  dplyr::pull(Intersection_Genes) %>%
  paste(collapse = ";") %>%
  stringr::str_split(";") %>%
  unlist() %>%
  trimws() %>%
  .[. != "" & !is.na(.)] %>%
  unique()

hiv1_target_key <- toupper(hiv1_target_genes)

cat("\nHIV-1共有靶基因总数：", length(hiv1_target_genes), "\n")

# 6. 分别对皮肤和PBMC样本进行差异表达分析 -------------------------------------------
# 6.1 皮肤样本差异分析
skin_samples <- sample_group %>% dplyr::filter(Type == "Skin")
skin_rpm <- rpm_matrix_filtered[, skin_samples$Sample]

skin_group <- factor(skin_samples$Group, levels = c("Healthy", "HIV"))

# 设计矩阵
skin_design <- model.matrix(~ 0 + skin_group)
colnames(skin_design) <- levels(skin_group)

# 线性模型拟合
skin_fit <- lmFit(skin_rpm, skin_design)
skin_contrast <- makeContrasts(HIV_vs_Healthy = HIV - Healthy, levels = skin_design)
skin_fit2 <- contrasts.fit(skin_fit, skin_contrast)
skin_fit2 <- eBayes(skin_fit2)

# 提取结果
skin_res <- topTable(skin_fit2, number = Inf, adjust.method = "BH")

skin_res$Gene <- rownames(skin_res)
skin_res <- skin_res %>%
  dplyr::mutate(
    GeneKey = toupper(Gene),
    is_HIV1_target = GeneKey %in% hiv1_target_key,
    DEG_status = dplyr::case_when(
      adj.P.Val < 0.05 & logFC > 1 ~ "Up in HIV",
      adj.P.Val < 0.05 & logFC < -1 ~ "Down in HIV",
      TRUE ~ "Not significant"
    )
  )


# 7. 导出完整分析结果 ------------------------------------------------------------
# 皮肤样本结果
write.csv(skin_res, "GSE184320_Skin_HIV_vs_Healthy_limma_results.csv", row.names = FALSE)


# 8. 提取HIV-1靶基因的差异表达结果 ------------------------------------------------
# 皮肤样本中的HIV-1靶基因
skin_hiv1_target_res <- skin_res %>% dplyr::filter(is_HIV1_target)
skin_hiv1_target_deg <- skin_hiv1_target_res %>%
  dplyr::filter(adj.P.Val < 0.05, abs(logFC) > 1)


# 导出
write.csv(skin_hiv1_target_res, "GSE184320_Skin_HIV1_target_genes_all.csv", row.names = FALSE)
write.csv(skin_hiv1_target_deg, "GSE184320_Skin_HIV1_target_genes_DEG.csv", row.names = FALSE)
write.csv(pbmc_hiv1_target_res, "GSE184320_PBMC_HIV1_target_genes_all.csv", row.names = FALSE)
write.csv(pbmc_hiv1_target_deg, "GSE184320_PBMC_HIV1_target_genes_DEG.csv", row.names = FALSE)

# 9. 统计信息 ------------------------------------------------------------
cat("\n========== 皮肤样本分析统计 ==========\n")
cat("总基因数：", nrow(skin_res), "\n")
cat("差异表达基因数(FDR<0.05, |logFC|>1)：", 
    sum(skin_res$adj.P.Val < 0.05 & abs(skin_res$logFC) > 1, na.rm = TRUE), "\n")
cat("上调基因数：", 
    sum(skin_res$adj.P.Val < 0.05 & skin_res$logFC > 1, na.rm = TRUE), "\n")
cat("下调基因数：", 
    sum(skin_res$adj.P.Val < 0.05 & skin_res$logFC < -1, na.rm = TRUE), "\n")
cat("表达的HIV-1靶基因数：", nrow(skin_hiv1_target_res), "\n")
cat("差异表达的HIV-1靶基因数：", nrow(skin_hiv1_target_deg), "\n")

# 10. 火山图 - 皮肤样本 ---------------------------------------------------------
min_nonzero_fdr_skin <- min(skin_res$adj.P.Val[skin_res$adj.P.Val > 0], na.rm = TRUE)

skin_plot <- skin_res %>%
  dplyr::mutate(
    adj.P.Val.plot = ifelse(adj.P.Val == 0, min_nonzero_fdr_skin / 10, adj.P.Val),
    neg_log10_FDR = -log10(adj.P.Val.plot),
    plot_group = dplyr::case_when(
      is_HIV1_target & adj.P.Val < 0.05 & logFC > 1 ~ "Up-regulated target genes",
      is_HIV1_target & adj.P.Val < 0.05 & logFC < -1 ~ "Down-regulated target genes",
      !is_HIV1_target & adj.P.Val < 0.05 & logFC > 1 ~ "Up-regulated other genes",
      !is_HIV1_target & adj.P.Val < 0.05 & logFC < -1 ~ "Down-regulated other genes",
      TRUE ~ "Not significant"
    )
  )

# 标注显著差异的HIV-1靶基因
label_genes_skin <- skin_plot %>%
  dplyr::filter(is_HIV1_target, adj.P.Val < 0.05, abs(logFC) > 1) %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = 20)

p_skin <- ggplot(skin_plot, aes(x = logFC, y = neg_log10_FDR)) +
  geom_point(
    data = skin_plot %>% dplyr::filter(plot_group == "Not significant"),
    color = "grey80",
    size = 0.8,
    alpha = 0.5
  ) +
  geom_point(
    data = skin_plot %>% dplyr::filter(plot_group %in% c("Up-regulated other genes", "Down-regulated other genes")),
    aes(color = plot_group),
    size = 1.5,
    alpha = 0.6
  ) +
  geom_point(
    data = skin_plot %>% dplyr::filter(plot_group %in% c("Up-regulated target genes", "Down-regulated target genes")),
    aes(color = plot_group),
    size = 2,
    alpha = 0.95
  ) +
  scale_color_manual(
    values = c(
      "Up-regulated other genes" = "#F4A6A6",
      "Down-regulated other genes" = "#A6C8E8",
      "Up-regulated target genes" = "#E21C21",
      "Down-regulated target genes" = "#3A7CB5"
    )
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4, color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "black") +
  geom_text_repel(
    data = label_genes_skin,
    aes(label = Gene),
    size = 3,
    max.overlaps = 100,
    box.padding = 0.4,
    point.padding = 0.2
  ) +
  labs(
    title = "GSE184320 Skin: HIV-1 vs Healthy",
    x = "log2 fold change",
    y = "-log10(FDR)",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.position = "right"
  )

print(p_skin)
#ggsave("GSE184320_Skin_volcano_highlight_HIV1_targets.pdf", p_skin, width = 8, height = 4.5)





# 12. GO和KEGG富集分析 - 皮肤样本差异表达的HIV-1靶基因 ----------------------------
# 提取皮肤样本中差异表达的HIV-1靶基因列表
skin_target_deg_genes <- skin_hiv1_target_deg$Gene

if(length(skin_target_deg_genes) >= 3) {
  # 转换为Entrez ID
  skin_entrez <- bitr(
    skin_target_deg_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  skin_entrez <- distinct(skin_entrez, SYMBOL, .keep_all = TRUE)
  skin_entrez_ids <- unique(skin_entrez$ENTREZID)
  
  cat("\n皮肤样本差异表达HIV-1靶基因数：", length(skin_target_deg_genes), "\n")
  cat("成功转换Entrez ID数：", length(skin_entrez_ids), "\n")
  
  # GO BP富集分析
  skin_ego_bp <- enrichGO(
    gene = skin_entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # GO CC富集分析
  skin_ego_cc <- enrichGO(
    gene = skin_entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # GO MF富集分析
  skin_ego_mf <- enrichGO(
    gene = skin_entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # KEGG富集分析
  skin_kegg <- enrichKEGG(
    gene = skin_entrez_ids,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  # 保存结果
  if(!is.null(skin_ego_bp) && nrow(as.data.frame(skin_ego_bp)) > 0) {
    write.csv(as.data.frame(skin_ego_bp), "Skin_HIV1_target_DEG_GO_BP.csv", row.names = FALSE)
  }
  if(!is.null(skin_ego_cc) && nrow(as.data.frame(skin_ego_cc)) > 0) {
    write.csv(as.data.frame(skin_ego_cc), "Skin_HIV1_target_DEG_GO_CC.csv", row.names = FALSE)
  }
  if(!is.null(skin_ego_mf) && nrow(as.data.frame(skin_ego_mf)) > 0) {
    write.csv(as.data.frame(skin_ego_mf), "Skin_HIV1_target_DEG_GO_MF.csv", row.names = FALSE)
  }
  if(!is.null(skin_kegg) && nrow(as.data.frame(skin_kegg)) > 0) {
    write.csv(as.data.frame(skin_kegg), "Skin_HIV1_target_DEG_KEGG.csv", row.names = FALSE)
  }
  
  # 绘制富集分析图
  if(!is.null(skin_ego_bp) && nrow(as.data.frame(skin_ego_bp)) > 0) {
    p_skin_bp <- dotplot(skin_ego_bp, showCategory = 15, title = "Skin: GO BP Enrichment of HIV-1 Target DEGs") +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11)
      )
    ggsave("Skin_HIV1_target_DEG_GO_BP_dotplot.pdf", p_skin_bp, width = 10, height = 8)
    print(p_skin_bp)
  }
  
  if(!is.null(skin_kegg) && nrow(as.data.frame(skin_kegg)) > 0) {
    p_skin_kegg <- dotplot(skin_kegg, showCategory = 15, title = "Skin: KEGG Enrichment of HIV-1 Target DEGs") +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11)
      )
    ggsave("Skin_HIV1_target_DEG_KEGG_dotplot.pdf", p_skin_kegg, width = 10, height = 8)
    print(p_skin_kegg)
  }
} else {
  cat("\n皮肤样本差异表达的HIV-1靶基因不足3个，跳过富集分析\n")
}

# 13. GO和KEGG富集分析 - PBMC样本差异表达的HIV-1靶基因 ----------------------------
pbmc_target_deg_genes <- pbmc_hiv1_target_deg$Gene

if(length(pbmc_target_deg_genes) >= 3) {
  pbmc_entrez <- bitr(
    pbmc_target_deg_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  pbmc_entrez <- distinct(pbmc_entrez, SYMBOL, .keep_all = TRUE)
  pbmc_entrez_ids <- unique(pbmc_entrez$ENTREZID)
  
  cat("\nPBMC样本差异表达HIV-1靶基因数：", length(pbmc_target_deg_genes), "\n")
  cat("成功转换Entrez ID数：", length(pbmc_entrez_ids), "\n")
  
  # GO BP富集分析
  pbmc_ego_bp <- enrichGO(
    gene = pbmc_entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # GO CC富集分析
  pbmc_ego_cc <- enrichGO(
    gene = pbmc_entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # GO MF富集分析
  pbmc_ego_mf <- enrichGO(
    gene = pbmc_entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # KEGG富集分析
  pbmc_kegg <- enrichKEGG(
    gene = pbmc_entrez_ids,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  # 保存结果
  if(!is.null(pbmc_ego_bp) && nrow(as.data.frame(pbmc_ego_bp)) > 0) {
    write.csv(as.data.frame(pbmc_ego_bp), "PBMC_HIV1_target_DEG_GO_BP.csv", row.names = FALSE)
  }
  if(!is.null(pbmc_ego_cc) && nrow(as.data.frame(pbmc_ego_cc)) > 0) {
    write.csv(as.data.frame(pbmc_ego_cc), "PBMC_HIV1_target_DEG_GO_CC.csv", row.names = FALSE)
  }
  if(!is.null(pbmc_ego_mf) && nrow(as.data.frame(pbmc_ego_mf)) > 0) {
    write.csv(as.data.frame(pbmc_ego_mf), "PBMC_HIV1_target_DEG_GO_MF.csv", row.names = FALSE)
  }
  if(!is.null(pbmc_kegg) && nrow(as.data.frame(pbmc_kegg)) > 0) {
    write.csv(as.data.frame(pbmc_kegg), "PBMC_HIV1_target_DEG_KEGG.csv", row.names = FALSE)
  }
  
  # 绘制富集分析图
  if(!is.null(pbmc_ego_bp) && nrow(as.data.frame(pbmc_ego_bp)) > 0) {
    p_pbmc_bp <- dotplot(pbmc_ego_bp, showCategory = 15, title = "PBMC: GO BP Enrichment of HIV-1 Target DEGs") +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11)
      )
    ggsave("PBMC_HIV1_target_DEG_GO_BP_dotplot.pdf", p_pbmc_bp, width = 10, height = 8)
    print(p_pbmc_bp)
  }
  
  if(!is.null(pbmc_kegg) && nrow(as.data.frame(pbmc_kegg)) > 0) {
    p_pbmc_kegg <- dotplot(pbmc_kegg, showCategory = 15, title = "PBMC: KEGG Enrichment of HIV-1 Target DEGs") +
      theme_bw(base_size = 11) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text = element_text(size = 10, color = "black"),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 12, hjust = 0.5),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 11)
      )
    ggsave("PBMC_HIV1_target_DEG_KEGG_dotplot.pdf", p_pbmc_kegg, width = 10, height = 8)
    print(p_pbmc_kegg)
  }
} else {
  cat("\nPBMC样本差异表达的HIV-1靶基因不足3个，跳过富集分析\n")
}

# 14. 热图 - 差异表达的HIV-1靶基因表达模式 ----------------------------------------
# 合并皮肤和PBMC中所有差异表达的HIV-1靶基因
all_target_deg <- unique(c(skin_target_deg_genes, pbmc_target_deg_genes))

if(length(all_target_deg) > 1) {
  # 提取这些基因在所有样本中的表达
  target_deg_expr <- rpm_matrix_log[rownames(rpm_matrix_log) %in% all_target_deg, ]
  
  # 标准化表达值（Z-score）
  target_deg_expr_z <- t(scale(t(target_deg_expr)))
  
  # 创建注释信息
  annotation_col <- data.frame(
    Sample = colnames(target_deg_expr_z),
    Group = ifelse(grepl("^HC_", colnames(target_deg_expr_z)), "Healthy", "HIV"),
    Tissue = ifelse(grepl("SKIN$", colnames(target_deg_expr_z)), "Skin", "PBMC")
  )
  rownames(annotation_col) <- colnames(target_deg_expr_z)
  
  # 保存表达数据用于外部热图绘制
  write.csv(target_deg_expr_z, "HIV1_target_DEG_expression_matrix_zscore.csv", row.names = TRUE)
  
  cat("\n差异表达HIV-1靶基因数：", length(all_target_deg), "\n")
  cat("基因列表：", paste(all_target_deg, collapse = ", "), "\n")
}
cat("\n========== 分析完成 ==========\n")









































######################## GSE184320 HIV-1 使用DESeq2进行差异表达分析 ########################

######################## GSE156072 HIV-1 CD4和MDM细胞分析 ########################

######################## GSE156072 HIV-1 CD4和MDM细胞分析（修正版）########################

# 设置工作目录
setwd("D:/0work/0wholetransctiptome/2sRNAminic/17HIV1/GEO_HIV_bulk")

# 加载必要的包
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pheatmap)

# 1. 读取原始counts数据 ------------------------------------------------------------
counts_data <- read.delim(
  "GSE156072_MDM_CD4_HIV_data.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("原始数据维度：", dim(counts_data), "\n")
cat("列名：", colnames(counts_data), "\n")
cat("前5行数据：\n")
print(head(counts_data[, 1:5]))

# 2. 转换Ensembl ID为Gene Symbol ------------------------------------------------
# 提取Ensembl ID（第一列）
ensembl_ids <- counts_data$Gene

# 去掉版本号
ensembl_ids_clean <- gsub("\\..*$", "", ensembl_ids)

# 创建数据框，包含原始Ensembl ID和表达数据
expr_df <- data.frame(
  ENSEMBL_original = ensembl_ids,
  ENSEMBL_clean = ensembl_ids_clean,
  counts_data[, -1],
  stringsAsFactors = FALSE
)

# 使用bitr转换ID
cat("\n正在将Ensembl ID转换为Gene Symbol...\n")
gene_conversion <- bitr(
  unique(ensembl_ids_clean),
  fromType = "ENSEMBL",
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)

cat("成功转换的基因数：", nrow(gene_conversion), "\n")

# 将转换结果合并到表达数据中
expr_df_converted <- merge(
  expr_df,
  gene_conversion,
  by.x = "ENSEMBL_clean",
  by.y = "ENSEMBL",
  all = FALSE  # 只保留成功转换的基因
)

cat("转换后保留的基因数：", nrow(expr_df_converted), "\n")

# 提取表达矩阵（使用Gene Symbol作为行名）
expr_matrix <- as.matrix(expr_df_converted[, grepl("\\.\\d+$", colnames(expr_df_converted)) | 
                                             grepl("CD4\\.", colnames(expr_df_converted)) |
                                             grepl("MDM\\.", colnames(expr_df_converted))])
rownames(expr_matrix) <- expr_df_converted$SYMBOL

# 确保数据为整数
mode(expr_matrix) <- "integer"

cat("\n表达矩阵维度：", dim(expr_matrix), "\n")
cat("前20个基因名：", paste(head(rownames(expr_matrix), 20), collapse = ", "), "\n")

# 检查是否有重复的Gene Symbol
if(any(duplicated(rownames(expr_matrix)))) {
  cat("\n发现重复的Gene Symbol，正在合并（求和）...\n")
  # 对于重复的基因，将counts相加
  expr_df_temp <- as.data.frame(expr_matrix)
  expr_df_temp$Gene <- rownames(expr_df_temp)
  
  expr_aggregated <- aggregate(. ~ Gene, data = expr_df_temp, FUN = sum)
  rownames(expr_aggregated) <- expr_aggregated$Gene
  expr_aggregated$Gene <- NULL
  
  expr_matrix <- as.matrix(expr_aggregated)
  mode(expr_matrix) <- "integer"
  cat("合并后基因数：", nrow(expr_matrix), "\n")
  cat("合并后前20个基因名：", paste(head(rownames(expr_matrix), 20), collapse = ", "), "\n")
}

# 3. 构建样本分组信息 ------------------------------------------------------------
sample_names <- colnames(expr_matrix)

# 根据列名解析分组信息
sample_group <- data.frame(
  Sample = sample_names,
  CellType = factor(ifelse(grepl("^CD4", sample_names), "CD4", "MDM"),
                    levels = c("CD4", "MDM")),
  Infection = factor(ifelse(grepl("\\.u\\.", sample_names), "Uninfected", "HIV"),
                     levels = c("Uninfected", "HIV"))
)
rownames(sample_group) <- sample_names

cat("\n样本分组信息：\n")
print(sample_group)

# 4. 过滤低表达基因 ------------------------------------------------------------
# 至少在2个样本中count > 10
keep <- rowSums(expr_matrix > 1) >= 2
counts_filtered <- expr_matrix[keep, ]

cat("\n过滤前基因数：", nrow(expr_matrix), "\n")
cat("过滤后基因数：", nrow(counts_filtered), "\n")
cat("过滤后前20个基因名：", paste(head(rownames(counts_filtered), 20), collapse = ", "), "\n")

# 5. 读取HIV-1共有靶基因信息 ----------------------------------------------------
target_info <- read.csv(
  "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 提取HIV-1共有靶基因
hiv1_target_genes <- target_info %>%
  dplyr::filter(abbreviation == "HIV-1") %>%
  dplyr::pull(Intersection_Genes) %>%
  paste(collapse = ";") %>%
  stringr::str_split(";") %>%
  unlist() %>%
  trimws() %>%
  .[. != "" & !is.na(.)] %>%
  unique()

hiv1_target_key <- toupper(hiv1_target_genes)

cat("\nHIV-1共有靶基因总数：", length(hiv1_target_genes), "\n")
cat("前20个靶基因：", paste(head(hiv1_target_genes, 20), collapse = ", "), "\n")

# 检查过滤后的数据中有多少靶基因
filtered_genes <- rownames(counts_filtered)
filtered_target_genes <- intersect(toupper(filtered_genes), hiv1_target_key)
cat("\n过滤后数据中匹配到的靶基因数：", length(filtered_target_genes), "\n")

if(length(filtered_target_genes) > 0) {
  cat("匹配到的靶基因示例：", paste(head(filtered_target_genes, 20), collapse = ", "), "\n")
  
  # 显示这些靶基因在数据中的实际表达情况
  cat("\n匹配到的靶基因在数据中的名称：\n")
  matched_genes <- filtered_genes[toupper(filtered_genes) %in% hiv1_target_key]
  cat(paste(head(matched_genes, 20), collapse = ", "), "\n")
}

# 6. CD4+ T细胞差异分析 ------------------------------------------------------------
cd4_samples <- sample_group %>% dplyr::filter(CellType == "CD4")
cd4_counts <- counts_filtered[, rownames(cd4_samples), drop = FALSE]

cat("\n========== CD4+ T细胞分析 ==========\n")
cat("样本数：", ncol(cd4_counts), "\n")
cat("样本名：", paste(colnames(cd4_counts), collapse = ", "), "\n")
cat("未感染样本数：", sum(cd4_samples$Infection == "Uninfected"), "\n")
cat("感染样本数：", sum(cd4_samples$Infection == "HIV"), "\n")
cat("CD4数据维度：", dim(cd4_counts), "\n")

# 检查数据是否为整数
cat("数据是否为整数：", all(cd4_counts == floor(cd4_counts)), "\n")

# 创建DESeq2对象
dds_cd4 <- DESeqDataSetFromMatrix(
  countData = cd4_counts,
  colData = cd4_samples,
  design = ~ Infection
)

# 运行DESeq2
dds_cd4 <- DESeq(dds_cd4)

# 提取结果
cd4_res <- results(dds_cd4, 
                   contrast = c("Infection", "HIV", "Uninfected"),
                   alpha = 0.05)

# 转换为数据框并添加基因信息
cd4_res_df <- as.data.frame(cd4_res)
cd4_res_df$Gene <- rownames(cd4_res_df)
cd4_res_df <- cd4_res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    GeneKey = toupper(Gene),
    is_HIV1_target = GeneKey %in% hiv1_target_key,
    DEG_status = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up in HIV",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in HIV",
      TRUE ~ "Not significant"
    )
  )

# 7. MDM细胞差异分析 ------------------------------------------------------------
mdm_samples <- sample_group %>% dplyr::filter(CellType == "MDM")
mdm_counts <- counts_filtered[, rownames(mdm_samples), drop = FALSE]

cat("\n========== MDM细胞分析 ==========\n")
cat("样本数：", ncol(mdm_counts), "\n")
cat("样本名：", paste(colnames(mdm_counts), collapse = ", "), "\n")
cat("未感染样本数：", sum(mdm_samples$Infection == "Uninfected"), "\n")
cat("感染样本数：", sum(mdm_samples$Infection == "HIV"), "\n")
cat("MDM数据维度：", dim(mdm_counts), "\n")

# 创建DESeq2对象
dds_mdm <- DESeqDataSetFromMatrix(
  countData = mdm_counts,
  colData = mdm_samples,
  design = ~ Infection
)

# 运行DESeq2
dds_mdm <- DESeq(dds_mdm)

# 提取结果
mdm_res <- results(dds_mdm, 
                   contrast = c("Infection", "HIV", "Uninfected"),
                   alpha = 0.05)

# 转换为数据框并添加基因信息
mdm_res_df <- as.data.frame(mdm_res)
mdm_res_df$Gene <- rownames(mdm_res_df)
mdm_res_df <- mdm_res_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    GeneKey = toupper(Gene),
    is_HIV1_target = GeneKey %in% hiv1_target_key,
    DEG_status = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up in HIV",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in HIV",
      TRUE ~ "Not significant"
    )
  )

# 8. 保存结果 ------------------------------------------------------------
write.csv(cd4_res_df, "GSE156072_CD4_DESeq2_results.csv", row.names = FALSE)
write.csv(mdm_res_df, "GSE156072_MDM_DESeq2_results.csv", row.names = FALSE)

# 9. 提取HIV-1靶基因的差异表达结果 ------------------------------------------------
# CD4细胞
cd4_hiv1_target_res <- cd4_res_df %>% dplyr::filter(is_HIV1_target)
cd4_hiv1_target_deg <- cd4_hiv1_target_res %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1)

# MDM细胞
mdm_hiv1_target_res <- mdm_res_df %>% dplyr::filter(is_HIV1_target)
mdm_hiv1_target_deg <- mdm_hiv1_target_res %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1)

# 导出
write.csv(cd4_hiv1_target_res, "GSE156072_CD4_HIV1_target_genes_all.csv", row.names = FALSE)
write.csv(cd4_hiv1_target_deg, "GSE156072_CD4_HIV1_target_genes_DEG.csv", row.names = FALSE)
write.csv(mdm_hiv1_target_res, "GSE156072_MDM_HIV1_target_genes_all.csv", row.names = FALSE)
write.csv(mdm_hiv1_target_deg, "GSE156072_MDM_HIV1_target_genes_DEG.csv", row.names = FALSE)

# 10. 统计信息 ------------------------------------------------------------
cat("\n========== CD4+ T细胞分析统计 ==========\n")
cat("总基因数：", nrow(cd4_res_df), "\n")
cat("差异表达基因数(padj<0.05, |log2FC|>1)：", 
    sum(cd4_res_df$padj < 0.05 & abs(cd4_res_df$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("上调基因数：", 
    sum(cd4_res_df$padj < 0.05 & cd4_res_df$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("下调基因数：", 
    sum(cd4_res_df$padj < 0.05 & cd4_res_df$log2FoldChange < -1, na.rm = TRUE), "\n")
cat("表达的HIV-1靶基因数：", nrow(cd4_hiv1_target_res), "\n")
cat("差异表达的HIV-1靶基因数：", nrow(cd4_hiv1_target_deg), "\n")

if(nrow(cd4_hiv1_target_deg) > 0) {
  cat("\n差异表达靶基因列表：\n")
  print(cd4_hiv1_target_deg[, c("Gene", "log2FoldChange", "padj")])
}

cat("\n========== MDM细胞分析统计 ==========\n")
cat("总基因数：", nrow(mdm_res_df), "\n")
cat("差异表达基因数(padj<0.05, |log2FC|>1)：", 
    sum(mdm_res_df$padj < 0.05 & abs(mdm_res_df$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("上调基因数：", 
    sum(mdm_res_df$padj < 0.05 & mdm_res_df$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("下调基因数：", 
    sum(mdm_res_df$padj < 0.05 & mdm_res_df$log2FoldChange < -1, na.rm = TRUE), "\n")
cat("表达的HIV-1靶基因数：", nrow(mdm_hiv1_target_res), "\n")
cat("差异表达的HIV-1靶基因数：", nrow(mdm_hiv1_target_deg), "\n")

if(nrow(mdm_hiv1_target_deg) > 0) {
  cat("\n差异表达靶基因列表：\n")
  print(mdm_hiv1_target_deg[, c("Gene", "log2FoldChange", "padj")])
}

# 11. 火山图 - CD4细胞 ---------------------------------------------------------
if(nrow(cd4_res_df) > 0) {
  # 处理可能的padj为0的情况
  padj_values <- cd4_res_df$padj[cd4_res_df$padj > 0]
  if(length(padj_values) > 0) {
    min_padj_cd4 <- min(padj_values)
  } else {
    min_padj_cd4 <- 1e-10
  }
  
  cd4_plot <- cd4_res_df %>%
    dplyr::mutate(
      padj_plot = ifelse(padj == 0, min_padj_cd4, padj),
      neg_log10_FDR = -log10(padj_plot),
      plot_group = dplyr::case_when(
        is_HIV1_target & padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated target genes",
        is_HIV1_target & padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated target genes",
        !is_HIV1_target & padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated other genes",
        !is_HIV1_target & padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated other genes",
        TRUE ~ "Not significant"
      )
    )
  
  # 标注显著差异的HIV-1靶基因
  label_genes_cd4 <- cd4_plot %>%
    dplyr::filter(is_HIV1_target, padj < 0.05, abs(log2FoldChange) > 1) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 30)
  
  p_cd4 <- ggplot(cd4_plot, aes(x = log2FoldChange, y = neg_log10_FDR)) +
    geom_point(data = cd4_plot %>% dplyr::filter(plot_group == "Not significant"),
               color = "grey80", size = 0.8, alpha = 0.5) +
    geom_point(data = cd4_plot %>% dplyr::filter(plot_group %in% c("Up-regulated other genes", "Down-regulated other genes")),
               aes(color = plot_group), size = 1.5, alpha = 0.6) +
    geom_point(data = cd4_plot %>% dplyr::filter(plot_group %in% c("Up-regulated target genes", "Down-regulated target genes")),
               aes(color = plot_group), size = 2, alpha = 0.95) +
    scale_color_manual(values = c("Up-regulated other genes" = "#F4A6A6",
                                  "Down-regulated other genes" = "#A6C8E8",
                                  "Up-regulated target genes" = "#E21C21",
                                  "Down-regulated target genes" = "#3A7CB5")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_text_repel(data = label_genes_cd4, aes(label = Gene), size = 3, max.overlaps = 100) +
    labs(title = "GSE156072 CD4+ T cells: HIV vs Uninfected",
         x = "log2 fold change", y = "-log10(FDR)", color = NULL) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  
  print(p_cd4)
  ggsave("GSE156072_CD4_DESeq2_volcano.pdf", p_cd4, width = 9, height = 5)
}

# 12. 火山图 - MDM细胞 ---------------------------------------------------------
if(nrow(mdm_res_df) > 0) {
  padj_values <- mdm_res_df$padj[mdm_res_df$padj > 0]
  if(length(padj_values) > 0) {
    min_padj_mdm <- min(padj_values)
  } else {
    min_padj_mdm <- 1e-10
  }
  
  mdm_plot <- mdm_res_df %>%
    dplyr::mutate(
      padj_plot = ifelse(padj == 0, min_padj_mdm, padj),
      neg_log10_FDR = -log10(padj_plot),
      plot_group = dplyr::case_when(
        is_HIV1_target & padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated target genes",
        is_HIV1_target & padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated target genes",
        !is_HIV1_target & padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated other genes",
        !is_HIV1_target & padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated other genes",
        TRUE ~ "Not significant"
      )
    )
  
  label_genes_mdm <- mdm_plot %>%
    dplyr::filter(is_HIV1_target, padj < 0.05, abs(log2FoldChange) > 1) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 30)
  
  p_mdm <- ggplot(mdm_plot, aes(x = log2FoldChange, y = neg_log10_FDR)) +
    geom_point(data = mdm_plot %>% dplyr::filter(plot_group == "Not significant"),
               color = "grey80", size = 0.8, alpha = 0.5) +
    geom_point(data = mdm_plot %>% dplyr::filter(plot_group %in% c("Up-regulated other genes", "Down-regulated other genes")),
               aes(color = plot_group), size = 1.5, alpha = 0.6) +
    geom_point(data = mdm_plot %>% dplyr::filter(plot_group %in% c("Up-regulated target genes", "Down-regulated target genes")),
               aes(color = plot_group), size = 2, alpha = 0.95) +
    scale_color_manual(values = c("Up-regulated other genes" = "#F4A6A6",
                                  "Down-regulated other genes" = "#A6C8E8",
                                  "Up-regulated target genes" = "#E21C21",
                                  "Down-regulated target genes" = "#3A7CB5")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_text_repel(data = label_genes_mdm, aes(label = Gene), size = 3, max.overlaps = 100) +
    labs(title = "GSE156072 MDM cells: HIV vs Uninfected",
         x = "log2 fold change", y = "-log10(FDR)", color = NULL) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  
  print(p_mdm)
  ggsave("GSE156072_MDM_DESeq2_volcano.pdf", p_mdm, width = 9, height = 5)
}

cat("\n========== 分析完成 ==========\n")
















##############################################GSE289893########################################################


######################## GSE289893 HIV-1 感染CD4 T细胞分析 ########################

# 设置工作目录
setwd("D:/0work/0wholetransctiptome/2sRNAminic/17HIV1/GEO_HIV_bulk")

# 加载必要的包
library(DESeq2)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(clusterProfiler)
library(pheatmap)

# 1. 读取原始counts数据 ------------------------------------------------------------
counts_data <- read.delim(
  "GSE289893_Gene_Rsem_Counts.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("原始数据维度：", dim(counts_data), "\n")
cat("列名：", colnames(counts_data), "\n")
cat("前5行数据：\n")
print(head(counts_data[, 1:4]))

# 2. 提取基因名（从转录本ID中提取）------------------------------------------------
# Gene列包含多个转录本，格式如"ENST00000373020_TSPAN6-201,ENST00000494424_TSPAN6-202,..."
# 我们需要提取基因名（如TSPAN6）

extract_gene_name <- function(gene_string) {
  # 提取第一个转录本中的基因名
  first_transcript <- strsplit(gene_string, ",")[[1]][1]
  # 提取下划线后面的部分，然后去掉版本号
  gene_name <- str_extract(first_transcript, "_([A-Z0-9]+)-?\\d*")
  gene_name <- gsub("_", "", gene_name)
  gene_name <- gsub("-\\d+$", "", gene_name)
  return(gene_name)
}

# 应用提取函数
gene_names <- sapply(counts_data$Gene, extract_gene_name)
gene_names <- trimws(gene_names)

cat("\n提取的基因名示例：", paste(head(gene_names, 10), collapse = ", "), "\n")

# 检查是否有空值或NA
if(any(is.na(gene_names)) | any(gene_names == "")) {
  cat("警告：部分基因名提取失败，请检查格式\n")
  # 手动修正常见情况
  gene_names[gene_names == ""] <- paste0("Gene_", which(gene_names == ""))
}

# 3. 准备表达矩阵 ------------------------------------------------------------
# 提取表达数据（从第2列开始）
expr_data <- counts_data[, -1]  # 移除Gene列
expr_matrix <- as.matrix(expr_data)
mode(expr_matrix) <- "integer"

# 设置行名为基因名
rownames(expr_matrix) <- gene_names

cat("\n表达矩阵维度：", dim(expr_matrix), "\n")
cat("前20个基因名：", paste(head(rownames(expr_matrix), 20), collapse = ", "), "\n")

# 检查是否有重复的Gene Symbol
if(any(duplicated(rownames(expr_matrix)))) {
  cat("\n发现重复的Gene Symbol，正在合并（求和）...\n")
  # 对于重复的基因，将counts相加
  expr_df_temp <- as.data.frame(expr_matrix)
  expr_df_temp$Gene <- rownames(expr_df_temp)
  
  expr_aggregated <- aggregate(. ~ Gene, data = expr_df_temp, FUN = sum)
  rownames(expr_aggregated) <- expr_aggregated$Gene
  expr_aggregated$Gene <- NULL
  
  expr_matrix <- as.matrix(expr_aggregated)
  mode(expr_matrix) <- "integer"
  cat("合并后基因数：", nrow(expr_matrix), "\n")
  cat("合并后前20个基因名：", paste(head(rownames(expr_matrix), 20), collapse = ", "), "\n")
}

# 4. 构建样本分组信息 ------------------------------------------------------------
sample_names <- colnames(expr_matrix)

# 解析样本信息
# 样本命名规则：1_1_Control-0h_Donor-1_S19_Count 等
sample_group <- data.frame(
  Sample = sample_names,
  Condition = factor(ifelse(grepl("Control", sample_names), "Control",
                            ifelse(grepl("UI-CD44", sample_names), "UI-CD44", "HIV")),
                     levels = c("Control", "UI-CD44", "HIV")),
  TimePoint = factor(ifelse(grepl("0h", sample_names), "0h", "24h"),
                     levels = c("0h", "24h")),
  Donor = factor(gsub(".*Donor-([0-9]+).*", "\\1", sample_names)),
  stringsAsFactors = FALSE
)
rownames(sample_group) <- sample_names

cat("\n样本分组信息：\n")
print(sample_group)

# 5. 过滤低表达基因 ------------------------------------------------------------
# 至少在2个样本中count > 10
keep <- rowSums(expr_matrix > 10) >= 2
counts_filtered <- expr_matrix[keep, ]

cat("\n过滤前基因数：", nrow(expr_matrix), "\n")
cat("过滤后基因数：", nrow(counts_filtered), "\n")
cat("过滤后前20个基因名：", paste(head(rownames(counts_filtered), 20), collapse = ", "), "\n")

# 6. 读取HIV-1共有靶基因信息 ----------------------------------------------------
target_info <- read.csv(
  "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 提取HIV-1共有靶基因
hiv1_target_genes <- target_info %>%
  dplyr::filter(abbreviation == "HIV-1") %>%
  dplyr::pull(Intersection_Genes) %>%
  paste(collapse = ";") %>%
  stringr::str_split(";") %>%
  unlist() %>%
  trimws() %>%
  .[. != "" & !is.na(.)] %>%
  unique()

hiv1_target_key <- toupper(hiv1_target_genes)

cat("\nHIV-1共有靶基因总数：", length(hiv1_target_genes), "\n")
cat("前20个靶基因：", paste(head(hiv1_target_genes, 20), collapse = ", "), "\n")

# 检查过滤后的数据中有多少靶基因
filtered_genes <- rownames(counts_filtered)
filtered_target_genes <- intersect(toupper(filtered_genes), hiv1_target_key)
cat("\n过滤后数据中匹配到的靶基因数：", length(filtered_target_genes), "\n")

if(length(filtered_target_genes) > 0) {
  cat("匹配到的靶基因示例：", paste(head(filtered_target_genes, 20), collapse = ", "), "\n")
}

# 7. HIV感染 vs Control 差异分析（24h时间点）----------------------------------------
# 筛选24h的样本
samples_24h <- sample_group %>% dplyr::filter(Condition %in% c("Control", "HIV"))
counts_24h <- counts_filtered[, rownames(samples_24h), drop = FALSE]

# 更新分组信息
samples_24h$Comparison <- factor(samples_24h$Condition, levels = c("Control", "HIV"))

cat("\n========== HIV vs Control (24h) 分析 ==========\n")
cat("样本数：", ncol(counts_24h), "\n")
cat("Control样本：", sum(samples_24h$Condition == "Control"), "\n")
cat("HIV样本：", sum(samples_24h$Condition == "HIV"), "\n")
cat("样本名：", paste(colnames(counts_24h), collapse = ", "), "\n")

# 创建DESeq2对象
dds_24h <- DESeqDataSetFromMatrix(
  countData = counts_24h,
  colData = samples_24h,
  design = ~ Donor + Condition  # 考虑供体差异
)

# 运行DESeq2
dds_24h <- DESeq(dds_24h)

# 提取结果
res_24h <- results(dds_24h, 
                   contrast = c("Condition", "HIV", "Control"),
                   alpha = 0.05)

# 转换为数据框并添加基因信息
res_24h_df <- as.data.frame(res_24h)
res_24h_df$Gene <- rownames(res_24h_df)
res_24h_df <- res_24h_df %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(
    GeneKey = toupper(Gene),
    is_HIV1_target = GeneKey %in% hiv1_target_key,
    DEG_status = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 1 ~ "Up in HIV",
      padj < 0.05 & log2FoldChange < -1 ~ "Down in HIV",
      TRUE ~ "Not significant"
    )
  )


# 10. 保存主要结果（HIV vs Control 24h）-------------------------------------------
write.csv(res_24h_df, "GSE289893_HIV_vs_Control_24h_DESeq2_results.csv", row.names = FALSE)

# 11. 提取HIV-1靶基因的差异表达结果 ------------------------------------------------
hiv1_target_res <- res_24h_df %>% dplyr::filter(is_HIV1_target)
hiv1_target_deg <- hiv1_target_res %>%
  dplyr::filter(padj < 0.05, abs(log2FoldChange) > 1)

# 导出
write.csv(hiv1_target_res, "GSE289893_HIV1_target_genes_all.csv", row.names = FALSE)
write.csv(hiv1_target_deg, "GSE289893_HIV1_target_genes_DEG.csv", row.names = FALSE)

# 12. 统计信息 ------------------------------------------------------------
cat("\n========== HIV vs Control (24h) 分析统计 ==========\n")
cat("总基因数：", nrow(res_24h_df), "\n")
cat("差异表达基因数(padj<0.05, |log2FC|>1)：", 
    sum(res_24h_df$padj < 0.05 & abs(res_24h_df$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("上调基因数：", 
    sum(res_24h_df$padj < 0.05 & res_24h_df$log2FoldChange > 1, na.rm = TRUE), "\n")
cat("下调基因数：", 
    sum(res_24h_df$padj < 0.05 & res_24h_df$log2FoldChange < -1, na.rm = TRUE), "\n")
cat("表达的HIV-1靶基因数：", nrow(hiv1_target_res), "\n")
cat("差异表达的HIV-1靶基因数：", nrow(hiv1_target_deg), "\n")

if(nrow(hiv1_target_deg) > 0) {
  cat("\n差异表达靶基因列表：\n")
  print(hiv1_target_deg[, c("Gene", "log2FoldChange", "padj")])
}







# 13. 火山图 - HIV vs Control 24h -------------------------------------------------
if(nrow(res_24h_df) > 0) {
  # 处理可能的padj为0的情况
  padj_values <- res_24h_df$padj[res_24h_df$padj > 0]
  if(length(padj_values) > 0) {
    min_padj <- min(padj_values)
  } else {
    min_padj <- 1e-10
  }
  
  # 定义你要标注的基因列表
  my_genes <- c("CD226", "CTIF", "JUN", "KLF2", "MAF", 
                "NCKIPSD", "OTUB2", "PPP1R15A", "SIRT2", "SRC", "ZNF598")
  
  plot_df <- res_24h_df %>%
    dplyr::mutate(
      padj_plot = ifelse(padj == 0, min_padj, padj),
      neg_log10_FDR = -log10(padj_plot),
      plot_group = dplyr::case_when(
        is_HIV1_target & padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated target genes",
        is_HIV1_target & padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated target genes",
        !is_HIV1_target & padj < 0.05 & log2FoldChange > 1 ~ "Up-regulated other genes",
        !is_HIV1_target & padj < 0.05 & log2FoldChange < -1 ~ "Down-regulated other genes",
        TRUE ~ "Not significant"
      )
    )
  
  # 修改这里：标注你的特定基因（不管是不是HIV-1靶基因，只要在列表中且显著差异就标注）
  label_genes <- plot_df %>%
    dplyr::filter(Gene %in% my_genes, padj < 0.05, abs(log2FoldChange) > 1) %>%
    dplyr::arrange(padj)
  
  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_FDR)) +
    geom_point(data = plot_df %>% dplyr::filter(plot_group == "Not significant"),
               color = "grey80", size = 0.8, alpha = 0.5) +
    geom_point(data = plot_df %>% dplyr::filter(plot_group %in% c("Up-regulated other genes", "Down-regulated other genes")),
               aes(color = plot_group), size = 1.5, alpha = 0.6) +
    geom_point(data = plot_df %>% dplyr::filter(plot_group %in% c("Up-regulated target genes", "Down-regulated target genes")),
               aes(color = plot_group), size = 1.5, alpha = 0.95) +
    scale_color_manual(values = c("Up-regulated other genes" = "#F4A6A6",
                                  "Down-regulated other genes" = "#A6C8E8",
                                  "Up-regulated target genes" = "#E21C21",
                                  "Down-regulated target genes" = "#3A7CB5")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.4, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.4, color = "black") +
    # 取消注释这行，启用文字标注
    geom_text_repel(data = label_genes, aes(label = Gene), size = 3, max.overlaps = 100,
                    box.padding = 0.5, point.padding = 0.3, fontface = "bold") +
    labs(title = "GSE289893 CD4+ T cells: HIV vs Control (24h)",
         x = "log2(fold change)", 
         y ="-log10(FDR)", 
         color = NULL) +
    theme_bw(base_size = 13) +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right",
          legend.text = element_text(size = 10))
  
  print(p)
  ggsave("GSE289893_HIV_vs_Control_24h_volcano.pdf", p, width = 7, height = 4)
}

















# 14. 热图 - 差异表达的HIV-1靶基因 ------------------------------------------------
# 14. 热图 - 差异表达的HIV-1靶基因（无聚类，原始顺序）--------------------------------
# 14. 热图 - 差异表达的HIV-1靶基因（按Condition排序，添加上调下调标注）----------------
if(nrow(hiv1_target_deg) > 1) {
  cat("\n========== 绘制差异表达靶基因热图 ==========\n")
  cat("差异表达靶基因数：", nrow(hiv1_target_deg), "\n")
  
  # 获取标准化后的counts
  normalized_counts <- counts(dds_24h, normalized = TRUE)
  
  # 提取差异表达靶基因的表达值
  target_expr <- normalized_counts[rownames(normalized_counts) %in% hiv1_target_deg$Gene, ]
  
  # 按照差异表达基因列表的顺序排序
  target_expr <- target_expr[hiv1_target_deg$Gene, ]
  
  if(nrow(target_expr) > 1) {
    # 只进行log2转换
    target_expr_log <- log2(target_expr + 1)
    
    # 按Condition排序样本
    current_order <- colnames(target_expr_log)
    
    # 创建样本排序的数据框
    sample_order_df <- data.frame(
      Sample = current_order,
      Condition = ifelse(grepl("Control", current_order), "Control", "HIV")
    )
    
    # 按Condition排序（Control在前，HIV在后）
    sample_order_df <- sample_order_df %>%
      arrange(Condition)
    
    # 重新排列表达矩阵的列
    target_expr_log_sorted <- target_expr_log[, sample_order_df$Sample]
    
    # 更新列注释信息
    annotation_col <- data.frame(
      Condition = sample_order_df$Condition,
      Donor = gsub(".*Donor-([0-9]+).*", "\\1", sample_order_df$Sample)
    )
    rownames(annotation_col) <- sample_order_df$Sample
    
    # 创建行注释信息（标记上调或下调）
    # 提取每个基因的log2FC和上调/下调状态
    annotation_row <- data.frame(
      Regulation = hiv1_target_deg$DEG_status
    )
    rownames(annotation_row) <- hiv1_target_deg$Gene
    
    # 将Regulation转换为因子并设置级别
    annotation_row$Regulation <- factor(annotation_row$Regulation, 
                                        levels = c("Up in HIV", "Down in HIV"))
    
    # 设置行注释的颜色
    ann_colors <- list(
      Condition = c(Control = "#92c5de", HIV = "#F39B7F"),
      Donor = c(`1` = "#aecfd4", `2` = "#fff2cc"),
      Regulation = c("Up in HIV" = "#fbe5d6", "Down in HIV" =  "#bdd7ee")
    )
    cat("\n样本排序顺序：\n")
    print(sample_order_df)
    
    cat("\n基因调控方向统计：\n")
    print(table(annotation_row$Regulation))
    blue_gradient <- colorRampPalette(c("#FBF9FA", "steelblue"))(40)
    
    # 绘制热图（添加上下调行注释）
    pheatmap(target_expr_log_sorted,
             scale = "none",
             cluster_rows = FALSE,  # 不对基因聚类
             cluster_cols = FALSE,  # 不对样本聚类
             annotation_col = annotation_col,
             annotation_row = annotation_row,  # 添加行注释
             color = blue_gradient,
             annotation_colors = ann_colors,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize_row = 6,  
             fontsize_col = 6,
             border_color = "white",
             main = "",
             filename = "GSE289893_HIV1_target_DEG_heatmap.pdf",
             width = 3.8,  # 稍微增加宽度以容纳行注释
             height = max(6, nrow(target_expr_log_sorted) * 0.1))
    
    # 保存表达数据
    write.csv(target_expr_log_sorted, "GSE289893_HIV1_target_DEG_expression_log2.csv", row.names = TRUE)
    
    cat("\n热图已保存为: GSE289893_HIV1_target_DEG_heatmap.pdf\n")
    cat("表达数据已保存为: GSE289893_HIV1_target_DEG_expression_log2.csv\n")
  }
}





# 15. GO富集分析（如果差异表达靶基因足够）----------------------------------------
  cat("\n========== GO富集分析 ==========\n")
  # 转换为Entrez ID
  deg_entrez <- bitr(
    hiv1_target_deg$Gene,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  
    ego_bp <- enrichGO(
      gene = deg_entrez$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "MF",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.05,
      readable = TRUE
    )
    
    ego_bp






    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ##################################miRNA的差异表达分析（所有miRNA，不过滤）######################################
    library(limma)
    library(edgeR)
    library(dplyr) 
    # =========================
    # 1. Input files
    # =========================
    
    base_dir <- "D:/0work/0wholetransctiptome/2sRNAminic/13Infect_mimic_correlation/Time"
    
    control_files <- c(
      file.path(base_dir, "HIV-1_12hpi_control1.txt"),
      file.path(base_dir, "HIV-1_12hpi_control2.txt"),
      file.path(base_dir, "HIV-1_12hpi_control3.txt")
    )
    
    hiv1_files <- c(
      file.path(base_dir, "HIV-1_12hpi_1.txt"),
      file.path(base_dir, "HIV-1_12hpi_2.txt"),
      file.path(base_dir, "HIV-1_12hpi_3.txt")
    )
    
    out_file <- file.path(base_dir, "HIV-1_12hpi_infected_vs_control_limma_DEG_all_miRNAs.csv")
    
    mimic_file <- "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv"
    
    # =========================
    # 2. Function to read miRNA count file
    # No header
    # Column 1 = miRNA
    # Column 3 = count
    # 修改：读取后对所有count加1
    # =========================
    
    read_miRNA_count <- function(file, sample_name) {
      
      dat <- read.table(
        file,
        header = FALSE,
        sep = "\t",
        stringsAsFactors = FALSE,
        quote = "",
        comment.char = ""
      )
      
      dat <- dat[, c(1, 3)]
      colnames(dat) <- c("miRNA", sample_name)
      
      dat$miRNA <- trimws(dat$miRNA)
      dat[[sample_name]] <- as.numeric(dat[[sample_name]])
      
      # 将NA替换为0
      dat[[sample_name]][is.na(dat[[sample_name]])] <- 0
      
      # 关键修改：对所有count加1（包括control和HIV样本）
      dat[[sample_name]] <- dat[[sample_name]] + 1
      
      dat <- aggregate(
        dat[[sample_name]],
        by = list(miRNA = dat$miRNA),
        FUN = sum
      )
      
      colnames(dat) <- c("miRNA", sample_name)
      return(dat)
    }
    
    # =========================
    # 3. Read all samples
    # =========================
    
    sample_names <- c(
      "Control_1",
      "Control_2",
      "Control_3",
      "HIV1_1",
      "HIV1_2",
      "HIV1_3"
    )
    
    all_files <- c(control_files, hiv1_files)
    
    dat_list <- mapply(
      FUN = read_miRNA_count,
      file = all_files,
      sample_name = sample_names,
      SIMPLIFY = FALSE
    )
    
    # =========================
    # 4. Merge all samples
    # =========================
    
    dat <- Reduce(
      function(x, y) merge(x, y, by = "miRNA", all = TRUE),
      dat_list
    )
    
    # 合并后可能还有NA（如果某个miRNA只出现在部分样本中），将其替换为1
    dat[is.na(dat)] <- 1
    
    count_mat <- as.matrix(dat[, sample_names])
    rownames(count_mat) <- dat$miRNA
    
    # 验证：检查是否有0值
    cat("检查count矩阵中是否有0值：", any(count_mat == 0), "\n")
    cat("count矩阵最小值：", min(count_mat), "\n")
    cat("总miRNA数：", nrow(count_mat), "\n")
    
    # =========================
    # 5. Build group information
    # =========================
    
    group <- factor(
      c(rep("Control", 3), rep("HIV1", 3)),
      levels = c("Control", "HIV1")
    )
    
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    
    # =========================
    # 6. 计算所有miRNA的logFC（不进行过滤）
    # =========================
    
    # 方法1：直接计算logFC（基于加1后的counts）
    cat("\n========== 计算所有miRNA的logFC ==========\n")
    
    # 计算每个miRNA在control和HIV组中的平均count
    mean_control_all <- rowMeans(count_mat[, grep("^Control", colnames(count_mat)), drop = FALSE])
    mean_HIV1_all <- rowMeans(count_mat[, grep("^HIV1", colnames(count_mat)), drop = FALSE])
    
    # 计算log2 fold change
    log2FC_all <- log2(mean_HIV1_all / mean_control_all)
    
    # 处理无穷值（如果control均值为0，但我们已经加了1，所以不会出现无穷）
    log2FC_all[is.infinite(log2FC_all)] <- NA
    log2FC_all[is.na(log2FC_all)] <- 0
    
    # 创建包含所有miRNA的结果数据框
    all_mirnas_res <- data.frame(
      miRNA = rownames(count_mat),
      mean_count_control = mean_control_all,
      mean_count_HIV1 = mean_HIV1_all,
      logFC = log2FC_all,
      stringsAsFactors = FALSE
    )
    
    # =========================
    # 7. 使用limma-voom进行差异表达分析（保留所有miRNA）
    # =========================
    
    cat("\n========== 运行limma-voom分析（保留所有miRNA） ==========\n")
    
    # 创建DGEList对象
    dge_all <- DGEList(counts = count_mat, group = group)
    dge_all <- calcNormFactors(dge_all)
    
    # 重要：跳过filterByExpr过滤，直接使用所有miRNA
    # 但voom要求至少有一定数量的基因，如果miRNA太少会有问题
    # 这里我们保留所有miRNA
    
    # 添加一个小的伪计数来避免voom时的数值问题（虽然我们已经加了1）
    # voom会进行log2转换，所以没问题
    v_all <- voom(dge_all, design, plot = TRUE)
    
    # 线性拟合
    fit_all <- lmFit(v_all, design)
    
    # 对比矩阵
    contrast.matrix_all <- makeContrasts(
      HIV1_vs_Control = HIV1 - Control,
      levels = design
    )
    
    fit2_all <- contrasts.fit(fit_all, contrast.matrix_all)
    fit2_all <- eBayes(fit2_all, trend = TRUE)
    
    # 提取所有miRNA的结果
    res_all <- topTable(
      fit2_all,
      coef = "HIV1_vs_Control",
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
    
    # 添加miRNA列
    res_all$miRNA <- rownames(res_all)
    res_all$miRNA <- trimws(res_all$miRNA)
    
    # =========================
    # 8. 合并两种方法的结果并添加注释
    # =========================
    
    # 将limma结果与直接计算的logFC合并
    final_res <- merge(all_mirnas_res, res_all, by = "miRNA", all = TRUE, suffixes = c("_manual", "_limma"))
    
    # 使用limma的logFC（更可靠），如果limma没有结果则使用手动计算的
    final_res$logFC_final <- ifelse(is.na(final_res$logFC_limma), 
                                    final_res$logFC_manual, 
                                    final_res$logFC_limma)
    
    # 提取需要的列
    final_res <- final_res %>%
      dplyr::mutate(
        mean_count_control = mean_count_control,
        mean_count_HIV1 = mean_count_HIV1,
        logFC = logFC_final,
        AveExpr = ifelse(is.na(AveExpr), (log2(mean_count_control) + log2(mean_count_HIV1))/2, AveExpr),
        t = ifelse(is.na(t), NA, t),
        P.Value = ifelse(is.na(P.Value), 1, P.Value),
        adj.P.Val = ifelse(is.na(adj.P.Val), 1, adj.P.Val),
        B = ifelse(is.na(B), NA, B),
        # 标记是否被limma过滤（原本会被过滤的miRNA）
        was_filtered = is.na(logFC_limma)
      )
    
    # 添加DEG标签
    final_res$DEG_status <- "Not_DEG"
    final_res$DEG_status[final_res$logFC > 1 & final_res$adj.P.Val < 0.05 & !final_res$was_filtered] <- "Up_in_HIV1"
    final_res$DEG_status[final_res$logFC < -1 & final_res$adj.P.Val < 0.05 & !final_res$was_filtered] <- "Down_in_HIV1"
    final_res$DEG_status[final_res$was_filtered] <- "Low_expression_filtered"
    
    # =========================
    # 9. Annotate whether the human miRNA is mimicked by HIV-1 only
    # =========================
    
    mimic_dat <- read.csv(
      mimic_file,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    # Check required columns
    required_cols <- c("abbreviation", "miRNA")
    missing_cols <- setdiff(required_cols, colnames(mimic_dat))
    
    if (length(missing_cols) > 0) {
      warning(
        paste(
          "Missing required columns in mimic file:",
          paste(missing_cols, collapse = ", ")
        )
      )
      final_res$is_mimicked_by_HIV1 <- "No"
    } else {
      # abbreviation = virus name or viral miRNA annotation
      # miRNA = human miRNA mimicked by the virus
      mimic_dat <- mimic_dat[, c("abbreviation", "miRNA")]
      
      mimic_dat$abbreviation <- trimws(mimic_dat$abbreviation)
      mimic_dat$miRNA <- trimws(mimic_dat$miRNA)
      
      mimic_dat <- mimic_dat[
        !is.na(mimic_dat$abbreviation) & mimic_dat$abbreviation != "" &
          !is.na(mimic_dat$miRNA) & mimic_dat$miRNA != "",
      ]
      
      # Keep HIV-1 records only
      mimic_dat_HIV1 <- mimic_dat[
        mimic_dat$abbreviation %in% c(
          "HIV-1",
          "HIV1",
          "Human immunodeficiency virus 1"
        ) |
          grepl("^HIV1_", mimic_dat$abbreviation) |
          grepl("^HIV-1_", mimic_dat$abbreviation),
      ]
      
      # Extract human miRNAs mimicked by HIV-1
      HIV1_mimicked_hsa_miRNAs <- unique(mimic_dat_HIV1$miRNA)
      
      # Annotate differential expression result
      final_res$is_mimicked_by_HIV1 <- ifelse(
        final_res$miRNA %in% HIV1_mimicked_hsa_miRNAs,
        "Yes",
        "No"
      )
    }
    
    # =========================
    # 10. Reorder columns
    # =========================
    
    final_res <- final_res[, c(
      "miRNA",
      "mean_count_control",
      "mean_count_HIV1",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B",
      "DEG_status",
      "was_filtered",
      "is_mimicked_by_HIV1"
    )]
    
    # 按P值排序
    final_res <- final_res[order(final_res$P.Value), ]
    
    # =========================
    # 11. Export result
    # =========================
    
    write.table(
      final_res,
      file = out_file,
      sep = ",",
      quote = FALSE,
      row.names = FALSE
    )
    
    # =========================
    # 12. Summary
    # =========================
    
    cat("\n========== 分析完成 ==========\n")
    cat("总miRNA数：", nrow(final_res), "\n")
    cat("原本会被limma过滤的miRNA数：", sum(final_res$was_filtered), "\n")
    cat("Up in HIV-1 infected cells:", sum(final_res$DEG_status == "Up_in_HIV1"), "\n")
    cat("Down in HIV-1 infected cells:", sum(final_res$DEG_status == "Down_in_HIV1"), "\n")
    cat("Low expression filtered (but still reported):", sum(final_res$DEG_status == "Low_expression_filtered"), "\n")
    cat("Human miRNAs mimicked by HIV-1 in mimic file:", length(unique(mimic_dat_HIV1$miRNA)), "\n")
    cat("Tested miRNAs mimicked by HIV-1:", sum(final_res$is_mimicked_by_HIV1 == "Yes"), "\n")
    cat("Result saved to:", out_file, "\n")
    
    # 输出简化的结果表格
    simple_res <- final_res[, c("miRNA", "logFC", "adj.P.Val", "DEG_status", "was_filtered", "is_mimicked_by_HIV1")]
    cat("\n========== Top 30 miRNAs by P-value ==========\n")
    print(head(simple_res, 30))
    
    # 输出被过滤的miRNA示例
    filtered_mirnas <- final_res[final_res$was_filtered, "miRNA"]
    cat("\n========== 被过滤的miRNA示例（前20个） ==========\n")
    print(head(filtered_mirnas, 20))