######################## GSE289893 ########################

library(tidyverse)
library(data.table)
library(fgsea)
library(igraph)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(patchwork)


############################################################
## 1. File paths
############################################################

work_dir <- "D:/0work/0wholetransctiptome/2sRNAminic/17HIV1"

reactome_file <- file.path(work_dir, "Reactome.csv")
target_file <- file.path(work_dir, "intersect_symbols.csv")

target_deg_file <- file.path(
  work_dir,
  "GEO_HIV_bulk/GSE289893_HIV1_target_genes_all.csv"
)

all_deg_file <- file.path(
  work_dir,
  "GEO_HIV_bulk/GSE289893_HIV_vs_Control_24h_DESeq2_results.csv"
)

out_dir <- file.path(work_dir, "Reactome_405_target_analysis_R")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


############################################################
## 2. Helper functions
############################################################

clean_gene <- function(x) {
  x %>%
    as.character() %>%
    str_trim() %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("\\.$", "") %>%
    str_replace_all("^\"|\"$", "") %>%
    toupper()
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}


############################################################
## 3. Read input files
############################################################

reactome_raw <- fread(reactome_file, data.table = FALSE)
target_raw <- fread(target_file, data.table = FALSE)
target_deg_raw <- fread(target_deg_file, data.table = FALSE)
all_deg_raw <- fread(all_deg_file, data.table = FALSE)

cat("Reactome columns:\n")
print(colnames(reactome_raw))

cat("405 target gene columns:\n")
print(colnames(target_raw))

cat("Target DEG columns:\n")
print(colnames(target_deg_raw))

cat("All DEG columns:\n")
print(colnames(all_deg_raw))


############################################################
## 4. Prepare 405 target genes
############################################################

## intersect_symbols.csv 中基因列为 V2
target_genes <- target_raw %>%
  transmute(Gene = clean_gene(V2)) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct() %>%
  pull(Gene)

cat("Number of target genes:\n")
print(length(target_genes))

write.csv(
  data.frame(Gene = target_genes),
  file.path(out_dir, "00_405_target_genes_cleaned.csv"),
  row.names = FALSE
)


############################################################
## 5. Prepare Reactome pathway to gene table
############################################################

## Reactome.csv 中：
## 通路列：Pathway name
## 基因列：Submitted entities found

reactome_long <- reactome_raw %>%
  transmute(
    pathway_id = as.character(`Pathway identifier`),
    pathway = as.character(`Pathway name`),
    entities_found = as.character(`Submitted entities found`),
    entities_pvalue = safe_numeric(`Entities pValue`),
    entities_fdr = safe_numeric(`Entities FDR`),
    species = as.character(`Species name`)
  ) %>%
  filter(
    !is.na(pathway),
    pathway != "",
    !is.na(entities_found),
    entities_found != ""
  ) %>%
  separate_rows(entities_found, sep = ";") %>%
  transmute(
    pathway_id = pathway_id,
    pathway = pathway,
    Gene = clean_gene(entities_found),
    entities_pvalue = entities_pvalue,
    entities_fdr = entities_fdr,
    species = species
  ) %>%
  filter(!is.na(Gene), Gene != "") %>%
  distinct()

## Reactome gene sets
reactome_sets <- split(reactome_long$Gene, reactome_long$pathway)
reactome_sets <- lapply(reactome_sets, unique)

cat("Number of Reactome pathways:\n")
print(length(reactome_sets))

write.csv(
  reactome_long,
  file.path(out_dir, "00_Reactome_pathway_gene_long_table.csv"),
  row.names = FALSE
)


############################################################
## Part 1.
## 405 target genes Reactome pathway composition
############################################################

target_reactome <- reactome_long %>%
  filter(Gene %in% target_genes)

mapped_targets <- unique(target_reactome$Gene)
unmapped_targets <- setdiff(target_genes, mapped_targets)

cat("405 target genes mapped to Reactome:\n")
print(length(mapped_targets))

cat("405 target genes without Reactome annotation:\n")
print(length(unmapped_targets))

## 以全部 405 个 target genes 作为分母
total_target_genes <- length(unique(target_genes))

## 每个 immune / viral / virus Reactome 通路中包含多少 405 靶基因
pathway_composition <- target_reactome %>%
  filter(
    !is.na(pathway),
    str_detect(
      str_to_lower(pathway),
      "immu|viral|virus|SARS"
    )
  ) %>%
  group_by(pathway_id, pathway) %>%
  summarise(
    target_gene_count = n_distinct(Gene),
    target_genes = paste(sort(unique(Gene)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    total_target_genes = total_target_genes,
    percent_in_405 = target_gene_count / total_target_genes * 100
  ) %>%
  arrange(desc(target_gene_count), pathway)

pathway_composition

write.csv(
  pathway_composition,
  file.path(out_dir, "Part1_405_targets_Reactome_pathway_composition.csv"),
  row.names = FALSE
)

## 每个靶基因参与多少 Reactome 通路
gene_pathway_count <- target_reactome %>%
  group_by(Gene) %>%
  summarise(
    pathway_count = n_distinct(pathway),
    pathways = paste(sort(unique(pathway)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(pathway_count), Gene)

write.csv(
  gene_pathway_count,
  file.path(out_dir, "Part1_405_targets_gene_pathway_count.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Gene = unmapped_targets),
  file.path(out_dir, "Part1_405_targets_without_Reactome_annotation.csv"),
  row.names = FALSE
)

## Top 30 Reactome pathway composition plot
top_n <- 30

p1 <- pathway_composition %>%
  slice_max(order_by = target_gene_count, n = top_n) %>%
  mutate(pathway = fct_reorder(pathway, target_gene_count)) %>%
  ggplot(aes(x = pathway, y = percent_in_405)) +
  geom_col(
    width = 0.75,
    fill = "#4292c6"
  )  +
  geom_text(
    aes(label = paste0(target_gene_count, "")),
    hjust = -0.05,
    size = 3.2
  ) +
  coord_flip() +
  labs(
    title = "",
    x = "Reactome pathway",
    y = "Percentage of all target genes (%)"
  ) +
  theme_bw(base_size = 8) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  ) +
  expand_limits(
    y = max(pathway_composition$percent_in_405, na.rm = TRUE) * 1.15
  )

ggsave(
  file.path(out_dir, "Part1_405_targets_Reactome_pathway_composition_top30.pdf"),
  p1,
  width = 7,
  height = 5
)


############################################################
## Part 2.
## HIV infection induced up or down regulation of 405 targets
############################################################

## 关键修改：提前创建 all_deg
## 后续 Part 2 密度图和 Part 3 GSEA 都共用这个对象
all_deg <- all_deg_raw %>%
  mutate(
    Gene_clean = clean_gene(Gene),
    log2FoldChange = safe_numeric(log2FoldChange),
    stat = safe_numeric(stat),
    padj = safe_numeric(padj),
    is_HIV1_target_clean = Gene_clean %in% target_genes
  ) %>%
  filter(!is.na(Gene_clean), Gene_clean != "")

write.csv(
  all_deg,
  file.path(out_dir, "Part2_all_genes_DEG_cleaned.csv"),
  row.names = FALSE
)

target_deg <- target_deg_raw %>%
  mutate(
    Gene_clean = clean_gene(Gene),
    log2FoldChange = safe_numeric(log2FoldChange),
    FoldChange = 2 ^ log2FoldChange,
    FC_direction = case_when(
      is.na(log2FoldChange) ~ "NA",
      log2FoldChange > 0 ~ "Up",
      log2FoldChange < 0 ~ "Down",
      log2FoldChange == 0 ~ "No change"
    )
  ) %>%
  filter(Gene_clean %in% target_genes)

## 上下调比例
fc_summary <- target_deg %>%
  filter(FC_direction != "NA") %>%
  count(FC_direction, name = "gene_count") %>%
  mutate(
    total_gene_count = sum(gene_count),
    percent = gene_count / total_gene_count * 100
  ) %>%
  arrange(desc(gene_count))

write.csv(
  fc_summary,
  file.path(out_dir, "Part2_405_targets_FC_direction_summary.csv"),
  row.names = FALSE
)

write.csv(
  target_deg,
  file.path(out_dir, "Part2_405_targets_DEG_with_FC_direction.csv"),
  row.names = FALSE
)

## 上下调比例柱状图
p2_bar <- fc_summary %>%
  mutate(
    FC_direction = factor(
      FC_direction,
      levels = c("Up", "Down", "No change")
    )
  ) %>%
  ggplot(aes(x = FC_direction, y = percent)) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0(gene_count, "\n", sprintf("%.1f%%", percent))),
    vjust = -0.25,
    size = 4
  ) +
  labs(
    title = "Expression direction of 405 HIV1 target genes after HIV infection",
    x = "Direction based on log2FoldChange",
    y = "Percentage (%)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  expand_limits(y = max(fc_summary$percent, na.rm = TRUE) * 1.15)

ggsave(
  file.path(out_dir, "Part2_405_targets_up_down_percentage.pdf"),
  p2_bar,
  width = 6,
  height = 5
)



############################################################
## Part 2B.
## log2FoldChange density distribution
## All genes as grey background and 405 HIV1 target genes in light blue
############################################################
density_df <- bind_rows(
  all_deg %>%
    filter(
      !is.na(log2FoldChange),
      log2FoldChange != 0
    ) %>%
    transmute(
      log2FoldChange = log2FoldChange,
      Group = "All genes"
    ),
  target_deg %>%
    filter(
      !is.na(log2FoldChange),
      log2FoldChange != 0
    ) %>%
    transmute(
      log2FoldChange = log2FoldChange,
      Group = "HIV1 target genes"
    )
)

density_df <- density_df %>%
  mutate(
    Group = factor(
      Group,
      levels = c("All genes", "HIV1 target genes")
    )
  )

## 提取两个分布的 peak 位置
p_tmp <- ggplot(
  density_df,
  aes(x = log2FoldChange, fill = Group, color = Group)
) +
  geom_density(alpha = 0.45, linewidth = 1.1)

density_build <- ggplot_build(p_tmp)$data[[1]]

peak_df <- density_build %>%
  group_by(group) %>%
  slice_max(order_by = density, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    Group = levels(density_df$Group)[group]
  ) %>%
  select(Group, peak_x = x, peak_density = density)

write.csv(
  density_df,
  file.path(out_dir, "Part2_all_genes_vs_405_targets_log2FC_density_input.csv"),
  row.names = FALSE
)

write.csv(
  peak_df,
  file.path(out_dir, "Part2_all_genes_vs_405_targets_log2FC_density_peak.csv"),
  row.names = FALSE
)

print(peak_df)

p2_density <- ggplot() +
  geom_density(
    data = density_df %>% filter(Group == "All genes"),
    aes(x = log2FoldChange, fill = Group, color = Group),
    alpha = 0.45,
    linewidth = 1.1
  ) +
  geom_density(
    data = density_df %>% filter(Group == "HIV1 target genes"),
    aes(x = log2FoldChange, fill = Group, color = Group),
    alpha = 0.55,
    linewidth = 1.1
  ) +
  geom_vline(
    data = peak_df,
    aes(xintercept = peak_x, color = Group),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.7,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "All genes" = "grey80",
      "HIV1 target genes" = "#9ecae1"
    )
  ) +
  scale_color_manual(
    values = c(
      "All genes" = "grey50",
      "HIV1 target genes" = "#3182bd"
    )
  ) +
  labs(
    title = "Density distribution of log2FoldChange",
    subtitle = "All genes versus 405 HIV1 target genes",
    x = "log2FoldChange, HIV infection versus control",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(
  file.path(out_dir, "Part2_all_genes_vs_405_targets_log2FC_density.pdf"),
  p2_density,
  width = 7,
  height = 5
)



## log2FC 分布统计：405 targets
fc_distribution_summary <- target_deg %>%
  summarise(
    n = sum(!is.na(log2FoldChange)),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    median_log2FC = median(log2FoldChange, na.rm = TRUE),
    min_log2FC = min(log2FoldChange, na.rm = TRUE),
    max_log2FC = max(log2FoldChange, na.rm = TRUE),
    q25_log2FC = quantile(log2FoldChange, 0.25, na.rm = TRUE),
    q75_log2FC = quantile(log2FoldChange, 0.75, na.rm = TRUE)
  )

write.csv(
  fc_distribution_summary,
  file.path(out_dir, "Part2_405_targets_log2FC_distribution_summary.csv"),
  row.names = FALSE
)

## log2FC 分布统计：全基因和405 targets对比
fc_distribution_summary_compare <- density_df %>%
  group_by(Group) %>%
  summarise(
    n = sum(!is.na(log2FoldChange)),
    mean_log2FC = mean(log2FoldChange, na.rm = TRUE),
    median_log2FC = median(log2FoldChange, na.rm = TRUE),
    min_log2FC = min(log2FoldChange, na.rm = TRUE),
    max_log2FC = max(log2FoldChange, na.rm = TRUE),
    q25_log2FC = quantile(log2FoldChange, 0.25, na.rm = TRUE),
    q75_log2FC = quantile(log2FoldChange, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  fc_distribution_summary_compare,
  file.path(out_dir, "Part2_all_genes_vs_405_targets_log2FC_distribution_summary.csv"),
  row.names = FALSE
)


############################################################
## Part 3.
## Reactome GSEA based on whole transcriptome expression profile
############################################################

## GSEA 推荐使用 DESeq2 stat 排序
## 如果 stat 大量缺失，则使用 log2FoldChange
if (sum(!is.na(all_deg$stat)) >= 100) {
  rank_df_all <- all_deg %>%
    filter(!is.na(stat)) %>%
    group_by(Gene_clean) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(Gene = Gene_clean, rank_score = stat)
} else {
  rank_df_all <- all_deg %>%
    filter(!is.na(log2FoldChange)) %>%
    group_by(Gene_clean) %>%
    slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(Gene = Gene_clean, rank_score = log2FoldChange)
}

rank_all <- rank_df_all$rank_score
names(rank_all) <- rank_df_all$Gene
rank_all <- sort(rank_all, decreasing = TRUE)

## Reactome gene sets 需要和全基因表达谱取交集
reactome_sets_all <- lapply(
  reactome_sets,
  function(g) intersect(g, names(rank_all))
)

reactome_sets_all <- reactome_sets_all[
  lengths(reactome_sets_all) >= 5 &
    lengths(reactome_sets_all) <= 500
]

set.seed(123)

fgsea_all <- fgsea(
  pathways = reactome_sets_all,
  stats = rank_all,
  minSize = 5,
  maxSize = 500,
  eps = 0
)

fgsea_all_res <- fgsea_all %>%
  as.data.frame() %>%
  arrange(padj, desc(abs(NES))) %>%
  mutate(
    leadingEdge = sapply(leadingEdge, paste, collapse = ";"),
    target_gene_count_in_pathway = sapply(pathway, function(pw) {
      length(intersect(reactome_sets[[pw]], target_genes))
    }),
    target_genes_in_pathway = sapply(pathway, function(pw) {
      paste(sort(intersect(reactome_sets[[pw]], target_genes)), collapse = ";")
    }),
    target_ratio_in_pathway = target_gene_count_in_pathway / size
  )

write.csv(
  fgsea_all_res,
  file.path(out_dir, "Part3_all_genes_Reactome_GSEA_results.csv"),
  row.names = FALSE
)


############################################################
## Function:
## Build Reactome pathway network
##
## Node = Reactome pathway
## Node color = NES from GSEA
## Node size = number of 405 target genes in this pathway
## Edge = shared genes between two pathways
## Edge width = Jaccard similarity
############################################################

build_target_ratio_pathway_network <- function(
    gsea_res,
    reactome_sets,
    target_genes,
    pval_cutoff = 0.05,
    padj_cutoff = NULL,
    top_n = 30,
    jaccard_cutoff = 0.01,
    min_overlap = 1,
    min_target_count = 1,
    out_prefix
) {
  
  ##########################################################
  ## 0. Check required columns
  ##########################################################
  
  required_cols <- c("pathway", "NES", "pval", "padj", "size")
  
  missing_cols <- setdiff(required_cols, colnames(gsea_res))
  
  if (length(missing_cols) > 0) {
    stop(
      paste0(
        "Missing required columns in gsea_res: ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  
  ##########################################################
  ## 1. Build node table
  ##########################################################
  
  nodes <- gsea_res %>%
    dplyr::filter(!is.na(pval)) %>%
    dplyr::mutate(
      pathway = as.character(pathway),
      
      ## 该通路中有多少 405 靶基因
      target_gene_count_in_pathway = sapply(pathway, function(pw) {
        length(intersect(reactome_sets[[pw]], target_genes))
      }),
      
      ## 该通路中的 405 靶基因列表
      target_genes_in_pathway = sapply(pathway, function(pw) {
        paste(sort(intersect(reactome_sets[[pw]], target_genes)), collapse = ";")
      }),
      
      ## 通路总基因数
      pathway_gene_count = sapply(pathway, function(pw) {
        length(unique(reactome_sets[[pw]]))
      }),
      
      ## 405 靶基因占该通路的比例
      target_ratio_in_pathway = target_gene_count_in_pathway / pathway_gene_count,
      
      ## 百分比形式，便于后续结果解释
      target_percent_in_pathway = target_ratio_in_pathway * 100,
      
      ## 显著性
      neg_log10_pval = -log10(pval),
      
      ## 标签换行，避免通路名称过长
      pathway_label = stringr::str_wrap(pathway, width = 30)
    )
  
  ##########################################################
  ## 2. Filter enriched pathways
  ##########################################################
  
  ## 优先使用 padj_cutoff
  ## 如果 padj_cutoff = NULL，则使用 nominal pvalue
  if (!is.null(padj_cutoff)) {
    nodes <- nodes %>%
      dplyr::filter(!is.na(padj), padj < padj_cutoff)
  } else {
    nodes <- nodes %>%
      dplyr::filter(pval < pval_cutoff)
  }
  
  ## 只展示包含 405 靶基因的通路
  nodes <- nodes %>%
    dplyr::filter(
      target_gene_count_in_pathway >= min_target_count,
      target_ratio_in_pathway > 0
    ) %>%
    dplyr::arrange(pval, dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::select(
      pathway,
      pathway_label,
      NES,
      pval,
      padj,
      neg_log10_pval,
      size,
      pathway_gene_count,
      target_gene_count_in_pathway,
      target_ratio_in_pathway,
      target_percent_in_pathway,
      target_genes_in_pathway
    )
  
  write.csv(
    nodes,
    paste0(out_prefix, "_nodes.csv"),
    row.names = FALSE
  )
  
  if (nrow(nodes) < 2) {
    warning("Fewer than two pathways contain 405 target genes. Network was not generated.")
    return(list(nodes = nodes, edges = NULL, plot = NULL))
  }
  
  ##########################################################
  ## 3. Build edge table
  ##########################################################
  
  pathway_vec <- nodes$pathway
  
  edges <- combn(pathway_vec, 2, simplify = FALSE) %>%
    purrr::map_dfr(function(x) {
      
      g1 <- unique(reactome_sets[[x[1]]])
      g2 <- unique(reactome_sets[[x[2]]])
      
      overlap_genes <- intersect(g1, g2)
      union_genes <- union(g1, g2)
      
      tibble::tibble(
        from = x[1],
        to = x[2],
        overlap_gene_count = length(overlap_genes),
        overlap_genes = paste(sort(overlap_genes), collapse = ";"),
        jaccard = length(overlap_genes) / length(union_genes)
      )
    }) %>%
    dplyr::filter(
      overlap_gene_count >= min_overlap,
      jaccard >= jaccard_cutoff
    )
  
  write.csv(
    edges,
    paste0(out_prefix, "_edges.csv"),
    row.names = FALSE
  )
  
  if (nrow(edges) == 0) {
    warning("No edges passed the cutoff. Try lowering jaccard_cutoff or min_overlap.")
    return(list(nodes = nodes, edges = edges, plot = NULL))
  }
  
  ##########################################################
  ## 4. Draw network
  ##########################################################
  
  graph_obj <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = FALSE
  )
  
  n_nodes <- nrow(nodes)
  
  label_size <- dplyr::case_when(
    n_nodes <= 30 ~ 3.2,
    n_nodes <= 60 ~ 2.8,
    TRUE ~ 2.4
  )
  
  point_size_max <- dplyr::case_when(
    n_nodes <= 30 ~ 10,
    n_nodes <= 60 ~ 8,
    TRUE ~ 6.5
  )
  
  p_net <- ggraph::ggraph(graph_obj, layout = "fr") +
    
    ## 边：两个通路之间的共享基因
    ggraph::geom_edge_link(
      ggplot2::aes(width = jaccard),
      alpha = 0.20,
      color = "grey65"
    ) +
    
    ## 点：Reactome 通路
    ## 点颜色：NES
    ## 点大小：该通路中 405 靶基因数量
    ggraph::geom_node_point(
      ggplot2::aes(
        size = target_gene_count_in_pathway,
        color = NES
      ),
      alpha = 0.92
    ) +
    
    ## 通路名称标签
    ggraph::geom_node_text(
      ggplot2::aes(label = pathway_label),
      repel = TRUE,
      size = label_size,
      max.overlaps = Inf,
      lineheight = 0.9
    ) +
    
    ## 线宽表示 Jaccard
    ggraph::scale_edge_width(
      range = c(0.2, 2.2)
    ) +
    
    ## 点大小表示 405 靶基因数量
    ggplot2::scale_size_continuous(
      range = c(2.5, point_size_max)
    ) +
    
    ## 点颜色表示 NES
    ## 红色：NES > 0，HIV 感染后排序上端富集
    ## 蓝色：NES < 0，HIV 感染后排序下端富集
    ggplot2::scale_color_gradient2(
      low = "#2166AC",
      mid = "grey90",
      high = "#B2182B",
      midpoint = 0
    ) +
    
    ggplot2::labs(
      title = "Reactome pathway network after HIV infection",
      subtitle = paste0(
        "Nodes are GSEA-enriched Reactome pathways with at least one 405 target gene; ",
        "edges indicate shared genes between pathways"
      ),
      color = "NES",
      size = "405 target\ngene count",
      edge_width = "Jaccard\nsimilarity"
    ) +
    
    ggplot2::coord_cartesian(clip = "off") +
    
    ggplot2::theme_void(base_size = 12) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0),
      legend.position = "right",
      plot.margin = ggplot2::margin(20, 50, 20, 50)
    )
  
  return(list(
    nodes = nodes,
    edges = edges,
    plot = p_net
  ))
}




############################################################
## Part 3 network
## Whole transcriptome Reactome GSEA network
## Node color = proportion of 405 target genes in each pathway
############################################################
library(ggplot2)

network_all_ratio <- build_target_ratio_pathway_network(
  gsea_res = fgsea_all_res,
  reactome_sets = reactome_sets,
  target_genes = target_genes,
  
  ## 如果 GSEA 结果 FDR 不显著，可以先用 nominal pvalue 展示趋势网络
  pval_cutoff = 1,
  padj_cutoff = NULL,
  
  ## 最多展示多少个通路
  top_n = 30,
  
  ## 边的筛选条件
  jaccard_cutoff = 0.01,
  min_overlap = 1,
  
  ## 只展示至少含有 1 个 405 靶基因的通路
  min_target_count = 1,
  
  out_prefix = file.path(
    out_dir,
    "Part3_all_genes_GSEA_Reactome_network_target_ratio"
  )
)
network_all_ratio
if (!is.null(network_all_ratio$plot)) {
  
  ggsave(
    file.path(out_dir, "Part3_all_genes_GSEA_Reactome_network_target_ratio.pdf"),
    network_all_ratio$plot,
    width = 8,
    height = 8.5
  )
  
}




















