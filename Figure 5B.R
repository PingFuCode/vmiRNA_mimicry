
setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\9cancer")
# 加载必要的库
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

# 读取miRNA靶基因数据
target_file <- "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv"
target_data <- read_csv(target_file)
target_data <- read_csv(target_file) %>%
  filter(abbreviation != "Herpesvirus")
#View(target_data)
# 读取癌基因与抑癌基因数据
cancer_gene_file <- "D:/0work/0wholetransctiptome/2sRNAminic/9cancer/CancerGene.xlsx"
tumor_suppressor_genes <- read_excel(cancer_gene_file, sheet = "TSGene") %>% distinct(GeneSymbol) %>% pull(GeneSymbol)
oncogenes <- read_excel(cancer_gene_file, sheet = "Oncogene") %>% distinct(GeneSymbol) %>% pull(GeneSymbol)

# 计算每个 miRNA 共有靶基因的抑癌和促癌基因比例
calculate_gene_proportion <- function(genes) {
  gene_list <- unlist(strsplit(genes, ";"))  # 拆分基因列表
  total_genes <- length(gene_list)  # 总基因数
  
  if (total_genes == 0) {
    return(c(0, 0))  # 避免除零错误
  }
  
  tumor_suppressor_count <- sum(gene_list %in% tumor_suppressor_genes)  # 抑癌基因数
  oncogene_count <- sum(gene_list %in% oncogenes)  # 促癌基因数
  
  tumor_suppressor_ratio <- tumor_suppressor_count / total_genes
  oncogene_ratio <- oncogene_count / total_genes
  
  return(c(tumor_suppressor_ratio, oncogene_ratio))
}

# 计算每一行的比例
target_data <- target_data %>%
  rowwise() %>%
  mutate(
    Tumor_Suppressor_Ratio = round(calculate_gene_proportion(Intersection_Genes)[1], 2),
    Oncogene_Ratio = round(calculate_gene_proportion(Intersection_Genes)[2], 2)
  ) %>%
  ungroup()


out_file<- "0006mer_target_with_hsaID_GENE.csv"
write_csv(target_data,out_file)

# 转换数据格式：从宽格式转换为长格式
target_long <- target_data %>%
  select(abbreviation, Tumor_Suppressor_Ratio, Oncogene_Ratio) %>%
  pivot_longer(cols = c(Tumor_Suppressor_Ratio, Oncogene_Ratio),
               names_to = "Gene_Type",
               values_to = "Ratio")

# 修改分类变量名称
target_long$Gene_Type <- factor(target_long$Gene_Type, levels = c("Tumor_Suppressor_Ratio", "Oncogene_Ratio"),
                                labels = c("Tumor Suppressor Gene", "Oncogene"))
#View(target_long)
# 绘制箱型图
ggplot(target_long, aes(x = abbreviation, y = Ratio, fill = Gene_Type)) +
  geom_boxplot() + # outlier.shape = NA
  facet_wrap(~ abbreviation, nrow = 4, ncol = 6, scales = "free_x") +  # 4行6列布局
  labs(title = "", y = "Gene Ratio", x = "Virus") +
  scale_fill_manual(values = c("#00BFC4", "#F8766D")) +  # 设定颜色("Mimic" = "#E41A1C", "UnMimic" = "#377EB8")
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # 移除 x 轴文本
    axis.ticks.x = element_blank(),  # 移除 x 轴刻度
    axis.text.y = element_text(size = 8),  # 调整 y 轴字体大小，避免重叠
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    panel.border = element_rect(colour = "black", fill = NA),  # 增加边框
    plot.background = element_blank(),  # 移除背景
    panel.background = element_blank()   # 移除面板背景
  )




##########################################考虑背景的影响##########################################
library(readr)
library(dplyr)
library(ggplot2)
library(readxl)
# 1. 读取目标基因数据
target_file <- "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv"
target_data <- read_csv(target_file) %>%
  filter(abbreviation %in% c("HIV-1", "HHV-4", "HHV-5", "HHV-8", "HSV-1"))  # 过滤掉没有显著促进疾病的病毒
# 提取病毒模拟的 miRNA（共 600 个）
mirna_list <- unique(target_data$miRNA)
length(mirna_list)
# 2. 读取所有 miRNA 的靶基因数据库
database_file <- "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/database.csv"
database_data <- read_csv(database_file)
# 3. 过滤病毒模拟的 miRNA，并提取其靶基因
mirna_targets <- database_data %>%
  filter(miRNA %in% mirna_list) %>%
  group_by(miRNA) %>%
  summarise(TargetGenes = paste(unique(TargetGene), collapse = ";"), .groups = "drop")

# 4. 替换 `Intersection_Genes` 列
target_data <- target_data %>%
  left_join(mirna_targets, by = "miRNA") %>%  # 关联 miRNA 靶基因数据
  mutate(Intersection_Genes = ifelse(is.na(TargetGenes), Intersection_Genes, TargetGenes)) %>%
  select(-TargetGenes)  # 移除额外的列


cancer_gene_file <- "D:/0work/0wholetransctiptome/2sRNAminic/9cancer/CancerGene.xlsx"
tumor_suppressor_genes <- read_excel(cancer_gene_file, sheet = "TSGene") %>%
  distinct(GeneSymbol) %>%
  pull(GeneSymbol)

oncogenes <- read_excel(cancer_gene_file, sheet = "Oncogene") %>%
  distinct(GeneSymbol) %>%
  pull(GeneSymbol)

# 计算背景比例
total_genes <- 21306  # 总基因数
total_tumor_suppressor_genes <- length(tumor_suppressor_genes)
total_oncogenes <- length(oncogenes)

background_tumor_suppressor_ratio <- total_tumor_suppressor_genes / total_genes
background_oncogene_ratio <- total_oncogenes / total_genes

# 定义计算基因比例的函数
calculate_gene_proportion <- function(genes) {
  gene_list <- unlist(strsplit(genes, ";"))  # 拆分基因列表
  total_genes <- length(gene_list)  # 总基因数
  
  if (total_genes == 0) {
    return(c(0, 0))  # 避免除零错误
  }
  
  tumor_suppressor_count <- sum(gene_list %in% tumor_suppressor_genes)  # 抑癌基因数
  oncogene_count <- sum(gene_list %in% oncogenes)  # 促癌基因数
  
  tumor_suppressor_ratio <- tumor_suppressor_count / total_genes
  oncogene_ratio <- oncogene_count / total_genes
  
  return(c(tumor_suppressor_ratio, oncogene_ratio))
}

target_data <- target_data %>%
  rowwise() %>%
  mutate(
    Tumor_Suppressor_Ratio = round(calculate_gene_proportion(Intersection_Genes)[1], 2),
    Oncogene_Ratio = round(calculate_gene_proportion(Intersection_Genes)[2], 2),
    Tumor_Suppressor_Ratio_vs_Background = round(Tumor_Suppressor_Ratio / background_tumor_suppressor_ratio, 2),
    Oncogene_Ratio_vs_Background = round(Oncogene_Ratio / background_oncogene_ratio, 2),
    Significant_Tumor_Suppressor = Tumor_Suppressor_Ratio_vs_Background > 1,
    Significant_Oncogene = Oncogene_Ratio_vs_Background > 1
  ) %>%
  ungroup()





# 5. 保存新的 target_data
out_file <- "0006mer_target_with_hsaID_GENE.csv"
write_csv(target_data, out_file)





target_data <- read.csv("0006mer_target_with_hsaID_GENE.csv", stringsAsFactors = FALSE)
filtered_data <- target_data[target_data$Tumor_Suppressor_Ratio_vs_Background > 0, ]

min(filtered_data$Tumor_Suppressor_Ratio_vs_Background)
max(filtered_data$Tumor_Suppressor_Ratio_vs_Background)
median(filtered_data$Tumor_Suppressor_Ratio_vs_Background)

count <- sum(filtered_data$Tumor_Suppressor_Ratio_vs_Background > 1)
total <- nrow(filtered_data)
percentage <- count / total
cat("占比为：", paste(round(percentage * 100, 2), "%", sep = ""))





####只画抑癌基因
ggplot(filtered_data, aes(x = Tumor_Suppressor_Ratio_vs_Background, fill = "Tumor Suppressor Ratio vs Background")) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 1) +
  labs(
    title = "",
    x = "Tumor Suppressor Ratio vs Backgroun",
    y = "Density",
    fill = "Legend"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Tumor Suppressor Ratio vs Background" = "#DB3124")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 12, color = "black", face = "bold"),
    axis.title.x = element_text(vjust = 2, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(vjust = 2, size = 12, face = "bold"),
    axis.line = element_line(color = 'black')
  )


####################################
# 筛选 Tumor_Suppressor_Ratio_vs_Background 大于 0 的数据
filtered_data <- target_data[target_data$Tumor_Suppressor_Ratio_vs_Background > 0, ]
head(filtered_data)
# 绘制箱型图
ggplot(filtered_data, aes(x="",y = Tumor_Suppressor_Ratio_vs_Background, fill = "Tumor Suppressor Ratio vs Background")) +
  geom_boxplot() +
  labs(
    title = "Tumor Suppressor Ratio vs Background",
    x = "",
    y = "Ratio",
    fill = "Legend"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Tumor Suppressor Ratio vs Background" = "#DB3124")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    axis.title.x = element_text(vjust = 2, size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    axis.title.y = element_text(vjust = 2, size = 12, face = "bold"),
    axis.line = element_line(color = 'black')
  )



filtered_data <- target_data[
  target_data$Oncogene_Ratio_vs_Background > 0 & 
    target_data$Tumor_Suppressor_Ratio_vs_Background > 0, 
]

# 8. 可视化分析--癌症和抑癌基因
ggplot(filtered_data, aes(x = Oncogene_Ratio_vs_Background, fill = "Oncogene Ratio vs Background")) +
  geom_density(alpha = 0.5) +
  geom_density(aes(x = Tumor_Suppressor_Ratio_vs_Background, fill = "Tumor Suppressor Ratio vs Background"), alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", size = 1)+
  labs(
    title = "",
    x = "",
    y = "Density",
    fill = "Legend"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Oncogene Ratio vs Background" = "#DB3124", 
                               "Tumor Suppressor Ratio vs Background" = "#4B74B2"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold"),axis.title.x=element_text(vjust=2, size=12, face = "bold"))+
  theme(axis.text.y = element_text(size = 12, color = "black", face = "bold"),axis.title.y=element_text(vjust=2, size=12, face = "bold"))+
  theme(axis.line = element_line(color = 'black'))


dim(target_data)

target_data <- target_data %>%
  filter(Tumor_Suppressor_Ratio_vs_Background > 1 | Oncogene_Ratio_vs_Background > 1)


target_data_long <- target_data%>%
  pivot_longer(
    cols = c(Oncogene_Ratio_vs_Background, Tumor_Suppressor_Ratio_vs_Background),
    names_to = "Ratio_Type",
    values_to = "Ratio_Value"
  ) %>%
  drop_na(Ratio_Value)

#View(target_data_long)

p_values <- target_data_long %>%
  group_by(abbreviation) %>%
  summarize(
    p_value = if (length(unique(Ratio_Type)) == 2) {
      wilcox.test(Ratio_Value ~ Ratio_Type, data = cur_data())$p.value
    } else {
      0  # 如果只有一组数据，则 p 值设置为 0
    }
  ) %>%
  mutate(p_value_text = ifelse(p_value == 0, "p = N/A", paste("p =", round(p_value, 3))))

print(p_values,n=21)

colnames(target_data_long)
# 绘制箱型图并按abbreviation分面
ggplot(target_data_long, aes(x = Ratio_Type, y = Ratio_Value, fill = Ratio_Type)) +
  geom_boxplot() +
  facet_wrap(~ abbreviation, scales = "free_y", ncol = 5) +  # 按abbreviation分面，自由y轴范围，每行4列
  labs(
    title = "",
    x = "Virus",
    y = "Ratio vs Background",
  ) +
  scale_fill_manual(values = c("#DB3124", "#4B74B2")) +  # 设定颜色
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8),  # 调整 y 轴字体大小，避免重叠
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    panel.border = element_rect(colour = "black", fill = NA),  # 增加边框
    plot.background = element_blank(),  # 移除背景
    panel.background = element_blank()   # 移除面板背景
  )

#("#DB3124", "#FC8C5A", "#4B74B2","#8DD3C7","#FFFFB3")















#####################################交集基因#######################
target_file <- "D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv"
target_data <- read_csv(target_file) %>%
  filter(abbreviation != "Herpesvirus")






# 读取癌基因与抑癌基因数据
cancer_gene_file <- "D:/0work/0wholetransctiptome/2sRNAminic/9cancer/CancerGene.xlsx"
tumor_suppressor_genes <- read_excel(cancer_gene_file, sheet = "TSGene") %>%
  distinct(GeneSymbol) %>%
  pull(GeneSymbol)

oncogenes <- read_excel(cancer_gene_file, sheet = "Oncogene") %>%
  distinct(GeneSymbol) %>%
  pull(GeneSymbol)

# 计算背景比例
total_genes <- 21306  # 总基因数
total_tumor_suppressor_genes <- length(tumor_suppressor_genes)
total_oncogenes <- length(oncogenes)

background_tumor_suppressor_ratio <- total_tumor_suppressor_genes / total_genes
background_oncogene_ratio <- total_oncogenes / total_genes

# 定义计算基因比例的函数
calculate_gene_proportion <- function(genes) {
  gene_list <- unlist(strsplit(genes, ";"))  # 拆分基因列表
  total_genes <- length(gene_list)  # 总基因数
  
  if (total_genes == 0) {
    return(c(0, 0))  # 避免除零错误
  }
  
  tumor_suppressor_count <- sum(gene_list %in% tumor_suppressor_genes)  # 抑癌基因数
  oncogene_count <- sum(gene_list %in% oncogenes)  # 促癌基因数
  
  tumor_suppressor_ratio <- tumor_suppressor_count / total_genes
  oncogene_ratio <- oncogene_count / total_genes
  
  return(c(tumor_suppressor_ratio, oncogene_ratio))
}

# 计算每一行的比例并添加到数据框中
target_data <- target_data %>%
  rowwise() %>%
  mutate(
    Tumor_Suppressor_Ratio = round(calculate_gene_proportion(Intersection_Genes)[1], 2),
    Oncogene_Ratio = round(calculate_gene_proportion(Intersection_Genes)[2], 2),
    Tumor_Suppressor_Ratio_vs_Background = round(Tumor_Suppressor_Ratio / background_tumor_suppressor_ratio, 2),
    Oncogene_Ratio_vs_Background = round(Oncogene_Ratio / background_oncogene_ratio, 2),
    Significant_Tumor_Suppressor = Tumor_Suppressor_Ratio_vs_Background > 1,
    Significant_Oncogene = Oncogene_Ratio_vs_Background > 1
  ) %>%
  ungroup()


out_file<- "0006mer_target_with_hsaID_GENE.csv"
write_csv(target_data,out_file)


ggplot(target_data, aes(x = Oncogene_Ratio_vs_Background, fill = "Oncogene Ratio vs Background")) +
  geom_density(alpha = 0.5) +
  geom_density(aes(x = Tumor_Suppressor_Ratio_vs_Background, fill = "Tumor Suppressor Ratio vs Background"), alpha = 0.5) +
  labs(
    title = "",
    x = "",
    y = "Density",
    fill = "Legend"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("Oncogene Ratio vs Background" = "#DB3124", 
                               "Tumor Suppressor Ratio vs Background" = "#4B74B2"))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())+
  theme(axis.text.x = element_text(size = 12, color = "black", face = "bold"),axis.title.x=element_text(vjust=2, size=12, face = "bold"))+
  theme(axis.text.y = element_text(size = 12, color = "black", face = "bold"),axis.title.y=element_text(vjust=2, size=12, face = "bold"))+
  theme(axis.line = element_line(color = 'black'))


dim(target_data)

target_data <- target_data %>%
  filter(Tumor_Suppressor_Ratio_vs_Background > 1 | Oncogene_Ratio_vs_Background > 1)


target_data_long <- target_data%>%
  pivot_longer(
    cols = c(Oncogene_Ratio_vs_Background, Tumor_Suppressor_Ratio_vs_Background),
    names_to = "Ratio_Type",
    values_to = "Ratio_Value"
  ) %>%
  drop_na(Ratio_Value)

#View(target_data_long)

p_values <- target_data_long %>%
  group_by(abbreviation) %>%
  summarize(
    p_value = if (length(unique(Ratio_Type)) == 2) {
      wilcox.test(Ratio_Value ~ Ratio_Type, data = cur_data())$p.value
    } else {
      0  # 如果只有一组数据，则 p 值设置为 0
    }
  ) %>%
  mutate(p_value_text = ifelse(p_value == 0, "p = N/A", paste("p =", round(p_value, 3))))

print(p_values,n=21)

colnames(target_data_long)
# 绘制箱型图并按abbreviation分面
ggplot(target_data_long, aes(x = Ratio_Type, y = Ratio_Value, fill = Ratio_Type)) +
  geom_boxplot() +
  facet_wrap(~ abbreviation, scales = "free_y", ncol = 5) +  # 按abbreviation分面，自由y轴范围，每行4列
  labs(
    title = "",
    x = "Virus",
    y = "Ratio vs Background",
  ) +
  scale_fill_manual(values = c("#DB3124", "#4B74B2")) +  # 设定颜色
  theme_minimal() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8),  # 调整 y 轴字体大小，避免重叠
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    panel.border = element_rect(colour = "black", fill = NA),  # 增加边框
    plot.background = element_blank(),  # 移除背景
    panel.background = element_blank()   # 移除面板背景
  )