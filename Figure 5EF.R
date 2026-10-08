
# 加载必要的库
library(GEOquery)
library(limma)
packageVersion("limma")
library(affy)
library(ggplot2)

##############################################sRNA做差异表达######################################
# 加载必要的R包
library(limma)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggsci)
# 设置工作目录
setwd("D:/0work/0wholetransctiptome/2sRNAminic/9cancer/EBV_Disease/GC/Expression/human")

# 1. 合并所有GSM文件
#file_list <- list.files(pattern = "GSM.*\\.gz", full.names = TRUE)
#combined_data <- data.frame(miRNA_ID = character())
#for (i in seq_along(file_list)) {
#  df <- read.delim(file_list[i], sep = "\t", stringsAsFactors = FALSE)
#  df <- df %>% 
#    select(miRNA_ID, Normalized_Value) %>%
#    rename(!!sample_names[i] := Normalized_Value)
#    combined_data <- full_join(combined_data, df, by = "miRNA_ID")
#}

#View(combined_data)
#write.csv(combined_data,"sRNA.csv")

# 保留Normalized_Value列

#final_data<-read.csv("sRNA.csv",row.names = 1,  header = TRUE)
final_data<-read.csv("sRNA1.csv",row.names = 1,  header = TRUE)
head(final_data)
final_data <- log2(final_data + 1)

# 2. 添加分组信息
# 分组信息
#group_info <- data.frame(
#  Library_Name = c("GSM7090061", "GSM7090062", "GSM7090063", "GSM7090064", 
#                   "GSM7090065", "GSM7090066", "GSM7090067", "GSM7090068", 
#                   "GSM7090069", "GSM7090070", "GSM7090071", "GSM7090072"),
#  Group = c(rep("EBV_negative", 6), rep("EBV_positive", 6)))

group_info <- data.frame(
  Library_Name = c("GSM7090061", "GSM7090062", "GSM7090063", "GSM7090067", "GSM7090068", "GSM7090071"),
  Group = c(rep("EBV_negative", 3), rep("EBV_positive", 3)))

sample_order <- colnames(final_data)[-1]
group_info <- group_info[match(sample_order, group_info$Library_Name), ]

final_data <- as.matrix(final_data[, -1])

groups <- factor(group_info$Group)
design <- model.matrix(~0 + groups)

colnames(design) <- levels(groups)
contrast_matrix <- makeContrasts(GroupBvsA = EBV_positive - EBV_negative, levels = design)

# 使用limma进行分析
fit <- lmFit(final_data, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

results <- topTable(fit2, coef = "GroupBvsA", number = Inf, adjust.method = "BH")
results$GeneSymbol <- rownames(results)  # 添加基因名列
#View(results)
#write.csv(results,"limma.csv",quote = TRUE)




results$change <- factor(
  ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
         ifelse(results$logFC > 1, "Up", "Down"), "Not Sig"),
  levels = c("Down", "Not Sig", "Up")
)

ggplot(results, aes(logFC, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999",linewidth=1) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "#999999",linewidth=1) +
  geom_point(aes(color = change), size = 1, alpha = 0.6) +
  theme_bw(base_size = 12) +
  ggsci::scale_color_npg() +  # npg期刊配色
  theme(panel.grid = element_blank(), legend.position = 'right') +
  geom_text_repel(
    data = subset(results, abs(logFC) > 1 & -log10(adj.P.Val) > 4),
    aes(label = GeneSymbol, color = change),
    size = 2,
    max.overlaps = 20
  ) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 (FDR Adjusted P-value)")+scale_x_continuous(limits = c(-10, 10))


#################################mimic miRNAs 查看比例和抑癌基因的比例#######################################


simulated_miRNA <- read.csv("D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv", stringsAsFactors = FALSE)
simulated_miRNA <- simulated_miRNA %>% filter(abbreviation == "HHV-4")

simulated_miRNA <- unique(simulated_miRNA$miRNA)


logFC_threshold <- 1  # 对数变化阈值
adjP_threshold <- 0.05  # 校正后的P值阈值

#head(results)
results1 <- results %>%
  mutate(miRNA_ID = rownames(.)) %>%  # 将行名添加为miRNA_ID列
  mutate(Status = case_when(
    logFC > logFC_threshold & adj.P.Val < adjP_threshold ~ "Up",
    logFC < -logFC_threshold & adj.P.Val < adjP_threshold ~ "Down",
    TRUE ~ "Not Sig"
  ))

simulated_results <- results1 %>%
  filter(miRNA_ID %in% simulated_miRNA)

simulated_summary <- simulated_results %>%
  count(Status) %>%
  rename(Simulated_n = n)

# 计算所有miRNA在每个状态的数量
all_summary <- results1 %>%
  count(Status) %>%
  rename(All_n = n)

summary <- simulated_summary %>%
  left_join(all_summary, by = "Status") %>%
  mutate(Proportion = Simulated_n / All_n * 100)  # 计算百分比

# 打印结果
print(summary)




#write.csv(simulated_results, file = "Mimic.csv",quote = FALSE)



highlight_miRNAs <-simulated_miRNA
geneList <- results[highlight_miRNAs,]


results$change <- factor(
  ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
         ifelse(results$logFC > 1, "Up", "Down"), "Not Sig"),
  levels = c("Down", "Not Sig", "Up")
)
head(results)

#library(RColorBrewer)
#brewer.pal(6, "Blues")

p<-ggplot(results, aes(logFC, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999",linewidth=1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999",linewidth=1) +
  geom_point(aes(color = change), size = 2.5, alpha = 0.5) +
  geom_point(data=geneList,aes(x = logFC, y = -log10(adj.P.Val)),colour="yellow",size=2)+
  theme_bw(base_size = 12) +
  #ggsci::scale_color_npg() +  # npg期刊配色
  scale_color_manual(values = c("Down" = "#6BAED6", "Up" = "#FB6A4A", "Not Sig" = "grey")) +
  theme(panel.grid = element_blank(), legend.position = 'right') +
  xlab("Log2 Fold Change") +
  ylab("-Log10 (FDR Adjusted P-value)")+scale_x_continuous(limits = c(-10, 10))
p



index <- which(simulated_results$Status %in% c("Down", "Up"))
index <- which(simulated_results$Status %in% c( "Down"))

highlight_DEmiRNAs <- row.names(simulated_results)[index]
highlight_DEmiRNAs
#highlight_DEmiRNAs <-c("hsa-miR-378c","hsa-miR-378d","hsa-miR-654-3p","hsa-miR-597-5p","hsa-miR-4525","hsa-miR-127-5p")
highlight_DEmiRNAs <-c("hsa-miR-429")


geneList1 <- results[rownames(results) %in% highlight_DEmiRNAs,]
geneList1 <- subset(geneList1, select = -change)
geneList1$label <- rownames(geneList1)
dim(geneList1)


geneList1$change <- ifelse(geneList1$logFC > 0, "Up", "Down")
psRNA<-p + geom_label_repel(data = geneList1, 
                     aes(x = logFC, y = -log10(adj.P.Val), label = label, color = change),  # 添加 color = change
                     size = 4,
                     box.padding = unit(0.4, "lines"), 
                     segment.color = "black",   # 连线的颜色
                     segment.size = 0.4, 
                     max.overlaps = 25) +
  scale_color_manual(values = c("Up" = "#DE2D26", "Down" = "#3182BD"))


psRNA

#hsa-miR-429靶向HOXD10（显著）
#hsa-miR-6864-3p、hsa-miR-6894-5p 靶向TCEA3（不显著）

############################################################

#gset <- getGEO('GSE51575', 
#               destdir = ".",
#               AnnotGPL = TRUE,   # 获取注释信息
#               getGPL = TRUE)  
#eset <- gset[[1]]
#dim(exprs(eset))
#exp<-exprs(gset[[1]])
#cli<-pData(gset[[1]])	## 获取临床信息
##group<-c(rep("control",3),rep("hht",3))	## 查看分组信息
#GPL<-fData(gset[[1]])
#gpl<-GPL[,c(1,6)]
#exp<-as.data.frame(exp)
#exp$ID<-rownames(exp)	# 增加新的一列（最后一列），存放基因ID信息
#exp_symbol<-merge(exp,gpl,by="ID")
#write.csv(exp_symbol,"GSE51575.csv")




################################################
setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\9cancer\\EBV_Disease")


final_data<-read.csv("GSE51575.csv", header = TRUE)
head(final_data)
library(dplyr)
final_data_median <- final_data %>%
  group_by(GeneName) %>%
  summarise(across(everything(), median))

# 将 GeneName 作为行名
rownames(final_data_median) <- final_data_median$GeneName  
final_data_median <- as.data.frame(final_data_median)
rownames(final_data_median) <- final_data_median$GeneName
final_data_median$GeneName <- NULL


group_info <- data.frame(
  Library_Name = c("GSM1248612","GSM1248614","GSM1248616","GSM1248618","GSM1248620","GSM1248622","GSM1248624","GSM1248626","GSM1248628","GSM1248630","GSM1248632","GSM1248634","GSM1248636","GSM1248638","GSM1248640","GSM1248642","GSM1248644","GSM1248646","GSM1248648","GSM1248650","GSM1248652","GSM1248654","GSM1248656","GSM1248660","GSM1248662","GSM1248664","GSM1248615","GSM1248617","GSM1248621","GSM1248635","GSM1248637","GSM1248639","GSM1248643","GSM1248645","GSM1248647","GSM1248649","GSM1248651","GSM1248653","GSM1248657","GSM1248663","GSM1248613","GSM1248619","GSM1248623","GSM1248625","GSM1248627","GSM1248629","GSM1248631","GSM1248633","GSM1248641","GSM1248655","GSM1248661","GSM1248665"),
  Group = c(rep("Gastric_normal_EBV_negative", 26), rep("Gastric_tumor_EBV_neagtive", 14), rep("Gastric_tumor_EBV_positive", 12)))



sample_order <- colnames(final_data_median)[-1]
group_info <- group_info[match(sample_order, group_info$Library_Name), ]

final_data_median <- as.matrix(final_data_median[, -1])

groups <- factor(group_info$Group)
design <- model.matrix(~0 + groups)


colnames(design) <- levels(groups)
contrast_matrix <- makeContrasts(GroupBvsA = Gastric_tumor_EBV_positive - Gastric_tumor_EBV_neagtive, levels = design)


# 使用limma进行分析
fit <- lmFit(final_data_median, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)


results <- topTable(fit2, coef = "GroupBvsA", number = Inf, adjust.method = "BH")
results$GeneSymbol <- rownames(results)  # 添加基因名列
#View(results)
#write.csv(results,"GSE51575_limma.csv",quote = TRUE)


results$change <- factor(
  ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
         ifelse(results$logFC > 1, "Up", "Down"), "Not Sig"),
  levels = c("Down", "Not Sig", "Up")
)

head(results)

#simulated_miRNA <- read.csv("D:/0work/0wholetransctiptome/2sRNAminic/9cancer/0006mer_target_with_hsaID_GENE_allTarget.csv", stringsAsFactors = FALSE)
simulated_miRNA <- read.csv("D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target_with_hsaID.csv", stringsAsFactors = FALSE)
simulated_miRNA <- simulated_miRNA %>% filter(abbreviation == "HHV-4")
genes_split <- strsplit(simulated_miRNA$Intersection_Genes, ";", fixed = TRUE)
all_genes <- unlist(genes_split)
unique_genes <- unique(all_genes)
print(unique_genes)


logFC_threshold <- 1  # 对数变化阈值
adjP_threshold <- 0.05  # 校正后的P值阈值
results1 <- results %>%
  mutate(miRNA_ID = rownames(.)) %>%  # 将行名添加为miRNA_ID列
  mutate(Status = case_when(
    logFC > logFC_threshold & adj.P.Val < adjP_threshold ~ "Up",
    logFC < -logFC_threshold & adj.P.Val < adjP_threshold ~ "Down",
    TRUE ~ "Not Sig"
  ))

# 读取抑癌基因列表
# 读取胃癌相关的抑癌基因列表
cancer_gene_file <- "D:/0work/0wholetransctiptome/2sRNAminic/9cancer/CancerGene.xlsx"
tumor_suppressor_genes <- readxl::read_excel(cancer_gene_file, sheet = "GC_TSGene") %>% 
  distinct(GeneSymbol) %>% 
  pull(GeneSymbol)


missing_suppressor_genes <- setdiff(tumor_suppressor_genes, results1$miRNA_ID)

# 打印缺失基因的数量和列表
if (length(missing_suppressor_genes) > 0) {
  cat("以下", length(missing_suppressor_genes), "个抑癌基因未在results1中找到:\n")
  print(missing_suppressor_genes)
} else {
  cat("所有抑癌基因均在results1中找到。\n")
}


# 筛选模拟基因中的抑癌基因
simulated_suppressor <- results1 %>%
  filter(miRNA_ID %in% unique_genes) %>%  # 先筛选模拟基因
  filter(miRNA_ID %in% tumor_suppressor_genes) %>%  # 再筛选抑癌基因
  mutate(Is_Suppressor = TRUE)

# 统计模拟基因中的抑癌基因Status
suppressor_in_simulated_summary <- simulated_suppressor %>%
  count(Status) %>%
  rename(Suppressor_in_simulated_n = n)
View(simulated_suppressor)

all_suppressor_summary <- results1 %>%
  filter(miRNA_ID %in% tumor_suppressor_genes) %>%
  count(Status, name = "All_suppressor_n") %>%
  # 添加总计行
  bind_rows(
    tibble(
      Status = "Total",
      All_suppressor_n = sum(.$All_suppressor_n)
    ))

suppressor_summary <- suppressor_summary %>%
  bind_rows(
    tibble(
      Status = "Total",
      Count = sum(suppressor_summary$Count),
      Proportion = 100,
      Category = "All Tumor Suppressors"
    )
  )

# 模拟基因整体统计（保持原有）
simulated_summary <- results1 %>%
  filter(miRNA_ID %in% unique_genes) %>%
  count(Status) %>%
  rename(Simulated_n = n)

# 合并统计结果并计算比例
summary <- simulated_summary %>%
  left_join(suppressor_in_simulated_summary, by = "Status") %>%
  # 替换NA为0
  mutate(Suppressor_in_simulated_n = replace_na(Suppressor_in_simulated_n, 0)) %>%
  # 计算比例
  mutate(
    # 抑癌基因占该状态模拟基因的比例
    Suppressor_prop_in_simulated = Suppressor_in_simulated_n / Simulated_n * 100,
    # 抑癌基因占所有模拟基因的比例
    Suppressor_overall_prop = Suppressor_in_simulated_n / sum(Suppressor_in_simulated_n) * 100
  )

# 添加总计行
total_simulated <- nrow(results1 %>% filter(miRNA_ID %in% unique_genes))
total_suppressor_in_simulated <- nrow(simulated_suppressor)

summary <- summary %>%
  bind_rows(
    tibble(
      Status = "Total",
      Simulated_n = total_simulated,
      Suppressor_in_simulated_n = total_suppressor_in_simulated,
      Suppressor_prop_in_simulated = total_suppressor_in_simulated / total_simulated * 100,
      Suppressor_overall_prop = 100
    )
  )

# 格式化输出
formatted_summary <- summary %>%
  mutate(across(where(is.numeric), ~round(., 2))) %>%
  rename(
    "Status" = Status,
    "模拟基因数" = Simulated_n,
    "抑癌基因数" = Suppressor_in_simulated_n,
    "该状态中抑癌基因占比(%)" = Suppressor_prop_in_simulated,
    "抑癌基因分布(%)" = Suppressor_overall_prop
  )

final_summary <- formatted_summary %>%
  left_join(all_suppressor_summary, by = "Status") %>%
  # 添加所有抑癌基因的比例列
  mutate(
    "所有抑癌基因占比(%)" = round(All_suppressor_n / sum(All_suppressor_n[Status == "Total"]) * 100, 2)
  ) %>%
  # 重新排列列顺序
  select(
    "Status", 
    "模拟基因数", 
    "抑癌基因数", 
    "所有抑癌基因数" = All_suppressor_n,
    "该状态中抑癌基因占比(%)",
    "抑癌基因分布(%)",
    "所有抑癌基因占比(%)"
  )
# 打印美观结果
print(final_summary)

#write_csv(formatted_summary, "D:/0work/0wholetransctiptome/2sRNAminic/9cancer/EBV_Disease/tumor_mimic_summary.csv")














highlight_miRNAs <-unique_genes
geneList <- results[highlight_miRNAs,]

results$change <- factor(
  ifelse(results$adj.P.Val < 0.05 & abs(results$logFC) > 1,
         ifelse(results$logFC > 1, "Up", "Down"), "Not Sig"),
  levels = c("Down", "Not Sig", "Up")
)
head(results)


p1<-ggplot(results, aes(logFC, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999",linewidth=1) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999",linewidth=1) +
  geom_point(aes(color = change), size = 2.5, alpha = 0.5) +
  geom_point(data=geneList,aes(x = logFC, y = -log10(adj.P.Val)),colour="yellow",size=2)+
  theme_bw(base_size = 12) +
  #ggsci::scale_color_npg() +  # npg期刊配色
  scale_color_manual(values = c("Down" = "#6BAED6", "Up" = "#FB6A4A", "Not Sig" = "grey")) +
  theme(panel.grid = element_blank(), legend.position = 'right') +
  xlab("Log2 Fold Change") +
  ylab("-Log10 (FDR Adjusted P-value)")+scale_x_continuous(limits = c(-10, 10))
p1



index <- which(simulated_results$Status %in% c("Down", "Up"))

highlight_DEmiRNAs <- row.names(simulated_results)[index]
#highlight_DEmiRNAs <- c("PTPRC","RASSF5","TNFAIP8L2","BASP1","CNTNAP2","LEFTY1","OLFM4")
highlight_DEmiRNAs <- c("HOXD10")


geneList1 <- results[rownames(results) %in% highlight_DEmiRNAs,]
geneList1 <- subset(geneList1, select = -change)
geneList1$label <- rownames(geneList1)
head(geneList1)
geneList1$change <- ifelse(geneList1$logFC > 0, "Up", "Down")
pRNA<-p1 + geom_label_repel(data = geneList1, 
                     aes(x = logFC, y = -log10(adj.P.Val), label = label, color = change),  # 添加 color = change
                     size = 3.5,
                     box.padding = unit(0.4, "lines"), 
                     segment.color = "black",   # 连线的颜色
                     segment.size = 0.4, 
                     max.overlaps = 20) +
  scale_color_manual(values = c("Up" = "#DE2D26", "Down" = "#3182BD"))

pRNA



library(ggpubr)
ggarrange(psRNA,pRNA,ncol = 2, nrow = 1,legend="right",align = "h",common.legend = FALSE, labels = c("A","B")) 

