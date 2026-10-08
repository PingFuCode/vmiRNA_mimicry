library(readxl)
library(dplyr)
library(tidyr)
setwd("D:/0work/0wholetransctiptome/2sRNAminic/4vRNA-Disease/RNAdiseaseResults/AllVirus")
# 读取所有的 .xlsx 文件（排除 "All.xlsx"）
files <- list.files("D:/0work/0wholetransctiptome/2sRNAminic/4vRNA-Disease/RNAdiseaseResults/AllVirus/", 
                    pattern = ".xlsx", full.names = TRUE)
files <- files[!grepl("All.xlsx", files)]  # 过滤掉 "All.xlsx"

all_data <- list()

for (file in files) {
  # 读取 Excel 文件
  virus_data <- read_excel(file)
  
  # 检查是否为空
  if (nrow(virus_data) == 0) {
    cat("警告: 文件", file, "为空，跳过该文件。\n")
    next
  }
  
  # 提取病毒名称（去除 .xlsx 后缀）
  virus_name <- sub(".xlsx$", "", basename(file))
  
  # 只保留 "Disease Name" 和 "FDR" 列
  virus_data <- virus_data %>% select(`Disease Name`, `FDR`)
  
  # 重命名 P Value 列（添加病毒名称）
  colnames(virus_data)[colnames(virus_data) == "FDR"] <- virus_name
  
  # 存入列表
  all_data[[virus_name]] <- virus_data
}

# 合并所有数据
combined_data <- Reduce(function(x, y) full_join(x, y, by = "Disease Name"), all_data)

# 填充缺失值（如果某些疾病在某些病毒中没有 P 值，则填充 NA）
combined_data[is.na(combined_data)] <- 1  # 设为1表示无显著性

# 保存合并后的数据
write.csv(combined_data, file = "merged_virus_disease.csv", row.names = FALSE)

# 查看数据
View(combined_data)

combined_data <- read.csv("merged_virus_disease.csv", header = TRUE, row.names = 1)

# Step 3: 筛选至少在两个病毒中出现的疾病
pathway_counts <- combined_data %>%
  mutate(VirusCount = rowSums(combined_data[-1] < 1, na.rm = TRUE)) %>%
  #filter(VirusCount >= 2) %>%
  filter(VirusCount >= 1) %>%
  select(-VirusCount)  # 移除临时列

# Step 4: 形成热图矩阵
heatmap_data <- pathway_counts %>%
  pivot_longer(cols = -`Disease Name`, names_to = "Virus", values_to = "FDR") %>%
  pivot_wider(names_from = "Virus", values_from = "FDR", values_fill = 1)

write.csv(heatmap_data, file = "merged_virus_disease_virusall.csv", row.names = FALSE)








df <- read.csv("merged_virus_disease_virusall.csv", row.names = 1, check.names = FALSE)
df_matrix <- as.matrix(df)
df_matrix
# 定义病毒和病毒家族的映射
virus_family_mapping <- data.frame(
  VirusFamily = c("Adenoviridae", "Anelloviridae", "Filoviridae", "Flaviviridae", "Flaviviridae", "Flaviviridae", 
                  "Flaviviridae", "Herpesviridae", "Herpesviridae", "Herpesviridae", "Herpesviridae", "Herpesviridae", 
                  "Herpesviridae", "Herpesviridae", "Orthomyxoviridae", "Papillomaviridae", "Picornaviridae", 
                  "Pneumoviridae", "Pneumoviridae", "Polyomaviridae", "Polyomaviridae", "Polyomaviridae", 
                  "Poxviridae", "Retroviridae", "Togaviridae"),
  Virus = c("HAdV-C", "TTV", "ZEBOV", "DENV", "HCV", "WNV", "ZIKV", "HSV-1", "HHV-5", "HHV-6", "HHV-4", 
            "HHV-8", "SaHV-2", "Herpesvirus", "IAV", "HPV", "EVA", "HMPV", "HRSV", "MCPyV", "HPyV-1", "HPyV-2", 
            "VACV", "HIV-1", "SINV")
)

# 假设df_matrix的列名是病毒的名称
virus_names <- colnames(df_matrix)

# 使用merge来匹配VirusFamily
merged_data <- merge(data.frame(Virus = virus_names), virus_family_mapping, by = "Virus", all.x = TRUE)

# 提取匹配的VirusFamily列
VirusFamily <- merged_data$VirusFamily

ordered_indices <- order(VirusFamily)
df_matrix_sorted <- df_matrix[, ordered_indices]  # 按照病毒家族排序

library(ComplexHeatmap)
# 创建行注释
ha_column <- HeatmapAnnotation(df = data.frame(Type = c(rep("High Pathogenicity Virus", 1), 
                                                        rep("Low Pathogenicity Virus", 4),
                                                        rep("Conditionally high pathogenicity", 2))),
                               col = list(Type = c("High Pathogenicity Virus" =  "#92C5DE", 
                                                   "Low Pathogenicity Virus" = "#FC8D62",
                                                   "Conditionally high pathogenicity" = "#999999"), 
                                          width = unit(0.5, "cm")))



disease_type_map <- data.frame(
  Disease = c(
    "Bladder Cancer", "Nasopharynx Carcinoma", "Amyotrophic Lateral Sclerosis",
    "Oral Squamous Cell Carcinoma", "Malignant Glioma", "Cervical Cancer",
    "Cardiovascular Disease", "Leukemia", "Breast Neoplasm", "Melanoma",
    "Stomach Cancer", "NonSmall Cell Lung Cancer", "Ovarian Cancer",
    "Multiple Sclerosis", "Osteosarcoma", "Lung Adenocarcinoma",
    "Renal Cell Carcinoma", "Lung Cancer", "Glioma", "Prostate Cancer",
    "Hepatocellular Carcinoma", "Breast Cancer", "Alzheimer Disease",
    "Glioblastoma", "Pancreatic Cancer", "Colon Cancer", "Colorectal Cancer",
    "Esophageal Cancer", "Gastric Cancer", "Parkinson Disease"
  ),
  Type = c(
    "Cancer",              # Bladder Cancer
    "Cancer",              # Nasopharynx Carcinoma
    "Neurological",        # ALS
    "Cancer",              # Oral SCC
    "Cancer",              # Malignant Glioma
    "Cancer",              # Cervical Cancer
    "Cardiovascular",      # Cardiovascular Disease
    "Cancer",              # Leukemia
    "Cancer",              # Breast Neoplasm
    "Cancer",              # Melanoma
    "Cancer",              # Stomach Cancer
    "Cancer",              # NSCLC
    "Cancer",              # Ovarian Cancer
    "Autoimmune/Inflammatory",  # Multiple Sclerosis
    "Cancer",              # Osteosarcoma
    "Cancer",              # Lung Adenocarcinoma
    "Cancer",              # Renal Cell Carcinoma
    "Cancer",              # Lung Cancer
    "Cancer",              # Glioma
    "Cancer",              # Prostate Cancer
    "Cancer",              # Hepatocellular Carcinoma
    "Cancer",              # Breast Cancer
    "Neurological",        # Alzheimer Disease
    "Cancer",              # Glioblastoma
    "Cancer",              # Pancreatic Cancer
    "Cancer",              # Colon Cancer
    "Cancer",              # Colorectal Cancer
    "Cancer",              # Esophageal Cancer
    "Cancer",              # Gastric Cancer
    "Neurological"         # Parkinson Disease
  ),
  stringsAsFactors = FALSE
)


disease_type_map <- disease_type_map[disease_type_map$Disease %in% rownames(df_matrix_sorted), ]
disease_type_map <- disease_type_map[order(disease_type_map$Type), ]
df_matrix_sorted <- df_matrix_sorted[disease_type_map$Disease, ]
#View(df_matrix_sorted)
# 4. 创建颜色映射
type_colors <- c(
  "Neurological" = "#92C5DE",
  "Cancer"  = "#FC8D62",
  "Cardiovascular" = "#91D1C2B2",
  "Autoimmune/Inflammatory"= "#F39B7FB2",
  "Infectious" = "#A6D854",
  "Metabolic/Renal" = "#E78AC3",
  "Cardiovascular/Metabolic" = "#B3B3B3",
  "Psychiatric/Addiction" = "#E5C494",
  "Reproductive/Other" = "#A6761D",
  "Congenital/Structural" = "#999999",
  "Metabolic/Pre-cancerous" = "#B15928"
)

disease_type_map$Malignancy <- c(
  "Malignant or acute",     # Bladder Cancer
  "Malignant or acute",     # Nasopharynx Carcinoma
  "Benign or chronic",        # ALS – 非肿瘤性，退行性神经疾病
  "Malignant or acute",     # Oral Squamous Cell Carcinoma
  "Malignant or acute",     # Malignant Glioma
  "Malignant or acute",     # Cervical Cancer
  "Benign or chronic",        # Cardiovascular Disease – 不是肿瘤
  "Malignant or acute",     # Leukemia
  "Unsure",     # Breast Neoplasm – 默认为恶性（若为 benign 则需特指）
  "Malignant or acute",     # Melanoma
  "Malignant or acute",     # Stomach Cancer
  "Malignant or acute",     # Non−Small Cell Lung Cancer
  "Malignant or acute",     # Ovarian Cancer
  "Benign or chronic",        # Multiple Sclerosis – 自免疾病，非恶性
  "Malignant or acute",     # Osteosarcoma
  "Malignant or acute",     # Lung Adenocarcinoma
  "Malignant or acute",     # Renal Cell Carcinoma
  "Malignant or acute",     # Lung Cancer
  "Malignant or acute",     # Glioma – 若未指明为 benign glioma，默认为 malignant
  "Malignant or acute",     # Prostate Cancer
  "Malignant or acute",     # Hepatocellular Carcinoma
  "Malignant or acute",     # Breast Cancer
  "Benign or chronic",        # Alzheimer Disease – 非肿瘤性
  "Malignant or acute",     # Glioblastoma
  "Malignant or acute",     # Pancreatic Cancer
  "Malignant or acute",     # Colon Cancer
  "Malignant or acute",     # Colorectal Cancer
  "Malignant or acute",     # Esophageal Cancer
  "Malignant or acute",     # Gastric Cancer
  "Benign or chronic"         # Parkinson Disease – 非恶性
)


malignancy_colors <- c(
  "Malignant or acute" = "#F39B7FB2",
  "Benign or chronic" = "#91D1C2B2",
  "Unsure" = "#bdbdbd"
)


# 5. 创建 row annotation
ha_row <- rowAnnotation(
  Type = disease_type_map$Type,
  Malignancy = disease_type_map$Malignancy,
  col = list(
    Type = type_colors,
    Malignancy = malignancy_colors
  ),
  width = unit(1, "cm")
)



# 在绘图之前，将行名中的点号替换为空格
rownames(df_matrix_sorted) <- gsub("\\.", " ", rownames(df_matrix_sorted))
rownames(df_matrix_sorted) <- tools::toTitleCase(rownames(df_matrix_sorted))


target_order <- c("HIV-1", "HHV-4", "HHV-5", "HHV-8", "DENV", "SINV", "HSV-1")
df_matrix_sorted <- df_matrix_sorted[, target_order]


df_matrix_sorted_log10 <- as.data.frame(df_matrix_sorted) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), NA, -log10(.))))
#View(df_matrix_sorted_log10)

# 创建并保存热图
pdf("heatmap_DISEASE_virusall.pdf", width =7.5, height = 6)  # 设置PDF尺寸
library(circlize)
#col_fun <- colorRamp2(c(0, 0.05,0.06,0.5, 1), c("#E41A1C", "#FFDFDF", "white", "white", "white"))
library(RColorBrewer)
reds <- colorRampPalette(brewer.pal(9, "Reds"))(8) 
#"#F7FBFF" "#DAE8F5" "#BAD6EB" "#88BEDC" "#539ECC" "#2A7AB9" "#0B559F" "#08306B"
#"#FFF5F0" "#FDDACB" "#FCAF93" "#FB8060" "#F44F38" "#D52221" "#AA1016" "#67000D"
col_fun <- colorRamp2(c(1, 10,  20, 30), c("#DAE8F5", "#88BEDC","#539ECC", "#2A7AB9"))
#col_fun <- colorRamp2(c(1, 10,  20, 30), c("#FDDACB", "#FCAF93","#FB8060", "#AA1016"))

ht <- Heatmap(df_matrix_sorted_log10, 
              name = "-log10 Transformed Adjusted p-values",
              na_col = "white",  # 设置NA的颜色
              rect_gp = gpar(col = "white", lwd = 0.5),  # 矩形边框设置
              border = TRUE,  # 显示外框线
              border_gp = gpar(col = "black", lwd = 1),  # 设置外框线颜色和宽度
              col = col_fun,  # 使用自定义颜色映射
              cluster_columns = FALSE,  # 进行列聚类
              cluster_rows = FALSE,  # 进行行聚类
              #bottom_annotation = ha_column,
              right_annotation = ha_row)
              #,
              #cell_fun = function(j, i, x, y, width, height, fill) {
              #  if (!is.na(df_matrix_sorted[i, j]) && df_matrix_sorted[i, j] <= 0.06) {
              #    grid.text(sprintf("%.2f", df_matrix_sorted[i, j]), x, y, gp = gpar(fontsize = 10))
              #  }
              #})  # 行注释


# 绘制热图并保存
draw(ht)
dev.off()





library(tidyverse)

df <- read_csv("D:/0work/0wholetransctiptome/2sRNAminic/3mRNATarget/0006mer_target.csv")

gene_cnt <- df %>% 
  separate_rows(Intersection_Genes, sep = ";") %>%   # 拆成行
  group_by(abbreviation) %>% 
  summarise(gene_n = n_distinct(Intersection_Genes)) # 去重后计数

print(max(gene_cnt$gene_n))
print(min(gene_cnt$gene_n))






