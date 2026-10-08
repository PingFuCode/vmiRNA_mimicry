library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
#setwd("E:\\修改成数据库格式\\VirusCircBase2.0\\heatmap")
setwd('D:\\0work\\0wholetransctiptome\\2sRNAminic\\5tissue_specific\\R\\')

display.brewer.pal(5,"Set3")
brewer.pal(5,"Set3")
getPalette = colorRampPalette(brewer.pal(6,"Reds"))
color=getPalette(12)
color

#host species 去掉Herpesvirus
df = read.delim("VirusTissue_mean.txt", row.names = 1,sep = '\t', stringsAsFactors = FALSE, check.names = FALSE )
virus_order <- c("HAdV-C", "HHV-4", "HHV-5", "HHV-6", "HHV-8", "HSV-1", 
                 "SaHV-2", "HPyV-1", "HPyV-2", "MCPyV", "TTV", "IAV", "HMPV", 
                 "HRSV", "DENV", "HCV", "WNV", "ZIKV", "EVA", "SINV", "HIV-1")

System_category <- c(
  "NCI_60_cancer_panel" , "blood" , "heart" , "artery" , "pericardium" ,  
  "vein" , "connective_tissue" , "adipose" , "bowel" , "liver" , "oral_cavity" , 
  "stomach" , "saliva" , "pancreas" , "feces" ,  "salivary_glands" , "ascites" , 
  "esophagus" ,  "submandibular_gland" ,  "gallbladder" , "glandular_breast_tissue" , 
  "adrenal_gland" , "thyroid" , "kidney" ,  "urine" , "bladder" , "urethra" ,  "bone_marrow" ,  
  "stem_cells" ,  "immune_system" , "tonsil" , "spleen" ,  "lymph_node" ,  "skin" , "sweat" ,
  "limb_muscle" , "paraspinal_muscle" , "diaphragm" , "smooth_muscle" , "brain" ,  "nerve" , 
  "spinal_cord" , "milk" , "placenta" , "cervix" ,  "testis" , "uterus" ,  "prostate" , "ovary" ,
  "umbilical_cord" , "vaginal_tissue" ,  "bronchus" ,  "lung" , "pleurae" ,  "airway" , 
  "trachea" , "otic_vesicle" , "retina" , "cornea" ,  "conjunctiva" , "aqueous_humor" , 
  "tears" , "tongue" , "sclera" , "bone" 
)
df$abbreviation <- factor(rownames(df), levels = virus_order)
# 按照因子的顺序排序
df_sorted <- df[order(df$abbreviation), ]
df <- subset(df_sorted, select = -c(virus_species, viral_family, viral_group))
dim(df)
length(System_category)
df <- df[, System_category]
dim(df)

df_log = log2(df + 1)
colnames(df_log) <- sapply(colnames(df_log), function(x) paste(toupper(substring(x, 1, 1)), substring(x, 2), sep = ""))
#df[is.na(df)] = 0
df_matrix <- as.matrix(df_log)
View(df_matrix)
df_matrix_rounded <- round(df_matrix, 2)
write.csv(df_matrix_rounded, "da_matrix.csv", row.names = TRUE)













################################################################################################
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)
setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\5tissue_specific\\")

# 读取病毒感染信息
virus_infection_table <- read.table("Humanvirus_organVirusData_1", header = TRUE, sep = "\t", fill = TRUE)

# 统一组织名称
virus_infection_table <- virus_infection_table %>%
  mutate(Organ = case_when(
    Organ == "urine" ~ "Urethra",
    Organ == "trachea_bronchi" ~ "Trachea",
    Organ == "saliva" ~ "Submandibular_gland",
    Organ == "skeletal muscle tissue" ~ "Limb_muscle",
    Organ == "muscle of leg" ~ "Limb_muscle",
    Organ == "thyroid gland" ~ "Thyroid",
    Organ == "salivary gland" ~ "Salivary_glands",
    Organ == "adrenal gland" ~ "Adrenal_gland",
    Organ == "blood vessel" ~ "Blood",
    Organ == "vagina" ~ "Vaginal_tissue",
    TRUE ~ Organ  # 其他值保持不变
  ))

virus_infection_table$Organ <- gsub(" ", "_", virus_infection_table$Organ)  # 替换空格为 "_"


virus_abbreviation <- data.frame(
  abbreviation = c("DENV", "EVA", "HAdV-C", "HCV", "HHV-4", "HHV-5", "HHV-6", "HHV-8", 
                   "HIV-1", "HMPV", "HPyV-1", "HPyV-2", "HRSV", "HSV-1", "Herpesvirus", 
                   "IAV", "MCPyV", "SINV", "SaHV-2", "TTV", "WNV", "ZIKV"),
  virus_species = c("Dengue virus", "Enterovirus A", "Human mastadenovirus C", 
                    "Hepacivirus C", "Human gammaherpesvirus 4", "Human betaherpesvirus 5", 
                    "Human betaherpesvirus 6", "Human gammaherpesvirus 8", 
                    "Human immunodeficiency virus 1", "Human metapneumovirus", 
                    "Betapolyomavirus hominis", "Betapolyomavirus secuhominis", 
                    "Human orthopneumovirus", "Human alphaherpesvirus 1", 
                    "unKnown", "Influenza A virus", "Alphapolyomavirus quintihominis", 
                    "Sindbis virus", "Saimiriine gammaherpesvirus 2", 
                    "Torque teno virus", "West Nile virus", "Zika virus"),
  stringsAsFactors = FALSE
)

# 将 virus_infection_table 中的 Species 替换为缩写
virus_infection_table <- virus_infection_table %>%
  left_join(virus_abbreviation, by = c("Species" = "virus_species")) %>%
  select(abbreviation, Organ) %>%
  rename(Virus = abbreviation) %>%
  mutate(Organ = tolower(Organ))



#View(virus_infection_table)

# 读取表达量矩阵
df <- read.csv("R\\da_matrix.csv", row.names = 1, check.names = FALSE)
df_matrix <- as.matrix(df)
colnames(df_matrix) <- tolower(colnames(df_matrix))

infection_matrix <- matrix(
  FALSE, 
  nrow = nrow(df_matrix), 
  ncol = ncol(df_matrix),
  dimnames = list(rownames(df_matrix), colnames(df_matrix))
)

# 填充感染矩阵
virus_organs <- virus_infection_table %>%
  select(Virus, Organ) %>%
  distinct()

infection_matrix <- matrix(
  "",  # 初始化为空字符串，而不是 FALSE
  nrow = nrow(df_matrix), 
  ncol = ncol(df_matrix),
  dimnames = list(rownames(df_matrix), colnames(df_matrix))
)

for (k in 1:nrow(virus_organs)) {
  virus <- virus_organs$Virus[k]
  organ <- virus_organs$Organ[k]
  if (virus %in% rownames(infection_matrix) && organ %in% colnames(infection_matrix)) {
    infection_matrix[virus, organ] <- "*"  # 将 TRUE 替换为 "*"
  }
}

#View(infection_matrix)

print(infection_matrix[1:5, 1:5])

VirusGroup <- c(
  "dsDNA", "dsDNA", "dsDNA", "dsDNA", "dsDNA", "dsDNA", "dsDNA", "dsDNA", "dsDNA", 
  "dsDNA", "ssDNA(-)", "ssRNA(-)", "ssRNA(-)", "ssRNA(-)", "ssRNA(+)", "ssRNA(+)", "ssRNA(+)", 
  "ssRNA(+)", "ssRNA(+)", "ssRNA(+)", "ssRNA-RT"
)
VirusFamily <- c(
  "Adenoviridae", "Herpesviridae", "Herpesviridae", "Herpesviridae", "Herpesviridae", 
  "Herpesviridae", "Herpesviridae", "Polyomaviridae", "Polyomaviridae", "Polyomaviridae", "Anelloviridae",
  "Orthomyxoviridae", "Pneumoviridae", "Pneumoviridae", "Flaviviridae", "Flaviviridae", "Flaviviridae", 
  "Flaviviridae", "Picornaviridae", "Togaviridae", "Retroviridae"
)

System_category<-c(
  "Cancer Research", "Circulatory System", "Circulatory System", "Circulatory System", 
  "Circulatory System", "Circulatory System", "Connective Tissue", "Connective Tissue", 
  "Digestive System", "Digestive System", "Digestive System", "Digestive System", 
  "Digestive System", "Digestive System", "Digestive System", "Digestive System", 
  "Digestive System", "Digestive System", "Digestive System", "Digestive System", 
  "Endocrine System", "Endocrine System", "Endocrine System", "Excretory System", 
  "Excretory System", "Excretory System", "Excretory System", "Hematopoietic System", 
  "Immune System", "Immune System", "Immune System", "Immune System", "Immune System", 
  "Integumentary System", "Integumentary System", "Muscular System", "Muscular System",
  "Muscular System", "Muscular System", "Nervous System", "Nervous System", "Nervous System", 
  "Reproductive System", "Reproductive System", "Reproductive System", "Reproductive System",
  "Reproductive System", "Reproductive System", "Reproductive System", "Reproductive System", 
  "Reproductive System", "Respiratory System", "Respiratory System", "Respiratory System",
  "Respiratory System", "Respiratory System", "Sensory Organs", "Sensory Organs", 
  "Sensory Organs", "Sensory Organs", "Sensory Organs", "Sensory Organs", "Sensory Organs", 
  "Sensory Organs", "Skeletal System"
)
length((System_category))
# 创建行注释
ha_row <- rowAnnotation(
  df = data.frame(VirusGroup = VirusGroup,VirusFamily = VirusFamily ),
  col = list(
    VirusGroup = c(
      "dsDNA" = "#8A5F93", "ssDNA(-)" = "#4AA956", "ssRNA(-)" = "#B65C73", "ssRNA(+)" = "#FF9D0C", 
      "ssRNA-RT" = "#E41A1C"
    ),
    VirusFamily = c(
      "Adenoviridae" = "#D6604D", "Herpesviridae" = "#F4A582", "Polyomaviridae" = "#FDDBC7",
      "Anelloviridae" = "#D1E5F0", "Orthomyxoviridae" = "#92C5DE", "Pneumoviridae" = "#3A85A8",
      "Flaviviridae" = "#66C2A5", "Picornaviridae" = "#FC8D62", "Togaviridae" = "#8DA0CB",
      "Retroviridae" = "#FF5733"
    )
  ),
  width = unit(1, "npc")
)


System_category_col = c(
  "Circulatory System" = "#D32F2F",   # Red (High Contrast)
  "Digestive System" = "#FFEB3B",     # Yellow
  "Nervous System" = "#2196F3",       # Blue
  "Respiratory System" = "#4CAF50",   # Green
  "Excretory System" = "#FF9800",     # Orange
  "Endocrine System" = "#9C27B0",     # Purple
  "Integumentary System" = "#795548", # Brown
  "Reproductive System" = "#607D8B",  # Greyish Blue
  "Sensory Organs" = "#9E9E9E",       # Light Grey
  "Muscular System" = "#673AB7",      # Deep Purple
  "Connective Tissue" = "#CDDC39",    # Lime
  "Hematopoietic System" = "#3F51B5", # Indigo Blue
  "Cancer Research" = "#FF5722",      # Red-Orange
  "Immune System" = "#8BC34A",        # Light Green
  "Skeletal System" = "#F4A582"       # Light Salmon
)

ha_column <- HeatmapAnnotation(
  df = data.frame(System_category = System_category),
  col = list(System_category = System_category_col),
  width = unit(1, "npc")  # 自动调整宽度
)

library(tools)
colnames(df_matrix) <- toTitleCase(colnames(df_matrix))
ht <- Heatmap(
  df_matrix, 
  name = "log2(RPM)",
  na_col = "white",
  rect_gp = gpar(col = "white", lwd = 0.5),
  col = rev(brewer.pal(n = 7, name = "RdYlBu")),
  cluster_columns = FALSE, 
  cluster_rows = FALSE, 
  row_order = 1:nrow(df_matrix),
  column_order = 1:ncol(df_matrix),
  bottom_annotation = ha_column,
  left_annotation = ha_row,
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 如果该 (病毒, 组织) 需要标记 "*"
    if (infection_matrix[i, j] == "*") {
      grid.text(
        "*", 
        x = x, 
        y = y,  # 使用原始 y 位置
        just = c("center", "center"),  # 水平和垂直居中
        gp = gpar(fontsize = 14, col = "black")
      )
    }
  }
)

# 绘制热图
draw(ht)


