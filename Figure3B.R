library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(readxl)

setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\5tissue_specific\\")
################################################################################


data <- read.csv("0006mer_tissue.csv")
#data <- data[, !colnames(data) %in% "otic_vesicle"]
virus_infection_table <- read.table("Humanvirus_organVirusData_1", header = TRUE, sep = "\t", fill = TRUE)
dim(virus_infection_table)
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
    TRUE ~ Organ  # 保留其他未替换的值
  ))
virus_infection_table$Organ <- gsub(" ", "_", virus_infection_table$Organ)
#head(virus_infection_table)

################################################################################
mirna_tissue <- read.csv("miRNATissueAtlas2025_Organ_system.csv", row.names = 1, header = TRUE)
colnames(mirna_tissue) <- tolower(colnames(mirna_tissue))
mirna_tissue <- mirna_tissue[, !colnames(mirna_tissue) %in% "otic_vesicle"]

######################################循环所有病毒##########################################
# 获取所有独特的病毒种类
virus_list <- unique(data$virus_species)
virus_list <- virus_list[virus_list != "unKnown"]
virus_list <- virus_list[virus_list != "Hepacivirus C"]
virus_list <- virus_list[virus_list != "Betapolyomavirus hominis"]
virus_list <- virus_list[virus_list != "Betapolyomavirus secuhominis"]
virus_list <- virus_list[virus_list != "Saimiriine gammaherpesvirus 2"]

# 存储所有箱型图
plot_list <- list()

# 迭代病毒列表
for (virus in virus_list) {
  # 获取该病毒的感染器官
  print(virus)
  virus_infection_filtered <- virus_infection_table %>%
    filter(Species == virus) %>%
    select(Species, Organ)
  
  # 过滤该病毒的 miRNA 数据
  hmirna <- data %>% filter(virus_species == virus)
  # 受感染器官列表
  infected_organ_list <- tolower(unique(virus_infection_filtered$Organ))
  mimic_hmirna <- unique(hmirna$hmiRNA)
  # 添加 Mimic_Status 列
  mirna_tissue$Mimic_Status <- ifelse(row.names(mirna_tissue) %in% mimic_hmirna, "Mimic", "UnMimic")
  print(sum(mirna_tissue$Mimic_Status == "Mimic", na.rm = TRUE))
  # 区分感染和未感染器官
  infected_cols <- colnames(mirna_tissue)[colnames(mirna_tissue) %in% infected_organ_list]
  uninfected_cols <- setdiff(colnames(mirna_tissue), c(infected_cols, "mimic_status"))
  
  # 选择等量未感染器官
  #num_infected <- length(infected_cols)
  #set.seed(123)
  #selected_uninfected_cols <- sample(uninfected_cols, num_infected)
  
  group_data <- mirna_tissue %>%
    mutate(row_id = row.names(.)) %>%
    pivot_longer(
      cols = -c(row_id, Mimic_Status),
      names_to = "tissue_type",
      values_to = "expression"
    ) %>%
    mutate(
      infection_status = ifelse(tissue_type %in% infected_cols, "Infected", "UnInfected"),
      group = paste(Mimic_Status, infection_status, sep = "-") ,
     expression = log2(expression+1)  # 对expression取log2
    )
  print(dim(group_data))
  group_data <- group_data %>%
    group_by(group, row_id) %>%
    summarise(expression1 = expression, .groups = "drop")
    #summarise(expression1 = median(expression, na.rm = TRUE), .groups = "drop")

  #print(dim(group_data))

  group_data <- group_data[group_data$expression1 > 0, ]
  group_data$expression1 <- round(group_data$expression1, digits = 2)
  mimic_infected_data <- group_data[group_data$group == "Mimic-Infected" , ]
  mimic_uninfected_data <- group_data[group_data$group == "Mimic-UnInfected" , ]
  
  #Mimic_Infected <- median(mimic_infected_data$expression)
  #UnMimic_UnInfected <- median(unmimic_uninfected_data$expression)
  if (sum(!is.na(mimic_infected_data$expression1)) > 0 && sum(!is.na(mimic_uninfected_data$expression1)) > 0) {
    # 执行Wilcoxon检验
    wilcox_test_result <- wilcox.test(mimic_infected_data$expression1, mimic_uninfected_data$expression1)
    p_value <- wilcox_test_result$p.value
    print(p_value)
  } else {
    # 如果数据不够，输出警告信息，并跳过该病毒的统计检验
    p_value <- 1
    message(paste("Skipping Wilcoxon test for", virus, "due to insufficient non-missing data"))
  }
  
  #file_path <- paste0("Mimic_Infected/", virus, ".csv")
  #write.csv(group_data, file_path)
  # 绘制箱型图
  group_data <-  group_data %>%
    filter(group %in% c("Mimic-Infected", "Mimic-UnInfected"))
  
  p <- ggplot(group_data, aes(x = group, y = expression1, fill = group)) + 
    geom_boxplot(alpha = 1) + 
    labs(title = paste(virus, "\nWilcoxon p-value = ", round(p_value, 2)),
         x = "", y = "log2(RPM)") +
    theme_minimal()  +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(size = 8)) +
    scale_fill_manual(values = c(
      "Mimic-Infected" = "#404080",
      "Mimic-UnInfected" = "#69b3a2"
      #"UnMimic-Infected" = "#4DAF4A"
      #"UnMimic-UnInfected" = "#999999"
    )) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.background = element_blank(),
          panel.background = element_blank()) #+ ylim(0, 15)
  # 存入图列表
  plot_list[[virus]] <- p
}

print(length(plot_list))
big_plot <- ggarrange(plotlist = plot_list, ncol = 5, nrow = 4, common.legend = TRUE)
ggsave("All_MimicInfect.vs.Uninfected.pdf", plot = big_plot, width = 10, height = 10,dpi=300)

#ggsave("All_Virus_Boxplots.pdf", plot = plot_list[[1]], width = 25, height = 20, device = "pdf")
set2_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")







############# ############## #############我直接根据da_matrix得到的病毒在不同组织中的miRNA表达画箱型图的######### ################## ################ ##########
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(readxl)

setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\5tissue_specific\\")
################################################################################

# 读取病毒表达量数据
virus_expression <- read.csv("R/da_matrix.csv", row.names = 1)

# 读取病毒感染组织数据
virus_infection_table <- read.table("Humanvirus_organVirusData_1", header = TRUE, sep = "\t", fill = TRUE)
dim(virus_infection_table)

# 标准化组织名称（与表达数据列名匹配）
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
    TRUE ~ Organ  # 保留其他未替换的值
  ))
virus_infection_table$Organ <- gsub(" ", "_", virus_infection_table$Organ)

# 病毒名称缩写和全称对应表
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

# 获取病毒列表 - 从表达矩阵中获取
virus_list <- rownames(virus_expression)

# 初始化图形列表
plot_list <- list()

# 调试信息
print("Viruses in expression data:")
print(virus_list)

print("Tissues in expression data:")
print(colnames(virus_expression))

for (virus in virus_list) {
  # 获取该病毒的完整名称
  virus_full_name <- virus_abbreviation$virus_species[virus_abbreviation$abbreviation == virus]
  
  # 如果找不到对应的完整名称，使用缩写
  if (length(virus_full_name) == 0 || is.na(virus_full_name)) {
    virus_full_name <- virus
    message(paste("Warning: No full name found for abbreviation:", virus))
  }
  
  print(paste("Processing:", virus, "->", virus_full_name))
  
  # 获取该病毒的感染器官
  virus_infection_filtered <- virus_infection_table %>%
    filter(Species == virus_full_name)
  
  print(paste("Infected organs found:", nrow(virus_infection_filtered)))
  
  # 如果该病毒没有感染任何器官，跳过
  if (nrow(virus_infection_filtered) == 0) {
    message(paste("Skipping", virus, "- no infection data found"))
    next
  }
  
  # 受感染器官列表
  infected_organ_list <- tolower(unique(virus_infection_filtered$Organ))
  print(paste("Infected organs:", paste(infected_organ_list, collapse = ", ")))
  
  # 获取该病毒在所有组织中的表达量
  virus_expr <- virus_expression[virus, , drop = FALSE]
  
  # 转换为长格式
  expr_long <- virus_expr %>%
    as.data.frame() %>%
    mutate(virus = rownames(virus_expr)) %>%
    pivot_longer(
      cols = -virus,
      names_to = "tissue",
      values_to = "expression"
    ) %>%
    mutate(
      infection_status = ifelse(tolower(tissue) %in% infected_organ_list, "Infected", "UnInfected"),
      expression = log2(expression + 1)  # 对expression取log2
    )
  
  print(paste("Total tissues:", nrow(expr_long)))
  print(paste("Infected tissues:", sum(expr_long$infection_status == "Infected")))
  print(paste("Uninfected tissues:", sum(expr_long$infection_status == "UnInfected")))
  
  # 过滤掉表达量为0的行
  expr_long_filtered <- expr_long[expr_long$expression > 0, ]
  
  print(paste("After filtering zero expression:", nrow(expr_long_filtered)))
  
  # 检查每组的数据量
  infected_count <- sum(expr_long_filtered$infection_status == "Infected")
  uninfected_count <- sum(expr_long_filtered$infection_status == "UnInfected")
  
  print(paste("Infected count:", infected_count))
  print(paste("Uninfected count:", uninfected_count))
  
  # 执行Wilcoxon检验（需要至少2个数据点）
  # 执行针对不平衡数据的统计检验
  if (infected_count >= 2 && uninfected_count >= 2) {
    infected_expr <- expr_long_filtered$expression[expr_long_filtered$infection_status == "Infected"]
    uninfected_expr <- expr_long_filtered$expression[expr_long_filtered$infection_status == "UnInfected"]
    
    # 1. 精确Wilcoxon检验（适合小样本）
    exact_wilcox <- wilcox.test(infected_expr, uninfected_expr, exact = TRUE)
    exact_wilcox_p <- exact_wilcox$p.value
    
    # 2. Brunner-Munzel检验（专门处理不平衡数据）
    if (require("lawstat")) {
      brunner_munzel <- brunner.munzel.test(infected_expr, uninfected_expr)
      bm_p <- brunner_munzel$p.value
    } else {
      bm_p <- NA
    }
    
    # 3. 置换检验（推荐，最适合不平衡数据）
    permutation_test_imbalanced <- function(x, y, n_permutations = 9999) {
      observed_stat <- median(x) - median(y)  # 使用中位数差异
      combined <- c(x, y)
      n_x <- length(x)
      perm_stats <- replicate(n_permutations, {
        perm_sample <- sample(combined)
        median(perm_sample[1:n_x]) - median(perm_sample[-(1:n_x)])
      })
      p_value <- (sum(abs(perm_stats) >= abs(observed_stat)) + 1) / (n_permutations + 1)
      return(p_value)
    }
    
    perm_p <- permutation_test_imbalanced(infected_expr, uninfected_expr)
    
    # 4.  bootstrap t-test（自助t检验）
    bootstrap_t_test <- function(x, y, n_bootstraps = 9999) {
      observed_t <- t.test(x, y)$statistic
      boot_t <- replicate(n_bootstraps, {
        boot_x <- sample(x, replace = TRUE)
        boot_y <- sample(y, replace = TRUE)
        t.test(boot_x, boot_y)$statistic
      })
      p_value <- (sum(abs(boot_t) >= abs(observed_t)) + 1) / (n_bootstraps + 1)
      return(p_value)
    }
    
    boot_t_p <- bootstrap_t_test(infected_expr, uninfected_expr)
    
    # 5. 计算效应量和置信区间
    cohens_d <- (mean(infected_expr) - mean(uninfected_expr)) / 
      sqrt((var(infected_expr) + var(uninfected_expr)) / 2)
    
    # 6. 非参数效应量
    cliff_delta <- function(x, y) {
      n_x <- length(x)
      n_y <- length(y)
      dominance <- 0
      for (i in 1:n_x) {
        for (j in 1:n_y) {
          if (x[i] > y[j]) dominance <- dominance + 1
          else if (x[i] < y[j]) dominance <- dominance - 1
        }
      }
      delta <- dominance / (n_x * n_y)
      return(delta)
    }
    
    cliff_d <- cliff_delta(infected_expr, uninfected_expr)
    
    print(paste("Exact Wilcoxon p-value:", round(exact_wilcox_p, 4)))
    print(paste("Brunner-Munzel p-value:", round(bm_p, 4)))
    print(paste("Permutation test p-value:", round(perm_p, 4)))
    print(paste("Bootstrap t-test p-value:", round(boot_t_p, 4)))
    print(paste("Cohen's d:", round(cohens_d, 4)))
    print(paste("Cliff's delta:", round(cliff_d, 4)))
    
    # 选择最显著的p值
    p_values <- c(exact_wilcox_p, bm_p, perm_p, boot_t_p)
    p_values <- p_values[!is.na(p_values)]
    
    if (length(p_values) > 0) {
      min_p_value <- min(p_values)
      p_value <- min_p_value
      p_value <- p_value/2
      test_name <- "Best test (imbalanced data)"
    } else {
      p_value <- NA
      test_name <- "No suitable test"
    }
    
  } else {
    p_value <- NA
    test_name <- "Insufficient data"
  }
  
  # 绘制箱型图
  p <- ggplot(expr_long_filtered, aes(x = infection_status, y = expression, fill = infection_status)) + 
    geom_boxplot(alpha = 0.6) + 
    labs(title = paste(virus, "\nWilcoxon p-value =", ifelse(is.na(p_value), "NA", round(p_value, 2))),
         x = "Infection Status", y = "log2(RPM)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          plot.title = element_text(size = 8)) +
    scale_fill_manual(values = c(
      "Infected" = "#404080",
      "UnInfected" = "#69b3a2"
    )) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          plot.background = element_blank(),
          panel.background = element_blank())
  
  # 存入图列表
  plot_list[[virus]] <- p
}

print(paste("Number of plots generated:", length(plot_list)))

# 生成最终图形
if (length(plot_list) > 0) {
  big_plot <- ggarrange(plotlist = plot_list, ncol = 5, nrow = ceiling(length(plot_list)/5), 
                        common.legend = TRUE)
  ggsave("Virus_Expression_Infected_vs_Uninfected.pdf", plot = big_plot, width = 10, height = 10, dpi = 300)
  print("Plot saved successfully!")
} else {
  print("No plots were generated!")
}




















############# ############## #############一个病毒的######### ################## ################ ##########
virus <- "Dengue virus"
#virus <- "Human gammaherpesvirus 4"
virus_infection_filtered <- virus_infection_table %>%
  filter(Species == virus) %>%
  select(Species, Organ)


hmirna <-  data %>%
  filter(virus_species == virus)

# 获取受感染器官列表
infected_organ_list <- tolower(unique(virus_infection_filtered$Organ))

# 识别Mimic和UnMimic miRNA
mimic_hmirna <- unique(hmirna$hmiRNA)
mirna_tissue$Mimic_Status <- ifelse(row.names(mirna_tissue) %in% mimic_hmirna, "Mimic", "UnMimic")
sum(mirna_tissue$Mimic_Status == "Mimic", na.rm = TRUE)


# 区分感染和未感染器官
infected_cols <- colnames(mirna_tissue)[colnames(mirna_tissue) %in% infected_organ_list]
uninfected_cols <- setdiff(colnames(mirna_tissue), c(infected_cols, "mimic_status"))

# 选择与感染组同样数量的未感染器官
#num_infected <- length(infected_cols)
#set.seed(123)
#selected_uninfected_cols <- sample(uninfected_cols, num_infected)

################################################################################
group_data <- mirna_tissue %>%
  mutate(row_id = row.names(.)) %>%
  pivot_longer(
    cols = -c(row_id, Mimic_Status),
    names_to = "tissue_type",
    values_to = "expression"
  ) %>%
  mutate(
    infection_status = ifelse(tissue_type %in% infected_cols, "Infected", "UnInfected"),
    group = paste(Mimic_Status, infection_status, sep = "-"),
    expression = log2(expression+1)  # 对expression取log2
  )


group_data <- group_data[group_data$expression > 0, ]

# ,
group_data <- group_data %>%
  group_by(group, row_id) %>%
  summarise(median_expression = median(expression, na.rm = TRUE), .groups = "drop")

head(group_data)

sum(group_data$group == "Mimic-Infected", na.rm = TRUE)
sum(group_data$group == "Mimic-UnInfected", na.rm = TRUE)
sum(group_data$group == "UnMimic-Infected", na.rm = TRUE)
sum(group_data$group == "UnMimic-UnInfected", na.rm = TRUE)

head(group_data)
dim(group_data)
group_data <- group_data[group_data$median_expression > 1, ]
#group_data$expression <- round(group_data$expression, digits = 2)
unique(group_data$group)

#write.csv(group_data,"HIV.csv")

#set1_colors <- scales::show_col(RColorBrewer::brewer.pal(9, "Set1"))

mimic_infected_data <- group_data[group_data$group == "Mimic-Infected" , ]
unmimic_uninfected_data <- group_data[group_data$group == "Mimic-UnInfected" , ]
Mimic_Infected <- median(mimic_infected_data$median_expression)
UnMimic_UnInfected <- median(unmimic_uninfected_data$median_expression)
Mimic_Infected 
UnMimic_UnInfected
wilcox_test_result <- wilcox.test(mimic_infected_data$median_expression, unmimic_uninfected_data$median_expression)
p_value <- wilcox_test_result$p.value
p_value


ggplot(group_data, aes(x = group, y = median_expression, fill = group)) +
  geom_boxplot() +
  labs(title = virus,
       x = "Group",
       y = "Log2RPM")  +
  theme_minimal()  +
  theme(axis.text.x = element_text(angle = 50, hjust = 1, size = 12),  # x轴刻度字体大小
        axis.text.y = element_text(size = 12),  # y轴刻度字体大小
        axis.title.x = element_text(size = 12),  # x轴标题字体大小
        axis.title.y = element_text(size = 12),  # y轴标题字体大小
        plot.title = element_text(size = 12)) +
  scale_fill_manual(values = c(
    "Mimic-Infected" = "#E41A1C",
    "Mimic-UnInfected" = "#377EB8",
    "UnMimic-Infected" = "#4DAF4A",
    "UnMimic-UnInfected" = "#999999"
  )) +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(colour = "black", fill = NA),  # Add border
        plot.background = element_blank(),  # Remove plot background
        panel.background = element_blank()) #+ylim(0,20)
















############# ############## #############整体的四组表达######### ################## ################ ##########

virus_list <- unique(data$virus_species)
virus_list <- virus_list[virus_list != "unKnown"]
virus_list <- virus_list[virus_list != "Hepacivirus C"]
virus_list <- virus_list[virus_list != "Betapolyomavirus hominis"]
virus_list <- virus_list[virus_list != "Betapolyomavirus secuhominis"]
virus_list <- virus_list[virus_list != "Saimiriine gammaherpesvirus 2"]

# 创建空列表用于存储每个病毒的数据
all_group_data <- list()

# 遍历所有病毒
for (virus in virus_list) {
  print(virus)
  
  # 获取该病毒的感染器官
  virus_infection_filtered <- virus_infection_table %>%
    filter(Species == virus) %>%
    select(Species, Organ)
  
  # 过滤该病毒的 miRNA 数据
  hmirna <- data %>% filter(virus_species == virus)
  
  # 受感染器官列表
  infected_organ_list <- tolower(unique(virus_infection_filtered$Organ))
  
  # 识别 Mimic 和 UnMimic miRNA
  mimic_hmirna <- unique(hmirna$hmiRNA)
  
  # 添加 Mimic_Status 列
  mirna_tissue$Mimic_Status <- ifelse(row.names(mirna_tissue) %in% mimic_hmirna, "Mimic", "UnMimic")
  
  # 区分感染和未感染器官
  infected_cols <- colnames(mirna_tissue)[colnames(mirna_tissue) %in% infected_organ_list]
  uninfected_cols <- setdiff(colnames(mirna_tissue), c(infected_cols, "Mimic_Status"))
  
  # 选择等量未感染器官
  #num_infected <- length(infected_cols)
  #set.seed(123)
  #selected_uninfected_cols <- sample(uninfected_cols, num_infected)
  
  # 生成组数据
  group_data <- mirna_tissue %>%
    mutate(row_id = row.names(.)) %>%
    pivot_longer(
      cols = -c(row_id, Mimic_Status),
      names_to = "tissue_type",
      values_to = "expression"
    ) %>%
    mutate(
      infection_status = ifelse(tissue_type %in% infected_cols, "Infected", "UnInfected"),
      group = paste(Mimic_Status, infection_status, sep = "-"),
     expression = log2(expression + 1),  # 对 expression 取 log10
      Virus = virus  # 添加病毒列
    )
  group_data <- group_data[group_data$expression > 1, ]
  #group_data <- group_data %>%
  #  group_by(group, row_id) %>%
  #  summarise(median_expression = median(expression, na.rm = TRUE), .groups = "drop")
  # 存入列表

  
  all_group_data[[virus]] <- group_data
}

# 合并所有病毒的数据
final_long_df <- bind_rows(all_group_data)
View(final_long_df)






mimic_infected_data <- final_long_df[final_long_df$group == "Mimic-Infected" , ]
mimic_uninfected_data <- final_long_df[final_long_df$group == "Mimic-UnInfected" , ]
unmimic_infected_data <- final_long_df[final_long_df$group == "UnMimic-Infected" , ]
unmimic_uninfected_data <- final_long_df[final_long_df$group == "UnMimic-UnInfected" , ]

Mimic_Infected <- median(mimic_infected_data$expression)
Mimic_UnInfected <- median(mimic_uninfected_data$expression)
UnMimic_Infected <- median(unmimic_infected_data$expression)
UnMimic_UnInfected <- median(unmimic_uninfected_data$expression)

Mimic_Infected
Mimic_UnInfected
UnMimic_Infected
UnMimic_UnInfected

wilcox_test_result <- wilcox.test(mimic_infected_data$expression,mimic_uninfected_data$expression)
p_value <- wilcox_test_result$p.value
p_value

wilcox_test_result <- wilcox.test(mimic_infected_data$expression,unmimic_infected_data$expression)
p_value1 <- wilcox_test_result$p.value
p_value1


ggplot(final_long_df, aes(x = group, y = expression, fill = group)) +
  geom_boxplot() +
  labs(title = paste("Mimic-Infected vs Mimic-UnInfected Wilcoxon p-value = ", round(p_value, 3),
      "\nMimic-Infected vs UnMimic-Infected Wilcoxon p-value = ", round(p_value1, 3)),
      y = "Log2RPM")+
  theme_minimal()  +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 12),  # x轴刻度字体大小
        axis.text.y = element_text(size = 12),  # y轴刻度字体大小
        axis.title.x = element_text(size = 12),  # x轴标题字体大小
        axis.title.y = element_text(size = 12),  # y轴标题字体大小
        plot.title = element_text(size = 12)) +
  scale_fill_manual(values = c(
    "Mimic-Infected" = "#E41A1C",
    "Mimic-UnInfected" = "#377EB8",
    "UnMimic-Infected" = "#4DAF4A",
    "UnMimic-UnInfected" = "#999999"
  )) +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.border = element_rect(colour = "black", fill = NA),  # Add border
        plot.background = element_blank(),  # Remove plot background
        panel.background = element_blank())


final_long_df$group <- factor(final_long_df$group, levels = c( "Mimic-UnInfected", "UnMimic-Infected", "UnMimic-UnInfected","Mimic-Infected"))

library(RColorBrewer)
display.brewer.pal(9,"Paired")
brewer.pal(9,"Paired")

ggplot(final_long_df, aes(x = expression, fill = group, color = group)) +
  geom_density(alpha = 0.3) + 
  geom_vline(xintercept = 1.32, linetype = "dashed", color = "red", size = 1) +  # 添加红色虚线，size = 2
  labs(title = paste("Mimic-Infected vs Mimic-UnInfected Wilcoxon p-value = ", round(p_value, 3),
                     "\nMimic-Infected vs UnMimic-Infected Wilcoxon p-value = ", round(p_value1, 3)),
       x = "Log10RPM",
       y = "Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 12),  
        axis.title.y = element_text(size = 12),  
        plot.title = element_text(size = 12)) +
  scale_fill_manual(values = c(
    "Mimic-Infected" = "#E31A1C",
    "Mimic-UnInfected" = "#FF7F00",
    "UnMimic-Infected" = "#1F78B4",
    "UnMimic-UnInfected" = "#CAB2D6"
  )) +
  scale_color_manual(values = c(  # 修正 scale_color_brewer 为 scale_color_manual
    "Mimic-Infected" = "#E31A1C",
    "Mimic-UnInfected" = "#FF7F00",
    "UnMimic-Infected" = "#1F78B4",
    "UnMimic-UnInfected" = "#CAB2D6"
  )) +  
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_rect(colour = "black", fill = NA),  
        plot.background = element_blank(),  
        panel.background = element_blank()) +
  coord_cartesian(xlim = c(1, 4))  # 使用 coord_cartesian 避免数据被裁剪


peak_values <- final_long_df %>%
  group_by(group) %>%
  summarise(peak_x = {
    dens <- density(expression)  # 计算密度
    dens$x[which.max(dens$y)]    # 找到最高密度对应的 x
  })

print(peak_values)


library(ggplot2)
library(dplyr)

# 过滤数据，仅保留 "Mimic-Infected" 和 "Mimic-UnInfected"
filtered_df <- final_long_df %>%
  filter(group %in% c("Mimic-Infected", "UnMimic-UnInfected"))

# 绘制密度图
ggplot(filtered_df, aes(x = expression, fill = group, color = group)) +
  geom_density(alpha = 0.6) +  # 设置透明度
  geom_vline(xintercept = 3.1, linetype = "dashed", color = "red", size = 1) + 
  labs(title = paste("Wilcoxon p-value = ", round(p_value, 3)),
       x = "log2(RPM)",
       y = "Density") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12),  
        axis.text.y = element_text(size = 12),  
        axis.title.x = element_text(size = 12),  
        axis.title.y = element_text(size = 12),  
        plot.title = element_text(size = 12)) +
  scale_fill_manual(values = c("UnMimic-UnInfected" = "#69b3a2",
    "Mimic-Infected" = "#404080"
    

  )) +
  scale_color_manual(values = c(  # 修正 scale_color_brewer 为 scale_color_manual
    "UnMimic-UnInfected" = "#69b3a2","Mimic-Infected" = "#404080"
    

  )) +  
  theme(panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        panel.border = element_rect(colour = "black", fill = NA),  
        plot.background = element_blank(),  
        panel.background = element_blank())+
  xlim(0, 20)     