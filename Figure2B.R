
setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\1seedregionInmiRNA\\Conservation_R")
library(grid)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
################################################共有的miRNA的靶向基因做富集分析########################################

RNA <- read.csv("..\\HostSpeciesFA\\0jaccard.csv",row.names = 1)
View(RNA)
df_RNA <- as.matrix(RNA)
pdf("0Jaccard4.pdf", width = 6, height = 5)  # 设置PDF尺寸

ht <- Heatmap(df_RNA, 
              name = "Jaccard index",
              na_col = "white",  # 设置NA的颜色
              rect_gp = gpar(col = "white", lwd = 0.5),  # 矩形边框设置
              border = TRUE,  # 显示外框线
              border_gp = gpar(col = "black", lwd = 1),  # 设置外框线颜色和宽度
              col = brewer.pal(n = 9, name = "Reds"),  # 使用颜色Blues Reds
              cluster_columns = FALSE, 
              cluster_rows = FALSE,
              row_order = 1:nrow(df_RNA),
              column_order = 1:ncol(df_RNA),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (!is.na(df_RNA[i, j])) {
                  grid.text(sprintf("%.2f", df_RNA[i, j]), x, y, gp = gpar(fontsize = 10))
                }
              })

draw(ht)
dev.off()





upper_data <- RNA[upper.tri(RNA, diag = FALSE)]
lower_data <- RNA[lower.tri(RNA, diag = FALSE)]
data <- data.frame(
  Jaccard = c(upper_data, lower_data),
  Triangle = c(rep("Virus Mimicked hmiRNA", length(upper_data)),
               rep("All hmiRNA", length(lower_data)))
)
print(data)
print(str(data$Jaccard))
p2 <- ggplot(data, aes(x = Jaccard, fill = Triangle)) +
  geom_density(alpha = 0.6) +
  labs(x = "Jaccard Index", y = "Density") +
  scale_fill_manual(values = c("#4E79A7", "#AAAAAA")) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # 去掉背景网格线
    panel.grid.minor = element_blank(),  # 去掉背景网格线
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(color = "black")  # 调整坐标轴颜色
  ) 
p2




############################改成画条形图############################

RNA <- read.csv("..\\HostSpeciesFA\\0jaccard-bar.csv",header = 1)
View(data_long)

data_long <- RNA %>%
  pivot_longer(
    cols = c(Virus_like_hmiRNA, hmiRNA),
    names_to = "Group",
    values_to = "Value"
  ) %>%
  mutate(Group = factor(Group, levels = c("Virus_like_hmiRNA", "hmiRNA"))) %>%
           mutate(Species = factor(Species, levels = rev(RNA$Species))) 


pdf("0Jaccard4.pdf", width = 6, height = 5) 
ggplot(data_long, aes(x = Value, y = Species, fill = Group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(
    values = c("Virus_like_hmiRNA" = "#4E79A7", "hmiRNA" = "#AAAAAA"),
    labels = c("Virus-like hmiRNA", "hmiRNA")
  ) +
  labs(
    x = "Value",
    y = "Species",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black"),
    legend.position = "top",
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.7)) 

dev.off()












