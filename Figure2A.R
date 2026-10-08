library(ggplot2)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)

# 设置数据
all_hmiRNA <- c(Unique = 320, NonUnique = 206)
virus_hmiRNA <- c(Unique = 203, NonUnique = 144)

# 设置颜色（统一配色方案）
colors <- c("#4E79A7", "#70A5D9") #灰色
colors1 <- c("#666666","#AAAAAA" ) #蓝色


# 设置数据
all_data <- data.frame(
  Category = factor(c("Unique", "NonUnique"), levels = c("NonUnique", "Unique")), # 反转顺序
  Count = c(320, 206),
  Group = "All vmiRNAs (n=526)"
)

virus_data <- data.frame(
  Category = factor(c("Unique", "NonUnique"), levels = c("NonUnique", "Unique")), # 反转顺序
  Count = c(203, 144),
  Group = "human-mimicry vmiRNAs (n=526)"
)

# 设置统一的配色方案
colors_all <- c("Unique" = "#666666", "NonUnique" = "#AAAAAA")  # 灰色系
colors_virus <- c("Unique" = "#4E79A7", "NonUnique" = "#70A5D9")  # 蓝色系

# 绘制条形图 - 修改后
p1 <- ggplot(all_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +  # 在条形上方添加数值标签
  scale_fill_manual(values = colors_all) +
  labs(title = "All vmiRNAs (n=526)", y = "Count", x = "") +
  ylim(0, 350) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.grid = element_blank(),  # 去除所有网格线
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 添加外框线
    axis.line = element_line(color = "black"),  # 添加坐标轴线
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )
p1

p2 <- ggplot(virus_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_fill_manual(values = colors_virus) +
  labs(title = "human-mimicry vmiRNAs (n=347)", y = "", x = "") +
  ylim(0, 250) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.grid = element_blank(),  # 去除所有网格线
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),  # 添加外框线
    axis.line = element_line(color = "black"),  # 添加坐标轴线
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

  
# 合并图形并添加总标题
library(ggpubr)
ggarrange(p1,p2,ncol = 2, nrow = 1,legend="right",align = "h",common.legend = FALSE) 












