####
setwd("D:\\0work\\0wholetransctiptome\\2sRNAminic\\2mismatchnumber\\barplot\\")
library(ggplot2)
library(scales)
library(ggpubr)
library(readxl)
library(ggpubr)
library(rstatix)

data <- read.table("Figure1.txt", sep = '\t', header = TRUE)
View(data)
#data <- data1 %>%
#  filter(Version != "Human")


data$Version<-factor(data$Version,
                     levels = c('Virus'),
                     labels = c("human-mimicry vmiRNAs"))
data$Balitimore = factor(data$Balitimore,levels=c("dsDNA", "ssRNA(+)", "ssRNA(-)", "ssRNA-RT", "ssDNA(-)"))
#因为不是按照mRNA的数量排序的，因此自己手动排序
data$Virus = factor(data$Virus,levels=c("HHV-4","HHV-5","HIV-1","HHV-8","HSV-1","DENV","SINV",
                                        "Herpesvirus","IAV","EVA","ZIKV","HHV-6","HPV","ZEBOV",
                                        "HAdV-C","TTV","SaHV-2","VACV","HMPV","HRSV","MCPyV",
                                        "HCV","WNV","HPyV-1","HPyV-2"))

p1<-ggplot(data,aes(x=Virus,y=Number,fill=Balitimore))+
  geom_bar(stat = 'identity',color = 'black',position = position_dodge(0.9))+ #position使柱子并排放置
  theme_bw(base_size = 18)+ 
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = c("#F4A99B","#019092","#80B1D3","#c55645","#CCD6D1"))+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5)) + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab("Numbers of human-mimicry vmiRNAs")+xlab("") + 
  #coord_cartesian(ylim = c(0, 10000))+
  facet_grid(.~Balitimore, scales="free_x", space = "free") +
  theme(strip.background = element_rect(color = "white", fill = "white"),
        panel.grid = element_blank())
p1


library("RColorBrewer")
display.brewer.pal(8,"Set3")
brewer.pal(8,"Set3")

p3<-ggplot(data,aes(x=Virus, y=nosense,fill=VirusFamily))+
  geom_tile(aes(fill=VirusFamily),linewidth=0)+ #color和size分别指定方块边线的颜色和粗细
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title = element_blank() #不显示图例title
  ) +xlab("")+ylab("")+#"#D9D9D9" ##FFFFB3  #888888 ##FB8072
  scale_fill_manual(values = c("Adenoviridae" = "#66C2A5", 
                               "Anelloviridae"="#FCCDE5",
                               "Filoviridae"="#AAAAAA",
                               "Picornaviridae"="#888888",
                               "Polyomaviridae"="#A6D854",
                               "Flaviviridae" = "#80B1D3",
                               "Herpesviridae" = "#FC8D62",
                               "Orthomyxoviridae"="#FFD92F",
                               "Papillomaviridae"= "#CCD6D1",
                               "Pneumoviridae"="#BEBADA",
                               "Poxviridae" = "#DDDDDD",
                               "Retroviridae" = "#E5C494",
                               "Togaviridae"= "#8DD3C7"
                               )) +
  facet_grid(.~Balitimore, scales="free_x", space = "free") + 
  coord_cartesian(ylim = c(0, 1))+
  theme(strip.background = element_rect(
    color = "white", fill = "white"),
    panel.grid = element_blank())

p3


p2 <- ggplot(data, aes(x = Virus, y = nosense, fill = Percentage)) +
  geom_tile(color = "black", linewidth = 0) + 
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title = element_blank() #不显示图例title
  ) +xlab("")+ylab("") +
  scale_fill_gradient(low = "#FFCCCC", high = "#990000", 
                      name = "Percentage") +
  facet_grid(.~Balitimore, scales="free_x", space = "free") + 
  coord_cartesian(ylim = c(0, 1))+
  theme(strip.background = element_rect(
    color = "white", fill = "white"),
    panel.grid = element_blank())
p2

all_1<-ggarrange(p1,p2,p3,heights=c(0.6,0.2,0.2),ncol = 1, nrow = 3,legend="right",align = "v",font.label = list(size = 18, face = "bold")) 
all_1


#############################################virus and human#############################################

data$Version<-factor(data$Version,
                     levels = c('Virus','Human'),
                     labels = c("human-mimicry vmiRNAs","virus-like hmiRNAs"))
data$Balitimore = factor(data$Balitimore,levels=c("dsDNA", "ssRNA(+)", "ssRNA(-)", "ssRNA-RT", "ssDNA(-)"))
#因为不是按照mRNA的数量排序的，因此自己手动排序
data$Virus = factor(data$Virus,levels=c("HHV-4","HHV-5","HHV-8","HSV-1","Herpesvirus",
                                        "HPV","HHV-6","VACV","SaHV-2","HAdV-C","MCPyV",
                                        "HPyV-1","HPyV-2","HIV-1","EVA","SINV","DENV",
                                        "ZIKV","HCV","WNV","IAV","ZEBOV","TTV","HMPV","HRSV"))

#p1<-ggplot(data,aes(x=reverse(Virus,-Number),y=Number,fill=Version))+
p1<-ggplot(data,aes(x=Virus,y=Number,fill=Version))+
  geom_bar(stat = 'identity',color = 'black',position = position_dodge(0.9))+ #position使柱子并排放置
  theme_bw(base_size = 18)+ 
  theme(axis.text = element_text(colour = 'black'))+
  scale_fill_manual(values = c("#DB3124","#4B74B2"))+
  theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5)) + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  ylab("Numbers of Viral transcripts")+xlab("") + 
  #coord_cartesian(ylim = c(0, 10000))+
  facet_grid(.~Balitimore, scales="free_x", space = "free") +
  theme(strip.background = element_rect(color = "white", fill = "white"),
        panel.grid = element_blank())
p1


p2<-ggplot(data,aes(x=Virus, y=nosense,fill=Balitimore))+
  geom_tile(aes(fill=Balitimore),linewidth=0)+ #color和size分别指定方块边线的颜色和粗细
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(), #不显示坐标轴刻度
        panel.grid.minor = element_blank(), panel.background = element_blank(),
        legend.title = element_blank())  +xlab("")+ylab("")+
  scale_fill_manual(values = c("dsDNA" = "#08519C", 
                               "ssDNA(-)"="#2171B5",
                               "ssRNA-RT" = "#4292C6", 
                               "ssRNA(-)" = "#9ECAE1",
                               "ssRNA(+)" = "#DEEBF7"
                               
  )) +
  facet_grid(.~Balitimore, scales="free_x", space = "free") + 
  coord_cartesian(ylim = c(0, 1))+
  theme(strip.background = element_rect(
    color = "white", fill = "white"),
    panel.grid = element_blank())

p2