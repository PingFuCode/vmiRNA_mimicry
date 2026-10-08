setwd("E:\\")

R.version.string
library(TCGAbiolinks)
library(data.table)
library(readxl)
library(data.table)
library(GenomicRanges)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(RIdeogram)
library(ChIPseeker)
library(AnnotationDbi)
library(clusterProfiler)
library(GenomeInfoDb)
############################################################
## 3 把这些区间里的蛋白编码基因捞出来
############################################################
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#BiocManager::install("ChIPseeker")

cfs.dir   <- "E:/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg"
## 0.2 缺失片段（你已经做好）
delFrag.GR <- seg.all[Segment_Mean < -0.3,   # 关键修改：用数值阈值代替 isDel
                      .(chr   = paste0("chr", Chromosome),
                        start = as.numeric(Start),
                        end   = as.numeric(End),
                        Sample,Segment_Mean)]
dim(delFrag.GR)

head(delFrag.GR)

delFrag.GR <- as(delFrag.GR, "GRanges")

annotated_peaks <- annotatePeak(
  peak = delFrag.GR,        # 确保这是GRanges对象
  tssRegion = c(0, 0),
  TxDb = txdb,              # 确保这是有效的TxDb对象
 annoDb = "org.Hs.eg.db"   # 确保已安装并加载org.Hs.eg.db
)

colnames(as.data.frame(annotated_peaks))

gene.df <- as.data.table(annotated_peaks)[
  , .(gene = SYMBOL, chr = seqnames, start, end)]

View(gene.df)



# 1. 读取HHV-4靶基因列表
target_genes <- unique(unlist(strsplit(
  fread(r"(E:\2sRNAminic\0006mer_target_with_hsaID.csv)")
  [abbreviation!=" "]$Intersection_Genes, ";")))
length(target_genes)
# 2. 转换基因符号为ENTREZID
gene.info <- bitr(target_genes,
                  fromType = "SYMBOL",
                  toType   = c("ENTREZID", "MAP"),
                  OrgDb    = org.Hs.eg.db)

# 将gene.info转换为data.table
gene.info <- as.data.table(gene.info)

# 将ENTREZID转换为字符型
gene.info[, ENTREZID := as.character(ENTREZID)]

## 3 用 TxDb 把 ENTREZID → 基因坐标（start/end 用基因全长）
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# 提取所有基因范围
# 提取所有基因范围
g_rng_list <- genes(txdb, columns = "GENEID", single.strand.genes.only = FALSE)

# 清理染色体名称
main_chroms <- paste0("chr", c(1:22, "X", "Y", "M"))
g_rng_list <- keepSeqlevels(g_rng_list, main_chroms, pruning.mode = "coarse")

# 获取基因ID
gene_ids <- mcols(g_rng_list)$GENEID

# 获取每个基因的范围数量
range_counts <- elementNROWS(g_rng_list)

# 展开GRangesList
g_rng_flat <- unlist(g_rng_list)

# 创建展开后的gene_id向量
expanded_gene_ids <- rep(gene_ids, range_counts)

# 将gene_id添加到元数据
mcols(g_rng_flat)$GENEID <- expanded_gene_ids

# 现在转换为data.table
g.df <- as.data.table(g_rng_flat)[, .(
  chr = seqnames,
  start = start,
  end = end,
  ENTREZID = as.character(GENEID)
)]

# 3. 合并基因信息与坐标信息
target_coord <- merge(
  gene.info[, .(SYMBOL, ENTREZID)],  # 只保留需要的列
  g.df[, .(ENTREZID, chr, start, end)],
  by = "ENTREZID",
  all.x = TRUE  # 保留所有靶基因，即使没有坐标信息
)

# 过滤掉没有坐标信息的基因
target_coord <- target_coord[!is.na(chr)]

# 4. 去重复（一个基因多个坐标记录时取并集）
target_coord <- target_coord[
  , .(start = min(start),
      end = max(end),
      ENTREZID = first(ENTREZID)),
  by = .(gene = SYMBOL, chr)
]

# 调整列顺序，移除不需要的ENTREZID列（可选）
target_coord <- target_coord[, .(gene, chr, start, end)]

# 5. 保存结果
fwrite(target_coord, "E:/HHV4_target_genes_hg38_coord.txt", sep = "\t")
cat(sprintf("已获取 %d 个 HHV-4 靶基因的 hg38 坐标\n", nrow(target_coord)))

# 可选：显示统计信息
cat(sprintf("原始靶基因数: %d\n", length(target_genes)))
cat(sprintf("成功转换的基因数: %d\n", nrow(gene.info)))
cat(sprintf("成功获取坐标的基因数: %d\n", nrow(target_coord)))

# 检查哪些基因没有成功获取坐标
#if (nrow(gene.info) > nrow(target_coord)) {
#  missing_genes <- setdiff(gene.info$SYMBOL, target_coord$gene)
 # cat(sprintf("有 %d 个基因无法获取坐标信息\n", length(missing_genes)))
  
 # # 保存这些基因以供参考
 # fwrite(data.table(gene = missing_genes), 
 #        "E:/HHV4_missing_coordinates.txt", 
 #        sep = "\t")}


library(GenomicRanges)

# 首先确保 target_coord 已经正确生成
head(target_coord)
dim(target_coord)
head(gene.df)
# 1. 构建两个 GRanges 对象
gene_gr <- makeGRangesFromDataFrame(gene.df, keep.extra.columns = TRUE)
target_gr <- makeGRangesFromDataFrame(target_coord, keep.extra.columns = TRUE)

# 2. 找出重叠的区域
# 方法1: 使用 findOverlaps 获取详细的对应关系
overlaps <- findOverlaps(gene_gr, target_gr)

#seqlevels(gene_gr)
#seqlevels(target_gr)

# 3. 从重叠中提取信息
# 获取重叠的 gene_gr 中的区域
common_gene_regions <- gene_gr[queryHits(overlaps)]
# 获取重叠的 target_gr 中的区域
common_target_regions <- target_gr[subjectHits(overlaps)]

# 5. 统计信息
cat(sprintf("gene.df 中的区域数: %d\n", length(gene_gr)))
cat(sprintf("target_coord 中的靶基因数: %d\n", length(target_gr)))
cat(sprintf("发现重叠区域数: %d\n", length(overlaps)))
# 7. 保存结果
# 保存所有重叠区域
dt_out <- as.data.table(common_target_regions)
# 如果还想保留基因名等元数据
dt_out[, gene := mcols(common_target_regions)$gene]

fwrite(dt_out, "E:/overlapping_HHV4_target_regions.txt", sep = "\t")




library(GenomicRanges)
library(data.table)
library(karyoploteR)

############################################################
## 1 把 GRanges → 核型条形数据
############################################################
# 1.1 处理脆弱区段数据
#cfs.df <- as.data.table(cfs.gr)
#bar_cfs <- cfs.df[, .(
#  chr   = as.character(seqnames),  # 保留原始染色体名
#  start = start,
#  end   = end
#)]
#View(bar_cfs)
# 1.2 处理缺失区段数据delFrag.GR
delFrag.dt <- as.data.table(delFrag.GR)
bar_del <- delFrag.dt[, .(
  chr   = as.character(seqnames),  # 保留原始染色体名
  start = start,
  end   = end
)]

head(common_target_regions)

#head(delFrag.GR)

############################################################
## 2 绘图
############################################################
# 1. 首先读取Excel文件中的基因列表
library(readxl)

# 读取GC_TSGene sheet的GeneSymbol列
cancer_genes <- read_excel("E:/2sRNAminic/CancerGene.xlsx", sheet = "TSGene")
cancer_gene_list <- cancer_genes$GeneSymbol
length(cancer_gene_list)
# 2. 筛选common_target_regions中属于癌症基因的条目
# 使用 %in% 操作符筛选
cancer_target_regions <- common_target_regions[common_target_regions$gene %in% cancer_gene_list, ]
# 查看匹配到的癌症基因
cat("找到", length(cancer_target_regions), "个抑癌癌症相关基因在共同靶区段中\n")
#print(cancer_target_regions$gene)


library(karyoploteR)   # 画染色体
library(viridisLite)   # 蓝色渐变
# 1.2 目标基因
library(readxl)
## ========== 1 读数据 ==========
# 1.1 缺失片段
cfs.dir      <- fread("E:/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.seg")
target_genes <- fread("E:/HHV4_target_genes_hg38_coord.txt")
length(target_genes$gene)
# 读取CancerGene.xlsx文件的CancerGeneCensus工作表
cancer_genes <- read_excel("E:/2sRNAminic/CancerGene.xlsx", sheet = "TSGene") #TSGene

# 获取去重的GENE_SYMBOL
genes_to_remove <- unique(cancer_genes$GeneSymbol)
length(genes_to_remove)
# 从target_genes中去除这些基因
# 假设target_genes的gene列对应的是基因符号,一共67个靶基因是抑癌基因
target_genes_label_TS <- target_genes[gene %in% genes_to_remove]
head(target_genes_label_TS)
target_genes_label <- target_genes
dim(target_genes)
dim(target_genes_label)
dim(target_genes_label_TS)

## ========== 2 统一chr格式 ==========
clean_chr <- function(x) {
  x <- sub("^chr","",x)
  x[x=="X"] <- "23"
  x[x=="Y"] <- "24"
  x
}

#View(target_genes_label_TS)
#target_genes$chr <- clean_chr(target_genes$chr)
target_genes_label_TS$chr <- clean_chr(target_genes_label_TS$chr)
target_genes_label$chr <- clean_chr(target_genes_label$chr)
cfs.dir$Chromosome <- clean_chr(cfs.dir$Chromosome)
cfs.dir$Segment_Mean <- as.numeric(cfs.dir$Segment_Mean)
dim(cfs.dir)
## ========== 3 建GR对象 ==========
# 3.1 缺失片段（Segment_Mean < -0.3）
#keepID <- c("TCGA-96-A4JL","TCGA-AA-3488","TCGA-AA-3549","TCGA-AA-3860",
#            "TCGA-AA-3989","TCGA-AG-3593","TCGA-B7-5818","TCGA-BA-7269",
 #           "TCGA-BR-4253","TCGA-BR-6706","TCGA-BR-6707","TCGA-BR-7196",
  #          "TCGA-BR-7958","TCGA-BR-8366","TCGA-BR-8381","TCGA-BR-8589",
  #          "TCGA-BR-8676","TCGA-BR-8686","TCGA-BR-A4J4","TCGA-CD-5801",
  #          "TCGA-CG-4304","TCGA-CG-5722","TCGA-CG-5725","TCGA-D1-A1O8",
  #          "TCGA-D7-5577","TCGA-D7-8570","TCGA-D7-8573","TCGA-D7-A4YX",
   #         "TCGA-EI-6511","TCGA-F7-A61S","TCGA-FI-A2CX","TCGA-FP-7916",
   #         "TCGA-FP-7998","TCGA-HU-8608","TCGA-HU-A4G2","TCGA-HU-A4G6",
  #          "TCGA-HU-A4GF","TCGA-HU-A4H0","TCGA-K7-A5RG")

#pat <- paste0(keepID, collapse = "|")   # 正则：TCGA-96-A4JL|TCGA-AA-3488|...

#tmp <- cfs.dir[grepl(pat, Sample)]
#dim(tmp)
library(dplyr)
#length(unique(cfs.dir$Sample))

###如果筛选在5%的病人中都出现，那么这些区域只集中在9、10 13等部分染色体上
#cfs.tmp <- cfs.dir %>%
  # 首先筛选符合条件的行
#  filter(
 #   Segment_Mean < -0.3,
#    as.numeric(substr(Sample, 14, 15)) == 01,  # 癌组织样本
 #   substr(Sample, 16, 16) == "A"  # 冰冻组织
#  ) %>%
  # 计算每个区域在多少样本中出现
 # group_by(Chromosome, Start, End) %>%
  # 添加样本计数
  #mutate(sample_count = n_distinct(Sample)) %>%
  # 筛选出在多个样本中出现的区域
 # filter(sample_count >= 1093) %>%
  # 移除辅助列
 # select(-sample_count) %>%
 # ungroup() %>%
  # 排序
 # arrange(Chromosome, Start)


# 1. 应用严格阈值
cfs.strict <- cfs.dir %>%
  filter(
    abs(as.numeric(Segment_Mean)) > 0.3,  # 强信号
    as.numeric(substr(Sample, 14, 15)) == 01,
    substr(Sample, 16, 16) == "A"
  )


setDT(cfs.strict)
cfs.strict[, region_length := End - Start]

# 3. 快速筛选：先去掉极长和极短的区域
#    保留1000bp到10Mb之间的区域
cfs.tmp <- cfs.strict[region_length >= 1000 & region_length <= 1000000]

dim(cfs.tmp)


####fragile区域：
library(dplyr)
fragile_dir <- "E:/2sRNAminic/fragile_site_bed/fragile_site_bed/"
bed_files <- list.files(path = fragile_dir, 
                        pattern = "chr.*_fragile_site\\.bed$",
                        full.names = TRUE)
# 使用lapply读取所有文件
fragile_list <- lapply(bed_files, function(bed_file) {
  df <- read.table(bed_file, header = FALSE, sep = "\t", 
                   stringsAsFactors = FALSE)
  colnames(df) <- c("chr", "start", "end", "name", "score", "strand")
  df$chr <- gsub("chr", "", df$chr)
  return(df)})
fragile <- bind_rows(fragile_list)
View(fragile)

# 处理染色体列，去掉"chr"前缀
fragile$chr <- gsub("chr", "", fragile$chr)
# 将X/Y染色体转换为数字表示（X=23, Y=24）
fragile$chr <- ifelse(fragile$chr == "x", 23, 
                      ifelse(fragile$chr == "Y", 24, fragile$chr))

# 将染色体列转换为整数
fragile$chr <- as.integer(fragile$chr)
head(fragile)
#View(fragile)
####染色体缺失和脆弱区域取交集
cfs_gr <- GRanges(
  seqnames = paste0("chr", cfs.tmp$Chromosome),
  ranges = IRanges(start = cfs.tmp$Start, end = cfs.tmp$End),
  Sample = cfs.tmp$Sample,
  Segment_Mean = cfs.tmp$Segment_Mean
)

# 创建脆弱位点的GenomicRanges对象
fragile_gr <- GRanges(
  seqnames = paste0("chr", fragile$chr),
  ranges = IRanges(start = fragile$start, end = fragile$end),
  name = fragile$name
)

seqlevels(cfs_gr)
seqlevels(fragile_gr)
seqlevels(cfs_gr)[seqlevels(cfs_gr) == "chrX"] <- "chr23"

# 查找重叠区间
overlaps <- findOverlaps(cfs_gr, fragile_gr)
# 获取重叠的结果
cfs.filt <- cfs.tmp[queryHits(overlaps), ]
fragile_overlap <- fragile[subjectHits(overlaps), ]
#合并结果，添加脆弱位点信息
cfs.filt[, fragile_name := fragile_overlap$name]
cfs.filt[, fragile_start := fragile_overlap$start]
cfs.filt[, fragile_end := fragile_overlap$end]

# 查看结果
head(cfs.filt)
dim(cfs.filt)
#unique(cfs.filt$Chromosome)

## 2. 建 GRanges（带 Segment_Mean）
delGR.raw <- makeGRangesFromDataFrame(
  cfs.filt,
  seqnames.field = "Chromosome",
  start.field    = "Start",
  end.field      = "End",
  keep.extra.columns = TRUE
)

delGR <- reduce(delGR.raw, min.gapwidth = 1000) 
unique(delGR$seqnames)



## 3. 用重叠原始片段给新区段补回 Segment_Mean（取均值）
olap <- findOverlaps(delGR.raw, delGR)          # query=原始, subject=合并后
segMean <- tapply(mcols(delGR.raw)$Segment_Mean[queryHits(olap)],
                  subjectHits(olap),
                  mean, na.rm = TRUE)

## 4. 把算好的值写回 mcols
mcols(delGR) <- DataFrame(Segment_Mean = segMean)

## 现在检查
head(delGR)
mcols(delGR)

# 3.2 目标基因
geneGR <- makeGRangesFromDataFrame(target_genes,
                                   seqnames.field="chr",
                                   start.field="start",
                                   end.field="end",
                                   keep.extra.columns=TRUE)
geneGR_label <- makeGRangesFromDataFrame(target_genes_label,
                                   seqnames.field="chr",
                                   start.field="start",
                                   end.field="end",
                                   keep.extra.columns=TRUE)

geneGR_label_TS <- makeGRangesFromDataFrame(target_genes_label_TS,
                                         seqnames.field="chr",
                                         start.field="start",
                                         end.field="end",
                                         keep.extra.columns=TRUE)

head(geneGR)
head(geneGR_label)
head(geneGR_label_TS)
length(unique(geneGR_label_TS$gene))
## ========== 4 交集：缺失区段里的基因 ==========

##TS

ol_TS <- findOverlaps(geneGR_label_TS,delGR)  #缺失区段里的基因
hitGR_TS <- geneGR_label_TS[queryHits(ol_TS)]          # 真正落在缺失区段的基因
hitGR_TS$gene <- target_genes_label_TS$gene[queryHits(ol_TS)]
length(hitGR_TS$gene)
head(hitGR_TS)

ol <- geneGR  #全部基因
hitGR <-geneGR    #全部基因
hitGR$gene <- target_genes$gene
hitGR <- hitGR[!duplicated(hitGR$gene)]
length(hitGR$gene)
seqlevels(hitGR) <- sub("^chr", "", seqlevels(hitGR))
seqlevels(hitGR)[seqlevels(hitGR) == "X"] <- "23"


ol_label <- geneGR_label
hitGR_label <- ol_label  # 或直接用 geneGR_label
hitGR_label$gene <- target_genes_label$gene  # 直接赋值，无需 queryHits
hitGR_label <- hitGR_label[!duplicated(hitGR_label$gene)]
length(hitGR_label$gene)

target_genes$chr <- as.character(target_genes$chr)
target_genes_label$chr <- as.character(target_genes_label$chr)
target_genes$chr <- clean_chr(target_genes$chr)
target_genes_label$chr <- clean_chr(target_genes_label$chr)
cfs.dir$Chromosome <- clean_chr(cfs.dir$Chromosome)

fragile$chr <- clean_chr(fragile$chr)
seqlevels(delGR)[seqlevels(delGR) == "X"] <- "23"

unique(seqnames(delGR))
unique(seqnames(hitGR))
unique(seqnames(custom_genome))
target_genes$chr[target_genes$chr == "X"] <- "23"
unique(target_genes$chr) 


## ========== 5 画染色体骨架 ==========
custom_genome <- GRanges(
  seqnames = c(as.character(1:22), "23", "24"), 
  ranges = IRanges(start = 1, width = c(
    248956422, 242193529, 198295559, 190214555, 181538259,
    170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189,
    90338345, 83257441, 80373285, 58617616, 64444167,
    46709983, 50818468, 156040895, 57227415
  ))
)

# 设置绘图参数
params <- getDefaultPlotParams(plot.type = 1)
params$ideogram$fill <- "transparent"
params$ideogram$col <- "black"

# 绘图
kp <- plotKaryotype(
  genome = custom_genome,
  plot.type = 1,
  plot.params = params
)


## ========== 6 画缺失区段（蓝色渐变） ==========
# 把Segment_Mean映射到0-1，再用viridis蓝
vals  <- as.numeric(mcols(delGR)$Segment_Mean)

blues <- colorRampPalette(c("lightblue", "darkblue"))(100)
cols <- blues[cut(vals, breaks = 100, labels = FALSE)]

kpPlotRegions(kp, data = delGR, r0 = 0, r1 = 0.2,
              col = cols)

## ========== 7 用红色条带标出“所有靶基” ==========

#if(length(hitGR)){
#  kpPlotRegions(kp, data = hitGR, r0 = 0.2, r1 = 0.5,
#                col = "#E41A1C")
#}
if(length(hitGR)){
  kpPlotRegions(kp, data = hitGR, 
                r0 = 0.3, r1 = 0.5,
                col = "black", 
                border = "black",
                avoid.overlapping = FALSE)  # 关闭重叠避免
}

if(length(hitGR_TS)){
  kpPlotRegions(kp, data = hitGR_TS, 
                r0 = 0.3, r1 = 0.8,
                col = "#E41A1C", 
                border = "#E41A1C",
                avoid.overlapping = FALSE)  # 关闭重叠避免
}

write.csv(hitGR,"SharedTargetGene.csv")
write.csv(hitGR_TS,"SharedTargetGene_TS.csv")
## ========== 8 基因标签-既是TS由位于deletion区域 ==========


if(length(hitGR_label)){
  ## 1 计算中点并新建 GRanges
  mid_pos <- (start(hitGR_label) + end(hitGR_label)) %/% 2
  midGR   <- GRanges(seqnames = seqnames(hitGR_label),
                     ranges   = IRanges(start = mid_pos, width = 1),
                     gene     = hitGR_label$gene)   # 把基因名带过去
  
  ## 2 画竖线（连接线）
  #kpPlotMarkers(kp, data = midGR, r0 = 0.5, r1 = 0.8,
  #              labels = NA, line.col = "black", cex = 0.3)
  
  ## 3 写基因名
  kpPlotMarkers(kp, data = midGR, r0 = 0.5, r1 = 0.6,
                labels = midGR$gene,
                label.color = "black",
                cex = 0.5,
                text.orientation = "horizontal",
                avoid.overlapping = TRUE)
}

as.data.frame(hitGR_label) 


## ========== 9 图例 ==========
legend("bottomright",
       legend = c("Genome-wide copy-number deletion & fragile region\n","Shared Target Genes"),
       fill = c("#4682B4", "#E41A1C"), border = NA, bty = "n", cex = 0.8)



write.csv(as.data.frame(hitGR_label),
          file = "EBVTargetGene_TS.csv",
          row.names = FALSE)
##########################统计靶基因出现的概率###########################
write.csv(as.data.frame(hitGR),
          file = "EBVTargetGene.csv",
          row.names = FALSE)
# 创建自定义基因组范围，包含所有染色体
custom_genome <- GRanges(
  seqnames = paste0("", c(1:22, "X", "Y")),
  ranges = IRanges(start = 1, width = c(
    248956422, 242193529, 198295559, 190214555, 181538259,
    170805979, 159345973, 145138636, 138394717, 133797422,
    135086622, 133275309, 114364328, 107043718, 101991189,
    90338345, 83257441, 80373285, 58617616, 64444167,
    46709983, 50818468, 156040895, 57227415
  ))
)

# 将每条染色体分成四等分
create_quartiles <- function(gr) {
  quartile_list <- list()
  
  for(i in 1:length(gr)) {
    chr <- as.character(seqnames(gr[i]))
    chr_length <- width(gr[i])
    
    # 计算四个等分的边界
    q1_end <- floor(chr_length * 0.25)
    q2_end <- floor(chr_length * 0.50)
    q3_end <- floor(chr_length * 0.75)
    q4_end <- chr_length
    
    # 创建四个等分的GRanges
    q1 <- GRanges(seqnames = chr, ranges = IRanges(start = 1, end = q1_end))
    q2 <- GRanges(seqnames = chr, ranges = IRanges(start = q1_end + 1, end = q2_end))
    q3 <- GRanges(seqnames = chr, ranges = IRanges(start = q2_end + 1, end = q3_end))
    q4 <- GRanges(seqnames = chr, ranges = IRanges(start = q3_end + 1, end = q4_end))
    
    # 添加等分标签
    mcols(q1)$quartile <- "Q1"
    mcols(q2)$quartile <- "Q2"
    mcols(q3)$quartile <- "Q3"
    mcols(q4)$quartile <- "Q4"
    
    quartile_list <- c(quartile_list, list(q1, q2, q3, q4))
  }
  
  # 合并所有等分
  return(do.call(c, quartile_list))
}

# 创建染色体四等分
quartile_regions <- create_quartiles(custom_genome)

# 统计靶基因在各个等分的分布
count_hits_in_quartiles <- function(hitGR, quartile_regions) {
  # 找到重叠区域
  overlaps <- findOverlaps(hitGR, quartile_regions)
  
  # 获取重叠的等分信息
  hit_quartiles <- quartile_regions[subjectHits(overlaps)]$quartile
  
  # 统计每个等分的数量
  quartile_counts <- table(hit_quartiles)
  
  # 转换为数据框
  result_df <- data.frame(
    Quartile = names(quartile_counts),
    Count = as.numeric(quartile_counts),
    Proportion = as.numeric(quartile_counts) / length(hitGR)
  )
  
  # 按等分排序
  result_df <- result_df[order(result_df$Quartile), ]
  
  return(result_df)
}

# 执行统计
quartile_stats <- count_hits_in_quartiles(hitGR, quartile_regions)

# 查看结果
print(quartile_stats)

# 可视化结果
library(ggplot2)

ggplot(quartile_stats,
       aes(x = Quartile, y = Proportion, fill = Quartile)) +
  geom_bar(stat = "identity", colour = "black", size = 1) +
  geom_text(aes(label = paste0(round(Proportion * 100, 1), "%")),
            vjust = -0.5, size = 4) +
  labs(title = "",
       x = "Chromosomal Quartile",
       y = "Proportion of EBV Shared Target Genes") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),   # 去掉网格
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 外框线
    axis.line    = element_line(colour = "black", size = 1)              # 坐标轴线
  ) +
  scale_fill_manual(values = rep("steelblue", 4)) +
  ylim(0, max(quartile_stats$Proportion) * 1.2)




########################三者的交集#####################################
# 如缺少包先安装
# install.packages(c("readxl","data.table","ggVennDiagram"))

library(readxl)   
library(data.table)  
library(ggVennDiagram)
library(dplyr)

## 1. 抑癌基因 —— Excel 工作表 TSGene
ts_genes <- read_excel("E:/2sRNAminic/CancerGene.xlsx", sheet = "TSGene") %>%
  pull(GeneSymbol) %>%              # 把“Gene”列抽出来；列名按实际改
  unique()

## 2. 癌基因 —— Excel 工作表 CancerGeneCensus
cancer_genes <- read_excel("E:/2sRNAminic/CancerGene.xlsx", sheet = "CancerGeneCensus") %>%
  pull(GeneSymbol) %>%              # 列名按实际改
  unique()

## 3. 靶基因 —— csv 文件，列 Intersection_Genes 用分号分隔
target_genes <- fread(r"(E:\2sRNAminic\0006mer_target_with_hsaID.csv)")[
  , unlist(strsplit(Intersection_Genes, ";"))
] %>% unique()
length(target_genes)
## 4. 整理成 list（ggVennDiagram 的输入格式）
venn_list <- list(
  TS_genes    = ts_genes,
  Oncogenes   = cancer_genes,
  Target_genes = target_genes
)

library(ggplot2)
## 5. 画三集合韦恩图
ggVennDiagram(venn_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(legend.position = "none")


# 1. 提取癌基因与靶基因的交集
oncogene_target_intersect <- intersect(cancer_genes, target_genes)

# 2. 转换为数据框并写入 CSV 文件
write.csv(data.frame(Gene = oncogene_target_intersect), 
          "E:/2sRNAminic/oncogene_target_intersect.csv", 
          row.names = FALSE)

# 可选：提取抑癌基因与靶基因的交集
tsgene_target_intersect <- intersect(ts_genes, target_genes)
write.csv(data.frame(Gene = tsgene_target_intersect), 
          "E:/2sRNAminic/tsgene_target_intersect.csv", 
          row.names = FALSE)

# 可选：提取三者的共同交集
all_three_intersect <- Reduce(intersect, list(ts_genes, cancer_genes, target_genes))
write.csv(data.frame(Gene = all_three_intersect), 
          "E:/2sRNAminic/all_three_intersect.csv", 
          row.names = FALSE)



















contingency_table <- matrix(
  c(69, 1586, 694, 18957), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    c("Common Target Genes", "Non-common Target Genes"),
    c("CancerGeneCensus", "Non-CancerGeneCensus")
  )
)

length(unique(target_genes$gene))

seqlevels(delGR) <- sub("^chr", "", seqlevels(delGR))
head(delGR)
#seqlevels(delGR) <- paste0("chr", seqlevels(delGR))
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- genes(txdb)


library(org.Hs.eg.db)

# 获取所有重叠基因的Entrez ID
overlaps <- findOverlaps(delGR, genes)
gene_ids <- names(genes)[subjectHits(overlaps)]
unique_gene_ids <- unique(gene_ids)

# 获取基因类型
gene_types <- select(org.Hs.eg.db, 
                     keys = unique_gene_ids,
                     columns = c("SYMBOL", "GENETYPE"),
                     keytype = "ENTREZID")

# 筛选蛋白质编码基因
protein_coding_ids <- gene_types$ENTREZID[gene_types$GENETYPE == "protein-coding"]
protein_coding_genes <- genes[names(genes) %in% protein_coding_ids]

# 重新计算与蛋白质编码基因的重叠
protein_coding_counts <- countOverlaps(delGR, protein_coding_genes)
pc_overlaps <- findOverlaps(delGR, protein_coding_genes)
pc_unique_genes <- protein_coding_genes[unique(subjectHits(pc_overlaps))]

cat("deletion区域中的蛋白质编码基因数量（去重）:", length(pc_unique_genes), "\n")

##547个在deletion中 1109-547
#21306-7316
#与缺失区段重叠的基因总数: 10540 

7316-505
21306-6811
contingency_table <- matrix(
  c(505,1150, 6811, 12840), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    c("Shared Target Genes", "Non-Target Genes"),
    c("Gene in deletion and fragile regions", "Gene not in deletion and fragile regions")
  )
)
contingency_table <- matrix(
  c(264,498, 7052, 13492), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    c("Shared Target Genes", "Non-Target Genes"),
    c("Gene in deletion and fragile regions", "Gene not in deletion and fragile regions")
  )
)

print("列联表：")
print(contingency_table)
#View(contingency_table)
# Fisher精确检验
fisher_result <- fisher.test(contingency_table)

print("Fisher精确检验结果：")
print(fisher_result)

# 计算比值比(OR)和95%置信区间
or_value <- fisher_result$estimate
ci_lower <- fisher_result$conf.int[1]
ci_upper <- fisher_result$conf.int[2]
p_value <- fisher_result$p.value

cat("\n结果摘要：\n")
cat("比值比(OR):", round(or_value, 3), "\n")
cat("95%置信区间:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n")
cat("p值:", format(p_value, scientific = TRUE), "\n")


mat <- matrix(c(264, 499,
               7052, 13493),
              nrow = 2, byrow = TRUE,
              dimnames = list(GeneSet = c("EBV", "non-EBV"),
                              Region  = c("InDel", "NotInDel")))

## Fisher 精确检验
fisher.test(mat)









########################################fisher
contingency_table <- matrix(
  c(52, 1050, 711, 19493), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    c("HHV-4 Target Genes", "Non-HHV-4 Target Genes"),
    c("CancerGeneCensus", "Non-CancerGeneCensus")
  )
)

contingency_table <- matrix(
  c(34, 1068, 454, 19750), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    c("HHV-4 Target Genes", "Non-HHV-4 Target Genes"),
    c("CancerGeneCensus", "Non-CancerGeneCensus")
  )
)


contingency_table <- matrix(
  c( 28  , 48, 476, 1102), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    c("HHV-4 Target Genes", "Non-HHV-4 Target Genes"),
    c("CancerGeneCensus", "Non-CancerGeneCensus")
  )
)
contingency_table <- matrix(
  c( 19  , 36, 486, 1114), 
  nrow = 2, 
  byrow = TRUE,
  dimnames = list(
    c("HHV-4 Target Genes", "Non-HHV-4 Target Genes"),
    c("CancerGeneCensus", "Non-CancerGeneCensus")
  )
)
print("列联表：")
print(contingency_table)

# Fisher精确检验
fisher_result <- fisher.test(contingency_table)

print("Fisher精确检验结果：")
print(fisher_result)

# 计算比值比(OR)和95%置信区间
or_value <- fisher_result$estimate
ci_lower <- fisher_result$conf.int[1]
ci_upper <- fisher_result$conf.int[2]
p_value <- fisher_result$p.value

cat("\n结果摘要：\n")
cat("比值比(OR):", round(or_value, 3), "\n")
cat("95%置信区间:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n")
cat("p值:", format(p_value, scientific = TRUE), "\n")





















