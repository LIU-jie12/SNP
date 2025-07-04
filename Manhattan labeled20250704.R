library(ieugwasr)
library(VariantAnnotation)
library(remotes)
library(gwasglue)
library(data.table)
library(dplyr)
library(tidyr)
library("plinkbinr")
library(MendelianRandomization)
library("reshape")
library("plyr")
library("CMplot")
library("TwoSampleMR")
library("FastGWASR")
library(LDlinkR)
vignette("qqman")
library(FastTraitR)
library(biomaRt)
library(ggplot2)

library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)  #尝试加载这个特定的包
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(MungeSumstats)
library(tidyverse)


setwd("F:/GCST")

data1 <- fread("1GCST90018808.tsv") %>% 
  dplyr::mutate(ranges = str_c(chromosome, ":", base_pair_location))
# 加载 SNP 数据库
snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37

# 根据 chromosome 和 base_pair_location 创建 GRanges 对象
query_ranges <- GRanges(seqnames = data1$chromosome,
                        ranges = IRanges(start = data1$base_pair_location, 
                                         end = data1$base_pair_location))

# 使用精确位置检索 SNP 编号
snp.res <- snpsByOverlaps(snps, query_ranges)
snp_info <- as.data.frame(snp.res) %>%
  dplyr::mutate(ranges = str_c(seqnames, ":", pos)) %>%
  dplyr::select(ranges, rsid = RefSNP_id)  # 仅保留 ranges 和 rsid


# 使用 left_join 合并，确保保留所有原始数据
merged_data <- data1 %>%
  left_join(snp_info, by = "ranges")

write.table(merged_data, file = "1GCST90018808-snp.csv", sep = ",", row.names = FALSE)


GWAS_data <- read.table("1GCST90018808-snp.csv", header = TRUE, sep = ",")

two_sample_MR_data <- GWAS_data[, c("rsid", "chromosome","base_pair_location", "p_value")]
colnames(two_sample_MR_data) <- c("SNP", "CHR", "BP", "pvalue")


# 要标注的 SNP 
snp_to <- c("rs12022676","rs10049390","rs1063192","rs7068313","rs66830472","rs7398375",
"rs3118233","rs4986080")


# 提取需要标注的 SNP 的信息
highlight_snp <- two_sample_MR_data[two_sample_MR_data$SNP %in% snp_to, ]

# ---绘制输出（“m”线性曼哈顿图，“c”环形曼哈顿图）
CMplot(two_sample_MR_data, 
       plot.type = "m",
       LOG10 = TRUE, 
       threshold = 1e-06, 
       threshold.lwd = 2, 
       threshold.lty = 2, 
       signal.cex = 0.2,
       chr.den.col = NULL, 
       cex = 0.2, 
       bin.size = 1e5, 
       ylim = c(0, 24), 
       width = 20, 
       height = 9,
       file.output = T,   #导出图片后，标注的位置会变化，距离点比较近
       file = "png",
       highlight = highlight_snp$SNP,  # 需要标注的 SNP
       highlight.col = "green",  # 标注点的颜色
       highlight.text = highlight_snp$SNP, 
       highlight.text.cex = 1.4,  # 标注文本大小
       highlight.text.col = "black",  # 标注文本颜色   
       
       axis.cex=1.6,  #刻度标签字体的大小
       axis.lwd=2,  #坐标轴的线宽
       lab.cex=2,   #坐标轴名称的字体大小
       lab.font=2,    #字体样式，1：普通字体； 2：粗体
       verbose = TRUE)
gc()

