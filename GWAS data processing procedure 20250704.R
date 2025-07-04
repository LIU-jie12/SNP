
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

####处理1GCST90018808.tsv
setwd("F:/GCST")
Gwas-data <- read.table(file = "1GCST90018808.tsv", sep = "\t", header = TRUE)

Gwas-data1<-subset(Gwas-data,p_value<1e-06) 
write.table(Gwas-data1,"1GCST90018808.csv", sep = ",", row.names = FALSE) 


#无samplesize列应添加

data1 <- fread("1GCST90018808.csv") %>% 
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
gc()



exposure_data <- read_exposure_data(
  filename = "1GCST90018808-snp-1e6.csv",
  sep = ",",
  snp_col = "SNP",
  phenotype_col = "Phenotype", 
  beta_col = "beta", 
  se_col = "standard_error", 
  eaf_col = "effect_allele_frequency", 
  effect_allele_col = "effect_allele", 
  other_allele_col = "other_allele", 
  pval_col = "p_value", 
  ncase_col = "ncase",
  samplesize_col = "samplesize", 
  id_col = "id", 
  chr_col = "chromosome", 
  pos_col = "base_pair_location"
)

str(exposure_data)  #检查数据是否为数值型, 将非数值型的列的数据转换为数值型
exposure_data$samplesize.exposure <- as.numeric(exposure_data$samplesize.exposure)

# ---补齐 R² 和 F 列
exposure_data$R2 <- TwoSampleMR::get_r_from_bsen(exposure_data$beta.exposure, exposure_data$se.exposure, exposure_data$samplesize.exposure)^2
exposure_data$Fval <- (exposure_data$samplesize.exposure - 2) * exposure_data$R2 / (1 - exposure_data$R2)

# ---过滤保留 F>10 的工具变量
exp_date <- exposure_data[exposure_data$F > 10,]
print(paste("剩余", nrow(exp_date), "个 SNP"))

# ---去除连锁不平衡的 SNP 并保存
unclump_snps <- ieugwasr::ld_clump(
  dat = dplyr::tibble(
    rsid = exp_date$SNP, 
    other_allele.exposure= exp_date$other_allele.exposure,
    effect_allele.exposure= exp_date$effect_allele.exposure,
    pval.exposure= exp_date$pval.exposure,
    exposure= exp_date$exposure,
    pval_origin.exposure= exp_date$pval_origin.exposure,
    id.exposure= exp_date$id.exposure,
  ),
  clump_kb = 100000,
  clump_r2 = 0.1,
  plink_bin = NULL,
  bfile = NULL
)
write.csv(unclump_snps,"1GCST90018808-snp-1e6-clump.csv")


# 获取当前的工具变量表型并保存到文件
snp_with_trait <- FastTraitR::look_trait(
  rsids = unclump_snps$rsid,
  out_file = '1GCST90018808-snp-1e6-clump-trait.csv')

#####################结束









