# 示例ASV丰度表
asv_abundance <- data.frame(
  ASV = c("ASV1", "ASV2", "ASV3", "ASV4"),
  Sample1 = c(10, 5, 0, 3),
  Sample2 = c(0, 2, 7, 1)
)

# 示例ASV分类阶元
asv_taxonomy <- data.frame(
  ASV = c("ASV1", "ASV2", "ASV3", "ASV4"),
  Phylum = c("Firmicutes", "Proteobacteria", "Firmicutes", "Actinobacteria"),
  Class = c("Bacilli", "Gammaproteobacteria", "Clostridia", "Actinobacteria")
)
rm(list=ls())
setwd("C:\\Users\\dell\\Desktop\\convert")
asv_abundance <- read.table("microbetab.txt",header = T)
asv_taxonomy <- read.table("taxonomy_separated.txt",header = T)
library(tidyverse)
combined_data <- asv_abundance %>%
  left_join(asv_taxonomy, by = "ASV")
phylum_abundance <- combined_data %>%
  group_by(Phylum) %>%
  summarise(across(starts_with("G_"), sum, na.rm = TRUE))
class_abundance <- combined_data %>%
  group_by(Class) %>%
  summarise(across(starts_with("G_"), sum, na.rm = TRUE))
write.table(phylum_abundance,"phylum_abundance.txt",quote = F,row.names = F,sep="\t")
write.table(class_abundance,"class_abundance.txt",quote = F,row.names = F,sep="\t")