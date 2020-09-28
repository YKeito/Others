#"~/Nakano_RNAseq/network_analysis/script/allCYFPKM_log2FC.R"
#package----
library(stringr)
library(dplyr)
#input data----
CY.FPKM <- readRDS("~/Nakano_RNAseq/network_analysis/.RData/RDS/allCYFPKM.rds")
#processing data----
T.sample <- str_split(colnames(CY.FPKM), pattern = "_rep", simplify = T)[, 1] %>% unique()
T.data <- c()
j <- 1
for(j in j:length(T.sample)){
  T.data <- cbind(T.data, CY.FPKM %>% select(contains(T.sample[j])) %>% apply(MARGIN = 1, FUN = mean))
  colnames(T.data)[j] <- T.sample[j]
}
#save AverageExp----
saveRDS(T.data, "~/Nakano_RNAseq/network_analysis/.RData/RDS/AverageOfallCYFPKM.rds")
#processing data----
T.data <- T.data %>% as.data.frame()
T.time <- str_split(colnames(T.data), pattern = "_", simplify = T)[, 2] %>% unique()
T.FC <- c()
j <- 1
for(j in j:length(T.time)){
  T.FPKM <- T.data %>% select(contains(T.time[j]))
  k <- 2
  for(k in k:ncol(T.FPKM)){
    T.FC <- cbind(T.FC, log2(c(c(T.FPKM[, k]+0.01)/c(T.FPKM[, 1] + 0.01))))
  }
}
colnames(T.FC) <- paste0(T.sample[5:length(T.sample)], "/control +0.01")
#save log2FC----
saveRDS(T.FC, "~/Nakano_RNAseq/network_analysis/.RData/RDS/log2FoldChangeOfallCYFPKM.rds")
#remove object----
rm(list = ls())