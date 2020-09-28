"~/Nakano_RNAseq/network_analysis/script/allCYFPKM_diff.R"
before <- proc.time()
#package---------
library(dplyr)
library(stringr)
#processing data----------
CY.FPKM <- c()
condition <- list.files(path = "/home/yasue/bigdata/yasue/RNASeq/CopyOf035_nakano_19021/cuffdiff_results/")
i <- 1
for(i in i:length(condition)){
  title <- paste0("/home/yasue/bigdata/yasue/RNASeq/CopyOf035_nakano_19021/cuffdiff_results/", condition[i], "/gene_exp.diff")
  T.data <- read.table(file = title, header = T, stringsAsFactors = F, fill = T, sep = "\t", quote = "")
  T.data <- T.data %>% select(test_id, gene, value_1, value_2, q_value)
  DMSO.data <- T.data %>% select(value_1) %>% +1 %>% as.matrix()
  CY.data <- T.data %>% select(value_2) %>% +1
  log2.data <- CY.data %>% sweep(STATS = DMSO.data, MARGIN = 1, FUN = "/") %>% log2()
  T.data <- T.data %>% mutate(log2.data$value_2)
  colnames(T.data) <- c("AGI", "ShortName", 
                        paste0("DMSO", str_sub(condition[i], start = 5), "_FPKM"), 
                        paste0(condition[i], "_FPKM"), 
                        paste0(condition[i], "_FDR"), 
                        paste0(condition[i], "_log2FC"))
  if(i == 1 | i == 2 | i == 3 | i == 4){
    if(is.null(CY.FPKM)){
      CY.FPKM <- T.data
    }else{
      CY.FPKM <- full_join(CY.FPKM, T.data %>% select(-ShortName), by = "AGI")
    }
  }else{
    CY.FPKM <- full_join(CY.FPKM, T.data %>% select(-ShortName, -3), by = "AGI")
  }
  print(i)
  i <- i+1
}
condition <- list.files(path = "/home/yasue/bigdata/yasue/RNASeq/CY_flg22/cuffdiff/")
i <- 1
for(i in i:length(condition)){
  title <- paste0("/home/yasue/bigdata/yasue/RNASeq/CY_flg22/cuffdiff/", condition[i], "/gene_exp.diff")
  T.data <- read.table(file = title, header = T, stringsAsFactors = F, fill = T, sep = "\t", quote = "")
  T.data <- T.data %>% select(test_id, gene, value_1, value_2, q_value)
  DMSO.data <- T.data %>% select(value_1) %>% +1 %>% as.matrix()
  CY.data <- T.data %>% select(value_2) %>% +1
  log2.data <- CY.data %>% sweep(STATS = DMSO.data, MARGIN = 1, FUN = "/") %>% log2()
  T.data <- T.data %>% mutate(log2.data$value_2)
  colnames(T.data) <- c("AGI", "ShortName", "DMSO_48h_FPKM", 
                        paste0(str_sub(condition[i], end = 4), "_48h_FPKM"), 
                        paste0(str_sub(condition[i], end = 4), "_48h_FDR"), 
                        paste0(str_sub(condition[i], end = 4), "_48h_log2FC"))
  if(i == 1){
    CY.FPKM <- full_join(CY.FPKM, T.data %>% select(-ShortName), by = "AGI")
  }else{
    CY.FPKM <- full_join(CY.FPKM, T.data %>% select(-ShortName, -3), by = "AGI")
  }
  print(i)
  i <- i+1
}
#save data-------------
saveRDS(CY.FPKM, "~/Nakano_RNAseq/network_analysis/.RData/RDS/allCY_cuffdiff.rds")
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#15.946 sec
#remove object-------
rm(list = ls())