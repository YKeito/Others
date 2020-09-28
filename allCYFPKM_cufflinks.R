#"~/Nakano_RNAseq/network_analysis/script/allCYFPKM_cufflinks.R"
before <- proc.time()
#package---------
library(dplyr)
#processing data----------
CY.FPKM <- c()
DirectoryInfo <- read.table("~/bigdata/yasue/RNASeq/NakanoRNASeq_Info.txt", sep = "\t", header = T, stringsAsFactors = F)
condition <- c("DMSO", "CY15", "CY16", "CY20")
T.time <- c("1h", "3h", "12h", "24h")
n <- 1
for(n in n:length(condition)){
  m <- 1
  for(m in m:length(T.time)){
    T.sample <- DirectoryInfo$sample[grep(paste(condition[n], T.time[m], sep = "_"), DirectoryInfo$sample)]
    T.directory <- DirectoryInfo$directory[grep(paste(condition[n], T.time[m], sep = "_"), DirectoryInfo$sample)]
    i <- 1
    for(i in i:length(T.sample)){
      title <- paste0("~/bigdata/yasue/RNASeq/CopyOf035_nakano_19021/cufflinks_results/", T.directory[i], "_outputs/", "genes.fpkm_tracking")
      T.FPKM <- read.table(title, sep = "\t", fill = T, header = T, stringsAsFactors = F, quote = "")
      T.FPKM <- T.FPKM[!duplicated(T.FPKM$tracking_id), ]
      T.FPKM <- T.FPKM %>% select("tracking_id", "FPKM")
      colnames(T.FPKM)[2] <- T.sample[i]
      if(is.null(CY.FPKM)){
        CY.FPKM <- T.FPKM
      }else{
        CY.FPKM <- full_join(CY.FPKM, T.FPKM, by = "tracking_id")
      }
      i <- i+1
    }
    m <- m+1
  }
  print(n)
  n <- n+1
}
rownames(CY.FPKM) <- CY.FPKM$tracking_id
CY.FPKM <- CY.FPKM[, 2:ncol(CY.FPKM)]
#save data-------------
saveRDS(CY.FPKM, "~/Nakano_RNAseq/network_analysis/.RData/RDS/allCYFPKM.rds")
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#15.946 sec
#remove object-------
rm(list = ls())