#"~/Nakano_RNAseq/network_analysis/script/CYPCC_allFPKM.R"
before <- proc.time()
#package--------------
#install.packages("Hmisc")
#install.packages("Hmisc")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Hmisc))
#processing data----------
CY.FPKM <- readRDS("~/Nakano_RNAseq/network_analysis/.RData/RDS/allCYFPKM.rds")
T.CY <- c("CY15", "CY16", "CY20")
check <- c()
i <- 1
for(i in i:length(T.CY)){
  T.data <- CY.FPKM %>% select(c(starts_with("DMSO"), starts_with(T.CY[i])))
  T.data <- t(T.data)
  base <- rcorr(as.matrix(T.data), type = "pearson")
  #####Cytoscape_format####
  temp <- combn(colnames(T.data), 2)
  CY15.PCC <- data.frame(source_genes = temp[1, ],
                         interaction_value = base$r[lower.tri(base$r)],
                         target_genes = temp[2, ],
                         p_value = base$P[lower.tri(base$P)],
                         q_value = p.adjust(base$P[lower.tri(base$P)], method = "BH"),
                         stringsAsFactors = F
  )
  CY15.PCC$interaction_value[is.na(CY15.PCC$interaction_value)] <- 0
  CY15.PCC$p_value[is.na(CY15.PCC$p_value)] <- 1
  CY15.PCC$q_value[is.na(CY15.PCC$q_value)] <- 1
  title <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/CYPCC/", T.CY[i], "_PCC_dataframe.rds")
  saveRDS(object = CY15.PCC, file = title)
  CY15.PCC <- CY15.PCC[CY15.PCC$q_value < 1e-3, ]
  check <- c(check, nrow(CY15.PCC))
  #g <- ggplot(CY15.PCC, aes(x = interaction_value, y = ..density..))
  #g <- g + geom_histogram(position = "identity", alpha = 0.8)
  #g <- g + geom_density()
  #title <- paste0("~/bigdata/yasue/PCCOfPCC/CY/supple/", T.CY[i], "_PCC.png")
  #ggsave(file = title, plot = g)
  title <- paste0("~/bigdata/yasue/PCCOfPCC/RDS/CYPCC/", T.CY[i], "_PCC_th_FDR1e3.rds")
  saveRDS(object = CY15.PCC, file = title)
  title <- paste0("~/bigdata/yasue/PCCOfPCC/CY/Table/", T.CY[i], "PCC.txt")
  write.table(CY15.PCC, file = title, sep = "\t", quote = F, row.names = F)
  rm(base)
  rm(temp)
  rm(CY15.PCC)
  #rm(g)
  rm(T.data)
  print(i)
  i <- i+1
}
#elapsed time------------------------------------------
after <- proc.time()
print(after - before)#7057.942 sec