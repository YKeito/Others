#read data
#CY15_FDR005 <- allRNASeq[allRNASeq$CY15_1h_q_value < 0.05 | allRNASeq$CY15_3h_q_value < 0.05 | allRNASeq$CY15_12h_q_value < 0.05 | allRNASeq$CY15_24h_q_value < 0.05 | allRNASeq$CY15_48h_q_value < 0.05, 1:5]
#CY16_FDR005 <- allRNASeq[allRNASeq$CY16_1h_q_value < 0.05 | allRNASeq$CY16_3h_q_value < 0.05 | allRNASeq$CY16_12h_q_value < 0.05 | allRNASeq$CY16_24h_q_value < 0.05 | allRNASeq$CY16_48h_q_value < 0.05, 6:10]
#CY20_FDR005 <- allRNASeq[allRNASeq$CY20_1h_q_value < 0.05 | allRNASeq$CY20_3h_q_value < 0.05 | allRNASeq$CY20_12h_q_value < 0.05 | allRNASeq$CY20_24h_q_value < 0.05 | allRNASeq$CY20_48h_q_value < 0.05, 11:15]
i <- 1
allfile <- list(CY15_FDR005, CY16_FDR005, CY20_FDR005)
for(i in i:length(allfile)){
  CY <- c("CY15", "CY16", "CY20")
  temp <- allfile[[i]]
  n <- 1
  m <- 2
  base <- c()
  PCC <- c()
  PCC_all <- list()
  PCC_pvalue <- c()
  PCC_pvalue_all <- list()
  for(n in n:c(nrow(temp)-1)){
    for(m in m:nrow(temp)){
      base <- cor.test(as.numeric(temp[n, ]), as.numeric(temp[m, ]), method = "pearson")
      PCC <- c(PCC, base$estimate)
      PCC_pvalue <- c(PCC_pvalue, base$p.value)
      m <- m+1
    }
    PCC_all <- c(PCC_all, list(PCC))
    PCC <- c()
    PCC_pvalue_all <- c(PCC_pvalue_all, list(PCC_pvalue))
    PCC_pvalue <- c()
    print(n)
    n <- n+1
    m <- n+1
  }
  ####PCC q-value####
  PCC_qvalue_all <- p.adjust(unlist(PCC_pvalue_all), method = "BH")
  #####Cytoscape_format####
  source_genes <- combn(rownames(temp), 2)[1, ]
  target_genes <- combn(rownames(temp), 2)[2, ]
  cytoscape <- data.frame(source_genes = source_genes, 
                          interaction_value = unlist(PCC_all),
                          abs = abs(unlist(PCC_all))
                          target_genes = target_genes,
                          p_value = unlist(PCC_pvalue_all),
                          q_value = PCC_qvalue_all,
                          stringsAsFactors = F
                          )
  cytoscape <- cytoscape[cytoscape$q_value < 0.05, ]
  File_title <- paste0("~/Nakano_RNAseq/network_analysis/cytoscape/CY_pattern/", paste0(CY[i], "_FDR005"),".txt")
  write.table(cytoscape, File_title, append=F, quote = F, sep = "\t", row.names = F)
}
