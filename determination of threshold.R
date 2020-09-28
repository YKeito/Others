#CY15_cytoscape <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/CY_pattern/Sum_CY15_FDR005.txt", sep = "\t", header = T, stringsAsFactors = F)
#CY16_cytoscape <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/CY_pattern/Sum_CY16_FDR005.txt", sep = "\t", header = T, stringsAsFactors = F)
#CY20_cytoscape <- read.table("~/Nakano_RNAseq/network_analysis/cytoscape/CY_pattern/Sum_CY20_FDR005.txt", sep = "\t", header = T, stringsAsFactors = F)
#nrow(CY15_cytoscape)*0.05 #355228.2
#nrow(CY16_cytoscape)*0.05 #129789.1
#nrow(CY20_cytoscape)*0.05 #39910.8
#hist(CY15_cytoscape$interaction_value[CY15_cytoscape$p_value < 0.01])
#hist(CY16_cytoscape$interaction_value[CY16_cytoscape$p_value < 0.01])
#hist(CY20_cytoscape$interaction_value[CY20_cytoscape$p_value < 0.01])

#threshold p value
sample <- list(CY15_cytoscape, CY16_cytoscape, CY20_cytoscape)
samplename <- c("CY15", "CY16", "CY20")
value <- c(0.001, 0.005, 0.01, 0.05)
n <- 1
for(n in n:length(sample)){
  temp <- sample[[n]]
  Numedges <- c()
  i <- 1
  for(i in i:length(value)){
    Numedges <- c(Numedges, sum(temp$p_value < value[i]))
    i <- i+1
  }
  plot(value, Numedges,
       main=paste0(samplename[n], ":Number of edges in various thresholds of PCC"),
       xlab="p-value used for cut-off", ylab="Numbers of edges",
  )
  borderline <- nrow(temp)*0.05
  abline(h=borderline, lty=2)
  n <- n+1
}

#threshold PCC
sample <- list(CY15_cytoscape, CY16_cytoscape, CY20_cytoscape)
samplename <- c("CY15", "CY16", "CY20")
value <- c(0.9, 0.91, 0.92, 0.93, 0.95)                                     #c(0.92, 0.921, 0.922, 0.923, 0.924, 0.925, 0.926, 0.927, 0.928, 0.929, 0.93)                        
n <- 1
for(n in n:length(sample)){
  temp <- sample[[n]]
  Numedges <- c()
  i <- 1
  for(i in i:length(value)){
    Numedges <- c(Numedges, sum(abs(temp$interaction_value) > value[i]))
    i <- i+1
  }
  plot(value, Numedges,
       main=paste0(samplename[n], ":Number of edges in various thresholds of PCC"),
       xlab="PCC used for cut-off", ylab="Numbers of edges",
  )
  borderline <- nrow(temp)*0.05
  abline(h=borderline, lty=2)
  n <- n+1
}

