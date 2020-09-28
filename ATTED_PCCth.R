####閾値の設定####
#T_sumはエッジの数を表す
th <- 0.6 #任意の閾値

i <- 1
j <- 2
T_sum <- 0
total <- nrow(CoExp_PCC)

for (i in i:total){
  T_sum <- T_sum + sum(CoExp_PCC[i, j:ncol(CoExp_PCC)] > th)
  print(i)
  i <- i+1
  j <- i+1
}

#エッジ数の表示
print(T_sum)


####閾値でカットオフ(通常)####
i <- 1
j <- 2
k <- -1
total <- nrow(CoExp_PCC)
temp_A <- c()
temp_B <- c()
temp_C <- c()

for (i in i:total){
  for (j in j:total){
    k <- k+1
    if (CoExp_PCC[i,j] > th){
      temp_A <- c(temp_A, i)
      temp_B <- c(temp_B, j)
      temp_C <- c(temp_C,CoExp_PCC[i,j])
      j <- j+1
    } else { j <- j+1}
  }
  print(i)
  i <- i+1
  j <- i+1
}

#GeneIDとID_AGIを含むファイル(GeneID_AGI)を読み込む
th60_PCC <- cbind(matrix(rownames(CoExp_PCC)[temp_A], ncol=1),
                  matrix(temp_C),
                  matrix(rownames(CoExp_PCC)[temp_B], ncol=1))

colnames(th60_PCC) <- c("Source_Interaction","Interaction_Type","Target_Interaction")

#任意の閾値でカットオフできているか確認
nrow(th60_PCC) == T_sum


#kは全部のセルをサーチできてるかの確認用
#TRUEがでれば全部サーチできてる
c(c(nrow(CoExp_PCC)*ncol(CoExp_PCC)) - nrow(CoExp_PCC))/2 == k


#ノードの数の確認
length(unique(c(unique(th60_PCC[,1]),unique(th60_PCC[,3]))))


#Cytoscapeに読み込む用のデータを出力
write.table(th60_PCC, file="~/bigdata/yasue/ATTEDII/PCC/th60_PCC.txt", append=F, quote=F, sep="\t", row.names = F, col.names = T)