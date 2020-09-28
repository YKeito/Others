pwd=/home/yasue/Nakano_RNAseq/fromKengo/
a=0
for sample in DMSO CY15 CY16 CY3 CY20 CY22 CY55
do
for rep in 1 2 3
do
a=$((a+1))
cufflinks -p 36 -o /home/yasue/bigdata/yasue/RNASeq/CY_flg22/cufflinks/${sample}_rep${rep}_cufflinks -G /root/db/TAIR10/TAIR10.gtf ${pwd}/160113_KMN_${a}.bam
done
done

