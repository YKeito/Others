#48h:CY_flg22
for condition in `cat /export/fastq_bk/Nakano_CY/SampleName.txt`
do
#mkdir /home/yasue/bigdata/yasue/RNASeq/CY_flg22/${condition}/
R1=`cat /home/yasue/AllFastqName.txt | grep ${condition}_ | grep 'R1'`
R2=`cat /home/yasue/AllFastqName.txt | grep ${condition}_ | grep 'R2'`
tophat2 -g 1 -p 36 -o /home/yasue/bigdata/yasue/RNASeq/CY_flg22/${condition}/${condition}_mapped /root/db/TAIR10/TAIR10 /export/fastq_bk/Nakano_CY/$R1 /export/fastq_bk/Nakano_CY/$R2
done
