##################ATAC_cor_Gene_150kb_split_Chr.r
rm(list=ls());options(stringsAsFactors=FALSE)
wd = "/home/pencode/2024H3K27ac/co_activity"
setwd(wd);system(paste("cd",wd))
options(scipen=200)
options(digits=6)
library(psych)
argv <- commandArgs(TRUE)
# ----------Read data-----------
datCHIP <- read.table("enhancer_151sample_TPM.txt",header=T)
datRNA <- read.table("promoter_151sample_TPM.txt",header=T)
# ------- Cor ---------------------
ResultAll <- NULL
#Chr 1-20
j <- argv[1]
posChr <- which(datCHIP$CHR==j)
for(k in 1:length(posChr)){
enhancerStar <- datCHIP[posChr[k],5]-150000
enhancerEnd <- datCHIP[posChr[k],5]+150000
posRNA <- which(datRNA$CHR==j & datRNA$midpos>enhancerStar & datRNA$midpos<enhancerEnd)
Ngene <- length(posRNA)
if(Ngene>0){
        Result <- NULL
        if(Ngene==1){
        Result <- datCHIP[posChr[k],1:5]
        }else{
        Result <- datCHIP[posChr[k],1:5]
        Result[2:Ngene,] <- datCHIP[posChr[k],1:5]
        }
        Result <- cbind(Result,datRNA[posRNA[1:Ngene],1:5])
        Cor1 <- corr.test(as.numeric(t(datCHIP[posChr[k],6:ncol(datCHIP)])),(t(datRNA[posRNA[1:Ngene],6:ncol(datRNA)])),method = "spearman",adjust="none")
        Cor <- t(Cor1$r)
        Pvalue <- t(Cor1$p)
        Result <- cbind(Result,Cor)
        Result <- cbind(Result,Pvalue)
        ResultAll <- rbind(ResultAll,Result)
        }
        cat(k,"\n")
}
#colnames(ResultAll)[16] <- "Pvalue"
write.table(ResultAll,file=paste0("./",argv[1],".csv"),quote=F,sep=",",col.names = T,row.names = F)



####################################PeakcorGene.sh
#!/bin/bash
#PBS -N co_activity
#PBS -l nodes=1:ppn=1,mem=10gb
#PBS -e /home/pencode/qsublog
#PBS -o /home/pencode/qsublog
#PBS -M 714534452@qq.com
#PBS -m e
#PBS -q cu
#PBS -t 1-20

set -e
nprocs=`wc -l < $PBS_NODEFILE`

INPUT=/home/pencode/2024H3K27ac/co_activity
OUTPUT=/home/pencode/2024H3K27ac/co_activity

cd $INPUT

/opt/software/miniconda3/bin/Rscript ATAC_cor_Gene_150kb_split_Chr.r ${PBS_ARRAYID}



## network.r 
rm(list=ls());options(stringsAsFactors=FALSE,scipen=200,digits=6)
wd = "/work/jiangtao/network"
setwd(wd);system(paste("cd",wd))
library(psych)
argv <- commandArgs(TRUE)
library(reshape2)
library(parallel)
j <- argv[1]
j=as.numeric(as.character(j))
# ----------Read data-----------



enh_cor_pro_gene=read.table("/work/jiangtao/network/enhancer_cor_promoter_AllChr_gene.csv")
gene_list=read.table("/work/jiangtao/network/pro_conn_gene.txt")

tmp_gene=as.character(unique(gene_list$V2)[j])
tmp_gene=gsub("\\|", "\\\\\\|", tmp_gene)

tmp_enh <- enh_cor_pro_gene[grepl(paste0("\\b",tmp_gene,"\\b"), enh_cor_pro_gene$V14), ]
if(length(unique(tmp_enh$V1)) >1 ){
write.table(tmp_enh$V1,file=paste0("tmp_enh_",j,".txt"),col.names=F,row.names=F,quote=F)


input_file <- paste0("tmp_enh_", j, ".txt")
output_file <- paste0("tmp_enh_", j, ".txt2")
reference_file <- "/home/pencode/2024H3K27ac/co_activity/enhancer_151sample_TPM.txt"


cmd <- paste0("awk 'NR==FNR {keys[$1]; next} $1 in keys' ", 
              input_file, " ", 
              reference_file, " > ", 
              output_file)


system(cmd)

tmp_enh_density <- read.table(output_file)
rownames(tmp_enh_density)=tmp_enh_density$V1
num=nrow(tmp_enh_density)
df=t(tmp_enh_density[,7:ncol(tmp_enh_density)])
Cor1 <- corr.test(df,method = "spearman",adjust="none")

Cor <- Cor1$r


Cor[!upper.tri(Cor)] <- NA

Corr <- melt(Cor,na.rm=T)
rownames(Corr)=paste(Corr$Var1,Corr$Var2,sep="_")
Pvalue <- Cor1$p
Pvalue[!upper.tri(Pvalue)] <- NA

Pvalue2 <- melt(Pvalue,na.rm=T)
rownames(Pvalue2)=paste(Pvalue2$Var1,Pvalue2$Var2,sep="_")

common=intersect(rownames(Corr),rownames(Pvalue2))
Result <- data.frame(Corr[common,1:3],Pvalue2[common,3])
colnames(Result)=c("Enh1","Enh2","r","Pvalue")
Result$Gene=tmp_gene

Result2=Result[Result$Enh1!=Result$Enh2,]
write.table(Result2,file=paste0("./result3/AllResult_",j,".txt"),col.names=T,row.names=F,quote=F)
cat(tmp_gene,"\n")
}




library(ggplot2)

dat=read.table("/tmpdisk/total.txt",head=F)

df <- data.frame(Value = abs(dat$V3))
q95 <- as.numeric(quantile(df$Value, 0.95))

pdf("/tmpdisk/thread_r.pdf",height=4,width=4.5)
ggplot(df, aes(x = Value)) +
geom_density(color = "blue", size = 1) +
  geom_vline(xintercept = q95, color = "red", linetype = "dashed") +
  annotate("text", x = q95, y = 0.01, color = "red", hjust = -0.1,label = paste0("95%: ", round(q95, 2))) +
  labs( x = "Value",y = "Density") +
  theme_minimal()+theme_classic()+
  theme(  
legend.background=element_rect(fill="transparent",colour = NA),  
legend.text=element_text(colour="black",size=15),
axis.text.x = element_text(size=18),
axis.text.y = element_text(size=18),
axis.title.x = element_text(size=20),  
axis.title.y = element_text(size=20),
plot.margin = unit(c(1, 1, 1, 1), "lines"))
dev.off()





####
ls AllResult*txt|while read id
do
connect=$(awk 'NR > 1 && $3 > 0.8158592 {count++} END {print count}' "$id")
node=$(awk 'NR > 1 {print $1; print $2}' "$id" | sort | uniq | wc -l)
gene=$(awk 'NR > 1 {print $5}' "$id" | sort | uniq )
 if [ "$connect" -gt 0 ]; then
connectivity=$(echo "scale=2; $connect * 2 / $node" | bc)
echo -e "$gene\t$node\t$connect\t$connectivity" >> total_connectivity
else
    echo -e "$gene\t$node\t0\t0" >> total_connectivity
  fi
done


cat *total_connectivity >total_connectivity.txt

le only_one_enh_gene_list|awk '{print $1" "1" "0" "0}' >>total_connectivity.txt



