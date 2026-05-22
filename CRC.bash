#/home/pencode/genome

##安装python
wget https://www.python.org/ftp/python/2.7.18/Python-2.7.18.tgz

 tar -xvf Python-2.7.18.tgz
 cd Python-2.7.18
./configure
make

##写入全局环境变量
vi ~/.bashrc
source ~/.bashrc

##依赖软件
alias python='/home/pencode/genome/Python-2.7.18/python'
alias fimo='/home/pencode/genome/meme-5.5.4/src/fimo'
alias fasta-get-markov='/home/pencode/genome/meme-5.5.4/src/fasta-get-markov'
export PATH="/home/pencode/genome/Python-2.7.18/:$PATH"
export PATH="/home/pencode/genome/samtools-1.3/bin:$PATH"
export PATH="/home/pencode/genome/R-4.3.3/bin:$PATH"
export PATH="/home/pencode/genome/meme-5.5.4/src:$PATH"


##安装ROSE
git clone https://bitbucket.org/young_computation/rose.git
##制作注释文件
#wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
#chmod 775 gtfToGenePred
#./gtfToGenePred ./pig_gold/DRCv20_final5.1.gtf  ./ROSE-1.3.2/annotation/DRC20_refseq.ucsc
#chmod 775 DRC20_refseq.ucsc

##另一种方法
##提供ref_UCSC格式的gtf文件
/home/pencode/genome/Python-2.7.18/python -m pip install pandas --user
/home/pencode/genome/Python-2.7.18/python gtf_to_refseq_ucsc.py your_input_file.gtf output_refseq


#将注释文件名字添加到 ROSE_main.py中genomeDict
vi   ROSE_main.py
#'PIGV20':'%s/annotation/PIGV20_refseq.ucsc' % (cwd),

sed -i 's/python ROSE_bamToGFF.py/\/home\/pencode\/genome\/Python-2.7.18\/python ROSE_bamToGFF.py/g' /home/pencode/software/young_computation-rose-feb35cb1d955/ROSE_main.py


#鉴定超级增强子
python  /home/pencode/genome/rose/ROSE_main.py \
-g PIGV20 \
-i /home/pencode/genome/rose/result/${Sample2}.gff \ #less narrowPeak|awk '{print $1"\t"$4"\t"" ""\t"$2"\t"$3"\t"" ""\t"$6}' >/home/pencode/genome/rose/result/${Sample2}.gff
-r ${R1[$PBS_ARRAYID]} \
-o /home/pencode/genome/rose/result/${Sample2} \
-s 12500 \
-t 1000   ##距离TSS位点1kb作为enhancer


##安装CRCmapper依赖
/home/pencode/genome/Python-2.7.18/python -m ensurepip --upgrade --user
/home/pencode/genome/Python-2.7.18/python -m pip install numpy --user
/home/pencode/genome/Python-2.7.18/python -m pip install scipy --user
/home/pencode/genome/Python-2.7.18/python -m pip install networkx --user
/home/pencode/genome/Python-2.7.18/python -m pip install matplotlib --user


##制作TFlist_NMid_DRCV20.txt；
根据PIGV20 gtf （PIGV20_gene 转录本加gene name）从VertebratePWMs.txt选择可用的TF文件

le VertebratePWMs.txt|grep MOTIF|awk '{print $3}'|sed 's/ //g'|sed 's/\t//g'|sed -i 's/\r$//g'|sort|uniq|while read id;
 do 
 grep -iw $id PIGV20_gene >>TFlist_NMid_DRCV20.txt;
 done

 
##运行CRC
python CRCmapper.py  \
-e /home/pencode/genome/rose/small_AllEnhancers.table.txt \
-b /home/pencode/genome/rose/small_sample_sorted.bam \
-g PIGV20 \
-f ./DRCV20_fa/ \ ##分染色体fa文件目录
-o ./ \   ##存放结果的目录，记得末尾加/
-n annotation \ ##结果文件的basename
-s /home/pencode/genome/rose/small.gff \ ##narrowPeak改成bed结尾
-x 33 \    #TF是否表达
-l 500     #左右扩500bp检索TF



##超级增强子分析#dbSUPER和SEdb这两个超级增强子数据库都是使用h3K27ac组蛋白修饰作为mark来识别超级增强子
##重新作图
setwd("/home/pencode/genome/rose/result/HS4M42a1L")
enhancer = read.table("HS4M42a1L_AllEnhancers.table.txt",  sep="\t",head=T)
H3K27ac = sort(enhancer[,7])

pdf("HS4M42a1L_SuperEnh.pdf",width=5,height=5)
par(mai=c(0.5,0.5,0.3,0.2),mgp=c(1.3,0.3,0))
plot(H3K27ac, col=2,lwd=3.5,type="l",xlab="Enhancers ranked",ylab="Signal at enhancers")


numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}

inputVector <- H3K27ac
#set those regions with more control than ranking equal to zero
inputVector[inputVector<0]<-0 

# This is the slope of the line we want to slide. This is the diagonal.
slope <- (max(inputVector)-min(inputVector))/length(inputVector) 

# Find the x-axis point where a line passing through that point has the minimum number 
xPt <- floor(optimize(numPts_below_line, lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) 

y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.

b <- y_cutoff-(slope* xPt)
abline(v= xPt,h= y_cutoff,lty=2,col=8)
points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
abline(coef=c(b,slope),col=2,lwd=3.5)


SupEnh=enhancer[enhancer$V2>=y_cutoff,1]
#filename=paste0("./SuperEnh_peaks/",args[1],"_super_enhancer_cluster.bed")
#write.table(SupEnh,filename,row.names=F,col.names=F,quote=F, sep="\t")

title(paste0("Super enhancer identified = ",length(SupEnh)))
# title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=", 
        # signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",
        # signif(y_cutoff/mean(inputVector),3),"x",sep=""))

#Number of regions with zero signal
# axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") 
dev.off()


##举例展示
#/home/pencode/2024H3K27ac/co_activity/figure_Brain12

super_enh: chr9    53685793        53739765
gene:chr9         53680436        53717972 DCPS



 singularity exec -B /home/PersonalFiles/Liusiyi/cutTag/H3K27AC/HS4M42a1L_H3K27AC/mapping:/tmpdisk \
 -B /home/pencode/genome/rose/result/HS4M42a1L:/tmpdisk2\
 /home/Singusoft/143.deeptools.sif bamCoverage -p 10 --binSize 20 -b /tmpdisk/HS4M42a1L_H3K27AC_sorted.bam -o /tmpdisk2/HS4M42a1L_H3K27AC_sorted.bw;

 
singularity exec  -B /home/pencode/2024H3K27ac/co_activity/bw_file/tissue:/tmpd \
-B /home/pencode/genome/rose/result/HS4M42a1L:/tmpd2 \
-B /home/pencode/2024H3K27ac/co_activity/figure_Brain12:/tmpd3  \
/home/Singusoft/pygenometracks.sif pyGenomeTracks \
--tracks /tmpd3/DCPS_track.ini --region chr9:53682793-53742765 \
--width 25  --height 40 --fontSize 10 --outFileName \
/tmpd3/DCPS_track_SE.pdf



##长度分布、数量   /home/pencode/genome/rose/result
##组成超级增强子的普通增强子数量分布
 find -name *SuperEnhancers.table.txt|while read id
 do
 le $id|grep chr|awk '{print $5}' >>./total/SE_content_enh.txt
 le $id|grep chr|awk '{print $4-$3}' >>./total/SE_len.txt
 le $id|grep chr|wc -l >>./total/SE_num.txt
 done

 

cd /home/pencode/genome/rose/result/total
library(dplyr)
library(magrittr)
library(ggplot2)


numbers=read.table("SE_len.txt")

 numbers$V2=numbers$V1/1000
  median_value <- median(numbers$V2)
   mean_value <- mean(numbers$V2)
   hist_density <- hist(numbers$V2, plot = FALSE)
  max_density <- max(hist_density$density)  
   mean_pos_x <- mean_value
  median_pos_x <- median_value
   mean_pos_y <- 80  # 均值文本距离顶部稍微下方
  median_pos_y <- 100  # 中位数文本距离顶部稍微下方
pdf("super_Enhancer_len.pdf",width=4,height=4)
hist(numbers$V2, 
     main="Distribution of SuperEnhancer length", # 图表标题
     xlab="Length",                   # x轴标签
     ylab="Frequency",               # y轴标签
     col="lightblue",                # 直方图颜色
     border="black",                 # 直方图边框颜色
     breaks=50)                      # 设置箱数（分割区间数）

# 如果需要平滑曲线，可以添加密度曲线
#lines(density(numbers$V2), col="red", lwd=2)  # 红色平滑曲线
  # 添加均值和中位数的标注
text(mean_pos_x, mean_pos_y, 
       labels = paste("Mean:", round(mean_value, 2)), col = "blue", pos = 4, cex = 0.8)
  
  text(median_pos_x, median_pos_y, 
       labels = paste("Median:", round(median_value, 2)), col = "green", pos = 4, cex = 0.8)
  
dev.off()

 
 
numbers=read.table("SE_content_enh.txt")

  median_value <- median(numbers$V1)
   mean_value <- mean(numbers$V1)
   hist_density <- hist(numbers$V1, plot = FALSE)
  max_density <- max(hist_density$density)  
   mean_pos_x <- mean_value
  median_pos_x <- median_value
   mean_pos_y <- 80  # 均值文本距离顶部稍微下方
  median_pos_y <- 100  # 中位数文本距离顶部稍微下方
pdf("super_Enhancer_num.pdf",width=4,height=4)
hist(numbers$V1, 
     main="Distribution of SuperEnhancer length", # 图表标题
     xlab="Number",                   # x轴标签
     ylab="Frequency",               # y轴标签
     col="lightblue",                # 直方图颜色
     border="black",                 # 直方图边框颜色
     breaks=10)                      # 设置箱数（分割区间数）

# 如果需要平滑曲线，可以添加密度曲线
lines(density(numbers$V1), col="red", lwd=2)  # 红色平滑曲线
  # 添加均值和中位数的标注
text(mean_pos_x, mean_pos_y, 
       labels = paste("Mean:", round(mean_value, 2)), col = "blue", pos = 4, cex = 0.8)
  
  text(median_pos_x, median_pos_y, 
       labels = paste("Median:", round(median_value, 2)), col = "green", pos = 4, cex = 0.8)
  
dev.off() 
 
 
 

numbers=read.table("SE_num.txt")

  median_value <- median(numbers$V1)
   mean_value <- mean(numbers$V1)
   hist_density <- hist(numbers$V1, plot = FALSE)
  max_density <- max(hist_density$density)  
   mean_pos_x <- mean_value
  median_pos_x <- median_value
   mean_pos_y <- 80  # 均值文本距离顶部稍微下方
  median_pos_y <- 100  # 中位数文本距离顶部稍微下方
pdf("super_Enhancer_num.pdf",width=4,height=4)
hist(numbers$V1, 
     main="Distribution of SuperEnhancer length", # 图表标题
     xlab="Number",                   # x轴标签
     ylab="Frequency",               # y轴标签
     col="lightblue",                # 直方图颜色
     border="black",                 # 直方图边框颜色
     breaks=10)                      # 设置箱数（分割区间数）

# 如果需要平滑曲线，可以添加密度曲线
lines(density(numbers$V1), col="red", lwd=2)  # 红色平滑曲线
  # 添加均值和中位数的标注
text(mean_pos_x, mean_pos_y, 
       labels = paste("Mean:", round(mean_value, 2)), col = "blue", pos = 4, cex = 0.8)
  
  text(median_pos_x, median_pos_y, 
       labels = paste("Median:", round(median_value, 2)), col = "green", pos = 4, cex = 0.8)
  
dev.off() 
  
 
 
 
 
 
##/home/pencode/genome/rose/result/HS4M42a1L

 le HS4M42a1L_AllEnhancers.table.txt|awk '{if($9=="0")print $2"\t"$3"\t"$4}' >nonSE.bed
 bedtools intersect -a nonSE.bed -b ./gff/HS4M42a1L.bed -wa -wb|awk '{print $4"\t"$5"\t"$6}' >nonSE_TE.bed
bedtools intersect -a nonSE_TE.bed -b /work/jiangtao/network/transcript_TSS_1kb_gene_rmfigure0.bed -v >nonSE_TE_distal.bed

# run compute matrix to collect the data needed for plotting
 singularity exec  -B /home/pencode/genome/rose/result/HS4M42a1L:/tmpdisk /home/Singusoft/143.deeptools.sif computeMatrix scale-regions \
 -p 10  -b 2500 -a 2500  \
 -R  /tmpdisk/HS4M42a1L_Gateway_SuperEnhancers.bed /tmpdisk/nonSE_TE_distal.bed \
 -S /tmpdisk/HS4M42a1L_H3K27AC_sorted.bw   \
--skipZeros -o /tmpdisk/matrix.mat.gz

singularity exec  -B /home/pencode/genome/rose/result/HS4M42a1L:/tmpdisk /home/Singusoft/143.deeptools.sif  plotHeatmap -m /tmpdisk/matrix.mat.gz  -out /tmpdisk/SE_TE_region_signal.pdf

 
 
##构建一个同一的超级增强子矩阵、density /home/pencode/genome/rose/result
find -name "*SuperEnhancers.table.txt"|while read id; do le $id|grep chr|awk '{print $2"\t"$3"\t"$4}'|sort|uniq >./total/sample/$(basename $id); done

cd ./total/sample/
ls *SuperEnhancers.table.txt|while read id
> do
> mv $id ${id%%_SuperEnhancers.table.txt*}_Gateway_SuperEnhancers.bed
> done

find -name "*Gateway_SuperEnhancers.bed"|while read id
 do
 cat $id >>151_SE.bed
 done


le 151_SE.bed |awk '{print $1"\t"$2"\t"$3}'|sort -k1,1 -k2,2n|bedtools merge >../151_SE_merge.bed
##7448
rm 151_SE.bed


ls *Gateway_SuperEnhancers.bed|while read id
 do
 bedtools intersect -a ../151_SE_merge.bed -b $id -wa|sort|uniq >${id%%_Gateway_SuperEnhancers.bed*}_intersect_SE.bed
 done

 
##01文件
file = dir(,pattern="intersect_SE.bed")

for(i in 1:length(file)){
 tmp=read.table(as.character(file[i]))
 tmp$name=paste(tmp$V1,tmp$V2,tmp$V3,sep="_")
 tmp$value=1
tmp=tmp[,4:5] 
 if(!exists("result")){
  result=tmp
   }else{
  result <- merge(result, tmp, by="name", all=TRUE)
  }
  }
  
rownames(result)=result$name
result=result[,2:ncol(result)]
colnames(result)=gsub("_intersect_SE.bed","",file) 
colnames(result)=sapply(colnames(result),function(x) {unlist(strsplit(x,"_"))[1]})
colnames(result)=gsub("HS4M56b1L","HS4M56b1L_A",colnames(result))
result[is.na(result)] <- 0

info=read.table("/home/pencode/promoter/tissue_id.txt2")
rownames(info)=info$V1

common=intersect(colnames(result),info$V1)
result2=result[,common]
info2=info[common,]

colnames(result2)=info2$V4
write.table(result2,"151SE_01file.txt",col.names=T,row.names=T,quote=F)
result3=t(result2)


# 安装并加载 vegan 包
install.packages("vegan")  # 如果尚未安装
library(vegan)

# 计算 pairwise Jaccard 距离
jaccard_dist <- vegdist(result3, method = "jaccard")

# 显示距离矩阵
jaccard_matrix=as.matrix(jaccard_dist)

# 经典多维缩放（MDS）转换 Jaccard 距离
pca_res <- cmdscale(jaccard_matrix, eig = TRUE, k = 2)  # 取前2个主成分
aa=pca_res$points
aa2=as.data.frame(as.matrix(aa))
colnames(aa2)=c("PC1","PC2")
aa2$group=rownames(aa2)
aa2$group=sapply(aa2$group,function(x) {unlist(strsplit(x,"\\."))[1]})


# 可视化 PCA 结果
library(ggplot2)

pdf(file="SE_PCA1_PCA2.pdf",width=14,height=6)
ggplot(data = aa2,aes(x = PC1,y = PC2))+geom_point(aes(color=group),size = 3)+ 
 theme_set(theme_bw())+theme(panel.grid = element_blank())+
 theme(panel.border = element_blank(),axis.line = element_line(colour="black",size=1))+
 theme(panel.grid = element_blank(),panel.background=element_blank())+
 theme(legend.text=element_text(size=15),legend.title=element_blank())+
 scale_color_manual(values=c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5",
 "#7A4900", "#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#8CD0FF","#5A0007",
 "#809693","#6A3A4C","#1B4400", "#4FC601","#3B5DFF","#4A3B53", "#FF2F80","#61615A",
 "#BA0900","#6B7900","#00C2A0","#FFAA92","#04F757","#D16100","#DDEFFF","#000035", "#7B4F4B","#A1C299",
 "#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF",
 "#C2FF99", "#001E09","#00489C", "#6F0062","#0CBD66","#EEC3FF","#456D75",
 "#B77B68","#7A87A1","#788D66","#788D66", "#FAD09F", "#FF913F","#D790FF","#922329","#A4E804"),
 breaks=c("Adipose","Adrenal_Gland","Anus","Bladder","Blood_Vessel","Bone_marrow","Brain", "Brain_stem",
 "Bronchi","Bulbourethral_gland","Cartilage","Cerebellum","Choroidal_plexus","Diaphragm","Endocardium_1","Endocardium_2",
 "Esophagus_membrane","Fallopian_tube","Gallbladder","Heart","Heart_valve","Kidney","Large_intestine","Ligament",
 "Liver","Lung","Lymph_node","Mammary_gland","Mesentery","Oralcavity_membrane","Ovary","Palate","Pancreas",
 "Parathyroid_gland","Parotid_gland","Penis","Pituitary","Respiratory_mucosa","Retina","Sciatic_nerve","Seminal_vesicles",
 "Skeletal_muscle","Skin","Small_intestine","Spinal_cord","Spleen","Stomach","Sublingual_gland","Submandibular_gland",
 "Thalamus","Thymus","Thyroid_gland","Tongue","Ureter","Urethra","Uterus","Vagina","Vocal_cords"),
 labels=c("Adipose","Adrenal_Gland","Anus","Bladder","Blood_Vessel","Bone_marrow","Brain","Brain_stem",
 "Bronchi","Bulbourethral_gland","Cartilage","Cerebellum","Choroidal_plexus","Diaphragm","Endocardium_1","Endocardium_2",
 "Esophagus_membrane","Fallopian_tube","Gallbladder","Heart","Heart_valve","Kidney","Large_intestine","Ligament",
 "Liver","Lung","Lymph_node","Mammary_gland","Mesentery","Oralcavity_membrane","Ovary","Palate","Pancreas",
 "Parathyroid_gland","Parotid_gland","Penis","Pituitary","Respiratory_mucosa","Retina","Sciatic_nerve","Seminal_vesicles",
 "Skeletal_muscle","Skin","Small_intestine","Spinal_cord","Spleen","Stomach","Sublingual_gland","Submandibular_gland",
 "Thalamus","Thymus","Thyroid_gland","Tongue","Ureter","Urethra","Uterus","Vagina","Vocal_cords"))+
 labs(x = "PC1", y = "PC2")+
 theme(axis.title= element_text(size=16,color="black",vjust=0.5, hjust=0.5))+
 theme(axis.text.x=element_text(hjust = 1,vjust =0.5,color ="black",size=16),
 axis.text.y=element_text(size=16,color="black"))
dev.off()   


# 安装必要的包
install.packages("ggplot2")
install.packages("ggdendro")
install.packages("dendextend")
install.packages("vegan")
# 加载包
library(ggplot2)
library(ggdendro)
library(dendextend)
library(vegan)

setwd("G:\\2024H3K27ac流程\\增强子网络")
result2=read.table("151SE_01file.txt",head=T)
colnames(result2)=sapply(colnames(result2),function(x) {unlist(strsplit(x,"\\."))[1]})
result3=t(result2)
jaccard_dist <- vegdist(result3, method = "jaccard")




# 计算层次聚类
hclust_res <- hclust(jaccard_dist, method = "complete")

# 转换为 dendrogram 对象
dendro <- as.dendrogram(hclust_res)

# 提取树叶的标签名称
labels <- labels(dendro)
# 定义样本名称对应的颜色
sample_colors <- c(
  "Adipose"="#FFFF00", "Adrenal_Gland"="#1CE6FF", "Anus_1"="#FF34FF",
  "Bladder"="#FF4A46", "Blood_Vessel"="#008941", "Bone_marrow"="#006FA6",
  "Brain"="#A30059", "Brain_stem"="#FFDBE5", "Bronchi"="#7A4900",
  "Bulbourethral_gland"="#0000A6", "Cartilage"="#63FFAC", "Cerebellum"="#B79762",
  "Choroidal_plexus"="#004D43", "Diaphragm"="#8FB0FF", "Endocardium_1"="#8CD0FF",
  "Endocardium_2"="#5A0007", "Esophagus_membrane"="#809693", "Fallopian_tube"="#6A3A4C",
  "Gallbladder"="#1B4400", "Heart"="#4FC601", "Heart_valve"="#3B5DFF",
  "Kidney"="#4A3B53", "Large_intestine"="#FF2F80", "Ligament"="#61615A",
  "Liver"="#BA0900", "Lung"="#6B7900", "Lymph_node"="#00C2A0",
  "Mammary_gland"="#FFAA92", "Mesentery"="#04F757", "Oralcavity_membrane"="#D16100",
  "Ovary"="#DDEFFF", "Palate"="#000035", "Pancreas"="#7B4F4B",
  "Parathyroid_gland"="#A1C299", "Parotid_gland"="#300018", "Penis"="#0AA6D8",
  "Pituitary"="#013349", "Respiratory_mucosa"="#00846F", "Retina"="#372101",
  "Sciatic_nerve"="#FFB500", "Seminal_vesicles"="#C2FFED", "Skeletal_muscle"="#A079BF",
  "Skin"="#C2FF99", "Small_intestine"="#001E09", "Spinal_cord"="#00489C",
  "Spleen"="#6F0062", "Stomach"="#0CBD66", "Sublingual_gland"="#EEC3FF",
  "Submandibular_gland"="#456D75", "Thalamus"="#B77B68", "Thymus"="#7A87A1",
  "Thyroid_gland"="#788D66", "Tongue"="#788D66", "Ureter"="#FAD09F",
  "Urethra"="#FF913F", "Uterus"="#D790FF", "Vagina"="#922329", "Vocal_cords"="#A4E804"
)

# 检查是否有未匹配的样本
missing_samples <- setdiff(labels, names(sample_colors))
if(length(missing_samples) > 0) {
  cat("⚠️ Warning: These samples have no defined color and will be black:\n", missing_samples, "\n")
}

# 提取树的数据
dendro_data <- ggdendro::dendro_data(dendro)

# 提取树叶信息
leaf_nodes <- dendro_data$labels
leaf_nodes$color <- sample_colors[leaf_nodes$label]
leaf_nodes$color[is.na(leaf_nodes$color)] <- "black"  # 默认黑色

# 绘制带有颜色实心圆点的树
ggplot() +
  geom_segment(data = dendro_data$segments, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               color = "grey50") +  # 灰色线条表示树的结构
  geom_point(data = leaf_nodes, 
             aes(x = x, y = y, color = color), 
             size = 4, shape = 16) +  # 实心圆点
  geom_text(data = leaf_nodes, 
            aes(x = x, y = y - 0.1, label = label, color = color), 
            hjust = 1, size = 3) +  # 叶子标签
  scale_color_identity() +  # 颜色直接使用样本颜色
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  labs(title = "Hierarchical Clustering with Colored Solid Circles")


###SE在样本间的分布情况
result2=read.table("151SE_01file.txt",head=T)

result2$sum=apply(result2,1,sum)
 data=as.data.frame(result2$sum)
 colnames(data)="value"

pdf("super_Enhancer_01.pdf",width=4,height=4)

# 绘制密度分布图
ggplot(data, aes(x = value)) +
  geom_density(fill = "red", alpha = 0.5) +
  theme_minimal() +
  ggtitle("Density Plot") +
  xlab("Value") +
  ylab("Density")+ theme_bw()+ 
 theme(legend.text=element_text(colour="black", size=12),panel.border=element_blank(),legend.position=c(0.75,0.9),legend.title = element_blank(),
  panel.grid.major.x = element_blank(),
         panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
 panel.grid.major.y = element_blank(),
 axis.line=element_line(color="black"),
 axis.text.x = element_text(vjust=0.5,colour="black",size=12),
         axis.title.x=element_text(colour="black", size=13),
         axis.text.y = element_text(vjust=0.5,colour="black",size=12),
         axis.title.y=element_text(colour="black", size=13))
 dev.off()
 

###CRC中超级增强子的数量以及基因/TF数量   hub  enhancer数量
 #/home/pencode/genome/crcmapper/result3
find -name "*crc_AUTOREG.txt"|while read id ;
 do 
 name=$(basename $id); 
 lujing=${name%%_crc_AUTOREG.txt*};
 cat $id|while read idd; 
 do 
 grep -iw $idd ./$lujing/${lujing}_crc_CANIDATE_TF_AND_SUPER_TABLE.txt >>./$lujing/${lujing}_SE.txt; 
 done; 
 done 
 
 find -name "*_SE.txt"|while read id
> do
> le $id|awk '{print $3"\t"$4"\t"$5"\t""'$(basename $id)'"}'|sort|uniq >> ./total/151_SE.txt
> done

le 151_SE.txt |awk '{print $4}'|sort|uniq -c >151_SE_num.txt
 
 
#find -name "*CRC_SCORES.txt" |while read id; do wc -l  $id >>./total/151_CRC_num.txt; done 
find -name "*crc_AUTOREG.txt" |while read id; 
do 
le $id|awk '{print $0"\t""'$(basename $id)'"}' >>./total/151_TF.txt; 
done

cd ./total

le 151_TF.txt |awk '{print $2}'|sort|uniq -c|sed 's/_crc_AUTOREG.txt//g' >151_TF_num.txt
le 151_TF.txt|awk '{print $1}'|sort|uniq -c >151_TF_share.txt
#dt=read.table("151_TF_num.txt")
> range(dt$V1)
[1]  13 156
> mean(dt$V1)
[1] 58.2053
> median(dt$V1)
[1] 56

dt=read.table("151_SE_num.txt")
> range(dt$V1)
[1]  14 173
> mean(dt$V1)
[1] 74.10596
> median(dt$V1)
[1] 70

dt$V2=sapply(as.character(dt$V2),function(x) {unlist(strsplit(x,"_"))[1]})
dt$V2=gsub("HS4M56b1L","HS4M56b1L_A",dt$V2)
rownames(dt)=dt$V2

info=read.table("/home/pencode/promoter/tissue_id.txt2")
rownames(info)=info$V1

common=intersect(dt$V2,info$V1)
dt=dt[common,]
info2=info[common,]
dt$group=info$V4

library(ggplot2)

pdf("SE_num.pdf",width=20,height=24)
ggplot(dt, aes(x = V1, y = V2)) + 
  geom_col(aes(fill = factor(group)), width = 0.2) +  # 使用填充颜色
  geom_point(aes(size = 3, color = factor(group))) +  # 颜色和大小
  scale_fill_manual(values=c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5",
 "#7A4900", "#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#8CD0FF","#5A0007",
 "#809693","#6A3A4C","#1B4400", "#4FC601","#3B5DFF","#4A3B53", "#FF2F80","#61615A",
 "#BA0900","#6B7900","#00C2A0","#FFAA92","#04F757","#D16100","#DDEFFF","#000035", "#7B4F4B","#A1C299",
 "#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF",
 "#C2FF99", "#001E09","#00489C", "#6F0062","#0CBD66","#EEC3FF","#456D75",
 "#B77B68","#7A87A1","#788D66","#788D66", "#FAD09F", "#FF913F","#D790FF","#922329","#A4E804"),
 breaks=c("Adipose","Adrenal_Gland","Anus","Bladder","Blood_Vessel","Bone_marrow","Brain", "Brain_stem",
 "Bronchi","Bulbourethral_gland","Cartilage","Cerebellum","Choroidal_plexus","Diaphragm","Endocardium_1","Endocardium_2",
 "Esophagus_membrane","Fallopian_tube","Gallbladder","Heart","Heart_valve","Kidney","Large_intestine","Ligament",
 "Liver","Lung","Lymph_node","Mammary_gland","Mesentery","Oralcavity_membrane","Ovary","Palate","Pancreas",
 "Parathyroid_gland","Parotid_gland","Penis","Pituitary","Respiratory_mucosa","Retina","Sciatic_nerve","Seminal_vesicles",
 "Skeletal_muscle","Skin","Small_intestine","Spinal_cord","Spleen","Stomach","Sublingual_gland","Submandibular_gland",
 "Thalamus","Thymus","Thyroid_gland","Tongue","Ureter","Urethra","Uterus","Vagina","Vocal_cords"),
 labels=c("Adipose","Adrenal_Gland","Anus","Bladder","Blood_Vessel","Bone_marrow","Brain","Brain_stem",
 "Bronchi","Bulbourethral_gland","Cartilage","Cerebellum","Choroidal_plexus","Diaphragm","Endocardium_1","Endocardium_2",
 "Esophagus_membrane","Fallopian_tube","Gallbladder","Heart","Heart_valve","Kidney","Large_intestine","Ligament",
 "Liver","Lung","Lymph_node","Mammary_gland","Mesentery","Oralcavity_membrane","Ovary","Palate","Pancreas",
 "Parathyroid_gland","Parotid_gland","Penis","Pituitary","Respiratory_mucosa","Retina","Sciatic_nerve","Seminal_vesicles",
 "Skeletal_muscle","Skin","Small_intestine","Spinal_cord","Spleen","Stomach","Sublingual_gland","Submandibular_gland",
 "Thalamus","Thymus","Thyroid_gland","Tongue","Ureter","Urethra","Uterus","Vagina","Vocal_cords")) +  # 设定自定义颜色
  scale_color_manual(values=c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941","#006FA6","#A30059","#FFDBE5",
 "#7A4900", "#0000A6","#63FFAC","#B79762","#004D43","#8FB0FF","#8CD0FF","#5A0007",
 "#809693","#6A3A4C","#1B4400", "#4FC601","#3B5DFF","#4A3B53", "#FF2F80","#61615A",
 "#BA0900","#6B7900","#00C2A0","#FFAA92","#04F757","#D16100","#DDEFFF","#000035", "#7B4F4B","#A1C299",
 "#300018","#0AA6D8","#013349","#00846F","#372101","#FFB500","#C2FFED","#A079BF",
 "#C2FF99", "#001E09","#00489C", "#6F0062","#0CBD66","#EEC3FF","#456D75",
 "#B77B68","#7A87A1","#788D66","#788D66", "#FAD09F", "#FF913F","#D790FF","#922329","#A4E804"),
 breaks=c("Adipose","Adrenal_Gland","Anus","Bladder","Blood_Vessel","Bone_marrow","Brain", "Brain_stem",
 "Bronchi","Bulbourethral_gland","Cartilage","Cerebellum","Choroidal_plexus","Diaphragm","Endocardium_1","Endocardium_2",
 "Esophagus_membrane","Fallopian_tube","Gallbladder","Heart","Heart_valve","Kidney","Large_intestine","Ligament",
 "Liver","Lung","Lymph_node","Mammary_gland","Mesentery","Oralcavity_membrane","Ovary","Palate","Pancreas",
 "Parathyroid_gland","Parotid_gland","Penis","Pituitary","Respiratory_mucosa","Retina","Sciatic_nerve","Seminal_vesicles",
 "Skeletal_muscle","Skin","Small_intestine","Spinal_cord","Spleen","Stomach","Sublingual_gland","Submandibular_gland",
 "Thalamus","Thymus","Thyroid_gland","Tongue","Ureter","Urethra","Uterus","Vagina","Vocal_cords"),
 labels=c("Adipose","Adrenal_Gland","Anus","Bladder","Blood_Vessel","Bone_marrow","Brain","Brain_stem",
 "Bronchi","Bulbourethral_gland","Cartilage","Cerebellum","Choroidal_plexus","Diaphragm","Endocardium_1","Endocardium_2",
 "Esophagus_membrane","Fallopian_tube","Gallbladder","Heart","Heart_valve","Kidney","Large_intestine","Ligament",
 "Liver","Lung","Lymph_node","Mammary_gland","Mesentery","Oralcavity_membrane","Ovary","Palate","Pancreas",
 "Parathyroid_gland","Parotid_gland","Penis","Pituitary","Respiratory_mucosa","Retina","Sciatic_nerve","Seminal_vesicles",
 "Skeletal_muscle","Skin","Small_intestine","Spinal_cord","Spleen","Stomach","Sublingual_gland","Submandibular_gland",
 "Thalamus","Thymus","Thyroid_gland","Tongue","Ureter","Urethra","Uterus","Vagina","Vocal_cords")) +  
  theme_classic() +
  ylab('') +scale_x_continuous(expand=c(0,0)) +
 theme(axis.title= element_text(size=16,color="black",vjust=0.5, hjust=0.5))+
 theme(axis.text.x=element_text(hjust = 1,vjust =0.5,color ="black",size=16),
 axis.text.y=element_text(size=16,color="black"))
dev.off()   



#####核心TF在样本间share情况

library(ggplot2)
df=read.table("151_TF_share.txt")
df2=as.data.frame(table(df$V1))

m<-c(0,1,5,10,20,40,60,100)#设置一个区间范围
df2$Var1=as.numeric(as.character(df2$Var1))

# 计算百分比
#df2$Percentage <- df2$Freq / sum(df2$Freq) * 100


library(ggplot2)
library(dplyr)

df2$Group <- cut(df2$Var1, breaks = c(0,1,5,10,20,40,76,180), labels = c("sample1", "sample5", "sample10","sample20", "sample40", "sample76", "sample180"), include.lowest = TRUE)

# 按区间汇总频数
df_summary <- df2 %>%
  group_by(Group) %>%
  summarise(Freq = sum(Freq))

# 画饼图
library(ggplot2)
library(RColorBrewer)

pdf("TF_share.pdf",width=4,height=4)
ggplot(df_summary, aes(x = "", y = Freq, fill = Group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = brewer.pal(n = 7, name = "Set1"))
dev.off()   





##不同样本共享的TF结合的SE数量
find -name "*crc_CANIDATE_TF_AND_SUPER_TABLE.txt"|while read id
> do
> le $id|awk '{print $2"\t"$3"\t"$4"\t"$5}'|sed 1d|sort|uniq|awk '{print $1}'|sort|uniq -c >>./total/TF_binding_SE.txt
> done


library(ggplot2)
df=read.table("151_TF_share.txt")
SE=read.table("TF_binding_SE.txt")

N=NULL
for(i in 1:max(df$V1)){

sample1_TF=df[df$V1==i,2]
if(length(sample1_TF)!=0){
sample1_SE=SE[SE$V2 %in% sample1_TF, ,drop=F]
sample1_SE$group=paste0("sample_",as.character(i))
N=rbind(N,sample1_SE)
}
}

library(dplyr)
library(ggplot2)

# 按 group 计算 V1 的均值
N_mean <- N %>%
  group_by(group) %>%
  summarise(V1_mean = mean(V1, na.rm = TRUE))  # 计算均值，忽略 NA 值

# 确保 group 是因子，按原始顺序排列
N_mean$group <- factor(N_mean$group, levels = as.character(unique(N_mean$group)))  

# 生成折线图 + 拟合直线
pdf("TF_share_binding_SE_num_line_fit.pdf", width = 4.5, height = 4)

ggplot(N_mean, aes(x = as.numeric(group), y = V1_mean)) +  # X轴转换为数值型以进行回归
   geom_point(color = "blue", size = 3) +  # 数据点
  geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE) +  # 线性拟合
  ylab("Average TF binding SuperEnhancer number") +  
  xlab("Group (Ordered)") +
  theme_classic() +
 theme(axis.title= element_text(size=16,color="black",vjust=0.5, hjust=0.5))+
 theme(axis.text.x=element_text(hjust = 1,vjust =0.5,color ="black",size=16),
 axis.text.y=element_text(size=16,color="black"))
dev.off()   


POU1F1 垂体 HS4M08q1W

superEnh:chr13    183184616       183234781
gene:PITX1  chr13 183200130       183218073
region:chr13:183181616-183237781 POU1F1_SEenh_track.ini;

##记住bed文件一定。bed结尾
singularity exec  -B /home/pencode/2024H3K27ac/co_activity/bw_file/tissue:/tmpd \
-B /home/pencode/2024H3K27ac/co_activity/figure_Brain12:/tmpd3  \
/home/Singusoft/pygenometracks.sif pyGenomeTracks \
--tracks /tmpd3/POU1F1_SEenh_track.ini --region chr13:183181616-183237781 \
--width 64  --height 16 --fontSize 10 --outFileName \
/tmpd3/POU1F1_SEenh_track.pdf



##track图展示comon TF
find -name "*CANIDATE_TF_AND_SUPER_TABLE.txt"|while read id; do grep -iw fos $id|awk '{print $3"\t"$4"\t"$5}'>> ./fos_SE.bed; done

superEnh:chr7    107828070       108160945
gene:FOS    chr7 107961568       107965187

region:chr7:107825070-108163945 FOS_SEenh_track.ini;

##记住bed文件一定。bed结尾
singularity exec  -B /home/pencode/2024H3K27ac/co_activity/bw_file/tissue:/tmpd \
-B /home/pencode/2024H3K27ac/co_activity/figure_Brain12:/tmpd3  \
/home/Singusoft/pygenometracks.sif pyGenomeTracks \
--tracks /tmpd3/POU1F1_SEenh_track.ini --region chr7:107825070-108163945 \
--width 64  --height 16 --fontSize 10 --outFileName \
/tmpd3/FOS_SEenh_track.pdf


####核心TF、SE、功能/home/pencode/genome/crcmapper/result3/total

awk '{if ($3 > $2) print $0; else print $1 "\t" $3 "\t" $2 "\t" $4}' 151_SE.txt >151_SE_sort.txt

bedtools intersect -a ../../../rose/result/total/151_SE_merge.bed -b 151_SE_sort.txt -wa -wb|awk '{print $1"\t"$2"\t"$3"\t"$7}'|sort|uniq >151sample_mergeSE.txt

data=read.table("151sample_mergeSE.txt")
data$V4=sapply(as.character(data$V4),function(x) {unlist(strsplit(x,"_"))[1]})
data$V4=gsub("HS4M56b1L","HS4M56b1L_A",data$V4)

info=read.table("/home/pencode/promoter/tissue_id.txt2")
rownames(info)=info$V1

for(i in 1:nrow(data)){
data$tissue[i]=as.character(info[as.character(data$V4[i]),4])
}

dat=unique(data[,c(1:3,5)])

 for(i in unique(dat$tissue)){
 tissue=dat[dat$tissue==as.character(i),]
 tissue$name=paste0(tissue$V1,"_",tissue$V2,"_",tissue$V3)
 tissue$value=1
 tissue2=tissue[,(ncol(tissue)-1):ncol(tissue)]
 colnames(tissue2)[2]=as.character(i)
 
  if(!exists("result")){
  result=tissue2
   }else{
  result <- merge(result, tissue2, by="name", all=TRUE)
  }
  }

  
rownames(result)=result$name
result=result[,2:ncol(result)]
result[is.na(result)] <- 0

result=result[,order(colnames(result))]
write.table(result,"151CRC_core_SE.txt",col.names=T,row.names=T,quote=F)
#477

###热图展示
rm(list=ls())
library(ComplexHeatmap)
setwd("G:\\2024H3K27ac流程\\增强子网络")
data <- read.table("151CRC_core_SE.txt", header = T)

#组织器官特异性SE

# 创建一个空的结果数据框来存储结果
result <- data.frame(RowName = character(), MaxCol = character(),  stringsAsFactors = FALSE)

# 对每一行进行操作
for (i in 1:nrow(data)) {
  row_values <- data[i, ]  # 取出每一行的值并保留列名
  sorted_values <- sort(as.numeric(row_values), decreasing = TRUE)  # 将值转换为数值并排序
  
  # 获取最大值和第二大值
  max_value <- sorted_values[1]
  second_max_value <- sorted_values[2]
  
  # 判断最大值是否大于等于第二大值的2倍
  if (max_value > second_max_value) {
    max_col <- names(row_values)[which.max(as.numeric(row_values))]  # 获取最大值所在的列名 
    # 确保列名一致时进行绑定
    result <- rbind(result, data.frame(RowName = rownames(data)[i], MaxCol = max_col,  stringsAsFactors = FALSE))
  }
}

result=result[order(result$MaxCol),]
dat1=data[result$RowName,]
dat2=data[setdiff(rownames(data),result$RowName),]
datt=rbind(dat1,dat2)

#genes <- result$RowName
#genes <- as.data.frame(genes)


datt <- as.matrix(datt) #将表达矩阵转化为matrix

#df <- Heatmap(datt,col = colorRampPalette(c("navy","white","firebrick3"))(100),
        show_row_names = F,cluster_rows = FALSE,cluster_columns = FALSE)
		
#+ rowAnnotation(link = anno_mark(at = which(rownames(datt) %in% genes$genes), 
                                      labels = genes$genes, labels_gp = gpar(fontsize = 10)))

pdf(file="151CRC_core_SE.pdf",width=13,height=5)
df 
dev.off()



####核心TF/home/pencode/genome/crcmapper/result3/total

data=read.table("151_TF.txt")
data$V2=sapply(as.character(data$V2),function(x) {unlist(strsplit(x,"_"))[1]})
data$V2=gsub("HS4M56b1L","HS4M56b1L_A",data$V2)

info=read.table("/home/pencode/promoter/tissue_id.txt2")
rownames(info)=info$V1

for(i in 1:nrow(data)){
data$tissue[i]=as.character(info[as.character(data$V2[i]),4])
}
data=unique(data[,c(1,3)])

 for(i in unique(data$tissue)){
 tissue=data[data$tissue==as.character(i),]
 tissue$name=tissue$V1
 tissue$value=1
 tissue2=tissue[,(ncol(tissue)-1):ncol(tissue)]
 colnames(tissue2)[2]=as.character(i)
 
  if(!exists("result")){
  result=tissue2
   }else{
  result <- merge(result, tissue2, by="name", all=TRUE)
  }
  cat(i,"\n")
  }

rownames(result)=result$name
result=result[,2:ncol(result)]
result[is.na(result)] <- 0

result=result[,order(colnames(result))]
write.table(result,"151CRC_core_TF.txt",col.names=T,row.names=T,quote=F)
#368

###热图展示
rm(list=ls())
library(ComplexHeatmap)
setwd("G:\\2024H3K27ac流程\\增强子网络")
data <- read.table("151CRC_core_TF.txt", header = T)

#组织器官特异性SE

# 创建一个空的结果数据框来存储结果
result <- data.frame(RowName = character(), MaxCol = character(),  stringsAsFactors = FALSE)

# 对每一行进行操作
for (i in 1:nrow(data)) {
  row_values <- data[i, ]  # 取出每一行的值并保留列名
  sorted_values <- sort(as.numeric(row_values), decreasing = TRUE)  # 将值转换为数值并排序
  
  # 获取最大值和第二大值
  max_value <- sorted_values[1]
  second_max_value <- sorted_values[2]
  
  # 判断最大值是否大于等于第二大值的2倍
  if (max_value > second_max_value) {
    max_col <- names(row_values)[which.max(as.numeric(row_values))]  # 获取最大值所在的列名 
    # 确保列名一致时进行绑定
    result <- rbind(result, data.frame(RowName = rownames(data)[i], MaxCol = max_col,  stringsAsFactors = FALSE))
  }
}

result=result[order(result$MaxCol),]
dat1=data[result$RowName,]
dat2=data[setdiff(rownames(data),result$RowName),]
datt=rbind(dat1,dat2)

genes <- result$RowName
genes <- as.data.frame(genes)


datt <- as.matrix(datt) #将表达矩阵转化为matrix

df <- Heatmap(datt,col = colorRampPalette(c("navy","white","firebrick3"))(100),
        show_row_names = F,cluster_rows = FALSE,cluster_columns = FALSE)
		


pdf(file="151CRC_core_TF2.pdf",width=13,height=5)
df + rowAnnotation(link = anno_mark(at = which(rownames(datt) %in% genes$genes), 
                                      labels = genes$genes, labels_gp = gpar(fontsize = 10)))
dev.off()




####/home/pencode/genome/crcmapper/result3/total
 for(i in 1:ncol(dat1)){
 cat(colnames(dat1)[i],"\n")
 print(rownames(dat1[dat1[,i]!=0,]))
 }




data=read.table("151_TF.txt")
data$V2=sapply(as.character(data$V2),function(x) {unlist(strsplit(x,"_"))[1]})
data$V2=gsub("HS4M56b1L","HS4M56b1L_A",data$V2)
info=read.table("/home/pencode/promoter/tissue_id.txt2")
rownames(info)=info$V1
for(i in 1:nrow(data)){
data$tissue[i]=as.character(info[as.character(data$V2[i]),4])
}

for(i in unique(data$tissue)){
tmp=as.data.frame(unique(data[data$tissue==as.character(i),1]))
write.table(tmp,file=paste0(as.character(i),"_TF.txt"),col.names=F,row.names=F,quote=F)
} 
 
 
 
 
 
 #bedtools intersect -a ../../../rose/result/total/sample/HS4M08q1W_Gateway_SuperEnhancers.bed -b FOXO3_HS4M08q1W_crc_motifs.bed  -wa|sort|uniq|wc -l

 #cat HS4M08q1W_crc_CRC_SCORES.txt |grep -iw pitx1|awk -F "\t" '{print $1}'|sed 's/,/\n/g'|sed "s/'//g"|sed 's/\[//g'|sed 's/\]//g'|sort|uniq|wc -l

 
 ls H*|grep ":"|sed 's/://g'|while read id
do
cd 
le /home/pencode/genome/crcmapper/result3/$id/${id}_crc_AUTOREG.txt|while read idd
do
bedtools intersect -a /home/pencode/genome/rose/result/total/sample/${id}_Gateway_SuperEnhancers.bed -b /home/pencode/genome/crcmapper/result3/$id/${idd}_${id}_crc_motifs.bed  -wa|sort|uniq|wc -l >>/home/pencode/genome/crcmapper/result3/$id/out_degree
cat /home/pencode/genome/crcmapper/result3/$id/${id}_crc_CRC_SCORES.txt |grep -iw ${idd}|awk -F "\t" '{print $1}'|sed 's/,/\n/g'|sed "s/'//g"|sed 's/\[//g'|sed 's/\]//g'|sort|uniq|wc -l>>/home/pencode/genome/crcmapper/result3/$id/in_degree
done
done
 
 
###
find -name "*crc_AUTOREG.txt"|while read id
> do
> paste $id ./$(echo $id|awk -F "/" '{print $2}')/in_degree  ./$(echo $id|awk -F "/" '{print $2}')/out_degree >./$(echo $id|awk -F "/" '{print $2}')/total_degree.txt
> done

 
## /home/pencode/genome/crcmapper/result3/HS4M08q1W
paste HS4M08q1W_crc_AUTOREG.txt in_degree out_degree >total_degree.txt


library(ggplot2)
data=read.table("total_degree.txt")
data$label <- ""
for(i in 1:nrow(data)){
if(as.character(data$V1[i]) %in% c("PITX1","PITX2","GLI2","EGR1","FOXO3","HES1","INSM1","ISL2","NR4A2","POU1F1","SIX6","SMAD3","SOX4","TCF4","TCF7L1","THRB")){
data$label[i]=as.character(data$V1[i])
}
}
library(ggplot2)
library(ggrepel)

pdf(file = "pituitary_degree.pdf", width = 4, height = 4)
ggplot(data, aes(V3, V2)) +
  geom_point(aes(color = "blue"), size = 3) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.9, 0.15),
        legend.background = element_blank()) +
  geom_text_repel(aes(color = "red", label = label),  # Specify labels and their colors
                  size = 3,                        # Label size
                  min.segment.length = 0,          # Skip drawing segments shorter than this length
                  seed = 123,                      # Random seed for reproducibility
                  max.overlaps = Inf,              # Remove excessive overlap of text labels
                  box.padding = 0.8,               # Padding around labels
                  arrow = arrow(length = unit(0.015, "npc"), ends = "first"),  # Arrow settings
                  nudge_x = 0, nudge_y = 0,        # Adjust label position
                  show.legend = FALSE)
dev.off()


#/home/pencode/genome/crcmapper/result3/HS4M08q1W
le   HS4M08q1W_crc_CANIDATE_TF_AND_SUPER_TABLE.txt |grep -iE 'pit|GLi'|awk '{print $3"\t"$4"\t"$5"\t"$2}' >PITX_GLI.bed

bedtools intersect -a ../../../rose/result/total/151_SE_merge.bed -b ./PITX_GLI.bed -wa -wb|awk '{print $1"\t"$2"\t"$3"\t"$7}' >PITX_GLI_SE.bed



##
chr15   53435538        53807258  GLI2  SE
chr15   53458742        53725257     gene 
chr15   53492456        53492466 PITX1  motif 
chr15   53509655        53509665 PITX2  motif  
chr15   53510549        53510559  GLI2  motif

region: chr15:53432538-53810258


chr2    147990791       148184212 PITX1  SE
chr2   148013616       148025727
chr2    148008611       148008627 PITX1 motif
chr2    148012875       148012884 PITX2 motif
chr2    148010987       148010997 GLI2  motif

region: chr2:147987791-148187212


chr8    114635805       114711909 PITX2  SE
chr8 114662493       114682557
chr8    114659024       114659040 PITX1 motif
chr8    114665240       114665256 PITX2 motif
chr8    114662581       114662592 GLI2  motif

region: chr8:114632805-114714909



##记住bed文件一定。bed结尾PITX_GLI2_TFBS.bed
singularity exec  -B /home/pencode/2024H3K27ac/co_activity/bw_file/tissue:/tmpd \
-B /home/pencode/2024H3K27ac/co_activity/figure_Brain12:/tmpd3  \
/home/Singusoft/pygenometracks.sif pyGenomeTracks \
--tracks /tmpd3/PITX1_SEenh_track.ini --region chr15:53432538-53810258 \
--width 10  --height 10 --fontSize 10 --outFileName \
/tmpd3/GLI2_SEenh_track.pdf


##提取hub 增强子/home/pencode/genome/crcmapper/result3/total
le  151CRC_core_SE.txt|awk '{print $1}'|sed 1d|sed 's/_/\t/g' > 151CRC_core_SE.bed

bedtools intersect -a  /home/pencode/2024H3K27ac/1Ref_all_peaks_merge_enhancer.bed -b 151CRC_core_SE.bed -wa -wb|awk '{print $1"\t"$2"\t"$3"\t"$4"_"$5"_"$6}' >151CRC_core_SE_enh.bed

进行CTCF的富集 每个SE找到top1

bedtools getfasta -fi /home/pencode/genome/pig_gold/DRC.MotherHap.fa \
-bed 151CRC_core_SE_enh.bed -fo 151CRC_core_SE_enh.fa

/home/pencode/genome/meme-5.5.4/src/fimo \
--oc 151CRC_core_SE_enh_fimo_out \
--parse-genomic-coord /home/pencode/set_script/region/motif_databases/JASPAR/CTCF \
151CRC_core_SE_enh.fa

le  fimo.tsv|awk '{print $3"\t"$4"\t"$5}'|grep chr|sort|uniq >fimo.bed
bedtools intersect -a ../151CRC_core_SE_enh.bed -b fimo.bed -wa|sort|uniq -c  >SE_CTCF.bed
#473SE  535Enh

data=read.table("SE_CTCF.bed")
N=NULL
for(i in unique(data$V5)){
tmp=data[data$V5==as.character(i),]
tmp2=tmp[tmp$V1==max(tmp$V1), ,drop=F]
N=rbind(N,tmp2)
}
write.table(N,"SE_hub_enh.txt",col.names=F,row.names=F,quote=F)

le  SE_hub_enh.txt|awk '{print $2"\t"$3"\t"$4}'|sort|uniq >SE_hub_enh.bed
 bedtools intersect -a SE_hub_enh.bed -b /work/jiangtao/network/result_Brain/tiisue_specific_pro_complex_hub_enhancer.bed -wa |wc -l
 
setwd("G:\\2024H3K27ac流程\\增强子网络")
library(VennDiagram)
library(RColorBrewer)
pdf("hubEnh_overlap.pdf",height=4,width=4.5)
draw.pairwise.venn(  area1 = 1511,  area2 = 535,  cross.area = 8,  
category = c("First", "second"), 
 fill = brewer.pal(5, "Spectral")[1:2],  
 cat.pos = c(0, 0))
dev.off()
 
 
###原因
#SE与enh cluster是不同的定义
#SE内部的enh相关性分布/home/pencode/genome/rose/result/total
bedtools intersect -a /home/pencode/2024H3K27ac/1Ref_all_peaks_merge_enhancer.bed  -b  151_SE_merge.bed -wa -wb |sort|uniq >SE_overlap_enh.bed
#/home/pencode/genome/rose/result/total/network3.r network3.sh

le /work/jiangtao/network/result3/total.txt|awk '{if($3>0.8158592)print $1" "$2" "$3}'|sort|uniq >/work/jiangtao/network/result3/total_thread.txt


enh cluster中 enh是SE的组成占比
le /work/jiangtao/network/result3/total_thread.txt|awk '{print $1}' >/work/jiangtao/network/result3/a
le /work/jiangtao/network/result3/total_thread.txt|awk '{print $2}' >/work/jiangtao/network/result3/b
cat /work/jiangtao/network/result3/a /work/jiangtao/network/result3/b|sort|uniq >/work/jiangtao/network/result3/total_thread_enh

sed -i 's/_/\t/g' /work/jiangtao/network/result3/total_thread_enh
#108046
bedtools intersect -a /work/jiangtao/network/result3/total_thread_enh -b /home/pencode/genome/rose/result/total/151_SE_merge.bed -wa|sort|uniq|wc -l
#90703

1266

98  

  library(ggplot2)

 df <- data.frame(
  value = c(1511, 1266,98),
  group = c("Hub_enhancer", "Hub_enhancer within SEs", "Hub_enhancer_within_SEs_of_CRC")
)
df$group <- factor(df$group, levels = c("Hub_enhancer", "Hub_enhancer within SEs", "Hub_enhancer_within_SEs_of_CRC"))

pdf("three_type_hub_enhancer.pdf",height=4,width=4.5)
ggplot(df, aes(x=as.factor(group),y=value, fill=as.factor(group) )) +  
  geom_bar(stat = "identity",width=0.7 ) +
scale_fill_brewer(palette="Set1")+
ylab("Number of hub Enhancers") + xlab("")+
theme_classic()+
theme(legend.position= "none",   
legend.background=element_rect(fill="transparent",colour = NA),  
legend.text=element_text(colour="black",size=15),
axis.text.x = element_text(size=18),
axis.text.y = element_text(size=18),
axis.title.x = element_text(size=20),  
axis.title.y = element_text(size=20),
plot.margin = unit(c(1, 1, 1, 1), "lines"))
dev.off()  
  














 
 
/work/jiangtao/network/result_Brain/hub_gene ##629
/home/pencode/genome/crcmapper/result3/total/hub_TF   #368
cat  /work/jiangtao/network/result_Brain/hub_gene /home/pencode/genome/crcmapper/result3/total/hub_TF|sort|uniq -c|awk '{if($1>1)print $2}'|wc -l
#21

setwd("G:\\2024H3K27ac流程\\增强子网络")
library(VennDiagram)
library(RColorBrewer)
pdf("hubGene_overlap.pdf",height=4,width=4.5)
draw.pairwise.venn(  area1 = 629,  area2 = 368,  cross.area = 21,  
category = c("First", "second"), 
 fill = brewer.pal(5, "Spectral")[1:2],  
 cat.pos = c(0, 0))
dev.off()


###对这8个进行enh进行k27ac、CTCF密度绘制/home/pencode/genome/crcmapper/result3/total/7enh
chr13   134397319       134402103   chr13_134378981_134412809  21 HS4M66a1W/Oralcavity_membrane
chr13   20331411        20333161      chr13_20328271_20428200     4    HS4M03a1W/Thymus   #  
chr17   53964291        53967447      chr17_53940079_53982205     23   HS4M36b1W/Small_intestine#
chr17   64696693        64699410      chr17_64642406_64733074     23   HS4M23a1R/Adrenal_Gland
chr17   76964203        76967032      chr17_75850911_77188449     31   保守  HS4M08n1T#
chr18   55513381        55518470      chr18_55470314_55520027     29   HS4M39b1W/Large_intestine#
chr3    63943787        63946091       chr3_63855958_63978770       14   HS4M11a1W/Thyroid_gland#
chr3    77257434        77260804       chr3_77103265_77308811        17  HS4M49c1R  /Lung#




split -l 1 7enh.txt ./line_
ls line*|while read id
> do
> le $id|awk '{print $1"\t"$2"\t"$3}' >${id}.bed
> done

#mv line_ac.bed HS4M36b1W.bed



#/home/PersonalFiles/Liusiyi/cutTag/H3K27AC/HS4M03a1W_H3K27AC/visual
#/work/jiahong/2024CTCF/05bw



 #singularity exec -B /home/PersonalFiles/Liusiyi/cutTag/H3K27AC/HS4M66a1W_H3K27AC/mapping:/tmpdisk \
# -B /home/pencode/genome/crcmapper/result3/total/7enh:/tmpdisk2 \
 #/home/Singusoft/143.deeptools.sif bamCoverage -p 10 --binSize 20 -b /tmpdisk/HS4M66a1W_H3K27AC_sorted.bam -o /tmpdisk2/HS4M66a1W_H3K27AC_sorted.bw;



# run compute matrix to collect the data needed for plotting
ls *bed|while read id
do
 singularity exec   -B /home/pencode/genome/crcmapper/result3/total/7enh:/tmpdisk \
 -B /work/jiahong/2024CTCF/05bw:/tmpdisk2 \
 /home/Singusoft/143.deeptools.sif computeMatrix scale-regions \
 -p 10  -b 1000 -a 1000  \
 -R  /tmpdisk/$id \
 -S /tmpdisk/${id%%.bed*}_H3K27AC.bw  /tmpdisk/${id%%.bed*}_CTCF.bw \
--skipZeros -o /tmpdisk/${id%%.bed*}_matrix.mat.gz
done



ls *
singularity exec  -B /home/pencode/genome/rose/result/HS4M42a1L:/tmpdisk /home/Singusoft/143.deeptools.sif  \
plotHeatmap -m /tmpdisk/matrix.mat.gz  \
-out /tmpdisk/SE_TE_region_signal.pdf

 
##对chr13   20331411        20333161      chr13_20328271_20428200     4    HS4M03a1W/Thymus  
#进行深入讨论

 chr13  20101137        20204236 SATB1 gene
chr13_20328271_20428200                SE
chr13   20331411        20333161 hub enh 
chr13   20331547        20332102  RAD21
chr13   20331593        20332448  CTCF
chr13_20196822_20197771          pro
chr13   20197007        20197421 RAD21
chr13   20197004        20197590 CTCF
chr13   20196970        20197717 H3K27ac

region: chr13:20098137-20431200  



##记住bed文件一定。bed结尾SATB1_SEenh_TFBS.bed
singularity exec  -B /home/pencode/genome/crcmapper/result3/total/7enh:/tmpd \
-B /home/pencode/2024H3K27ac/co_activity/figure_Brain12:/tmpd3  \
/home/Singusoft/pygenometracks.sif pyGenomeTracks \
--tracks /tmpd3/SATB1_SEenh_track.ini --region chr13:20098137-20431200 \
--width 40  --height 10 --fontSize 10 --outFileName \
/tmpd3/SATB1_SEenh_track.pdf


##SE区域内的enh
chr13   20328218        20329808   17
chr13   20331411        20333161   18
chr13   20337627        20338060
chr13   20342306        20342590
chr13   20343068        20343581   19
chr13   20344924        20345470
chr13   20347965        20349576   20
chr13   20351795        20352521
chr13   20352578        20352829
chr13   20355446        20356161
chr13   20358711        20359110
chr13   20361687        20362118
chr13   20362569        20362951
chr13   20369451        20369744
chr13   20375936        20376386
chr13   20381040        20381584
chr13   20387580        20388825
chr13   20393221        20393601
chr13   20397013        20398181
chr13   20399868        20400320
chr13   20405054        20405565
chr13   20409383        20410123
chr13   20413146        20413825
chr13   20418507        20419091
chr13   20420782        20421168
chr13   20422425        20422937
chr13   20426561        20427074
chr13   20427637        20429337



##NET 
rm(list=ls())
setwd("G:\\2024H3K27ac流程\\增强子网络")
Allresult=read.table("total_AllResult.txt")
PCSK2=Allresult[Allresult$V5=="SATB1",]

enh=read.table("H3K27ac_58tissue_01file.txt",head=T)

tmp_enh_cluster1=as.character(PCSK2[,1])
tmp_enh_cluster2=as.character(PCSK2[,2])
tmp_enh_cluster=sort(unique(na.omit(c(tmp_enh_cluster1,tmp_enh_cluster2))))
tmp_enh=enh[as.character(tmp_enh_cluster), , drop = FALSE]

###给增强子进行编号
code=as.data.frame(tmp_enh_cluster)
code$code_num=seq(1,nrow(code))
rownames(code)=code$tmp_enh_cluster

target_tissue_enh=rownames(tmp_enh[tmp_enh[,"Thymus"]!=0,])

EEI1=unique(PCSK2[PCSK2$V1 %in% as.character(target_tissue_enh), , drop = FALSE])
EEI2=unique(EEI1[EEI1$V2 %in% as.character(target_tissue_enh), , drop = FALSE])

for(i in 1:nrow(EEI2)){
EEI2$code1[i]=code[as.character(EEI2$V1[i]),2]
EEI2$code2[i]=code[as.character(EEI2$V2[i]),2]
}


custom_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#1F78B4", "#6A3D9A", "#17BECF", "#A50F15"
)



# 绘制网络图 
library(ggraph)
library(tidygraph)

# 创建边列表
edges <- data.frame(
  from = EEI2$code1,
  to = EEI2$code2
)

# 创建节点列表，包含分组信息
nodes <- data.frame(
  name = code$code_num,
  group = paste0("Group",code$code_num)
)

# 转换为 tbl_graph 对象
graph <- tbl_graph(edges = edges, nodes = nodes, directed = FALSE)

# 绘制网络图
pdf("Thymus_SATB1_hub_enh.pdf", height = 6.7, width = 8)
ggraph(graph, layout = "circle") + # 使用圆形布局
  geom_edge_link(color = "gray") +
  geom_node_point(aes(color = group), size =14) + # 根据分组颜色绘制节点
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_color_manual(values = custom_colors)+
    theme_void()

dev.off()


#code
                                tmp_enh_cluster code_num
chr13_20187452_20188632 chr13_20187452_20188632        1
chr13_20189464_20190186 chr13_20189464_20190186        2
chr13_20194873_20195453 chr13_20194873_20195453        3
chr13_20225254_20225968 chr13_20225254_20225968        4
chr13_20227019_20228245 chr13_20227019_20228245        5
chr13_20238074_20240206 chr13_20238074_20240206        6
chr13_20248615_20249230 chr13_20248615_20249230        7
chr13_20251503_20251987 chr13_20251503_20251987        8
chr13_20257662_20258399 chr13_20257662_20258399        9
chr13_20260190_20260817 chr13_20260190_20260817       10
chr13_20268138_20268621 chr13_20268138_20268621       11
chr13_20270219_20270688 chr13_20270219_20270688       12
chr13_20271631_20272175 chr13_20271631_20272175       13
chr13_20278409_20280735 chr13_20278409_20280735       14
chr13_20288291_20289593 chr13_20288291_20289593       15
chr13_20298192_20299543 chr13_20298192_20299543       16
chr13_20328218_20329808 chr13_20328218_20329808       17
chr13_20331411_20333161 chr13_20331411_20333161       18
chr13_20343068_20343581 chr13_20343068_20343581       19
chr13_20347965_20349576 chr13_20347965_20349576       20



 chr13  20101137        20204236 SATB1 gene
chr13_20328271_20428200                SE
chr13   20331411        20333161 hub enh 
chr13   20331547        20332102  RAD21
chr13   20331593        20332448  CTCF
chr13_20196822_20197771          pro
chr13   20197007        20197421 RAD21
chr13   20197004        20197590 CTCF
chr13   20196970        20197717 H3K27ac



/home/pencode/genome/liftOver SATB1_DRCv20.bed /home/pencode/genome/pig_gold/DRCv20_Sus11.minimap2.chain DRCv20_sus11.bed unmap.bed

chr13   20101137        20204236
chr13   20328271        20428200
chr13   20331411        20333161
chr13   20196822        20197771


chr13   5314148 5417610
chr13   5550003 5649324
chr13   5553152 5554918
chr13   5410203 5411152



##mm10 UCSC
chr17	51736429	51833218	chr13:5314149-5417610	gene
chr17	51966021	52081336	chr13:5550004-5649324	SE
#chr4	103564240	103565246	chr13:5550004-5649324	2
chr17	51971115	51972784	chr13:5553153-5554918	hub enh
chr17	51826061	51827012	chr13:5410204-5411152	pro




##mm9  UCSC
chr17:51875754-51972543
chr17:52105346-52220661
#chr4:103236846-103237851
chr17:52110440-52112109
chr17:51965386-51966337

##mm9

singularity exec  -B /home/pencode/genome/crcmapper/result3/total/7enh:/tmpd \
-B /home/pencode/2024H3K27ac/co_activity/figure_Brain12:/tmpd3  \
/home/Singusoft/pygenometracks.sif pyGenomeTracks \
--tracks /tmpd3/SATB1_mm9_SEenh_track.ini --region chr17:51872512-52223661 \
--width 40  --height 10 --fontSize 10 --outFileName \
/tmpd/SATB1_mm9_SEenh_track.pdf

 
#SATB1 peak mm9
chr17   51874075        51875275        MACS_peak_10453 96.99
chr17   51879382        51880181        MACS_peak_10454 63.21
chr17   51884637        51886319        MACS_peak_10455 131.23
chr17   51898897        51900348        MACS_peak_10456 87.44
chr17   51916736        51918758        MACS_peak_10457 138.09
chr17   51936836        51938221        MACS_peak_10458 151.15
chr17   51949077        51957431        MACS_peak_10459 1512.68
chr17   51958857        51960282        MACS_peak_10460 99.22
chr17   51960618        51974272        MACS_peak_10461 1582.66
chr17   51990931        51992432        MACS_peak_10462 133.17
chr17   52015404        52016452        MACS_peak_10463 67.45
chr17   52018076        52023431        MACS_peak_10464 1177.02
chr17   52030383        52033014        MACS_peak_10465 572.69
chr17   52038459        52039869        MACS_peak_10466 253.85
chr17   52042904        52044016        MACS_peak_10467 85.00
chr17   52045122        52047085        MACS_peak_10468 490.39
chr17   52047651        52048532        MACS_peak_10469 51.28
chr17   52095530        52096398        MACS_peak_10470 57.82
chr17   52143835        52144558        MACS_peak_10471 51.31
chr17   52168893        52171440        MACS_peak_10472 279.95
chr17   52171472        52172736        MACS_peak_10473 91.19
chr17   52201032        52203112        MACS_peak_10474 215.00
chr17   52206583        52207735        MACS_peak_10475 80.58

##RAD21 peak mm9转mm10

#mm10                                mm9
chr17	51826217	51826834	chr17:51965543-51966159	1
chr17	51828273	51828597	chr17:51967599-51967922	1
chr17	51830581	51830885	chr17:51969907-51970210	1
chr17	51832367	51832785	chr17:51971693-51972110	1
chr17	51832866	51833094	chr17:51972192-51972419	1
chr17	51833104	51833923	chr17:51972430-51973248	1
chr17	51834035	51834505	chr17:51973361-51973830	1
chr17	51879149	51879686	chr17:52018475-52019011	1
chr17	51882719	51882973	chr17:52022045-52022298	1
chr17	51899453	51899803	chr17:52038779-52039128	1
chr17	51906210	51906867	chr17:52045536-52046192	1
chr17	51956279	51956556	chr17:52095605-52095881	1
chr17	51971506	51971804	chr17:52110832-52111129	1
chr17	52031106	52031370	chr17:52170432-52170695	1
chr17	52064819	52065149	chr17:52204145-52204474	1





##liftover mm10 SATB1 peak
chr17:51734750-51735950
chr17:51740057-51740856
chr17:51745312-51746994
chr17:51759572-51761023
chr17:51777411-51779433
chr17:51797511-51798896
chr17:51809752-51818106
chr17:51819532-51820957
chr17:51821293-51834947
chr17:51851606-51853107
chr17:51876079-51877127
chr17:51878751-51884106
chr17:51891058-51893689
chr17:51899134-51900544
chr17:51903579-51904691
chr17:51905797-51907760
chr17:51908326-51909207
chr17:51956205-51957073
chr17:52004510-52005233
chr17:52029568-52032115
chr17:52032147-52033411
chr17:52061707-52063787
chr17:52067258-52068410



##165
ls *hic|while read id
do
singularity exec -B /home/tmpdir/jiangtao/old_script:/tmpdisk \
/home/guilu/software/hicexplorer_3_7_2.sif hicConvertFormat \
-m /tmpdisk/$id --inputFormat hic \
--outputFormat cool -o /tmpdisk/${id%%.hic*}.cool \
--resolutions 5000
done

#mm10 H3K4me1
#bw file       https://www.encodeproject.org/files/ENCFF726SHR/
#narrowPeak    https://www.encodeproject.org/files/ENCFF695SKH/


#CTCF
https://www.encodeproject.org/files/ENCFF696UII/
https://www.encodeproject.org/files/ENCFF655VZA/

#H3K4me3
https://www.encodeproject.org/files/ENCFF708XKV/
https://www.encodeproject.org/files/ENCFF331TER/


##H3K27AC
https://www.encodeproject.org/files/ENCFF653GDA/
https://www.encodeproject.org/files/ENCFF019NFC/

GSE90635 SATB1
GSE48763 RAD21
CTCF

singularity exec  -B /home/pencode/genome/crcmapper/result3/total/7enh/ren:/tmpd \
-B /home/pencode/2024H3K27ac/co_activity/figure_Brain12:/tmpd3  \
/home/Singusoft/pygenometracks.sif pyGenomeTracks \
--tracks /tmpd3/SATB1_mm10_SEenh_track.ini --region chr17:51731750-52071410 \
--width 40  --height 10 --fontSize 10 --outFileName \
/tmpd/SATB1_mm10_SEenh_track.pdf
