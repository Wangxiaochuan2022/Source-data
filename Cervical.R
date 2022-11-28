

######################    02.差异####
setwd("02.差异分析/")
library('gplots')
library('limma')
library(dplyr)

mRNA=read.table("../00.data/mRNA.symbol.uniq.txt",
                header=TRUE,row.names=1,check.names = FALSE)
par(mfrow=c(1,2))
mRNA <- log2(mRNA+1)
mRNA[1:5,1:5]
mRNA.group <- c(rep("normal_A",3),rep('tomor_N',296)) %>% factor(.,levels = c("normal_A","tomor_N"),ordered = F)
#group <- group[,1] #定义比较组，按照癌症和正常样品数#
mRNA.group <- model.matrix(~factor(mRNA.group))#把group设置成一个model matrix#
mRNA.fit <- lmFit(mRNA,mRNA.group)
mRNA.fit <- eBayes(mRNA.fit)
tempOutput = topTable(mRNA.fit,coef=2,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
mRNA.diff <- nrDEG
write.csv(mRNA.diff, "limmaOut.csv")
foldChange=2
padj=0.05
mRNA.diffSig = mRNA.diff[(mRNA.diff$adj.P.Val < padj & (mRNA.diff$logFC>foldChange | mRNA.diff$logFC<(-foldChange))),]
write.csv(mRNA.diffSig, "diffSig.csv")
mRNA.diffUp = mRNA.diff[(mRNA.diff$adj.P.Val < padj & (mRNA.diff$logFC>foldChange)),]
write.csv(mRNA.diffUp, "diffUp.csv")
mRNA.diffDown = mRNA.diff[(mRNA.diff$adj.P.Val < padj & (mRNA.diff$logFC<(-foldChange))),]
write.csv(mRNA.diffDown, "diffDown.csv")
logFC <-mRNA.diff$logFC
deg.padj <- mRNA.diff$P.Value
data <- data.frame(logFC=logFC,padj=deg.padj)
data$mRNA.diffsig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < foldChange)& data$logFC -foldChange] <- "no"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC >= foldChange] <- "up"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC <= -foldChange] <- "down"
x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)
pdf(file = "volcano.pdf",width=8,height=8)
theme_set(theme_bw())
p<-ggplot(data, aes(logFC,-1*log10(padj), color=mRNA.diffsig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#546de5", "#d2dae2","#ff4757"))+
  geom_vline(xintercept=c(-2,2),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (FDR)")+
  theme_bw()+  theme(plot.title = element_text(hjust = 0.5),
                     legend.position="right",
                     legend.title = element_blank())

p
dev.off()
print(p)

##########   差异热图
library(RColorBrewer)
library(pheatmap)
gene.a <- mRNA[row.names(mRNA.diffSig)[1:100],]
annotation_col=read.table("anno.txt",header=TRUE,
                          row.names=1,check.names = FALSE,sep = '\t')
gene.b <- t(scale(t(gene.a)))
gene.b[which(gene.b > quantile(gene.b,0.95))] <- quantile(gene.b,0.95)
gene.b[which(gene.b < quantile(gene.b,0.05))] <- quantile(gene.b,0.05)
pheatmap(gene.b,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,clustering_method = "ward.D",
         annotation_col =annotation_col,
         show_rownames = F, #不显示行名
         show_colnames = F,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(rev(brewer.pal(n = 5, name ="RdYlBu")))(100),
         scale= "row",border_color = NA, cluster_cols = FALSE)




##################  04.一致性聚类####
setwd("04.一致性聚类/")
library(ConsensusClusterPlus)
TCGA_9=read.table("聚类input.txt",sep="\t",header=T,check.names=F)

TCGA_9=as.matrix(TCGA_9)
rownames(TCGA_9)=TCGA_9[,1]
tmp1 <- TCGA_9[,-1]
Consensusdata1 <- apply(tmp1,2,function(x){as.numeric(x)})
Consensusdata1 <- data.matrix(Consensusdata1)

set.seed(1)
for (hclusti in c('km','pam','hc')) {
  for (disti in c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) {
    title = paste0(disti,hclusti)
    results <- ConsensusClusterPlus(Consensusdata1, maxK = 5,
                                    reps = 100, pItem = 0.95,
                                    pFeature = 0.95,
                                    clusterAlg = hclusti,
                                    distance = disti,
                                    title = title,
                                    plot = "pdf",writeTable=T)
  }
}



###############################km OS####
library(survival)
library(survminer)
rt=read.table("一致性km.input绘图.txt",header=T,sep="\t")
rt$futime=rt$futime/365       #如果以月为单位，除以30；以年为单位，除以365
rt$futime <- ifelse(rt$futime>5,5,rt$futime)
rt$fustat <- ifelse(rt$futime>5,0,rt$fustat)
outTab=data.frame()
for(gene in colnames(rt[,4:ncol(rt)])){
  #a=rt[,gene]
  a=rt[,gene]
  diff=survdiff(Surv(futime, fustat) ~a,data = rt)###对分组进行检验
  pValue=1-pchisq(diff$chisq,df=1)
  outTab=rbind(outTab,cbind(gene=gene,pvalue=pValue))
  #pValue=round(pValue,3)
  pValue=signif(pValue,2)
  #pValue=format(pValue, scientific = TRUE) ##Pvalue科学计数法
  fit <- survfit(Surv(futime, fustat) ~ a, data = rt) ##描绘生存曲线
  summary(fit)
  pdf(file = paste(gene,".survival.pdf"),width = 9,height = 9)
  surplot = ggsurvplot(fit,
                       conf.int = TRUE,
                       risk.table = T,
                       risk.table.title="",
                       risk.table.height = .25,
                       palette = c("red","blue"),
                       xlab="Time (years)",
                       #ylab="Survival rate",
                       break.time.by = 1,
                       pval = pValue,
                       legend.labs =  c(paste(gene," cluster1",sep=""),
                                        paste(gene," cluster2",sep="")),
                       legend.title = gene)
  print(surplot)
  dev.off()
}
write.table(outTab,file="survival.OS.xls",sep="\t",row.names=F,quote=F)




####################    05.PCA+tsne+热图####
setwd("05.聚类效果验证/")
dat <- read.table("pca.input.txt",head = TRUE,row.names=1,sep="\t")
groups <- read.table("pca.input.group.txt",header = TRUE,sep="\t")
library(rgl)
pca <- prcomp(dat)
library(FactoMineR)
library(factoextra)

pca.plot = function(dat,col){
  df.pca <- PCA(dat, graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = TRUE,
               legend.title = "Groups"
  )
}
pdf("pca.pdf",height = 6,width=8)
pca.plot(dat,factor(groups$group))
dev.off()



library(Rtsne)
cluster <- read.csv("tsne.input.csv",header = T,check.names = F,row.names = 1)
#cluster$Cluster <- ifelse(cluster$Cluster==1,'Cluster1',ifelse(cluster$Cluster==2,'Cluster2','Cluster3'))
set.seed(0326)
tsne<- Rtsne(cluster, dims = 2, perplexity=20, verbose=TRUE, max_iter = 1000)
data <- data.frame(tsne$Y,cluster$Cluster)
colnames(data) <- c('Coordinate1','Coordinate2','Cluster')
#png('Fig8.GEO.tSNE.png',width = 600,height = 500)
pdf('tSNE.pdf',width = 8,height = 8)
ggplot(data,aes(Coordinate1,Coordinate2,colour=Cluster))+labs(title = 'TSNE')+
  geom_point(size=3.2)+scale_colour_manual(values = c('salmon','DodgerBlue1',"yellow"))+
  theme_bw()+theme(panel.grid =element_blank())+theme(axis.title = element_text(size=18),
                                                      axis.text = element_text(size=18),
                                                      legend.text=element_text(size=18),
                                                      legend.title = element_text(size=18),
                                                      plot.title = element_text(size = 20,hjust = 0.5))
dev.off()




#################     06.亚型间GSEA####
setwd("06.亚型间GSEA富集/")
tmp1=read.table("mRNA.symbol.uniq.cluster.txt",
                header=TRUE,check.names = FALSE, row.names = 1)
tmp2 <- data.frame(t(tmp1),check.names = F)
tmp3 <- tmp2[order(tmp2$group,decreasing = F),]
tmp4 <- data.frame(t(tmp3),check.names = F)
write.table(tmp4,file="mRNA.symbol.uniq.cluster.sort.tmp.txt",sep = '\t',
            row.names = T,quote = F)

library('gplots')
library('limma')
library(dplyr)

TCGA=read.table("mRNA.symbol.uniq.cluster.sort.txt",
                header=TRUE,row.names=1,check.names = FALSE)
par(mfrow=c(1,2))
#boxplot(data.frame(TCGA),col="blue")    ####画箱式图，比较数据分布情况，数据分布好，则不用进行log2转换
TCGA <- log2(TCGA+1)
TCGA[1:5,1:5]
TCGA.group <- c(rep("low",228),rep('high',68)) %>% factor(.,levels = c("low","high"),ordered = F)
#group <- group[,1] #定义比较组，按照癌症和正常样品数目修改#
TCGA.group <- model.matrix(~factor(TCGA.group))#把group设置成一个model matrix#
TCGA.fit <- lmFit(TCGA,TCGA.group)
TCGA.fit <- eBayes(TCGA.fit)
tempOutput = topTable(TCGA.fit,coef=2,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
TCGA.diff <- nrDEG
write.csv(TCGA.diff, "TCGA.cluster.limmaOut.csv")

library(clusterProfiler)
gene_df <-  read.csv("TCGA.cluster.limmaOut.csv",header=T,sep=",")
gene <- gene_df$SYMBOL

gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- merge(gene_df,gene,by="SYMBOL")
write.table(gene,file="gene对应id.txt",sep = '\t',row.names = F,quote = F)

geneList <- gene_df$logFC
names(geneList)=gene_df$ENTREZID
geneList=sort(geneList,decreasing = T)

GOgmt<-read.gmt("c5.all.v7.5.1.entrez.gmt")
GO <-GSEA(geneList,TERM2GENE = GOgmt)
write.table(GO,file="GO.txt",sep = '\t',row.names = F,quote = F)

library(enrichplot)
library(RColorBrewer)
pdf('GO.top10.pdf',height=10,width=14)
gseaplot2(GO, geneSetID = 1:10,pvalue_table = FALSE,
          rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))
dev.off()

KEGGgmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt")
KEGG <-GSEA(geneList,TERM2GENE = KEGGgmt)
write.table(KEGG,file="KEGG.txt",sep = '\t',row.names = F,quote = F)

pdf('KEGG.top10.pdf',height=10,width=14)
gseaplot2(KEGG, geneSetID = 1:10,pvalue_table = FALSE,
          rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))
dev.off()




###################     07.亚型间临床性状 ####
setwd("07.亚型与临床相关性/")
library("data.table")
clinicalall <- read.csv('TCGA.clinical.csv',header = T,sep = ',')
colnames(clinicalall)[1] <- 'id.small'
heatmap.tmp <- read.csv('热图input.txt',header = T,sep = '\t')
sample_long <- heatmap.tmp$id
heatmap.tmp$id.small <- substr(sample_long,1,12)
heatmap <- merge(heatmap.tmp,clinicalall,by="id.small")

library(pheatmap)
heatmap = heatmap[order(heatmap[,"cluster"]),]
table(heatmap$cluster)

df1 = heatmap[,4:216]
df1 = data.frame(t(df1))

annotation_col = data.frame(age = heatmap$age,
                            gender = heatmap$gender,
                            STAGE = heatmap$STAGE,
                            M = heatmap$M,
                            N = heatmap$N,
                            T = heatmap$T,
                            Risk_Level = factor(c(rep("cluster1",228),
                                                  rep("cluster2",68))))

table(annotation_col$STAGE)
table(annotation_col$M)
table(annotation_col$N)
table(annotation_col$T)

annotation_col$STAGE <- factor(annotation_col$STAGE ,
                               labels=c("Stage I","Stage I","Stage I","Stage II",'Stage II','Stage II',
                                        'Stage II','Stage III','Stage III','Stage III','Stage III',
                                        'Stage IV','Stage IV','Stage IV'),
                               levels = c('Stage I','Stage IA',"Stage IB",'Stage II','Stage IIA','Stage IIB',
                                          'Stage IIC','Stage III','Stage IIIA','Stage IIIB','Stage IIIC',
                                          'Stage IV','Stage IVA','Stage IVB'))

annotation_col$T <- factor(annotation_col$T ,
                           labels=c("T1","T1","T1","T2","T2","T2","T3",'T4',
                                    'T4','T4','TX','Tis'),
                           levels = c("T1","T1a","T1b","T2","T2a","T2b","T3",
                                      "T4","T4a","T4b",'TX','Tis'))

annotation_col$M <- factor(annotation_col$M ,
                           labels=c("M0","M1","M1",'M1','MX'),
                           levels = c("M0","M1","M1a","M1b","MX"))

annotation_col$N <- factor(annotation_col$N ,
                           labels=c("N0","N1","N1",'N1','N1',"N2",'N2','N2',
                                    "N3",'N3','N3','NX'),
                           levels = c("N0","N1","N1a",'N1b','N1c',"N2","N2a",
                                      'N2b',"N3",'N3a','N3b','NX'))

rownames(annotation_col) = colnames(df1)
gene.b <- t(scale(t(df1)))
gene.b[which(gene.b > quantile(gene.b,0.8))] <- quantile(gene.b,0.8)
gene.b[which(gene.b < quantile(gene.b,0.2))] <- quantile(gene.b,0.2)


pdf(file = 'clinical_heatmap.pdf',width = 12,height = 9)
pheatmap(gene.b,
         #color = greenred(75),
         #main = 'heatmap', # 图标题
         scale = 'row', #值集中的方向，“column”，“row” “none”
         annotation_col = annotation_col, #列注释
         #annotation_row = annotation_row, #行注释
         #legend_labels = NA,
         cluster_cols = F,          # 以列聚类
         #cluster_rows = FALSE,         # 以行聚类
         clustering_method = "complete", # 聚类方法 “complete” “average” “median”
         show_rownames = F, #不显示行名
         show_colnames = F, #不显示列名
         #gaps_row = 1169, # 分片
         fontsize = 10,
         angle_col=45) #列名的显示角度
dev.off()

library(table1)
annotation_col1 <- annotation_col
row.names(annotation_col1) <- heatmap$id.small
annotation_col1 <- na.omit(annotation_col1)
annotation_col1$Risk_Level <- factor(annotation_col1$Risk_Level,
                                     levels = c("cluster1","cluster2","P-value"),
                                     labels = c("cluster1","cluster2","P-value"))

## 分类变量标签调整
annotation_col1$age = ifelse(annotation_col1$age > 60,">60","<=60")
#annotation_col1$gender <- factor(annotation_col1$gender)
annotation_col1$STAGE <- factor(annotation_col1$STAGE)
annotation_col1$M <- factor(annotation_col1$M)
annotation_col1$N <- factor(annotation_col1$N)
annotation_col1$T <- factor(annotation_col1$T)


## 左侧标签名调整
labels <- list(
  variables=list(age="age (years)",#gender="gender",
                 STAGE="STAGE",M = "M",
                 N = "N",T = "T"),
  groups=list("", "Risk-Level",""))

## 设置短横线亚组
strata <- c(list(Total=annotation_col1), split(annotation_col1, annotation_col1$Risk_Level))

### 共有三处变量需要修改
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- annotation_col1[[name]]##修改annotation_col1
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.vaild(y ~ annotation_col1$Risk_Level)$p.value##修改annotation_col1
    } else {
      p <- chisq.test(table(y, droplevels(annotation_col1$Risk_Level)))$p.value###修改annotation_col1
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

## 绘制三线表
table1(strata, labels, groupspan=c(1,2,1),droplevels = F,
       render=rndr, render.strat=rndr.strat,
       topclass="Rtable1-zebra")





################    08.免疫相关####
################    08.1 estimate####
setwd("08.亚型间免疫/1.estimate/")
library(estimate)
filterCommonGenes(input.f='../../00.data/mRNA.symbol.uniq.txt',
                  output.f="commonGenes.gct",
                  id="GeneSymbol")
estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct",
              platform="affymetrix")
#输出肿瘤样本的打分表
scores=read.table("estimateScore.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.","\\-",rownames(scores))
out=rbind(ID=colnames(scores),scores[1:nrow(scores),])
write.table(out,file="estimate_scores.tumor.txt",sep="\t",quote=F,col.names=F)

library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
Data <- read.table('StromalScore.txt',sep = '\t',header = T,check.names = F)
colnames(Data)[1] <- 'Type'
#boxdat2 <- gather(boxdat1,gene,expr,2:ncol(boxdat1))
ggboxplot(Data, x= 'Type', y='StromalScore',
          ylab = "StromalScore", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),
                     label = 'p.signif',label.x = 1.5)


Data <- read.table('ImmuneScore.txt',sep = '\t',header = T,check.names = F)
colnames(Data)[1] <- 'Type'
#boxdat2 <- gather(boxdat1,gene,expr,2:ncol(boxdat1))
ggboxplot(Data, x= 'Type', y='ImmuneScore',
          ylab = "ImmuneScore", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),
                     label = 'p.signif',label.x = 1.5)


Data <- read.table('ESTIMATEScore.txt',sep = '\t',header = T,check.names = F)
colnames(Data)[1] <- 'Type'
#boxdat2 <- gather(boxdat1,gene,expr,2:ncol(boxdat1))
ggboxplot(Data, x= 'Type', y='ESTIMATEScore',
          ylab = "ESTIMATEScore", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),
                     label = 'p.signif',label.x = 1.5)


Data <- read.table('TumorPurity.txt',sep = '\t',header = T,check.names = F)
colnames(Data)[1] <- 'Type'
#boxdat2 <- gather(boxdat1,gene,expr,2:ncol(boxdat1))
ggboxplot(Data, x= 'Type', y='TumorPurity',
          ylab = "TumorPurity", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),
                     label = 'p.signif',label.x = 1.5)




####################   08.2 ssgsea   ####
setwd("08.亚型间免疫/2.免疫浸润/")
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
exprSet1 <- read.table("../../06.亚型间GSEA富集/mRNA.symbol.uniq.cluster.sort.txt",
                       header=T,sep="\t",row.names=1,check.names = F)
gene_set<- read.table('mmc3.txt',
                      sep = '\t',header = T)
##读取已经下载好的免疫细胞和对应基因列表，来源见文献附件
gene_set<-gene_set[, 1:2]#选取特异基因和对应的免疫细胞两行
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
gsva_matrix <- gsva(as.matrix(exprSet1), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.table(gsva_matrix,file = 'ssGSEA推算.csv',quote = F,sep = ',')

library(RColorBrewer)
library(pheatmap)
annotation_col=read.table("../../06.亚型间GSEA富集/anno.txt",header=TRUE,
                          row.names=1,check.names = FALSE,sep = '\t')

pheatmap(gsva_matrix,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,clustering_method = "ward.D",
         annotation_col =annotation_col,
         show_rownames = T, #不显示行名
         show_colnames = F,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
         scale= "row",border_color = NA, cluster_cols = FALSE)


library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
C <- read.table('box.txt',sep = '\t',header = T,check.names = F)
#C[,2:ncol(C)] <- apply(C[,2:ncol(C)],2,function(x){log2(x+1)})
colnames(C)[1] <- 'Type'
C1 <- gather(C,gene,expr,2:ncol(C))
ggboxplot(C1, x= 'gene', y='expr',
          ylab = "proportion", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),label = 'p.signif',label.x = 1.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))####  横坐标倾斜45度



################    08.3 GSVA####
library('gplots')
library('limma')
library(dplyr)
setwd("08.亚型间免疫/3.GSVA/")
library(msigdbr)
library(dplyr)
library(data.table)
library(GSVA)
library(limma)
library(stringr)
library(ggplot2)
library(GSEABase)
msigdbr_show_species()
gmt<-getGmt("c7.all.v2022.1.Hs.symbols.gmt")
gsym.expr <- read.table("../../06.亚型间GSEA富集/mRNA.symbol.uniq.cluster.sort.txt", header=T, row.names = 1)
gsva_es <- gsva(as.matrix(gsym.expr), gmt)
write.csv(gsva_es, "gsva.go_output.csv", quote = F)


group_list <- data.frame(sample = colnames(gsva_es), 
                         group = c(rep("cluster1", 228), rep("cluster2", 68)))
design <- model.matrix(~ 0 + factor(group_list$group))
colnames(design) <- levels(factor(group_list$group))
rownames(design) <- colnames(gsva_es)
design

contrast.matrix <- makeContrasts(cluster2-cluster1, levels = design)
# 差异分析，low vs. high
fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
x <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(x)
write.csv(x, "gsva.go_limma.csv", quote = F)

#输出t值，用做FigureYa39bar的输入数据
pathway <- row.names(x)
df <- data.frame(ID = pathway, score = x$t)
df <- df[1:30,]
write.csv(df, "go.easy_input2_for39bar.top30.csv", quote = F, row.names = F)

df <- read.csv("go.easy_input2_for39bar.top30.csv")
cutoff <- 2
df$group <- cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3))

#按照score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)

ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') +
  scale_fill_manual(values = c('palegreen3', 'snow3', 'dodgerblue4'), guide = FALSE) +
  
  #画2条虚线
  geom_hline(yintercept = c(-cutoff,cutoff),
             color="white",
             linetype = 2, #画虚线
             size = 0.6) + #线的粗细
  
  #写label
  geom_text(data = subset(df, score < 0),
            aes(x=ID, y= 0, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 2.5, #字的大小
            hjust = "bottom",angle=90) +  #字的对齐方式
  geom_text(data = subset(df, score > 0),
            aes(x=ID, y= -0.1, label=paste0(" ", ID), color = group),
            size = 2.5, hjust = "bottom",angle=270) +
  scale_colour_manual(values = c("black","snow3","black"), guide = FALSE) +
  
  xlab("") +ylab("t value of GSVA score \n cluster1 versus cluster2")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_rect(size = 1)) + #边框粗细
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_blank())+ #去除x轴
  theme(axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15))
dev.off()


library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
C <- read.table('box.txt',sep = '\t',header = T,check.names = F)
#C[,2:ncol(C)] <- apply(C[,2:ncol(C)],2,function(x){log2(x+1)})
colnames(C)[1] <- 'Type'
C1 <- gather(C,gene,expr,2:ncol(C))
ggboxplot(C1, x= 'gene', y='expr',
          ylab = "gsva score", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),label = 'p.signif',label.x = 1.5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))####  横坐标倾斜45度




##################      09.免疫治疗####
setwd("09.亚型间免疫治疗/")
library(pheatmap)
dat <- read.table("../06.亚型间GSEA富集/mRNA.symbol.uniq.cluster.sort.txt",
                  sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
ann <- read.table("anno.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)

TIDE <- dat
# 这里为了得到比较好的结果，采用two direction median centered
TIDE <- sweep(TIDE,2, apply(TIDE,2,median,na.rm=T))
TIDE <- sweep(TIDE,1, apply(TIDE,1,median,na.rm=T))
write.table(TIDE,"TIDE_input.self_subtract",sep = "\t",row.names = T,col.names = NA,quote = F)


############## 在网站上得到TIDE output之后，将TIDE分数画成箱线图
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
boxdat33 <- read.table('tide.box.txt',sep = '\t' ,header = T,check.names = F)
colnames(boxdat33)[1] <- 'Type'
pdf("TIDEscore.box.pdf",height=6,width=6)
ggboxplot(boxdat33, x= 'Type', y='TIDE',
          ylab = "TIDE value", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),label = 'p.signif',label.x = 1.5)
dev.off()



TIDE.res <- read.csv("tide.out.csv",header = T,row.names = 1,check.names = F,stringsAsFactors = F)
ann$TIDE <- TIDE.res[rownames(ann),"Responder"]
print(table(ann$TIDE,ann$risk))
print(fisher.test(table(ann$TIDE,ann$risk)))

generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct,
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}

# 创建submap需要的数据格式 (SKCM)
skcm.immunotherapy.logNC <- read.table("../../../学习资料/skcm.immunotherapy.47samples.log2CountsNorm.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- toupper(rownames(skcm.immunotherapy.logNC)) # 基因大写，因为我使用的数据是把基因名都大写的
skcm.immunotherapy.info <- read.table("../../../学习资料/skcm.immunotherapy.47sampleInfo.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

# 创建submap需要的数据格式 (TCGA)
tmp <- read.table("../06.亚型间GSEA富集/mRNA.symbol.uniq.cluster.sort.txt",
                  sep = "\t",check.names = F,stringsAsFactors = F,header = T,row.names = 1)
# submap不允许出现flat value, 因此最好选取过滤掉低表达的表达谱，这里使用的数据过滤了超过90%样本表达值均<1的基因
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 取交集后的基因列表

sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生输出数据的文件名
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

# 提出亚型的样本，顺序排列
samples.low <- rownames(ann[which(ann$risk == "cluster1"),])
samples.high <- rownames(ann[which(ann$risk == "cluster2"),])

sam_info <- data.frame("ImmClust"=c(samples.low,samples.high),row.names = c(samples.low,samples.high))
sam_info$rank <- rep(c(1,2),times=c(length(samples.low),length(samples.high))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

# 产生输出数据的文件名
gct_file <- "Immune2.for.SubMap.gct"
cls_file <- "Immune2.for.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) # 产生和示例数据类似的形式，log2转化的标准化count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

####   画图
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
cherry    <- "#700353"
lightgrey <- "#dcddde"

# 把submap结果/130452/SubMap_SubMapResult.txt文件中的值填入相应的位置
# 输入文件中的名义p值和校正p值绘制热图
tmp <- matrix(c(0.08191808, 0.6073926, 0.9820180, 0.5024975,
                0.96503497, 0.1288711, 0.2267732, 0.6513487,
                1,0.6553447,1,1,1,1,1,1), # Bonferroni校正p值
              nrow = 4,byrow = T,dimnames = list(c("low_p","high_p","low_b","high_b"),
                                                 c("CTAL4-noR","CTLA4-R","PD1-noR","PD1-R")))

pheatmap(tmp, cellwidth = 30, cellheight = 30,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[5:1],
         display_numbers = matrix(ifelse(tmp < 0.05,'p < 0.05',''),nrow(tmp)),number_format = "%.3f",
         gaps_row = 2,
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         filename = "heatmap_submap.pdf")









########################   循环.data ####
setwd("00.data/GEO宫颈癌/")
tmp=read.table("GSE44001_series_matrix.id.txt",header=T,sep="\t",check.names = F)
tmp_uniq <- aggregate(.~ID_REF,tmp,max)
write.table(tmp_uniq,file="GSE44001_series_matrix.id.uniq.txt",
            sep = '\t',row.names = F,quote = F)


setwd("00.data/GEO头颈鳞癌/")
tmp=read.table("GSE65858_series_matrix.id.txt",header=T,sep="\t",check.names = F)
tmp_uniq <- aggregate(.~id,tmp,max)
write.table(tmp_uniq,file="GSE65858_series_matrix.id.uniq.txt",
            sep = '\t',row.names = F,quote = F)







###############   确定种子####
######################    10.单因素####
############################  分群体 7:3
setwd("10.单因素")
all.cox <- read.table("TCGA.主.input.txt",header=T,row.names=1)
#all.cox[,3:ncol(all.cox)] <- apply(all.cox[,3:ncol(all.cox)],2,function(x){log2(x+1)})
all.cox$group <- "case"
set.seed(23923)
index <- caret::createDataPartition(all.cox[,"group"], p =0.7)
cox_train <- all.cox[index$Resample1,]
cox_train <- cox_train[,-ncol(cox_train)]
cox_test<- all.cox[-index$Resample1,]
cox_test <- cox_test[,-ncol(cox_test)]
write.table(cox_test,"test.clinical.fpkm.txt",sep="\t",row.names=T)

pFilter=0.05                                                     #定义单因素显著性
library(survival)                                                 #引用包
library(UpSetR)
outTab=data.frame()
for(i in colnames(cox_train[,3:ncol(cox_train)])){
  cox <- coxph(Surv(futime, fustat) ~ cox_train[,i], data = cox_train)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
  )
}
#输出所有单因素的结果
outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="uniCoxResult.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著的结果
sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<pFilter,]
write.table(sigTab,file="uniCoxResult.Sig.txt",sep="\t",row.names=F,quote=F)
#输出单因素显著AS的PSI值，用于后续建模
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=cox_train[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExp.txt",sep="\t",row.names=F,quote=F)

#####    单因素森林图
univariateCox_Result=read.table("uniCoxResult.Sig.txt",header=T,sep="\t",row.names=1,check.names=F)
p.value<-signif(univariateCox_Result[,5], digits=2)
#wald.test<-signif(x$wald["test"], digits=2)
Coefficient<-signif(univariateCox_Result[,1], digits=2);#coeficient beta
HR <-signif(univariateCox_Result[,2], digits=2);#exp(beta)
HR.confint.lower <- signif(univariateCox_Result[,3], digits=2)
HR.confint.upper <- signif(univariateCox_Result[,4], digits=2)
HR.combine <- paste0(HR, " (",
                     HR.confint.lower, "-", HR.confint.upper, ")")
rescox.temp.1<-cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
names(rescox.temp.1)<-c("HR", "HR.confint.lower", "HR.confint.upper",'HR.combine',
                        "p.value")
rownames(rescox.temp.1) <- rownames(univariateCox_Result)

univOut_sig.plot.univar <- rescox.temp.1
gene <- rownames(univOut_sig.plot.univar)
hr <- univOut_sig.plot.univar[,1]
hrLow <- univOut_sig.plot.univar[,2]
hrHigh <- univOut_sig.plot.univar[,3]
Hazard.ratio <- univOut_sig.plot.univar[,4]
pVal <- univOut_sig.plot.univar[,5]

pdf(file='单因素cox森林图.pdf', width = 8,height = 6)
n <- nrow(univOut_sig.plot.univar)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=1.3
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5+0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
highlim <- max(as.numeric(hrHigh))+0.1
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(max(c(0,min(as.numeric(hrLow))))-0.1,highlim)
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",cex.lab=2)
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)
dev.off()




#######################    11.lasso####
library("glmnet")
library(survivalROC)
library(survminer)
library(survival)
setwd('11.lasso/')
result <- list()
for (sed in c(23923)) {
  cox.data <- read.table("../10.单因素/uniSigExp.txt",header=T,row.names=1)
  colnames(cox.data)[1:2] <- c('Time','Status')
  x=as.matrix(cox.data[,colnames(cox.data)[-(1:2)]])
  y=data.matrix(Surv(cox.data$Time,cox.data$Status))
  set.seed(23923)
  cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
  fit.train <- cvfit$glmnet.fit
  coef <- coef(cvfit, s = cvfit$lambda.min)
  pdf(file = "lasso.Binomial.Deviance.pdf",height = 9,width = 9)
  plot(cvfit,xlab='Log Lambda',cex.lab = 2.4)+
    text(x = log(cvfit$lambda.min),y = 11,
         paste('Lambda.min\n',round(cvfit$lambda.min,3)),cex=2,adj=0.5)+
    text(x = log(cvfit$lambda.1se),y = 11.5,
         paste('Lambda.lse\n',round(cvfit$lambda.1se,3)),cex=2)
  dev.off()
  pdf(file = "lasso.coefficients.venalty.pdf",height = 9,width = 9)
  par(mgp = c(4,1,0),mai=c(2,2,1,1))
  plot(fit.train, xvar="lambda",cex.lab = 3)+
    abline(v = c(log(cvfit$lambda.min), log(cvfit$lambda.1se)),lty=2)+
    text(x = log(cvfit$lambda.min),y = 0.5,
         paste('Lambda.min\n',round(cvfit$lambda.min,4)),cex=1.6,adj=0.9)+
    text(x = log(cvfit$lambda.1se),y = 1,
         paste('Lambda.lse\n',round(cvfit$lambda.1se,4)),cex=1.6,adj=0.9)
  dev.off()
  index <- which(as.numeric(coef) != 0)
  actCoef <- coef[index]
  lassoGene=row.names(coef)[index]
  lassoGene=c("Status","Time",lassoGene)
  cox.data=cox.data[,lassoGene]
  if(ncol(cox.data) < 3) next
  cox.data.step <- na.omit(cox.data[,c('Time','Status',colnames(cox.data)[-(1:2)])])
  fml <- as.formula(paste0('Surv(Time,Status)~',paste0(colnames(cox.data)[-(1:2)],collapse = '+')))
  f <- coxph(fml, data=cox.data.step,id = rownames(cox.data.step))
  cox=f
  riskScore=predict(cox,type="risk",newdata=cox.data.step)
  cox.data.plot <- cbind(cox.data.step,riskScore)
  cox.data.plot$Status <- as.numeric(as.character(cox.data.plot$Status))
  cox.data.plot$riskScore <- as.numeric(cox.data.plot$riskScore)
  cox.data.plot <- cox.data.plot[order(cox.data.plot$riskScore),]
  write.table(cox.data.plot,file="risk.train.txt",sep="\t",quote=F,row.names=T,col.names=T)
  train <- c()
  for(i in c(1:1)){
    roc_test <- survivalROC(Stime = cox.data.plot$Time,
                            status = cox.data.plot$Status,
                            marker = cox.data.plot$riskScore,
                            predict.time = 365*i,method = 'KM')
    train <- c(train,roc_test$AUC)
    pdf(file="135.train.ROC.pdf",height = 6,width = 6)
    plot(roc_test$FP, roc_test$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
         xlab="False positive rate", ylab="True positive rate",
         main="ROC curve",
         lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
    abline(0,1)
    aucText=c()
    rocCol <- c('#FA8072','#63B8FF','#FFC1C1','#ADFF2F','#FFFF00')
    aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc_test$AUC),")"))
    j =0
    for (i in c(3,5)){
      roc1=survivalROC(Stime=cox.data.plot$Time, status=cox.data.plot$Status,
                       marker = cox.data.plot$riskScore,predict.time =365*i,method = 'KM')
      j=j+1
      aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
      lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
    }
    legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
    abline(0,1)
    dev.off()
  }
  ####  KM
  bestthreshold.surv <- surv_cutpoint(cox.data.plot,time = "Time",
                                      event = "Status",
                                      'riskScore',
                                      minprop = 0.4,
                                      progressbar = TRUE)
  bestthreshold <- bestthreshold.surv$cutpoint$cutpoint
  bestthreshold.train <- bestthreshold
  # 若只需要对高低风险作KM曲线，则执行以下代码：
  #survdiff(Surv(Time,Status) ~ risk,data = cox.data.plot) # 单因素变量分析计算P值
  group <- ifelse(cox.data.plot$riskScore > bestthreshold, 'high','low')
  cox.data.plot$Time=cox.data.plot$Time/365
  fit <- survfit(Surv(Time,Status) ~ group, data = cox.data.plot)
  survdf <- survdiff(Surv(Time,Status) ~ group, data = cox.data.plot)
  p.value <- 1 - pchisq(survdf$chisq, length(survdf$n) -1)
  train <- c(train,p.value)
  write.table(cbind(cox.data.plot,group),file="risk.train.group.txt",sep="\t",quote=F,row.names=T)
  
  pdf(file="survival.coefficient.train.pdf",height = 9,width = 9)
  ggsurvplot(fit,
             conf.int=TRUE,
             pval=TRUE,
             risk.table=TRUE,
             #legend.labs=legend.labs,
             legend.title="Risk",
             #linetype = lty,
             palette = c( "red", "blue"),
             title="Kaplan-Meier Curve for Survival",
             risk.table.height=.15) %>% print()
  dev.off()
  
  
  ###########   test   ####
  cox.data <- read.table("../10.单因素/test.clinical.fpkm.txt",header=T,row.names=1)
  colnames(cox.data)[1:2] <- c('Time','Status')
  cox.data=cox.data[,lassoGene]
  cox.data.step <- na.omit(cox.data[,c('Time','Status',colnames(cox.data)[-(1:2)])])
  riskScore=predict(cox,type="risk",newdata=cox.data.step)
  cox.data.plot <- cbind(cox.data.step,riskScore)
  cox.data.plot$Status <- as.numeric(as.character(cox.data.plot$Status))
  cox.data.plot$riskScore <- as.numeric(cox.data.plot$riskScore)
  cox.data.plot <- cox.data.plot[order(cox.data.plot$riskScore),]
  write.table(cox.data.plot,file="risk.test.txt",sep="\t",quote=F,row.names=T,col.names=T)
  test <- c()
  for(i in c(1)){
    roc_test <- survivalROC(Stime = cox.data.plot$Time,
                            status = cox.data.plot$Status,
                            marker = cox.data.plot$riskScore,
                            predict.time = 365*i,method = 'KM')
    test <- c(test,roc_test$AUC)
    pdf(file="135.test.ROC.pdf",height = 6,width = 6)
    plot(roc_test$FP, roc_test$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
         xlab="False positive rate", ylab="True positive rate",
         main="ROC curve",
         lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
    abline(0,1)
    aucText=c()
    rocCol <- c('#FA8072','#63B8FF','#FFC1C1','#ADFF2F','#FFFF00')
    aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc_test$AUC),")"))
    j =0
    for (i in c(3,5)){
      roc1=survivalROC(Stime=cox.data.plot$Time, status=cox.data.plot$Status,
                       marker = cox.data.plot$riskScore,predict.time =365*i,method = 'KM')
      j=j+1
      aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
      lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
    }
    legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
    abline(0,1)
    dev.off()
  }
  ####  KM
  bestthreshold.surv <- surv_cutpoint(cox.data.plot,time = "Time",
                                      event = "Status",
                                      'riskScore',
                                      minprop = 0.4,
                                      progressbar = TRUE)
  bestthreshold <- bestthreshold.surv$cutpoint$cutpoint
  bestthreshold.test <- bestthreshold
  # 若只需要对高低风险作KM曲线，则执行以下代码：
  #survdiff(Surv(Time,Status) ~ risk,data = cox.data.plot) # 单因素变量分析计算P值
  group <- ifelse(cox.data.plot$riskScore > bestthreshold, 'high','low')
  cox.data.plot$Time=cox.data.plot$Time/365
  fit <- survfit(Surv(Time,Status) ~ group, data = cox.data.plot)
  survdf <- survdiff(Surv(Time,Status) ~ group, data = cox.data.plot)
  p.value <- 1 - pchisq(survdf$chisq, length(survdf$n) -1)
  test <- c(test,p.value)
  write.table(cbind(cox.data.plot,group),file="risk.test.group.txt",sep="\t",quote=F,row.names=T)
  
  pdf(file="survival.coefficient.test.pdf",height = 9,width = 9)
  ggsurvplot(fit,
             conf.int=TRUE,
             pval=TRUE,
             risk.table=TRUE,
             #legend.labs=legend.labs,
             legend.title="Risk",
             #linetype = lty,
             palette = c( "red", "blue"),
             title="Kaplan-Meier Curve for Survival",
             risk.table.height=.15) %>% print()
  dev.off()
  
  ###########   valid1  ####
  cox.data <- read.table("GEO宫颈癌input.txt",header=T,row.names=1)
  cox.data[,3:ncol(cox.data)] <- apply(cox.data[,3:ncol(cox.data)],2,function(x){log2(x+1)})
  colnames(cox.data)[1:2] <- c('Time','Status')
  cox.data=cox.data[,lassoGene]
  cox.data.step <- na.omit(cox.data[,c('Time','Status',colnames(cox.data)[-(1:2)])])
  riskScore=predict(cox,type="risk",newdata=cox.data.step)
  cox.data.plot <- cbind(cox.data.step,riskScore)
  cox.data.plot$Status <- as.numeric(as.character(cox.data.plot$Status))
  cox.data.plot$riskScore <- as.numeric(cox.data.plot$riskScore)
  cox.data.plot <- cox.data.plot[order(cox.data.plot$riskScore),]
  write.table(cox.data.plot,file="risk.valid宫颈癌.txt",sep="\t",quote=F,row.names=T,col.names=T)
  valid <- c()
  for(i in c(1)){
    roc_valid <- survivalROC(Stime = cox.data.plot$Time,
                             status = cox.data.plot$Status,
                             marker = cox.data.plot$riskScore,
                             predict.time = 365*i,method = 'KM')
    valid <- c(valid,roc_valid$AUC)
    pdf(file="135.valid宫颈癌.ROC.pdf",height = 6,width = 6)
    plot(roc_valid$FP, roc_valid$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
         xlab="False positive rate", ylab="True positive rate",
         main="ROC curve",
         lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
    abline(0,1)
    aucText=c()
    rocCol <- c('#FA8072','#63B8FF','#FFC1C1','#ADFF2F','#FFFF00')
    aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc_valid$AUC),")"))
    j =0
    for (i in c(3,5)){
      roc1=survivalROC(Stime=cox.data.plot$Time, status=cox.data.plot$Status,
                       marker = cox.data.plot$riskScore,predict.time =365*i, method="KM")
      j=j+1
      aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
      lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
    }
    legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
    abline(0,1)
    dev.off()
  }
  ####  KM
  bestthreshold.surv <- surv_cutpoint(cox.data.plot,time = "Time",
                                      event = "Status",
                                      'riskScore',
                                      minprop = 0.4,
                                      progressbar = TRUE)
  bestthreshold <- bestthreshold.surv$cutpoint$cutpoint
  bestthreshold.vaild1 <- bestthreshold
  # 若只需要对高低风险作KM曲线，则执行以下代码：
  #survdiff(Surv(Time,Status) ~ risk,data = cox.data.plot) # 单因素变量分析计算P值
  group <- ifelse(cox.data.plot$riskScore > bestthreshold, 'high','low')
  cox.data.plot$Time=cox.data.plot$Time/365
  fit <- survfit(Surv(Time,Status) ~ group, data = cox.data.plot)
  survdf <- survdiff(Surv(Time,Status) ~ group, data = cox.data.plot)
  p.value <- 1 - pchisq(survdf$chisq, length(survdf$n) -1)
  test <- c(test,p.value)
  write.table(cbind(cox.data.plot,group),file="risk.valid.group宫颈癌.txt",sep="\t",quote=F,row.names=T)
  
  pdf(file="survival.coefficient.valid宫颈癌.pdf",height = 9,width = 9)
  ggsurvplot(fit,
             conf.int=TRUE,
             pval=TRUE,
             risk.table=TRUE,
             #legend.labs=legend.labs,
             legend.title="Risk",
             #linetype = lty,
             palette = c( "red", "blue"),
             title="Kaplan-Meier Curve for Survival",
             risk.table.height=.15) %>% print()
  dev.off()
  
  ###########   valid2  ####
  cox.data <- read.table("GEO头颈鳞癌input.txt",header=T,row.names=1)
  cox.data[,3:ncol(cox.data)] <- apply(cox.data[,3:ncol(cox.data)],2,function(x){log2(x+1)})
  colnames(cox.data)[1:2] <- c('Time','Status')
  cox.data=cox.data[,lassoGene]
  cox.data.step <- na.omit(cox.data[,c('Time','Status',colnames(cox.data)[-(1:2)])])
  riskScore=predict(cox,type="risk",newdata=cox.data.step)
  cox.data.plot <- cbind(cox.data.step,riskScore)
  cox.data.plot$Status <- as.numeric(as.character(cox.data.plot$Status))
  cox.data.plot$riskScore <- as.numeric(cox.data.plot$riskScore)
  cox.data.plot <- cox.data.plot[order(cox.data.plot$riskScore),]
  write.table(cox.data.plot,file="risk.valid头颈鳞癌.txt",sep="\t",quote=F,row.names=T,col.names=T)
  valid <- c()
  for(i in c(1)){
    roc_valid <- survivalROC(Stime = cox.data.plot$Time,
                             status = cox.data.plot$Status,
                             marker = cox.data.plot$riskScore,
                             predict.time = 365*i,method = 'KM')
    valid <- c(valid,roc_valid$AUC)
    pdf(file="135.valid头颈鳞癌.ROC.pdf",height = 6,width = 6)
    plot(roc_valid$FP, roc_valid$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
         xlab="False positive rate", ylab="True positive rate",
         main="ROC curve",
         lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
    abline(0,1)
    aucText=c()
    rocCol <- c('#FA8072','#63B8FF','#FFC1C1','#ADFF2F','#FFFF00')
    aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",roc_valid$AUC),")"))
    j =0
    for (i in c(3,5)){
      roc1=survivalROC(Stime=cox.data.plot$Time, status=cox.data.plot$Status,
                       marker = cox.data.plot$riskScore,predict.time =365*i, method="KM")
      j=j+1
      aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",roc1$AUC),")"))
      lines(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
    }
    legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
    abline(0,1)
    dev.off()
  }
  ####  KM
  bestthreshold.surv <- surv_cutpoint(cox.data.plot,time = "Time",
                                      event = "Status",
                                      'riskScore',
                                      minprop = 0.4,
                                      progressbar = TRUE)
  bestthreshold <- bestthreshold.surv$cutpoint$cutpoint
  bestthreshold.vaild2 <- bestthreshold
  # 若只需要对高低风险作KM曲线，则执行以下代码：
  #survdiff(Surv(Time,Status) ~ risk,data = cox.data.plot) # 单因素变量分析计算P值
  group <- ifelse(cox.data.plot$riskScore > bestthreshold, 'high','low')
  cox.data.plot$Time=cox.data.plot$Time/365
  fit <- survfit(Surv(Time,Status) ~ group, data = cox.data.plot)
  survdf <- survdiff(Surv(Time,Status) ~ group, data = cox.data.plot)
  p.value <- 1 - pchisq(survdf$chisq, length(survdf$n) -1)
  test <- c(test,p.value)
  write.table(cbind(cox.data.plot,group),file="risk.valid.group头颈鳞癌.txt",sep="\t",quote=F,row.names=T)
  
  pdf(file="survival.coefficient.valid头颈鳞癌.pdf",height = 9,width = 9)
  ggsurvplot(fit,
             conf.int=TRUE,
             pval=TRUE,
             risk.table=TRUE,
             #legend.labs=legend.labs,
             legend.title="Risk",
             #linetype = lty,
             palette = c( "red", "blue"),
             title="Kaplan-Meier Curve for Survival",
             risk.table.height=.15) %>% print()
  dev.off()
}

################   12.13.14 其他图####
setwd('12.模型评估/')
heatmap405_train <- read.table("../11.lasso/risk.train.txt",header=T,row.names=1)
#heatmap405_train$riskScore <- scale(heatmap405_train$riskScore,center = F)
loc_train = match(c('Time','Status','riskScore'),colnames(heatmap405_train))
datalast_train = heatmap405_train[,loc_train]
datalast_train = data.frame(datalast_train)
#datalast_train$riskScore[which(datalast_train$riskScore >10)] <- 10
#图一
pdf(file="个体生存情况train.pdf",height = 12,width = 12)
par(mfrow=c(2,1))
par(mar=c(2,4.5,1,1))
col=c()
#0.978997006是risk的中位数
col[sort(datalast_train$riskScore) <= bestthreshold.train]="blue"
col[sort(datalast_train$riskScore) > bestthreshold.train]="red"
plot(sort(datalast_train$riskScore),axes=F,xlab = NA,ylab = "Risk Score",col=col,
     mgp=c(2.5,1,0),cex.lab=2)
box(lwd=2)
abline(v = length(which(sort(datalast_train$riskScore) <=  bestthreshold.train))+0.5,lty="dashed")
abline(h =  bestthreshold.train,lty="dashed")
text(50, bestthreshold.train+0.5,paste0('threshold = ',round( bestthreshold.train,4)),cex = 1.5,font=2)
axis(2,seq(0,100,5),cex.axis=2)
axis(1,seq(0,nrow(datalast_train),50),cex.axis=2)
legend("topleft", c("High risk","Low risk"), pch=16:16, col=c("red","blue"),cex = 1.5)
#图2
datalastSORT = datalast_train[order(datalast_train[,"riskScore"]),]
datalastSORT[,"Time"] = datalastSORT[,"Time"]/12/30
col=c()
col[datalastSORT[,"Status"]==1]= "red"
col[datalastSORT[,"Status"]==0]= "blue"
par(mar=c(2,4.5,1,1))
plot(datalastSORT[,"Time"],col=col,pch=16,axes=F,xlab = NA,mgp=c(2.5,1,0),cex.lab=1.8,
     ylab = "Following up (years)")
legend("topright", c("Alive","Dead"), pch=16:16, col=c("blue","red"),cex = 1.5)
box(lwd=2)
abline(v = length(which(sort(datalast_train$riskScore) <= bestthreshold.train))+0.5,lty="dashed")
# abline(h = 0.978997006,lty="dashed")
axis(2,seq(0,max(datalastSORT[,"Time"]),2),cex.axis=2)
axis(1,seq(0,nrow(datalast_train),10),cex.axis=2)
dev.off()


####  test
heatmap405_test <- read.table("../11.lasso/risk.test.txt",header=T,row.names=1)


loc_test = match(c('Time','Status','riskScore'),colnames(heatmap405_test))
datalast_test = heatmap405_test[,loc_test]
datalast_test = data.frame(datalast_test)

#图一
pdf(file="个体生存情况test.pdf",height = 12,width = 12)
par(mfrow=c(2,1))
par(mar=c(2,4.5,1,1))
col=c()
#0.978997006是risk的中位数
col[sort(datalast_test$riskScore) <=  bestthreshold.test]="blue"
col[sort(datalast_test$riskScore) >  bestthreshold.test]="red"
plot(sort(datalast_test$riskScore),axes=F,xlab = NA,ylab = "Risk Score",col=col,
     mgp=c(2.5,1,0),cex.lab=2)
box(lwd=2)
abline(v = length(which(sort(datalast_test$riskScore) <=  bestthreshold.test))+0.5,lty="dashed")
abline(h =  bestthreshold.test,lty="dashed")
text(50, bestthreshold.test+0.5,paste0('threshold = ',round( bestthreshold.test,4)),cex = 1.5,font=2)
axis(2,seq(0,20,5),cex.axis=2)
axis(1,seq(0,nrow(datalast_test),50),cex.axis=2)
legend("topleft", c("High risk","Low risk"), pch=16:16, col=c("red","blue"),cex = 1.5)
#图2
datalastSORT = datalast_test[order(datalast_test[,"riskScore"]),]
datalastSORT[,"Time"] = datalastSORT[,"Time"]/12/30
col=c()
col[datalastSORT[,"Status"]==1]= "red"
col[datalastSORT[,"Status"]==0]= "blue"
par(mar=c(2,4.5,1,1))
plot(datalastSORT[,"Time"],col=col,pch=16,axes=F,xlab = NA,mgp=c(2.5,1,0),cex.lab=1.8,
     ylab = "Following up (years)")
legend("topright", c("Alive","Dead"), pch=16:16, col=c("blue","red"),cex = 1.5)
box(lwd=2)
abline(v = length(which(sort(datalast_test$riskScore) <=  bestthreshold.test))+0.5,lty="dashed")
# abline(h = 0.978997006,lty="dashed")
axis(2,seq(0,max(datalastSORT[,"Time"]),2),cex.axis=2)
axis(1,seq(0,nrow(datalast_test),10),cex.axis=2)
dev.off()


####  valid1
heatmap405_valid <- read.table("../11.lasso/risk.valid宫颈癌.txt",header=T,row.names=1)

loc_valid = match(c('Time','Status','riskScore'),colnames(heatmap405_valid))
datalast_valid = heatmap405_valid[,loc_valid]
datalast_valid = data.frame(datalast_valid)

#图一
pdf(file="个体生存情况valid宫颈癌.pdf",height = 12,width = 12)
par(mfrow=c(2,1))
par(mar=c(2,4.5,1,1))
col=c()
#0.978997006是risk的中位数
col[sort(datalast_valid$riskScore) <=  bestthreshold.vaild1]="blue"
col[sort(datalast_valid$riskScore) >  bestthreshold.vaild1]="red"
plot(sort(datalast_valid$riskScore),axes=F,xlab = NA,ylab = "Risk Score",col=col,
     mgp=c(2.5,1,0),cex.lab=2)
box(lwd=2)
abline(v = length(which(sort(datalast_valid$riskScore) <=  bestthreshold.vaild1))+0.5,lty="dashed")
abline(h =  bestthreshold.vaild1,lty="dashed")
text(50, bestthreshold.vaild1+0.5,paste0('threshold = ',round( bestthreshold.vaild1,4)),cex = 1.5,font=2)
axis(2,seq(0,200,10),cex.axis=2)
axis(1,seq(0,nrow(datalast_valid),50),cex.axis=2)
legend("topleft", c("High risk","Low risk"), pch=16:16, col=c("red","blue"),cex = 1.5)
#图2
datalastSORT = datalast_valid[order(datalast_valid[,"riskScore"]),]
datalastSORT[,"Time"] = datalastSORT[,"Time"]/12/30
col=c()
col[datalastSORT[,"Status"]==1]= "red"
col[datalastSORT[,"Status"]==0]= "blue"
par(mar=c(2,4.5,1,1))
plot(datalastSORT[,"Time"],col=col,pch=16,axes=F,xlab = NA,mgp=c(2.5,1,0),cex.lab=1.8,
     ylab = "Following up (years)")
legend("topright", c("Alive","Dead"), pch=16:16, col=c("blue","red"),cex = 1.5)
box(lwd=2)
abline(v = length(which(sort(datalast_valid$riskScore) <=  bestthreshold.vaild1))+0.5,lty="dashed")
# abline(h = 0.978997006,lty="dashed")
axis(2,seq(0,max(datalastSORT[,"Time"]),2),cex.axis=2)
axis(1,seq(0,nrow(datalast_valid),10),cex.axis=2)
dev.off()



####  valid2
heatmap405_valid <- read.table("../11.lasso/risk.valid头颈鳞癌.txt",header=T,row.names=1)

loc_valid = match(c('Time','Status','riskScore'),colnames(heatmap405_valid))
datalast_valid = heatmap405_valid[,loc_valid]
datalast_valid = data.frame(datalast_valid)

#图一
pdf(file="个体生存情况valid头颈鳞癌.pdf",height = 12,width = 12)
par(mfrow=c(2,1))
par(mar=c(2,4.5,1,1))
col=c()
#0.978997006是risk的中位数
col[sort(datalast_valid$riskScore) <=  bestthreshold.vaild2]="blue"
col[sort(datalast_valid$riskScore) >  bestthreshold.vaild2]="red"
plot(sort(datalast_valid$riskScore),axes=F,xlab = NA,ylab = "Risk Score",col=col,
     mgp=c(2.5,1,0),cex.lab=2)
box(lwd=2)
abline(v = length(which(sort(datalast_valid$riskScore) <=  bestthreshold.vaild2))+0.5,lty="dashed")
abline(h =  bestthreshold.vaild2,lty="dashed")
text(50, bestthreshold.vaild2+0.5,paste0('threshold = ',round( bestthreshold.vaild2,4)),cex = 1.5,font=2)
axis(2,seq(0,200,10),cex.axis=2)
axis(1,seq(0,nrow(datalast_valid),50),cex.axis=2)
legend("topleft", c("High risk","Low risk"), pch=16:16, col=c("red","blue"),cex = 1.5)
#图2
datalastSORT = datalast_valid[order(datalast_valid[,"riskScore"]),]
datalastSORT[,"Time"] = datalastSORT[,"Time"]/12/30
col=c()
col[datalastSORT[,"Status"]==1]= "red"
col[datalastSORT[,"Status"]==0]= "blue"
par(mar=c(2,4.5,1,1))
plot(datalastSORT[,"Time"],col=col,pch=16,axes=F,xlab = NA,mgp=c(2.5,1,0),cex.lab=1.8,
     ylab = "Following up (years)")
legend("topright", c("Alive","Dead"), pch=16:16, col=c("blue","red"),cex = 1.5)
box(lwd=2)
abline(v = length(which(sort(datalast_valid$riskScore) <=  bestthreshold.vaild2))+0.5,lty="dashed")
# abline(h = 0.978997006,lty="dashed")
axis(2,seq(0,max(datalastSORT[,"Time"]),2),cex.axis=2)
axis(1,seq(0,nrow(datalast_valid),10),cex.axis=2)
dev.off()



################   热图三线表####
###    train####
clinicalall <- read.csv('F:/zy/！！！TCGA数据集整理/宫颈癌/clinical.cart.2022-10-11/TCGA.clinical.csv',
                        header = T,sep = ',')
colnames(clinicalall)[1] <- 'id.small'
heatmap.tmp <- read.csv('risk.train.group.txt',
                        header = T,sep = '\t')
sample_long <- row.names(heatmap.tmp)
heatmap.tmp$id.small <- substr(sample_long,1,12)
heatmap <- merge(heatmap.tmp,clinicalall,by="id.small")
write.table(heatmap,file="../16.独立预后/tmp.txt",sep = '\t',row.names = F,quote = F)

library(pheatmap)
heatmap = heatmap[order(heatmap[,"riskScore"]),]
table(heatmap$group)

df1 = heatmap[,c("PDE2A","ILK","MCM5","IGSF9","RNASEH2A")]
df1 = data.frame(t(df1))

annotation_col = data.frame(age = heatmap$age,
                            gender = heatmap$gender,
                            STAGE = heatmap$STAGE,
                            M = heatmap$M,
                            N = heatmap$N,
                            T = heatmap$T,
                            Risk_Level = factor(c(rep("Low",84),rep("High",115))))

table(annotation_col$STAGE)
annotation_col$STAGE <- factor(annotation_col$STAGE ,
                               labels=c("Stage I","Stage I","Stage I","Stage I","Stage I","Stage I",
                                        "Stage II",'Stage II','Stage II',"Stage II",'Stage II','Stage II',
                                        'Stage III','Stage III','Stage III','Stage III',
                                        'Stage IV','Stage IV','Stage IV'),
                               levels = c('Stage I','Stage IA','Stage IA2',"Stage IB","Stage IB1","Stage IB2",
                                          'Stage II','Stage IIA','Stage IIA1','Stage IIA2','Stage IIB','Stage IIC',
                                          'Stage III','Stage IIIA','Stage IIIB','Stage IIIC',
                                          'Stage IV','Stage IVA','Stage IVB'))


table(annotation_col$M)
annotation_col$M <- factor(annotation_col$M ,
                           labels=c("M0","M0","M1","M1",'M1','MX'),
                           levels = c("M0","cM0 (i+)","M1","M1a","M1b",'MX'))


table(annotation_col$N)
annotation_col$N <- factor(annotation_col$N ,
                           labels=c("N0","N0","N0","N1","N1",'N1','N1','N1',"N2",'N2','N2',
                                    "N3",'N3','N3','NX'),
                           levels = c("N0","N0 (i-)","N0 (i+)","N1",'N1mi',"N1a",'N1b','N1c',"N2","N2a",
                                      'N2b',"N3",'N3a','N3b','NX'))


table(annotation_col$T)
annotation_col$T <- factor(annotation_col$T ,
                           labels=c("T1","T1","T1","T1","T2","T2","T2","T3","T3",'T4',
                                    'T4','T4','Tis',"T1","T1","T2","T2","T3","TX"),
                           levels = c("T1","T1a","T1b","T1c","T2","T2a","T2b","T3","T3a",
                                      "T4","T4a","T4b",'Tis',"T1b1","T1b2","T2a1","T2a2","T3b","TX"))



#################    临床性状种类持续更新中


rownames(annotation_col) = colnames(df1)
df1 <- t(scale(t(df1)))
df1[which(df1 > quantile(df1,0.95))] <- quantile(df1,0.95)
df1[which(df1 < quantile(df1,0.05))] <- quantile(df1,0.05)

pdf(file = '与临床特征相关性train.更新.pdf',width = 10,height = 8)
pheatmap(df1,
         #color = greenred(75),
         #main = 'heatmap', # 图标题
         scale = 'row', #值集中的方向，“column”，“row” “none”
         annotation_col = annotation_col, #列注释
         #annotation_row = annotation_row, #行注释
         #legend_labels = NA,
         cluster_cols = F,          # 以列聚类
         #cluster_rows = FALSE,         # 以行聚类
         clustering_method = "complete", # 聚类方法 “complete” “average” “median”
         show_rownames = T, #不显示行名
         show_colnames = F, #不显示列名
         #gaps_row = 1169, # 分片
         fontsize = 10,
         angle_col=45)
dev.off()


library(table1)
annotation_col1 <- annotation_col
row.names(annotation_col1) <- heatmap$id.small
annotation_col1 <- na.omit(annotation_col1)
annotation_col1$Risk_Level <- factor(annotation_col1$Risk_Level,
                                     levels = c("High","Low","P-value"),
                                     labels = c("High","Low","P-value"))

## 分类变量标签调整
annotation_col1$age = ifelse(annotation_col1$age > 60,">60","<=60")
annotation_col1$gender <- factor(annotation_col1$gender)
annotation_col1$STAGE <- factor(annotation_col1$STAGE)
annotation_col1$M <- factor(annotation_col1$M)
annotation_col1$N <- factor(annotation_col1$N)
annotation_col1$T <- factor(annotation_col1$T)


## 左侧标签名调整
labels <- list(
  variables=list(age="age (years)",gender="gender",
                 STAGE="STAGE",M = "M",
                 N = "N",T = "T"),
  groups=list("", "Risk-Level",""))

## 设置短横线亚组
strata <- c(list(Total=annotation_col1), split(annotation_col1, annotation_col1$Risk_Level))

### 共有三处变量需要修改
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- annotation_col1[[name]]##修改annotation_col1
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.vaild(y ~ annotation_col1$Risk_Level)$p.value##修改annotation_col1
    } else {
      p <- chisq.test(table(y, droplevels(annotation_col1$Risk_Level)))$p.value###修改annotation_col1
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

## 绘制三线表
table1(strata, labels, groupspan=c(1,2,1),droplevels = F,
       render=rndr, render.strat=rndr.strat,
       topclass="Rtable1-zebra")



###     test####
setwd('13.模型测试/')
clinicalall <- read.csv('F:/zy/！！！TCGA数据集整理/宫颈癌/clinical.cart.2022-10-11/TCGA.clinical.csv',
                        header = T,sep = ',')
colnames(clinicalall)[1] <- 'id.small'
heatmap.tmp <- read.csv('risk.test.group.txt',
                        header = T,sep = '\t')
sample_long <- row.names(heatmap.tmp)
heatmap.tmp$id.small <- substr(sample_long,1,12)
heatmap <- merge(heatmap.tmp,clinicalall,by="id.small")

library(pheatmap)
heatmap = heatmap[order(heatmap[,"riskScore"]),]
table(heatmap$group)

df1 = heatmap[,c("PDE2A","ILK","MCM5","IGSF9","RNASEH2A")]
df1 = data.frame(t(df1))

annotation_col = data.frame(age = heatmap$age,
                            gender = heatmap$gender,
                            STAGE = heatmap$STAGE,
                            M = heatmap$M,
                            N = heatmap$N,
                            T = heatmap$T,
                            Risk_Level = factor(c(rep("Low",41),rep("High",43))))

table(annotation_col$STAGE)
annotation_col$STAGE <- factor(annotation_col$STAGE ,
                               labels=c("Stage I",'Stage I',"Stage I","Stage I","Stage I","Stage I","Stage I",
                                        "Stage II",'Stage II','Stage II',"Stage II",'Stage II','Stage II',
                                        'Stage III','Stage III','Stage III','Stage III',
                                        'Stage IV','Stage IV','Stage IV'),
                               levels = c('Stage I','Stage IA','Stage IA1','Stage IA2',"Stage IB","Stage IB1","Stage IB2",
                                          'Stage II','Stage IIA','Stage IIA1','Stage IIA2','Stage IIB','Stage IIC',
                                          'Stage III','Stage IIIA','Stage IIIB','Stage IIIC',
                                          'Stage IV','Stage IVA','Stage IVB'))


table(annotation_col$M)
annotation_col$M <- factor(annotation_col$M ,
                           labels=c("M0","M0","M1","M1",'M1','MX'),
                           levels = c("M0","cM0 (i+)","M1","M1a","M1b",'MX'))


table(annotation_col$N)
annotation_col$N <- factor(annotation_col$N ,
                           labels=c("N0","N0","N0","N1","N1",'N1','N1','N1',"N2",'N2','N2',
                                    "N3",'N3','N3','NX'),
                           levels = c("N0","N0 (i-)","N0 (i+)","N1",'N1mi',"N1a",'N1b','N1c',"N2","N2a",
                                      'N2b',"N3",'N3a','N3b','NX'))


table(annotation_col$T)
annotation_col$T <- factor(annotation_col$T ,
                           labels=c("T1","T1","T1","T1","T1","T2","T2","T2","T3","T3",'T4',
                                    'T4','T4','Tis',"T1","T1","T2","T2","T3","TX"),
                           levels = c("T1a1","T1","T1a","T1b","T1c","T2","T2a","T2b","T3","T3a",
                                      "T4","T4a","T4b",'Tis',"T1b1","T1b2","T2a1","T2a2","T3b","TX"))



rownames(annotation_col) = colnames(df1)
df1 <- t(scale(t(df1)))
df1[which(df1 > quantile(df1,0.95))] <- quantile(df1,0.95)
df1[which(df1 < quantile(df1,0.05))] <- quantile(df1,0.05)


pdf(file = '与临床特征相关性test.更新.pdf',width = 10,height = 8)
pheatmap(df1,
         #color = greenred(75),
         #main = 'heatmap', # 图标题
         scale = 'row', #值集中的方向，“column”，“row” “none”
         annotation_col = annotation_col, #列注释
         #annotation_row = annotation_row, #行注释
         #legend_labels = NA,
         cluster_cols = F,          # 以列聚类
         #cluster_rows = FALSE,         # 以行聚类
         clustering_method = "complete", # 聚类方法 “complete” “average” “median”
         show_rownames = T, #不显示行名
         show_colnames = F, #不显示列名
         #gaps_row = 1169, # 分片
         fontsize = 10,
         angle_col=45)
dev.off()


library(table1)
annotation_col1 <- annotation_col
row.names(annotation_col1) <- heatmap$id.small
annotation_col1 <- na.omit(annotation_col1)
annotation_col1$Risk_Level <- factor(annotation_col1$Risk_Level,
                                     levels = c("High","Low","P-value"),
                                     labels = c("High","Low","P-value"))

## 分类变量标签调整
annotation_col1$age = ifelse(annotation_col1$age > 60,">60","<=60")
annotation_col1$gender <- factor(annotation_col1$gender)
annotation_col1$STAGE <- factor(annotation_col1$STAGE)
annotation_col1$M <- factor(annotation_col1$M)
annotation_col1$N <- factor(annotation_col1$N)
annotation_col1$T <- factor(annotation_col1$T)


## 左侧标签名调整
labels <- list(
  variables=list(age="age (years)",gender="gender",
                 STAGE="STAGE",M = "M",
                 N = "N",T = "T"),
  groups=list("", "Risk-Level",""))

## 设置短横线亚组
strata <- c(list(Total=annotation_col1), split(annotation_col1, annotation_col1$Risk_Level))

### 共有三处变量需要修改
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- annotation_col1[[name]]##修改annotation_col1
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.vaild(y ~ annotation_col1$Risk_Level)$p.value##修改annotation_col1
    } else {
      p <- chisq.test(table(y, droplevels(annotation_col1$Risk_Level)))$p.value###修改annotation_col1
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

## 绘制三线表
table1(strata, labels, groupspan=c(1,2,1),droplevels = F,
       render=rndr, render.strat=rndr.strat,
       topclass="Rtable1-zebra")



###     vaild1####
clinicalall <- read.csv('anno.vaild宫颈癌.txt',header = T,sep = '\t')
heatmap.tmp <- read.csv('../11.lasso/risk.valid.group宫颈癌.txt',
                        header = T,sep = '\t')
heatmap.tmp$id <- row.names(heatmap.tmp)
heatmap <- merge(heatmap.tmp,clinicalall,by="id")

library(pheatmap)
heatmap = heatmap[order(heatmap[,"riskScore"]),]
table(heatmap$group)

df1 = heatmap[,c("PDE2A","ILK","MCM5","IGSF9","RNASEH2A")]
df1 = data.frame(t(df1))

annotation_col = data.frame(histology = heatmap$histology,
                            stage = heatmap$stage,
                            Risk_Level = factor(c(rep("Low",26),rep("High",29))))

table(annotation_col$stage)
annotation_col$stage <- factor(annotation_col$stage ,
                               labels=c("I","I","II","II",'III','IV','IV'),
                               levels = c("IB1","IB2","IIA","IIB",'IIIB','IVA','IVB'))


rownames(annotation_col) = colnames(df1)

pdf(file = '与临床特征相关性vaild宫颈癌.pdf',width = 10,height = 8)
pheatmap(df1,
         #color = greenred(75),
         #main = 'heatmap', # 图标题
         scale = 'row', #值集中的方向，“column”，“row” “none”
         annotation_col = annotation_col, #列注释
         #annotation_row = annotation_row, #行注释
         #legend_labels = NA,
         cluster_cols = F,          # 以列聚类
         #cluster_rows = FALSE,         # 以行聚类
         clustering_method = "complete", # 聚类方法 “complete” “average” “median”
         show_rownames = T, #不显示行名
         show_colnames = F, #不显示列名
         #gaps_row = 1169, # 分片
         fontsize = 10,
         angle_col=45)
dev.off()


library(table1)
annotation_col1 <- annotation_col
row.names(annotation_col1) <- heatmap$id
annotation_col1 <- na.omit(annotation_col1)
annotation_col1$Risk_Level <- factor(annotation_col1$Risk_Level,
                                     levels = c("High","Low","P-value"),
                                     labels = c("High","Low","P-value"))

## 分类变量标签调整
annotation_col1$histology <- factor(annotation_col1$histology)
annotation_col1$stage <- factor(annotation_col1$stage)

## 左侧标签名调整
labels <- list(
  variables=list(histology="histology",
                 stage="stage"),
  groups=list("", "Risk-Level",""))

## 设置短横线亚组
strata <- c(list(Total=annotation_col1), split(annotation_col1, annotation_col1$Risk_Level))

### 共有三处变量需要修改
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- annotation_col1[[name]]##修改annotation_col1
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.vaild(y ~ annotation_col1$Risk_Level)$p.value##修改annotation_col1
    } else {
      p <- chisq.test(table(y, droplevels(annotation_col1$Risk_Level)))$p.value###修改annotation_col1
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

## 绘制三线表
table1(strata, labels, groupspan=c(1,2,1),droplevels = F,
       render=rndr, render.strat=rndr.strat,
       topclass="Rtable1-zebra")



###     vaild2####
setwd('15.模型验证（头颈鳞癌）/')
clinicalall <- read.csv('anno.vaild头颈鳞癌.txt',header = T,sep = '\t')
heatmap.tmp <- read.csv('risk.valid.group头颈鳞癌.txt',
                        header = T,sep = '\t')
heatmap.tmp$id <- row.names(heatmap.tmp)
heatmap <- merge(heatmap.tmp,clinicalall,by="id")

library(pheatmap)
heatmap = heatmap[order(heatmap[,"riskScore"]),]
table(heatmap$group)

df1 = heatmap[,c("PDE2A","ILK","MCM5","IGSF9","RNASEH2A")]
df1 = data.frame(t(df1))

annotation_col = data.frame(age = heatmap$age,
                            gender = heatmap$gender,
                            smoking = heatmap$smoking,
                            T = heatmap$T,
                            N = heatmap$N,
                            HPV = heatmap$HPV,
                            stage = heatmap$stage,
                            Risk_Level = factor(c(rep("Low",153),rep("High",117))))

table(annotation_col$smoking)
table(annotation_col$T)
annotation_col$T <- factor(annotation_col$T ,
                               labels=c("1","2","3","4",'4'),
                               levels = c("1","2","3","4a",'4b'))

table(annotation_col$N)
annotation_col$N <- factor(annotation_col$N ,
                           labels=c("0","1","2","2","2",'3'),
                           levels = c("0","1","2a","2b","2c",'3'))

table(annotation_col$HPV)
table(annotation_col$stage)
annotation_col$stage <- factor(annotation_col$stage ,
                               labels=c("I","II","III",'IV','IV','IV'),
                               levels = c("I","II","III",'IVA','IVB','IVC'))



rownames(annotation_col) = colnames(df1)

pdf(file = '与临床特征相关性vaild头颈鳞癌.更新.pdf',width = 10,height = 8)
pheatmap(df1,
         #color = greenred(75),
         #main = 'heatmap', # 图标题
         scale = 'row', #值集中的方向，“column”，“row” “none”
         annotation_col = annotation_col, #列注释
         #annotation_row = annotation_row, #行注释
         #legend_labels = NA,
         cluster_cols = F,          # 以列聚类
         #cluster_rows = FALSE,         # 以行聚类
         clustering_method = "complete", # 聚类方法 “complete” “average” “median”
         show_rownames = T, #不显示行名
         show_colnames = F, #不显示列名
         #gaps_row = 1169, # 分片
         fontsize = 10,
         angle_col=45)
dev.off()


library(table1)
annotation_col1 <- annotation_col
row.names(annotation_col1) <- heatmap$id
annotation_col1 <- na.omit(annotation_col1)
annotation_col1$Risk_Level <- factor(annotation_col1$Risk_Level,
                                     levels = c("High","Low","P-value"),
                                     labels = c("High","Low","P-value"))

## 分类变量标签调整
annotation_col1$age = ifelse(annotation_col1$age > 60,">60","<=60")
annotation_col1$gender <- factor(annotation_col1$gender)
annotation_col1$smoking <- factor(annotation_col1$smoking)
annotation_col1$T <- factor(annotation_col1$T)
annotation_col1$N <- factor(annotation_col1$N)
annotation_col1$HPV <- factor(annotation_col1$HPV)
annotation_col1$stage <- factor(annotation_col1$stage)


## 左侧标签名调整
labels <- list(
  variables=list(age="age (years)",gender="gender",
                 smoking="smoking",T="T",
                 N="N",
                 HPV="HPV",
                 stage="stage"),
  groups=list("", "Risk-Level",""))

## 设置短横线亚组
strata <- c(list(Total=annotation_col1), split(annotation_col1, annotation_col1$Risk_Level))

### 共有三处变量需要修改
rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- annotation_col1[[name]]##修改annotation_col1
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.vaild(y ~ annotation_col1$Risk_Level)$p.value##修改annotation_col1
    } else {
      p <- chisq.test(table(y, droplevels(annotation_col1$Risk_Level)))$p.value###修改annotation_col1
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

## 绘制三线表
table1(strata, labels, groupspan=c(1,2,1),droplevels = F,
       render=rndr, render.strat=rndr.strat,
       topclass="Rtable1-zebra")




############################    15.独立预后分析####
setwd('16.独立预后/')
library(survival)
outTab=data.frame()
##准备数据：clinical.txt，最后只剩下377个没有UNKNOWN信息的个体
univariate_rt=read.table("1.input.2.txt",header=T,sep="\t",row.names=1,check.names=F)
for(i in colnames(univariate_rt[,3:ncol(univariate_rt)])){
  univar_cox_duli_ <- coxph(Surv(futime, fustat) ~ univariate_rt[,i], data = univariate_rt)
  #univar_cox_duli_=step(univar_cox_duli_,direction = "both")
  univar_cox_duli_ = summary(univar_cox_duli_)
  outTab=rbind(outTab,cbind(variable=i,
                            coef=univar_cox_duli_$coefficients[,"coef"],
                            HR=univar_cox_duli_$conf.int[,"exp(coef)"],
                            HR.95L=univar_cox_duli_$conf.int[,"lower .95"],
                            HR.95H=univar_cox_duli_$conf.int[,"upper .95"],
                            pvalue=univar_cox_duli_$coefficients[,"Pr(>|z|)"]))
}
write.table(outTab,file="univariateCox.Result2.xls",sep="\t",row.names=F,quote=F)

outTab=data.frame()
univariate_rt=read.table("1.input.5.txt",header=T,sep="\t",row.names=1,check.names=F)
univar_cox_duli_ <- coxph(Surv(futime, fustat) ~ univariate_rt[,3], data = univariate_rt)
#univar_cox_duli_=step(univar_cox_duli_,direction = "both")
univar_cox_duli_ = summary(univar_cox_duli_)
outTab=rbind(outTab,cbind(variable=3,
                          coef=univar_cox_duli_$coefficients[,"coef"],
                          HR=univar_cox_duli_$conf.int[,"exp(coef)"],
                          HR.95L=univar_cox_duli_$conf.int[,"lower .95"],
                          HR.95H=univar_cox_duli_$conf.int[,"upper .95"],
                          pvalue=univar_cox_duli_$coefficients[,"Pr(>|z|)"]))
write.table(outTab,file="univariateCox.Result5.xls",sep="\t",row.names=F,quote=F)


###    单因素森林图
##  去掉age
univariateCox_Result=read.table("univariateCox.Result.all.txt",header=T,sep="\t",row.names=1,check.names=F)
p.value<-signif(univariateCox_Result[,5], digits=2)
#wald.test<-signif(x$wald["test"], digits=2)
Coefficient<-signif(univariateCox_Result[,1], digits=2);#coeficient beta
HR <-signif(univariateCox_Result[,2], digits=2);#exp(beta)
HR.confint.lower <- signif(univariateCox_Result[,3], digits=2)
HR.confint.upper <- signif(univariateCox_Result[,4], digits=2)
HR.combine <- paste0(HR, " (",
                     HR.confint.lower, "-", HR.confint.upper, ")")
rescox.temp.1<-cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
names(rescox.temp.1)<-c("HR", "HR.confint.lower", "HR.confint.upper",'HR.combine',
                        "p.value")
rownames(rescox.temp.1) <- rownames(univariateCox_Result)

univOut_sig.plot.univar <- rescox.temp.1
gene <- rownames(univOut_sig.plot.univar)
hr <- univOut_sig.plot.univar[,1]
hrLow <- univOut_sig.plot.univar[,2]
hrHigh <- univOut_sig.plot.univar[,3]
Hazard.ratio <- univOut_sig.plot.univar[,4]
pVal <- univOut_sig.plot.univar[,5]

pdf(file='单因素cox森林图.pdf', width = 8,height = 4)
n <- nrow(univOut_sig.plot.univar)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=1.3
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5+0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
highlim <- max(as.numeric(hrHigh))+0.1
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(max(c(0,min(as.numeric(hrLow))))-0.1,highlim)
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",cex.lab=2)
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)
dev.off()

###################       多因素cox分析####
###准备好文件:将上一步显著的因素提出来，做多因素
multi_rt=read.table("2.input.txt",header=T,sep="\t",row.names=1,check.names=F)
cox <- coxph(Surv(futime, fustat) ~ ., data = multi_rt)
cox=step(cox,direction = "both")
coxSum_multi=summary(cox)
outTab=data.frame()
outTab=cbind(
  coef=coxSum_multi$coefficients[,"coef"],
  HR=coxSum_multi$conf.int[,"exp(coef)"],
  HR.95L=coxSum_multi$conf.int[,"lower .95"],
  HR.95H=coxSum_multi$conf.int[,"upper .95"],
  pvalue=coxSum_multi$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.Result.xls",sep="\t",row.names=F,quote=F)

###########    森林图
p.value<-signif(coxSum_multi$coefficients[,5], digits=2)
#wald.test<-signif(x$wald["test"], digits=2)
Coefficient<-signif(coxSum_multi$coefficients[,1], digits=2);#coeficient beta
HR <-signif(coxSum_multi$coefficients[,2], digits=2);#exp(beta)
HR.confint.lower <- signif(coxSum_multi$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(coxSum_multi$conf.int[,"upper .95"],2)
z<-signif(coxSum_multi$coefficients[,4],2)
HR.combine <- paste0(HR, " (",
                     HR.confint.lower, "-", HR.confint.upper, ")")
rescox.temp<-cbind(HR, HR.confint.lower,HR.confint.upper,HR.combine,p.value)
names(rescox.temp)<-c("HR", "HR.confint.lower", "HR.confint.upper",'HR.combine',
                      "p.value")
rownames(rescox.temp) <- rownames(coxSum_multi$coefficients)

univOut_sig.plot <- rescox.temp
gene <- rownames(univOut_sig.plot)
hr <- univOut_sig.plot[,1]
hrLow <- univOut_sig.plot[,2]
hrHigh <- univOut_sig.plot[,3]
Hazard.ratio <- univOut_sig.plot[,4]
pVal <- univOut_sig.plot[,5]
pdf(file='多因素cox森林图.pdf', width = 8,height = 3)
n <- nrow(univOut_sig.plot)
nRow <- n+1
ylim <- c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2.5))
xlim = c(0,3)
par(mar=c(4,2.5,2,1))
layout(matrix(c(1,1,1,2,2,2), 1, 6, byrow = TRUE))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
text.cex=1.3
text(0,n:1,gene,adj=0,cex=text.cex)
text(1.5+0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5+0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
highlim <- max(as.numeric(hrHigh))+0.1
par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
xlim = c(max(c(0,min(as.numeric(hrLow))))-0.1,highlim)
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio",cex.lab=2)
arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
abline(v=1,col="black",lty=2,lwd=2)
boxcolor = ifelse(as.numeric(hr) > 1, 'red', 'green')
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
axis(1, lwd = 1.5,cex.axis=text.cex,cex.lab=1.5)
dev.off()


########################################    列线图####
##  用多因素显著的做
pbc<-read.table("2.input.txt",header=TRUE,row.names=1)
library(rms)
dd<-datadist(pbc)
# pbc$age <- pbc$age/365
#pbc$futime <- pbc$futime*365
options(datadist="dd")
options(na.action="na.delete")
summary(pbc$futime)
coxpbc<-cph(formula = Surv(futime,fustat) ~  N + riskScore ,
            data=pbc,x=T,y=T,surv = T,na.action=na.delete)
print(coxpbc)
surv<-Survival(coxpbc)
surv1<-function(x) surv(1*365,x)
surv3<-function(x) surv(3*365,x)
surv5<-function(x) surv(5*365,x)


x<-nomogram(coxpbc,fun = list(surv1,surv3,surv5),lp=T,
            funlabel = c('1-year survival Probability','3-year survival Probability','5-year survival Probability'),
            maxscale = 100,fun.at = c(0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))

pdf("135年列线图.pdf",width = 14,height = 10)
plot(x, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2,
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE,
     cap.labels=FALSE,cex.var = 1.6,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))
dev.off()

#########         校准曲线####
#############    五年的数据太少了，所以画123年的
set.seed(38)
f1<-cph(formula = Surv(futime,fustat) ~ N + riskScore ,data=pbc,x=T,y=T,surv = T,
        na.action=na.delete, time.inc = 365)
#参数m=30表示每组30个样本进行重复计算
cal1<-calibrate(f1, cmethod="KM",method="boot",u=365,m=40,B=1000)

f3<-cph(formula = Surv(futime,fustat) ~  N + riskScore ,data=pbc,x=T,y=T,surv = T,
        na.action=na.delete, time.inc = 1095)
cal3<-calibrate(f3, cmethod="KM",method="boot",u=1095,m=40,B=1000)

f5<-cph(formula = Surv(futime,fustat) ~   N + riskScore ,data=pbc,x=T,y=T,surv = T,
        na.action=na.delete,time.inc = 1825)
cal5<-calibrate(f5, cmethod="KM", method="boot",u=1825,m=40 ,B=1000)


pdf("1_3_5校正曲线.pdf",width = 8,height = 8)
#lty = 0  不加误差线
plot(cal1,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     bty = "l", #只画左边和下边框
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#2166AC"),
     cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#2166AC"), pch = 16)
mtext("")

plot(cal3,lwd = 2,lty = 1,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#B2182B"), pch = 16)
mtext("")

plot(cal5,lwd = 2,lty = 1,errbar.col = c("#00A000"),
     xlim = c(0,1),ylim= c(0,1),col = c("#00A000"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 1, col = c("#00A000"), pch = 16)
mtext("")


abline(0,1, lwd = 2, lty = 3, col = c("#224444"))

legend("topleft", #图例的位置
       legend = c("1-year","3-year","5-year"), #图例文字
       col =c("#2166AC","#B2182B","#00A000"), #图例线的颜色，与文字对应
       lwd = 2,#图例中线的粗细
       cex = 1.2,#图例字体大小
       bty = "n")#不显示图例边框
dev.off()

set.seed(1)
v <- validate(f1, dxy=TRUE, B=1000)
Dxy = v[rownames(v)=='Dxy', colnames(v)=='index.corrected']
orig_Dxy = v[rownames(v)=='Dxy', colnames(v)=='index.orig']
bias_corrected_c_index  <- abs(Dxy)/2+0.5  # 计算校正c-index
orig_c_index <- abs(orig_Dxy)/2+0.5  # 计算未校正c-index
bias_corrected_c_index
orig_c_index

########   计算斜率
library(stringr)
caldat <- data.frame(summary(cal1))
cal1rate <- lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -1) ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -1))[["coefficients"]][["(Intercept)"]]
# summary(lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -1) ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -1)))
caldat <- data.frame(summary(cal3))
cal3rate <- lm(str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -3)[1:6] ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -8, -3)[1:6])[["coefficients"]][["(Intercept)"]]
caldat <- data.frame(summary(cal5))
cal5rate <- lm( str_sub(caldat[caldat$Var2 == " KM.corrected","Freq"], -8, -3)[1:6] ~ str_sub(caldat[caldat$Var2 == "mean.predicted","Freq"], -10, -3)[1:6])[["coefficients"]][["(Intercept)"]]



########################################      列线图ROC
##  用多因素显著的做
pbc<-read.table("4.列线图input.txt",header=TRUE,row.names=1)
library(rms)
dd<-datadist(pbc)
nobs<-NROW(pbc)
cutoff1<-1
cutoff2<-3
cutoff3<-5
options(datadist="dd")
options(na.action="na.delete")
Srv=Surv(pbc$futime,pbc$fustat)
coxmod<-coxph(Srv ~ N + riskScore,data=pbc)
summary(coxmod)

Nn<- as.numeric(pbc$N)
pbc$Npoint<- ifelse(Nn== 0,0,0.90083)

riskScoren <-  as.numeric(pbc$riskScore)
pbc$riskScorepoint <- ifelse(riskScoren == 0,0,1.3429)

pbc$PI<-pbc$points<-rowSums(pbc[,c("Npoint","riskScorepoint")])
data<-pbc[which(pbc$fustat!="NA"),]
data$futime <- data$futime/365


valid <- c()
for(i in c(1)){
  SROC <- survivalROC(Stime = data$futime,
                      status =  data$fustat,
                      marker = data$PI,
                      predict.time = i,method = 'KM')
  valid <- c(valid,SROC$AUC)
  pdf(file="135.列线图ROC.pdf",height = 6,width = 6)
  plot(SROC$FP,SROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red',
       xlab="False positive rate", ylab="True positive rate",
       main="ROC curve",
       lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
  abline(0,1)
  aucText=c()
  rocCol <- c('#FA8072','#63B8FF','#FFC1C1','#ADFF2F','#FFFF00')
  aucText=c(aucText,paste0("1 years"," (AUC=",sprintf("%.3f",SROC$AUC),")"))
  j =0
  for (i in c(3,5)){
    SROC=survivalROC(Stime = data$futime,
                     status =  data$fustat,
                     marker = data$PI,
                     predict.time =i, method="KM")
    j=j+1
    aucText=c(aucText,paste0(i," years"," (AUC=",sprintf("%.3f",SROC$AUC),")"))
    lines(SROC$FP,SROC$TP, type="l", xlim=c(0,1), ylim=c(0,1),col=rocCol[j+1],lwd = 3)
  }
  legend("bottomright", aucText,lwd=2,bty="n",col=rocCol,cex = 1)
  abline(0,1)
  dev.off()
}



##############    决策曲线####
library(survival)
library(ggplotify)
library(magick)
library(prodlim)
library(cmprsk)
source("stdca.R")
pbc <- read.table("2.input.txt",header=TRUE,row.names=1)
#pbc <- pbc[complete.cases(pbc),] #删掉缺失数据
pbc$futime <- pbc$futime/365
head(pbc)

Srv = Surv(pbc$futime, pbc$fustat)
N = coxph(Srv ~ N, data=pbc)
pbc$N = c(1 - (summary(survfit(N,newdata=pbc), times=5)$surv))

riskScore = coxph(Srv ~ riskScore, data=pbc)
pbc$riskScore = c(1 - (summary(survfit(riskScore,newdata=pbc), times=5)$surv))

N_riskScore = coxph(Srv ~ N + riskScore, data=pbc)
pbc$N_riskScore = c(1 - (summary(survfit(N_riskScore,newdata=pbc), times=5)$surv))

mod1 <- stdca(data=pbc, outcome="fustat", ttoutcome="futime", timepoint=3,
              predictors="N", cmprsk=TRUE, smooth=TRUE, xstop=1,intervention="FALSE")
mod2 <- stdca(data=pbc,outcome="fustat", ttoutcome="futime", timepoint=3,
              predictors="riskScore", cmprsk=TRUE, smooth=TRUE, xstop=1,intervention="FALSE")
mod3 <- stdca(data=pbc,outcome="fustat", ttoutcome="futime", timepoint=3,
              predictors="N_riskScore", cmprsk=TRUE, smooth=TRUE, xstop=1,intervention="FALSE")
pdf("决策曲线.pdf",width = 6,height = 6)
stdca(data=pbc, outcome="fustat", ttoutcome="futime", timepoint=3,
      predictors=c("N","riskScore","N_riskScore"),
      cmprsk=TRUE, smooth=TRUE,
      xstop=1,intervention="FALSE")
dev.off()





#############################   16.高低分组GSEA####
library('gplots')
library('limma')
library(dplyr)
setwd('17.高低风险组富集/')
### 去样本
tmp1=read.table("mRNA.symbol.uniq.group.txt",
                header=TRUE,check.names = FALSE, row.names = 1)
tmp2 <- data.frame(t(tmp1),check.names = F)
tmp3 <- tmp2[order(tmp2$group,decreasing = T),]
tmp4 <- data.frame(t(tmp3),check.names = F)
write.table(tmp4,file="mRNA.symbol.uniq.group.tmp.txt",sep = '\t',
            row.names = T,quote = F)

TCGA=read.table("mRNA.symbol.uniq.group.sort.txt",
                header=TRUE,row.names=1,check.names = FALSE)
par(mfrow=c(1,2))
#boxplot(data.frame(TCGA),col="blue")    ####画箱式图，比较数据分布情况，数据分布好，则不用进行log2转换
TCGA <- log2(TCGA+1)
TCGA[1:5,1:5]
TCGA.group <- c(rep("low",84),rep('high',115)) %>% factor(.,levels = c("low","high"),ordered = F)
#group <- group[,1] #定义比较组，按照癌症和正常样品数目修改#
TCGA.group <- model.matrix(~factor(TCGA.group))#把group设置成一个model matrix#
TCGA.fit <- lmFit(TCGA,TCGA.group)
TCGA.fit <- eBayes(TCGA.fit)
tempOutput = topTable(TCGA.fit,coef=2,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
TCGA.diff <- nrDEG
write.csv(TCGA.diff, "TCGA.high.low.limmaOut.csv")

library(clusterProfiler)
gene_df <-  read.csv("TCGA.high.low.limmaOut.csv",header=T,sep=",")
gene <- gene_df$SYMBOL
##################   ENSEMBL格式：ENSG00000108551  后面没有.1 .12等等
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- merge(gene_df,gene,by="SYMBOL")
write.table(gene,file="gene对应id.txt",sep = '\t',row.names = F,quote = F)

geneList <- gene_df$logFC
names(geneList)=gene_df$ENTREZID
geneList=sort(geneList,decreasing = T)

GOgmt<-read.gmt("c5.all.v7.5.1.entrez.gmt")
GO <-GSEA(geneList,TERM2GENE = GOgmt)
write.table(GO,file="GO.GSEA.txt",sep = '\t',row.names = F,quote = F)
library(enrichplot)
library(RColorBrewer)
pdf('GO.top10.pdf',height=10,width=14)
gseaplot2(GO, geneSetID = 2:11,pvalue_table = FALSE,
          rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))
dev.off()

KEGGgmt<-read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt")
KEGG <-GSEA(geneList,TERM2GENE = KEGGgmt)
write.table(KEGG,file="KEGG.GSEA.txt",sep = '\t',row.names = F,quote = F)
pdf('KEGG.top10.pdf',height=10,width=14)
gseaplot2(KEGG, geneSetID = 1:10,pvalue_table = FALSE,
          rel_heights = c(1.5, 0.5, 0.5),color= brewer.pal(10,'Paired'))

dev.off()



############################     18.高低风险组ssgsea####
setwd("18.高低风险组ssgsea/")

library(RColorBrewer)
library(pheatmap)
gene.a <- read.table("heatmap.txt",header=T,sep="\t",check.names = FALSE,row.names = 1)
annotation_col=read.table("../17.高低风险组富集/anno.txt",header=TRUE,
                          row.names=1,check.names = FALSE,sep = '\t')

pheatmap(gene.a,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,clustering_method = "ward.D",
         annotation_col =annotation_col,
         show_rownames = T, #不显示行名
         show_colnames = F,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
         scale= "row",border_color = NA, cluster_cols = FALSE)


library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
C <- read.table('box.txt',sep = '\t',header = T,check.names = F)
#C[,2:ncol(C)] <- apply(C[,2:ncol(C)],2,function(x){log2(x+1)})
colnames(C)[1] <- 'Type'
C1 <- gather(C,gene,expr,2:ncol(C))
ggboxplot(C1, x= 'gene', y='expr',
          ylab = "proportion", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),label = 'p.signif',label.x = 1.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))####  横坐标倾斜45度





################    19.estimate####
setwd("19.高低风险组estimate/")

library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
Data <- read.table('StromalScore.txt',sep = '\t',header = T,check.names = F)
colnames(Data)[1] <- 'Type'
#boxdat2 <- gather(boxdat1,gene,expr,2:ncol(boxdat1))
ggboxplot(Data, x= 'Type', y='StromalScore',
          ylab = "StromalScore", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),
                     label = 'p.signif',label.x = 1.5)


Data <- read.table('ImmuneScore.txt',sep = '\t',header = T,check.names = F)
colnames(Data)[1] <- 'Type'
#boxdat2 <- gather(boxdat1,gene,expr,2:ncol(boxdat1))
ggboxplot(Data, x= 'Type', y='ImmuneScore',
          ylab = "ImmuneScore", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),
                     label = 'p.signif',label.x = 1.5)


Data <- read.table('ESTIMATEScore.txt',sep = '\t',header = T,check.names = F)
colnames(Data)[1] <- 'Type'
#boxdat2 <- gather(boxdat1,gene,expr,2:ncol(boxdat1))
ggboxplot(Data, x= 'Type', y='ESTIMATEScore',
          ylab = "ESTIMATEScore", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),
                     label = 'p.signif',label.x = 1.5)


Data <- read.table('TumorPurity.txt',sep = '\t',header = T,check.names = F)
colnames(Data)[1] <- 'Type'
#boxdat2 <- gather(boxdat1,gene,expr,2:ncol(boxdat1))
ggboxplot(Data, x= 'Type', y='TumorPurity',
          ylab = "TumorPurity", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),
                     label = 'p.signif',label.x = 1.5)


