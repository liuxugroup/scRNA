#如何挑选目的基因集的基因
getwd()
setwd("/Users/louliming/Desktop/原发复发数据")

#输入矩阵----
setwd("~/Desktop/原发复发数据/03评价免疫微环境/免疫浸润0815")

##输入表型信息，分组
clinical_data=read.table('/Users/louliming/Desktop/原发复发数据/00工作文件夹/表型.txt',header = T,sep = '\t')
## 导入基因ID转化信息(仅保留编码基因)
ft <- read.csv('/Users/louliming/Desktop/生信/基因集/feature.csv', header = T)
ft <- filter(ft, ft$gene_biotype=='protein_coding')
##输入count
rawcount =read.table('/Users/louliming/Desktop/原发复发数据/02矩阵合集/公司原始ensg\ count.txt',header = T,sep = '\t')  #'data.txt'是文件名
rawcount$id <- ft$Symbol[match(rawcount$id, ft$ensembl_gene_id)]#注释
rawcount=rawcount[!duplicated(rawcount$id),] %>% na.omit()#count去重去na
rownames(rawcount)=rawcount[,1]
rawcount=rawcount[,-1]
names(rawcount) <-clinical_data$sample_name#根据某一列修改列名
range(rawcount)
#rawcount <- rawcount[apply(rawcount,1,sd)>0.5,] #质量控制！

##输入tpm
rawtpm=read.table('/Users/louliming/Desktop/原发复发数据/02矩阵合集/公司原始ensg\ tpm.txt',header = T,sep = '\t')
rawtpm$id <- ft$Symbol[match(rawtpm$id, ft$ensembl_gene_id)]#注释
rawtpm=rawtpm[!duplicated(rawtpm$id),] %>% na.omit()#tpm去重去na
rownames(rawtpm)=rawtpm[,1]
rawtpm=rawtpm[,-1]
names(rawtpm) <-clinical_data$sample_name
range(rawtpm)
#rawtpm <- rawtpm[apply(rawtpm,1,sd)>0.5,] #质量控制！总表达大于50%

tpm = rawtpm
# 将行名作为第一列
tpm1 <- tpm %>% rownames_to_column(var = "Symbol")
# 将列名作为第一行
tpm1 <- rbind(colnames(tpm1), tpm1)
write.table(tpm1,file='规范tpm.txt',quote = F,sep = '\t',row.names = F)

#筛选样本
#1,筛选出配对
clinical_data=clinical_data[order(clinical_data$sample),] 
group_list=clinical_data[which(clinical_data$event=='relapse'),]
count=rawcount[,group_list$sample_name]#提取筛选后的的count矩阵
tpm=rawtpm[,group_list$sample_name]#提取筛选后的的tpm矩阵

#2, 分别筛选出原发组和复发组
group_list=c(clinical_data$sample_name[which(clinical_data$event=='relapse'&clinical_data$paired=='pri')],#筛选过程：根据某一列的数值提取另一列
             clinical_data$sample_name[which(clinical_data$event=='relapse'&clinical_data$paired=='case')])
group_list
count=rawcount[,group_list]#提取筛选后的的count矩阵
tpm=rawtpm[,group_list]#提取筛选后的的tpm矩阵



#筛选基因
#某几个
gene_names <- c("EGFL6","ELFN2","ITGA11","ADGRB1","MFAP2","ACTA1","BPIFB3","CKM","MYL2","MYH7","XIRP2","MB","HNRNPA0","RNF138")
genes_tpm <- as.data.frame(tpm [rownames(tpm ) %in% gene_names ,])
#某个数据框      
gene_names <-read.table('/Users/louliming/Desktop/原发复发数据/00工作文件夹/想对比的基因.txt',header = T,sep = '\t')  #'data.txt'是文件名
genes_tpm <- as.data.frame(tpm [rownames(tpm ) %in% gene_names$gene ,])
#输出
genes_tpm=as.data.frame(t(genes_tpm))
genes_tpm <- rownames_to_column(genes_tpm, var = "gene")
write.table(genes_tpm,file='/Users/louliming/Desktop/原发复发数据/00工作文件夹/想仔细看的基因.txt',quote = F,sep = '\t',row.names = F)

#普通差异分析



#配对差异基因 分析----
# DEseq2的配对差异分析跑的有点慢，43对我的电脑需要约10分钟
# 所以，这里我为了加快运行速度，选了5对5.

suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
library(ggplot2)
library(ggpubr)
data <- as.data.frame(lapply(count, as.integer))
rownames(data) =rownames (count)
data[1:3,1:4]

##       TCGA.CV.6933.01 TCGA.CV.6933.11 TCGA.CV.6934.01 TCGA.CV.6934.11
## AAAS             2993            1328            2816            1977
## AACS             5124            4997            3522            4108
## AADAC              98               1              16              63

#配对分组
coldata_paired <- group_list[,c('sample','paired')]
rownames(coldata_paired)=group_list[,c('sample_name')]
coldata_paired$sample <- as.factor(coldata_paired$sample)
coldata_paired$paired <- as.factor(coldata_paired$paired)

# DEseq2的配对差异分析
logFC_cutoff = 0.584962501
#核心代码！！！
dds_paired <- DESeqDataSetFromMatrix(countData = data,
                                     colData = coldata_paired,
                                     design = ~ sample + paired)


dds_paired$paired <- relevel(dds_paired$paired, ref = "pri") # 指定哪一组作为对照

dds_paired <- DESeq(dds_paired)

resultsNames(dds_paired)
nrDEG_paired <- as.data.frame(results(dds_paired, tidy = TRUE))

#输入change
nrDEG_paired = nrDEG_paired[order(nrDEG_paired$log2FoldChange),] # 按照logFC排序
nrDEG_paired[1:3,1:6]
nrDEG_paired$change <-  as.factor(ifelse(nrDEG_paired$pvalue < 0.05 & abs(nrDEG_paired$log2FoldChange) > logFC_cutoff,
                                         ifelse(nrDEG_paired$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))
table(nrDEG_paired$change)

DEG <- nrDEG_paired
DEG <- nrDEG_paired %>% subset(pvalue < 0.05 & abs(log2FoldChange) >= log2(logFC_cutoff))
DEG.padj <- nrDEG_paired %>% subset(padj < 0.05 & abs(log2FoldChange) >= log2(logFC_cutoff))
DEG <- rownames_to_column(DEG, var = "ID")
write.table(DEG,file='鼻咽原发差异 deseq.txt',quote = F,sep = '\t',row.names = F)


#可视化----

##绘制火山图----

if (!require("theme_base", quietly = TRUE))
  install.packages("theme_base")
library(theme_base)

Volcano_data <- na.omit(nrDEG_paired)
Volcano_data$logP <- -log10(Volcano_data$pvalue)

Volcano_paired <- ggscatter(Volcano_data, x = "log2FoldChange", y = "logP",
                            color = "change",  
                            palette = c("#2f5688", "#BBBBBB", "#CC0000"),
                            size = 1,
                            font.label = 10, 
                            repel = T, 
                            xlab = "log2 FoldChange", 
                            ylab = "-log10 (pvalue)",   
                            title="Volcano Plot") + 
  geom_hline(yintercept = 1.30, linetype="dashed") +
  geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff), linetype="dashed")+ 
  theme(legend.position = "right",plot.title = element_text(size = 14,color="black",hjust = 0.5))
Volcano_paired

#热图----
table = tpm
table=table[!duplicated(rownames(table)),] %>% na.omit()#tpm去重去na
#table2 <- rownames_to_column(table2, var = "gene")
#write.table(table2,file='table2.txt',quote = F,sep = '\t',row.names = F)
table=as.data.frame(table)

norm_lim <- function(x, th){#归一化函数
  if(sum(!is.na(x)) == 0)
    return(x)
  if(sum(!is.na(x)) < 2)
    return(x-mean(x, na.rm = TRUE))
  sdx <- sd(x, na.rm = TRUE)
  if(sdx == 0)
    x <- (x-mean(x, na.rm = TRUE))
  else
    x <- (x-mean(x, na.rm = TRUE))/sdx
  if(missing(th))
    return(x)
  else{
    x[x > th] <- th
    x[x < -th] <- -th
  }
  return(x)
}

heatmap=table
dim(heatmap)
hmset_scaled <- apply(heatmap, 1, norm_lim, 1.5) %>% t() %>% as.data.frame()
range(hmset_scaled)


# 列注释
group_list1 =ifelse(str_sub(colnames(hmset_scaled),0,1) =="P","pri","case")
group_list1=factor(ifelse(str_sub(colnames(hmset_scaled),0,1) =="P","pri","case"),levels = c('pri','case'))#将样品带上属性，分为两个level，前面是对照组
str(group_list1)
annotation_col <- data.frame(Group = group_list1)
rownames(annotation_col) <- colnames(hmset_scaled)

heatmap_PriVsCase <-pheatmap(hmset_scaled, 
                           scale = 'none',#归一化已经在上面完成
                           # scale = 'row',
                           annotation_col = annotation_col,
                           cluster_cols = T,  show_colnames = T,
                           cluster_rows = T, show_rownames = T,
                           treeheight_row = 10,
                           border_color = NA)
print(heatmap_PriVsCase)

dev.off()
ann_colors = list(
  Groups = c(1="#172A88", 2="#F4A016",3="E60012"))


#富集分析
#Gsea
##富集分析
```{r 富集分析}
#数据整理
#手动转一下格式
library(org.Hs.eg.db)
DEG=read.table('~/Desktop/deg.txt',header = T,sep = '\t')
genelist_input <- DEG
colnames(genelist_input)[1]='Gene'
genename <- as.character(genelist_input[,1]) #提取第一列基因名
gene_map <- select(org.Hs.eg.db, keys=genename, 
                   keytype="SYMBOL", columns=c("ENTREZID"))
colnames(gene_map)[1]<-"Gene"
aaa<-inner_join(gene_map,genelist_input,by = "Gene")
aaa<-aaa[,-1]
aaa<-na.omit(aaa)
aaa<-aaa[order(aaa$log2FoldChange,decreasing = TRUE),]

bbb <- aaa[aaa$change == "UP", ]
geneList = aaa$log2FoldChange
names(geneList) = as.character(aaa[,1])

#GO
GO_gseresult <- gseGO(geneList, 'org.Hs.eg.db', 
                      keyType = "ENTREZID", ont="all",
                      nPerm = 1000, minGSSize = 10, maxGSSize = 1000, 
                      pvalueCutoff=1)
GO=summary(GO_gseresult)
write.csv(GO,'GO.csv')

#KEGG
KEGG_gseresult <- gseKEGG(geneList, 
                          nPerm = 1000, minGSSize = 10, 
                          maxGSSize = 1000, pvalueCutoff=1)
KEGG=summary(KEGG_gseresult)
write.csv(KEGG,'KEGG.csv')

#gsea图
gseaplot2(GO_gseresult, 
          "GO:0005200", 
          color = "firebrick",
          rel_heights=c(1, .2, .6))



#Hallmark
genelist_input <- DEG
colnames(genelist_input)[1]='Gene'
genelist_input<-genelist_input[order(genelist_input$logFC,
                                     decreasing = TRUE),]
hallmark <- read.gmt("h.all.v2022.1.Hs.symbols.gmt")
genelist <- genelist_input$logFC
names(genelist) <- genelist_input$Gene
gsea.re1<- GSEA(genelist,  #待富集的基因列表
                TERM2GENE = hallmark,  #基因集
                pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                pAdjustMethod = 'BH') 
Hallmark=summary(gsea.re1)
write.csv(Hallmark,'Hallmark.csv')

#自己基因集
genelist_input <- DEG
colnames(genelist_input)[1]='Gene'
genelist_input<-genelist_input[order(genelist_input$logFC,
                                     decreasing = TRUE),]
metaaa <- read.gmt("NC代谢数据.gmt")
genelist <- genelist_input$logFC
names(genelist) <- genelist_input$Gene
gsea.re1<- GSEA(genelist,  #待富集的基因列表
                TERM2GENE = metaaa,  #基因集
                pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                pAdjustMethod = 'BH') 
Hallmark=summary(gsea.re1)
write.csv(Hallmark,'自己基因.csv')
```



