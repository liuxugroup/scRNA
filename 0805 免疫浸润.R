#https://blog.csdn.net/Ayue0616/article/details/131009335
#https://blog.csdn.net/Ayue0616/article/details/131009363

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
rawcount <- rawcount[apply(rawcount,1,sd)>0.5,] #质量控制！

##输入tpm
rawtpm=read.table('/Users/louliming/Desktop/原发复发数据/02矩阵合集/公司原始ensg\ tpm.txt',header = T,sep = '\t')
rawtpm$id <- ft$Symbol[match(rawtpm$id, ft$ensembl_gene_id)]#注释
rawtpm=rawtpm[!duplicated(rawtpm$id),] %>% na.omit()#tpm去重去na
rownames(rawtpm)=rawtpm[,1]
rawtpm=rawtpm[,-1]
names(rawtpm) <-clinical_data$sample_name
range(rawtpm)
rawtpm <- rawtpm[apply(rawtpm,1,sd)>0.5,] #质量控制！总表达大于50%

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

#加载R包----
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}

if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")

library(IOBR)



#代码免疫浸润----
# 数据处理
logtpm <- log2(tpm+0.1)
#logtpm <- logtpm[apply(logtpm,1,sd)>0.5,]
logtpm[1:4,1:4]
dim(logtpm)

# MCPcounter
im_mcpcounter <- deconvo_tme(eset = logtpm,
                             method = "mcpcounter")
# EPIC
im_epic <- deconvo_tme(eset = logtpm,
                       method = "epic",
                       arrays = F)

# xCell
im_xcell <- deconvo_tme(eset = logtpm,
                        method = "xcell",
                        arrays = F)

# CIBERSORT
im_cibersort <- deconvo_tme(eset = logtpm,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000)

# IPS
im_ips <- deconvo_tme(eset = logtpm,
                      method = "ips",
                      plot = F)

# quanTIseq
im_quantiseq <- deconvo_tme(eset = logtpm,
                            method = "quantiseq",
                            scale_mrna = T)

# ESTIMATE
im_estimate <- deconvo_tme(eset = logtpm,
                           method = "estimate")


# TIMER
im_timer <- deconvo_tme(eset = logtpm
                        ,method = "timer"
                        ,group_list = rep("coad",dim(logtpm)[2]))

##ssGSEA计算免疫浸润分数----
load("/Users/louliming/Desktop/原发复发数据/03评价免疫微环境/ssGSEA28.Rdata")
im_ssgsea <- calculate_sig_score(eset = logtpm
                                 , signature = cellMarker # 这个28种细胞的文件需要自己准备
                                 , method = "ssgsea" # 选这个就好了
)
im_ssgsea[1:4,1:4]
##immunoai----

immunoai=read.table('/Users/louliming/Desktop/原发复发数据/03评价免疫微环境/免疫浸润0815/immunoai/原发复发\ ImmuCellAI_abundance_result.txt',header = T,sep = '\t')  #'data.txt'是文件名
immunoai <- as.data.frame(immunoai)
rownames(immunoai)=immunoai[,1]
immunoai=immunoai[,-1]
immunoai<- immunoai %>% rownames_to_column(var = "ID")
immunoai<- rbind(colnames(immunoai), immunoai)

##自己的ssgsea----
#https://mp.weixin.qq.com/s/6Y44sz_faS9dsCznJX6SzA
library(GSVA)
library(limma)
library(GSEABase)


symbol=as.matrix(tpm)
exp_data=symbol
dimnames=list(rownames(exp_data),colnames(exp_data))
exp_data1=matrix(as.numeric(as.matrix(exp_data)),nrow=nrow(exp_data),dimnames=dimnames)
exp_data1=avereps(exp_data1)
exp_data1=exp_data1[rowMeans(exp_data1)>0,]
geneSet=getGmt("/Users/louliming/Desktop/原发复发数据/03评价免疫微环境/佳蔚\ ssgsea\ 评分/nm\ 免疫浸润基因集.gmt", geneIdType=SymbolIdentifier())#读取我们所提供的免疫细胞基因文件
ssGSEA_Score=gsva(exp_data1, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)#ssGSEA计算

normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}#定义ssGSEA_Score矫正函数
norm_ssGSEA_Score=normalize(ssGSEA_Score)#对ssGSEA_Score进行矫正
norm_ssGSEA_Score=as.data.frame(norm_ssGSEA_Score)

norm_ssGSEA_Score2 = as.data.frame(t(norm_ssGSEA_Score))
norm_ssGSEA_Score2 <- norm_ssGSEA_Score2 %>% rownames_to_column(var = "ID")
norm_ssGSEA_Score2 <- rbind(colnames(norm_ssGSEA_Score2), norm_ssGSEA_Score2)


#write.table(norm_ssGSEA_Score2,file='佳蔚 免疫浸润.txt',quote = F,sep = '\t',row.names = F)



#数据合并----
rm(tme_combine) 
tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
  inner_join(im_ips, by= "ID") %>% 
  inner_join(im_quantiseq, by="ID") %>% 
  inner_join(im_estimate, by= "ID") %>% 
  inner_join(im_timer, by= "ID") %>% 
  inner_join(im_ssgsea, by= "ID") %>% 
  inner_join(immunoai, by= "ID") %>% 
  inner_join(norm_ssGSEA_Score2, by = "ID")

tme_combine[1:4,1:4]
saveRDS(tme_combine, paste0("原发复发免疫浸润.rds"))
dim(tme_combine)


#算比例----
# 创建一个函数来计算细胞比例
calculate_cell_ratio <- function(matrix, cell1, cell2) {
  cell1_col <- match(cell1, colnames(matrix))
  cell2_col <- match(cell2, colnames(matrix))
  cell_ratio <- matrix[, cell1_col] / matrix[, cell2_col]
  matrix[, new_col_name] <- cell_ratio
  return(matrix)}


#给 ratio 赋值#输入想要计算的评分
ratio<- as.data.frame(norm_ssGSEA_Score2)#输入想要计算的评分
rownames(ratio)=ratio[,1]
ratio=ratio[,-1]
colnames(ratio)=ratio[1,]
ratio=ratio[-1,]
#检查每一列的数据类型
for (col in colnames(ratio)) {
  col_class <- class(ratio[, col])
  print(paste(col, ":", col_class))}
# 遍历每一列
for (col in colnames(ratio)) {
  if (!is.numeric(ratio[[col]])) {
    # 将非数值列转换为数值列
  ratio[[col]] <- as.numeric(as.character(ratio[[col]]))
  }}

result <- data.frame(matrix(0, nrow = nrow(ratio), ncol = 0))
#开始算
for (i in 1:(ncol(ratio)-1)) {
  for (j in 1:(ncol(ratio)-1)) {
    # 获取细胞1和细胞2的列名
    cell1 <- colnames(ratio)[i]
    cell2 <- colnames(ratio)[j]
    # 计算细胞比例
    cell_ratio <- ratio[, i] / ratio[, j]
    # 创建新的列名 
    new_col_name <- paste(cell1, cell2, sep = "XX")  
    # 将细胞比例添加为新的列到结果数据框
    result[, new_col_name] <- cell_ratio
  }
}

#输出比例result,result1
rownames(result) <- rownames(ratio)
result <- t(result)
result <- as.data.frame(result)
result1 <- result %>% rownames_to_column(var = "ID")
#write.table(result,file='原发复发细胞比例（倒置）.txt',quote = F,sep = '\t',row.names = F)
#saveRDS(result, paste0("原发复发免疫评分比例XX.rds"))

# 将结果数据框与原始矩阵合并
updated_matrix <- cbind(tme_combine, result)

#可视化----
library(tidyHeatmap)
library(tidyverse)
library(RColorBrewer)


#简便，看总体评分
res_tmp <- cell_bar_plot(im_ssgsea)

#先分组
tmp =im_ssgsea
tmp=as.data.frame(tmp)
rownames(tmp)=tmp[,1]


# 分组

group_list=c(clinical_data$tpm[which(clinical_data$event=='relapse'&clinical_data$paired=='pri')],#筛选过程：根据某一列的数值提取另一列
             clinical_data$tpm[which(clinical_data$event=='relapse'&clinical_data$paired=='case')])
group_list
selcted_tmp=tmp[group_list,]#提取筛选后的样本

#重命名
selcted_tmp$ID <- ifelse(seq_along(selcted_tmp$ID) <= 17,
                         paste("00_", selcted_tmp$ID),
                         paste("11_", selcted_tmp$ID))



# 首先变为长数据
cibersort_long <- selcted_tmp %>% 
  select(`P-value_CIBERSORT`,Correlation_CIBERSORT, RMSE_CIBERSORT,ID,everything()) %>% 
  pivot_longer(- c(1:4),names_to = "cell_type",values_to = "fraction") %>% 
  dplyr::mutate(cell_type = gsub("_CIBERSORT","",cell_type),
                cell_type = gsub("_"," ",cell_type))

head(cibersort_long[,4:6])
## # A tibble: 6 × 3
##   ID                           cell_type                  fraction
##   <chr>                        <chr>                         <dbl>
## 1 TCGA-D5-6540-01A-11R-1723-07 B cells naive                0.0504
## 2 TCGA-D5-6540-01A-11R-1723-07 B cells memory               0     
## 3 TCGA-D5-6540-01A-11R-1723-07 Plasma cells                 0     
## 4 TCGA-D5-6540-01A-11R-1723-07 T cells CD8                  0.119 
## 5 TCGA-D5-6540-01A-11R-1723-07 T cells CD4 naive            0     
## 6 TCGA-D5-6540-01A-11R-1723-07 T cells CD4 memory resting   0.0951
p1 <- cibersort_long %>% 
  ggplot(aes(ID,fraction))+
  geom_bar(stat = "identity",position = "stack",aes(fill=cell_type))+
  labs(x=NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette4,name=NULL)+ # iobr还给大家准备了几个色盘，贴心！
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
  )
p1

# 有顺序的箱线图
library(forcats)

p2 <- ggplot(cibersort_long,aes(fct_reorder(cell_type, fraction),fraction,fill = cell_type)) + 
  geom_boxplot() + 
  #geom_jitter(width = 0.2,aes(color=cell_type))+
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_fill_manual(values = palette4)
p2
library(ggpubr)
library(stringr)


cibersort_long$Group = ifelse(as.numeric(str_sub(cibersort_long$ID,1,2))>0,"case","pri")

p3 <- ggplot(cibersort_long,aes(fct_reorder(cell_type,fraction),fraction,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  scale_fill_manual(values = palette1[c(2,4)])+ 
  theme_bw() + 
  labs(x = NULL, y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,hjust = 1),
        axis.text = element_text(color = "black",size = 12))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),
                     method = "kruskal.test",label.y = 0.4)
p3

library(tidyHeatmap)

p4 <- heatmap(.data = cibersort_long
              ,.row = cell_type
              ,.column = ID
              ,.value = fraction
              ,scale = "column"
              ,palette_value = circlize::colorRamp2(
                seq(-2, 2, length.out = 11), 
                RColorBrewer::brewer.pal(11, "RdBu")
              )
              ,show_column_names=F
              ,row_names_gp = gpar(fontsize = 10),
              column_names_gp = gpar(fontsize = 7),
              column_title_gp = gpar(fontsize = 7),
              row_title_gp = gpar(fontsize = 7)
) %>% 
  add_tile(Group) # 新版本已经改了，注意
p4


ssgsea_long <- im_ssgsea %>% 
  pivot_longer(- ID,names_to = "cell_type",values_to = "Score")
head(ssgsea_long)
## # A tibble: 6 × 3
##   ID                           cell_type                         Score
##   <chr>                        <chr>                             <dbl>
## 1 TCGA-3L-AA1B-01A-11R-A37K-07 Activated B cell               -0.192  
## 2 TCGA-3L-AA1B-01A-11R-A37K-07 Activated CD4 T cell            0.0120 
## 3 TCGA-3L-AA1B-01A-11R-A37K-07 Activated CD8 T cell            0.170  
## 4 TCGA-3L-AA1B-01A-11R-A37K-07 Activated dendritic cell        0.00545
## 5 TCGA-3L-AA1B-01A-11R-A37K-07 CD56bright natural killer cell  0.199  
## 6 TCGA-3L-AA1B-01A-11R-A37K-07 CD56dim natural killer cell     0.351
ggplot(ssgsea_long, aes(cell_type, Score))+
  geom_violin(width=2.0,aes(color=cell_type))+
  geom_boxplot(width=0.2,fill="black") + 
  theme_bw() + 
  labs(x = NULL) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  #scale_fill_manual(values = palette4)+
  scale_color_manual(values = palette4,name=NULL)


#批量配对t检验----
##https://blog.csdn.net/qazplm12_3/article/details/116141619

#整理数据
library(dplyr)
library(reshape2)
view(result)

paired_data <- as.data.frame(tme_combine)#列名是评分，行名是样本
rownames(paired_data)=paired_data[,1]
paired_data=paired_data[,-1]
paired_data$Group <- ifelse(substr(rownames(paired_data), 9, 9) == "P","pri","case")
paired_data$Paired <- substr(rownames(paired_data), 7, 8)
paired_data <- paired_data %>% 
  dplyr::select(Group, Paired, everything())
paired_data[1:4,1:4]
#paired_data1  <- rownames_to_column(paired_data, var = "ID")
#write.table(bb,file=' 代码 免疫浸润 转置.txt',quote = F,sep = '\t',row.names = T)

bb <- paired_data#对列的名字非常的敏感，不能有空格,-,+等奇怪符号
names(bb) <- gsub(" ", "_", names(bb))
names(bb) <- gsub("\\(", "_", names(bb))
names(bb) <- gsub("\\+", "_", names(bb))

bb1 <- melt(bb)
bb1 <- melt(bb, id.vars = c("Group", "Paired"))
bb$Group <- factor(bb$Group,levels = c("pri","case"))
bb1$Group <- factor(bb1$Group,levels = c("pri","case"))

bb.sample <- colnames(bb)[3:ncol(bb)]
test.b <- c()
x <- c("a","b")
y <- c("a","a")

###如果是评分比例，输出 test.b是否显著，test.p 显示p 值
p_values <- c()
test.b <- c()
test.P <- c()
col_names <- c()  # 创建一个空向量来存储列名
for (i in bb.sample) {
  tryCatch({
    fit1 <- t.test(as.formula(sprintf("%s ~ Group", i)),
                   data = bb, paired = TRUE)
    if (fit1$p.value < 0.05) {
      res1 <- x
    } else {
      res1 <- y
    }
    test.b <- cbind(test.b, res1)
    p_values <- c(p_values, fit1$p.value)
    col_names <- c(col_names, i)  # 将当前变量名添加到列名向量中
    colnames(test.b) <- col_names  # 设置列名
  }, error = function(e) {
    res1 <- NA
    test.b <- cbind(test.b, res1)
    col_names <- c(col_names, i)
    colnames(test.b) <- col_names
    
  })
}

for (i in bb.sample) {
  tryCatch({
    fit1 <- t.test(as.formula(sprintf("%s ~ Group", i)),
                   data = bb, paired = TRUE)
    p_values <-  fit1$p.value
    test.P <- cbind(test.P, p_values)
    col_names <- c(col_names, i)  # 将当前变量名添加到列名向量中
    colnames(test.b) <- col_names  # 设置列名
  }, error = function(e) {
    res1 <- NA
    test.P <- cbind(test.P, p_values)
    col_names <- c(col_names, i)
    colnames(test.b) <- col_names
  })
}

#saveRDS(test.b, paste0("test.b.rds"))
rownames(test.b) <- levels(bb$Group)
test.b <- melt(test.b)
colnames(test.b) <- c("Group","variable","value")

######如果是普通评分，输出 test.b是否显著，test.p 显示p 值
p_values <- c()
test.b <- c()
for (i in bb.sample) {
  fit1 <- t.test(as.formula(sprintf("%s ~ Group", i)),
                 data = bb, paired = TRUE)       
  if (fit1$p.value < 0.05) {
    res1 <- "x"
  } else {
    res1 <- "y"
  }
  test.b <- cbind(test.b, res1)
}

#输出有无意义，test.b
colnames(test.b) <- colnames(bb)[3:ncol(bb)]
rownames(test.b) <- levels(bb$Group)
test.b <- melt(test.b)
colnames(test.b) <- c("Group","variable","value")

#输出p值，test.P
test.P = as.data.frame(p_values)
rownames(test.P) <- colnames(bb)[3:ncol(bb)]
test.P <-test.P %>% rownames_to_column(var = "Cell")

#输出
saveRDS(test.b, paste0("原发复发佳蔚评分配对,ttest.rds"))
write.table(test.b,file='原发复发佳蔚评分配对,ttest.txt',quote = F,sep = '\t',row.names = F)



#作图数据处理----
library(ggplot2)
library(forcats)
library(dplyr)
library(tidyverse)

#数据处理
test.b1 <- bb[,c(1,c(1,3:ncol(bb)))] %>% gather(variable,value,-Group) %>% group_by(variable,Group) %>% 
  summarise(Max = max(value))
test.b11 <- dcast(test.b1,Group~variable) 
for (i in 2:ncol(test.b11)) {
  test.b11[,i] <- as.numeric(test.b11[,i]) + max(as.numeric(test.b11[,i])) * 0.015
  }
test.b11 <- melt(test.b11)
test.b1 <- merge(test.b1,test.b11,by = c("variable","Group"))
test.b2 <- merge(test.b,test.b1,by = c("variable","Group"))

bb1 <- melt(bb)
bb1$Group <- factor(bb1$Group,levels = c("pri","case"))

#输入感兴趣的评分/基因
ganxingqu=read.table('/Users/louliming/Desktop/原发复发数据/03评价免疫微环境/免疫浸润0815/想看的评分.txt',header = T,sep = '\t') 
ganxingqu1 <- unlist(ganxingqu$ID)
ganxingqu1

bb1 <- melt(bb)
bb1$Group <- factor(bb1$Group,levels = c("pri","case"))

bb1$variable <- factor(bb1$variable,levels = ganxingqu1)
test.b2 <- merge(test.b,test.b1,by = c("variable","Group"))
test.b2$variable <- factor(test.b2$variable,levels = ganxingqu1)



####显示P值----
test.p.selected<- data.frame(Cell = character(), P_value = numeric(), stringsAsFactors = FALSE)
for (cell in ganxingqu1) {
  # 检查cell是否在test.P的Cell列中
  if (cell %in% test.P$Cell) {
    # 提取对应细胞的P值并添加到新的数据框中
    p_value <- test.P[test.P$Cell == cell, "p_values"]
    test.p.selected <- rbind(test.p.selected, data.frame(Cell = cell, P_value = p_value))
  }
}
view(test.p.selected)

####作图----
library(forcats)
library(ggplot2)
cbbPalette <- c("#B2182B","#56B4E9")
p <- ggplot(bb1,aes(Group,value)) + 
  geom_boxplot(aes(color = Group),outlier.size = 0,size = 0.8) +
  geom_line(aes(group = Paired),color = "grey70") +
  geom_point(aes(color = Group),size = 2.5,alpha = 0.5) +
  scale_color_manual(values = cbbPalette) +
  ylab("T test") +#标题名字
  facet_wrap(.~variable,ncol = 8,scales = "free_y") +
  theme_bw()+
  theme(axis.ticks.length = unit(0.4,"lines"), 
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"), 
        axis.title.x=element_blank(),
        axis.title.y=element_text(colour='black', size=18,face = "bold",
                                  vjust = 1.5),
        axis.text.y=element_text(colour='black',size=13,face = "bold"),
        axis.text.x=element_text(colour = "black",size = 14,face = "bold"),
        strip.text = element_text(colour = "black",size = 8,face = "bold"),#小标题大小
        legend.position = "none")

print(p)
dev.off()

#记得每次做完图要重启bb1
bb1 <- melt(bb)
bb1$Group <- factor(bb1$Group,levels = c("pri","case"))


# geom_text(data = test.b2,aes(x = Group,y = value.y,label = value.x),#用a、b 显示显著性
#          size = 6,color = "black",fontface = "bold") +



