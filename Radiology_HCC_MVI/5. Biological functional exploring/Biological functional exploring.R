#### Identifying DEGs ####

library(DESeq2)
library(tidyverse)

#The data can download from TCGA dataset.
counts_01A <- read.table("LIHC_counts_mRNA_01A.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp <- read.table("LIHC_fpkm_mRNA_01A.txt", sep = "\t",row.names = 1,check.names = F,header = T)

#In this file, 35 cases of TCIA data (row-gene expression; column-samples). We added 6 imaging variables in the last 6 rows.
exp2 <- read.csv("TCIA_counts.csv",row.names = 1) 
exp <- exp[,colnames(exp2)]
counts_01A <- counts_01A[,colnames(exp2)]
identical(colnames(counts_01A),colnames(exp))
counts_01A <- rbind(counts_01A,exp2[c(60489:60494),])
exp <- rbind(exp,exp2[c(60489:60494),])
#write.table(counts_01A, file = "LIHC_counts_mRNA_01A_TCIA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#write.table(exp, file = "LIHC_fpkm_mRNA_01A_TCIA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

rownames(exp)[19621:19626] #shape;tumor;peritumor;peritumor enhancement; TTPVI; pseudocapsule
gene <- "shape" #Enter the above variable name

med=median(as.numeric(exp[gene,])) #When Enter shape;tumor;peritumor 
conditions=data.frame(sample=colnames(exp),
                      group=factor(ifelse(exp[gene,]>=med,"high","low"),levels = c("low","high"))) %>% 
  column_to_rownames("sample")

# Enter others
# conditions=data.frame(sample=colnames(exp),
#                       group=factor(ifelse(exp[gene,]>=1,"high","low"),levels = c("low","high"))) %>% 
#   column_to_rownames("sample")

dds <- DESeqDataSetFromMatrix(
  countData = counts_01A[c(1:19620),],
  colData = conditions,
  design = ~ group)

dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds)
# save(res,file="res_deseq_PE.Rda")
res_deseq2 <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)

#### GO/KEGG analysis  ####

library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
library(readxl)

#In this file, We take intersection and union of DEGs.
gene <- read_excel("DEGs.xlsx")
gene <- gene$Gene[1:623] 

#'DEG' refers to the DEGs results of tumors and liver parenchyma.
DEG <- as.data.frame(res)%>% 
  arrange(padj) %>% 
  dplyr::filter(abs(log2FoldChange) > 0, padj < 0.05)
DEG <- DEG[gene,]

DEG <- DEG %>% rownames_to_column("Gene")

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))

#GO
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                readable = TRUE)

ego_res <- ego@result
#save(ego,ego_res,file = "GO_radiomics_DEG.Rdata")
dotplot(ego, showCategory = 20)

#KEGG
#install.packages("R.utils")
R.utils::setOption("clusterProfiler.download.method",'auto') 
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

kk_res <- kk@result
#save(kk,kk_res,file = "KEGG_radiomics_DEG.Rdata")
dotplot(kk, showCategory = 20)


####  ssGSEA  ####
library(tidyverse)
library(data.table)
library(GSVA)

#In this file, the genes sets of each immune cells were defined.
cellMarker <- data.table::fread("cellMarker.csv",data.table = F)
colnames(cellMarker)[2] <- "celltype"
type <- split(cellMarker,cellMarker$celltype)

cellMarker <- lapply(type, function(x){
  dd = x$Metagene
  unique(dd)
})

#save(cellMarker,file = "cellMarker_ssGSEA.Rdata")
#load("cellMarker_ssGSEA.Rdata")

expr <- data.table::fread("LIHC_fpkm_mRNA_01A_TCIA.txt",data.table = F)   
rownames(expr) <- expr[,1]   
expr <- expr[,-1]   
expr <- expr[c(1:19620),]
expr <- as.matrix(expr)   

gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
a <- gsva_data %>% t() %>% as.data.frame()

#In this file, the 6 variables divided into 2 classes (continuous variable-median).
group <- read.csv("group.csv",row.names = "ID")
group <- group["shape"] #Enter the above variable name
names(group) <- "group"
name <- intersect(rownames(group),rownames(a))
a <- a[name,]

identical(rownames(a),rownames(group))
a$group <- group$group
a <- a %>% rownames_to_column("sample")
#write.table(a,"ssGSEA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
b <- gather(a,key=ssGSEA,value = Expression,-c(group,sample))
write.csv(b,"ssGSEA_shape.csv") #The enrichment score of each immune cell is distributed in 'shape' groups. 
