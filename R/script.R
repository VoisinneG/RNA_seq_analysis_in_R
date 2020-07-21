
### Import raw data #####################################################################

seqdata <- read.delim("./data/S16082_allresV2_rawdata.csv", 
                      stringsAsFactors = FALSE)

sampleinfo <- read.delim("./data/Sample_Info.csv")

### Format input data ###################################################################

# Edit countdata columns
idx_remove <- which(! names(seqdata) %in% sampleinfo$sampleName)
countdata <- seqdata[,-idx_remove]
# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]

rownames(sampleinfo) <- sampleinfo[,1]
sampleinfo <- sampleinfo[, -1]
table(colnames(countdata)==rownames(sampleinfo))



### Get gene annotations #################################################################
library(org.Hs.eg.db)
columns(org.Hs.eg.db)

annot <- AnnotationDbi::select(org.Hs.eg.db,keys=keys(org.Hs.eg.db),
                columns=c("ENTREZID", "ENSEMBL", "SYMBOL"))

idx_match <- match(row.names(countdata), annot$ENSEMBL)
annot <- annot[idx_match, ]
rownames(annot) <- rownames(countdata)


#seqdata does contain a Gene_name column. Let's keep it in annot
annot$Gene_name <- seqdata$Gene_name[match(rownames(annot), seqdata$Ensembl.gene.id)]



### filter genes with low counts or no gene symbol ######################################

#filter out genes with no symbol
#idx_no_symbol <- which(is.na(annot$SYMBOL))
#countdata <- countdata[-idx_no_symbol, ]
#annot <- annot[-idx_no_symbol, ]

#filter out genes with low counts
length( which(rowMeans(countdata) < 1) )
keep <- rowMeans(countdata) >= 1
countdata <- countdata[keep, ]

### Plot heatmap #########################################################################

library("pheatmap")
library(dplyr)

select <- order(rowMeans(countdata), decreasing=TRUE)[1:20]
select <- 1:30

pheatmap(log10(countdata[select,] + 1),
         cluster_rows=FALSE,
         cluster_cols=FALSE)


pheatmap(countdata[select,], 
         show_rownames=TRUE, 
         cluster_cols=TRUE, 
         annotation_col=sampleinfo[ , c("Treatment","Status")] ,
         scale = "row",
         labels_row = annot$Gene_name[select]
)


library(heatmaply)

heatmaply(log10(countdata[select,]+1), scale = "none", Rowv = NULL, Colv = NULL)

### Data visualization with ggplot2 ######################################################

library(reshape2)
library(dplyr)
df <- countdata
df$GeneID <- rownames(countdata)
df_melt <- melt(df, id.vars = "GeneID")
df_melt <- rename(df_melt,  sample = variable)

plot <- ggplot(df_melt, aes(x=asinh(value), fill = sample)) + 
  geom_density(alpha = 0.5)
plot


### Normalisation ######################################################################## 

df_median <- 
  df_melt %>% 
  group_by(sample) %>% 
  summarise(median = median(value, na.rm = TRUE))

df_median

median_mean <- mean(df_median$median)

df_melt <- 
  df_melt %>% 
  group_by(sample) %>% 
  mutate(value_norm = value/median(value, na.rm = TRUE)*median_mean,
         value_norm_log = log10(value_norm + 1))

ggplot(df_melt, aes(x=sample, y = value_norm_log, fill = sample)) + 
  theme(axis.text.x = element_text(angle = 90)) +
  geom_boxplot(alpha = 0.5)

data.norm.log <- dcast(df_melt, GeneID~sample, value.var = "value_norm_log")


### QC ##################################################################################

df_melt <- 
  df_melt %>% 
  group_by(GeneID, sample) %>% 
  mutate(value_log = log10(value + 1),
         value_norm_log = log10(value_norm +1))

mat <- dcast(df_melt, GeneID ~ sample, value.var = "value_norm_log")

corr_mat = cor( mat[-1] )
heatmaply(corr_mat)

library(ggpointdensity)

plist <- lapply(rownames(sampleinfo), function(sample){
  ggplot(mat, aes_string(x = "X12_54_AM_TNF", y = sample)) + 
    geom_pointdensity(size = 0.3, alpha = 0.1 ) + 
    scale_color_viridis() + 
    annotate("segment", x =0, y =0, xend = 7, yend = 7)
  
})

plist[[1]]

library(gridExtra)
marrangeGrob(plist, nrow = floor(length(plist)/3), ncol = 3)

### PCA #################################################################################
data <- mat[-1]
#data <- mat[rownames(sampleinfo)[sampleinfo$Status %in% c("quiescent", "active") & 
#                                   sampleinfo$Treatment == "NT"]]
rownames(data) <- mat[,1]

pca_res <- prcomp( t(data) )
summary(pca_res)


max_loading <- apply(pca_res$rotation, 1, max)
pca_max <- apply(pca_res$rotation, 1, which.max)

idx_select<- c(order(pca_res$rotation[, 2], decreasing = TRUE)[1:20],
               order(pca_res$rotation[, 2], decreasing = FALSE)[1:20])
annotation_row <- data.frame(pc = c(rep("up", 20), rep("down", 20) ))
rownames(annotation_row) <- rownames(data)[idx_select]

idx_select <- order(max_loading, decreasing = TRUE)[1:70]
annotation_row <- data.frame(pc = factor(pca_max[idx_select]))
rownames(annotation_row) <- rownames(data)[idx_select]
  
pheatmap(data[idx_select, ],
         show_rownames=TRUE,
         #cluster_rows = FALSE,
         #cluster_cols=TRUE,
         annotation_col=sampleinfo[ , c("Treatment","Status")],
         annotation_row = annotation_row,
         #scale = "row",
         labels_row = annot$Gene_name[ match(rownames(data)[idx_select], rownames(annot))]
)

heatmaply(data[idx_select, ], 
          scale = "row",
          labRow = annot$Gene_name[ match(rownames(data)[idx_select], rownames(annot))])

### plot pca ###########################################################################

df_pca <- as.data.frame(pca_res$x)
df_pca$sample <- rownames(df_pca)
names(df_pca)
idx_match <- match(rownames(df_pca), rownames(sampleinfo))
df_pca$Treatment <- sampleinfo$Treatment[idx_match]
df_pca$Status <- sampleinfo$Status[idx_match]

library(ggplot2)
library(ggrepel)

pca_plot <- ggplot(df_pca, aes(x=PC1, y=PC2, color = Status, label=sample)) + 
  geom_point(show.legend = FALSE) + 
  geom_text_repel(show.legend = FALSE)
pca_plot

### DE analysis #######################################################################

samples <- rownames(sampleinfo)[sampleinfo$Status %in% c("quiescent", "active") & 
                                  sampleinfo$Treatment == "NT"]

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata[,samples],
                              colData = sampleinfo[rownames(sampleinfo) %in% samples, ],
                              design = ~ Status)

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, alpha = 0.05, 
               pAdjustMethod="BH", name="Status_quiescent_vs_active")
head(res)

enrich <- data.frame(ID=rownames(res), res@listData)

enrich  <- enrich  %>% 
  filter(padj < 0.05, log2FoldChange > 2) %>%
  arrange(desc(log2FoldChange))

enrich %>% head(10)


d <- plotCounts(dds, gene=as.character(enrich$ID[1]), 
                intgroup=c("Treatment","Status"), 
                returnData=TRUE)

library("ggplot2")
p <- ggplot(d, aes(x=Status, y=asinh(count))) + 
  geom_point(position=position_jitter(w=0.1,h=0))  +
  facet_wrap(~Treatment)
p


df <- as.data.frame(res@listData)
df$ID <- res@rownames
df$names <- annot$SYMBOL[match(df$ID, annot$ENSEMBL)]
df$EntrezID <- annot$ENTREZID[match(df$ID, annot$ENSEMBL)]
select <- abs(df$log2FoldChange) > 2 & df$padj < 0.05
ggplot(df, aes(x=log2FoldChange, y=-log10(padj), label = names)) + 
  geom_point(color = "gray") +
  geom_point(data = df[select, ], color = "red") +
  geom_text_repel(data = df[select, ], color = "black", 
                  min.segment.length = 0, segment.colour = "gray")


### Annotation enrichment ##############################################################

library(clusterProfiler)

enrichgo_result_over_MF = enrichGO( gene = df$ID[select], 
                                    keyType = "ENSEMBL",
                                    universe = df$ID,
                                    OrgDb = org.Hs.eg.db,
                                    ont = "MF",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.2,
                                    readable = TRUE)

enrichgo_result_over_KEGG = enrichKEGG( gene = df$EntrezID[select], organism = "hsa",
                                    keyType = "kegg",
                                    universe = df$EntrezID,
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.2)

library(DT)
datatable( enrichgo_result_over_MF@result )
datatable( enrichgo_result_over_KEGG@result )
