BiocManager::install("dittoSeq")
BiocManager::install("SingleR")
BiocManager::install("slingshot")
BiocManager::install("celldex")
BiocManager::install("harmony")
BiocManager::install("SeuratData")
install.packages("svglite")

library("data.table")
library("rlang")
library("Seurat")
library('dplyr')
library('dittoSeq')
library(SingleR)
library(slingshot)
library(celldex)
library(dplyr)
library(ggplot2)
library(scales)
library(svglite)

##### Set pathway #####
setwd('/home/sbm/kyungmin/HYIBB/scRNA/Data/GSE130157/Data/')

#### Read dataset ####

GI_count_1 <- fread("GSE130157_14546_14554_14562_RawCounts.txt", sep = 'auto', header=T, , na.strings = c("",NA), check.names = F)
GI_count_2 <- fread("GSE130157_15424R_15435R_RawCounts.txt", sep = 'auto', header=T, na.strings = c("",NA), check.names = F)

GI_count_1_new <- na.omit(GI_count_1)
GI_count_2_new <- na.omit(GI_count_2)

GI_count_1_new <- GI_count_1_new[!duplicated(GI_count_1_new$`Gene Symbol`)]
GI_count_2_new <- GI_count_2_new[!duplicated(GI_count_2_new$`Gene Symbol`)]

GI_count_1_new <- as.data.frame(GI_count_1_new)
GI_count_2_new <- as.data.frame(GI_count_2_new)

rownames(GI_count_1_new) <- GI_count_1_new$`Gene Symbol`
rownames(GI_count_2_new) <- GI_count_2_new$`Gene Symbol`

GI_count_1_new <- GI_count_1_new[, -1:-3]
GI_count_2_new <- GI_count_2_new[, -1:-3]


# Read metadata
metadata <- fread("GSE130157.cell_annotations.txt", sep = 'auto', header=T,  check.names = F ,na.strings = c("",NA),) 
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$`Cell ID`


Genes_intersect_1 <- intersect(colnames(GI_count_1_new) , metadata$`Cell ID`)
Genes_intersect_2 <- intersect(colnames(GI_count_2_new) , metadata$`Cell ID`)

# Retrive intersecting cells
GI_count_1_new <- GI_count_1_new[,Genes_intersect_1]
GI_count_2_new <- GI_count_2_new[,Genes_intersect_2]

GI_count <- cbind(GI_count_1_new,GI_count_2_new)



##### Create seurat object #####
GI_seurat <- CreateSeuratObject(counts = GI_count , meta.data = metadata , project = "GI" )
rm(GI_count_1, GI_count_2, GI_count_1_new, GI_count_2_new, Genes_intersect_1, Genes_intersect_2, GI_count)


#### Quality control ####

GI_seurat[["percent.mito"]] <- PercentageFeatureSet(GI_seurat, pattern = "^MT-")
GI_seurat[["percent.ribo"]] <- PercentageFeatureSet(GI_seurat, pattern = "^RPS|^RPL")


# Check barcode rank plot
barcode <- colnames(GI_seurat)
umis <- colSums(GI_seurat@assays$RNA)
barcode_umi_df <- data.frame(Barcode = barcode, UMI_Count = umis)

sorted_barcodes <- names(sort(umis, decreasing = TRUE))
barcode_rank_df <- data.frame(Barcode = sorted_barcodes,  UMI_Count = umis[sorted_barcodes])
barcode_rank_df$rank <- 1:nrow(barcode_rank_df)

key <- c("nCount_RNA_min" , "nCount_RNA_max" , "nFeature_RNA_min", "nFeature_RNA_max")
value <- c(min(GI_seurat$nCount_RNA), max(GI_seurat$nCount_RNA), min(GI_seurat$nFeature_RNA), max(GI_seurat$nFeature_RNA))
RNA_value <- data.frame(key, value)

barcode_rank_plot <-ggplot(barcode_rank_df , aes(x=rank, y = UMI_Count  )) +
  geom_point() +
  scale_x_log10(breaks = 10^(0:5), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "Waterfall plot of read counts (log)", x = "Rank", y = "UMI counts") +
  theme_light()+
  geom_hline(yintercept = c(min(GI_seurat$nCount_RNA),max(GI_seurat$nCount_RNA)) , linetype='dashed', color='red', size=0.3)+
  theme(axis.text=element_text(size=16,face = "bold"), axis.title=element_text(size=20,face="bold"))+
  theme(plot.title = element_text(size = 30))

barcode_rank_plot


QC_plot <- VlnPlot(GI_seurat, features = c("nFeature_RNA", "nCount_RNA" , "percent.mito", "percent.ribo"), ncol = 2, pt.size=0 , group.by ="Patient.ID" )

# pdf(file = "barcode_plot.pdf", width = 4.33, height = 6.75)
barcode_rank_plot+QC_plot
# dev.off()

GI_seurat_QC <- subset(GI_seurat, subset = percent.mito<10 & nFeature_RNA>450)

dittoScatterPlot(x.var = "nCount_RNA",  y.var = "nFeature_RNA",  color.var = "percent.mito" , 
                 object = GI_seurat_QC, cells.use = colnames(GI_seurat), show.others = TRUE)

barcode_QC <- colnames(GI_seurat_QC)
umis_QC <- colSums(GI_seurat_QC@assays$RNA)
barcode_umi_df_QC <- data.frame(Barcode = barcode_QC, UMI_Count = umis_QC)

sorted_barcodes <- names(sort(umis_QC, decreasing = TRUE))
barcode_rank_df_QC <- data.frame(Barcode = sorted_barcodes,  UMI_Count = umis[sorted_barcodes])
barcode_rank_df_QC$rank <- 1:nrow(barcode_rank_df_QC)

barcode_rank_plot_QC <-ggplot(barcode_rank_df_QC , aes(x=rank, y = UMI_Count  )) +
  geom_point() +
  scale_x_log10(breaks = 10^(0:5), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  labs(title = "Waterfall plot of read counts (log)", x = "Rank", y = "UMI counts") +
  theme_light()+
  geom_hline(yintercept = c(min(GI_seurat_QC$nCount_RNA), max(GI_seurat_QC$nCount_RNA)) , linetype='dashed', color='red', size=0.3)+
  theme(axis.text=element_text(size=16,face = "bold"), axis.title=element_text(size=20,face="bold"))+
  theme(plot.title = element_text(size = 30))

QC_plot <- VlnPlot(GI_seurat_QC, features = c("nFeature_RNA", "nCount_RNA" , "percent.mito", "percent.ribo"), ncol = 2, pt.size=0 , group.by ="Patient.ID" )

barcode_rank_plot_QC+QC_plot


#### Analysis ####
GI_seurat_QC <- NormalizeData(GI_seurat_QC, verbose = FALSE)
GI_seurat_QC <- CellCycleScoring(object = GI_seurat_QC, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
GI_seurat_QC <- FindVariableFeatures(GI_seurat_QC, selection.method = "vst", nfeatures = 2000, verbose = TRUE)

features <- GI_seurat_QC@assays$RNA@var.features
GI_seurat_QC <- ScaleData(GI_seurat_QC, features = features,  vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))

##### Perform linear dimensional reduction #####
GI_seurat_QC = RunPCA(GI_seurat_QC, npcs = 40, verbose = TRUE , features = features)
DimPlot(GI_seurat_QC, reduction = "pca", na.value = "grey50" , group.by = "Patient.ID")
ElbowPlot(GI_seurat_QC,ndims=40)
VizDimLoadings(GI_seurat_QC, dims = 1:2, reduction = "pca")

##### Run Umap #####
GI_seurat_QC <- RunUMAP(GI_seurat_QC, reduction = "pca",dims = 1:20) 
GI_seurat_QC <- FindNeighbors(GI_seurat_QC, reduction = "pca", dims = 1:20)
GI_seurat_QC <- FindClusters(GI_seurat_QC, resolution = 0.4)
DimPlot(GI_seurat_QC, reduction = "umap"  , group.by = "DICE")
DimPlot(GI_seurat_QC, reduction = "umap" , group.by = c("Cancer", "Patient.ID" , "Responder", "Time.Point")  , ncol = 2 , pt.size = 0.05) # Batch effect low level

FeaturePlot(GI_seurat_QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "S.Score","G2M.Score") , coord.fixed = T, ncol= 3)

##### Cell type annotation (DICE reference) #####
DICE_ref <- readRDS('/home/sbm/kyungmin/HYIBB/scRNA/Data/Annotation_ref/DICE_ref.rds')
DICE <- SingleR(test = as.SingleCellExperiment(GI_seurat_QC),  ref = DICE_ref, labels = DICE_ref$label.fine)

# Grab labels, but change pruned NAs to character "NA"s
GI_seurat_QC@meta.data$DICE <- DICE$labels
dittoDimPlot(GI_seurat_QC, "DICE",
             reduction.use = "umap",
             do.letter = FALSE,
             cells.use = !is.na(DICE$pruned.labels),
             do.label = F  ,
             split.by = "DICE")


dittoDimPlot(GI_seurat_QC, "Sub.Cluster",
             reduction.use = "umap",
             do.letter = FALSE,
             do.label = F  ,
             split.by = "Sub.Cluster")

setwd('/home/sbm/kyungmin/HYIBB/scRNA/Data/GSE130157/Result/R/')
# saveRDS(GI_seurat_QC, "GI_seurat_QC.rds")
GI_seurat_QC <- readRDS('/home/sbm/kyungmin/HYIBB/scRNA/Data/GSE130157/Result/R/GI_seurat_QC.rds')

##### Subset only naive CD4+ T cells (annotation from metadata) #####
GI_naive_CD4_T <- subset(GI_seurat_QC, subset = DICE == "T cells, CD4+, naive")

# Analysis (performed batch correction)

##### Identifies variable features #####
GI_naive_CD4_T.list <- SplitObject(GI_naive_CD4_T, split.by = "Patient.ID")
for (i in 1:length(GI_naive_CD4_T.list)) {
  GI_naive_CD4_T.list[[i]] <- FindVariableFeatures(GI_naive_CD4_T.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = TRUE)
}
features_naive <- SelectIntegrationFeatures(GI_naive_CD4_T.list)

GI_naive_CD4_T <- ScaleData(GI_naive_CD4_T, features = features_naive,  vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score"))

GI_naive_CD4_T = RunPCA(GI_naive_CD4_T, npcs = 40, verbose = TRUE , features = features)
DimPlot(GI_naive_CD4_T, reduction = "pca", na.value = "grey50" , group.by = "Patient.ID")
ElbowPlot(GI_naive_CD4_T,ndims=40)
VizDimLoadings(GI_naive_CD4_T, dims = 1:2, reduction = "pca")

##### Run Umap #####
GI_naive_CD4_T <- RunUMAP(GI_naive_CD4_T, reduction = "pca",dims = 1:20) 
GI_naive_CD4_T <- FindNeighbors(GI_naive_CD4_T, reduction = "pca", dims = 1:20)
GI_naive_CD4_T <- FindClusters(GI_naive_CD4_T, resolution = 0.4)
DimPlot(GI_naive_CD4_T, reduction = "umap")
DimPlot(GI_naive_CD4_T, reduction = "umap" , group.by = c("Cancer", "Patient.ID" , "Responder", "Time.Point")  , ncol = 2 , pt.size = 0.1)

FeaturePlot(GI_naive_CD4_T, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.ribo", "S.Score","G2M.Score") , coord.fixed = T, ncol= 3)


##### Perform Label Transfer #####
#### 1. Ortholog Mapping (mouse data) ####

##Cell type classification using an integrated reference (project the PCA structure of a reference onto the query)
Tcell.query <- GI_naive_CD4_T
Tcell.anchors <- FindTransferAnchors(reference = Naive_ortholog_final_mouse, query = GI_naive_CD4_T,
                                     dims = 1:30, reference.reduction = "pca" )

predictions <- TransferData(anchorset = Tcell.anchors, refdata = Naive_ortholog_final_mouse$RNA_snn_res.0.2,
                            dims = 1:30)

Tcell.query <- AddMetaData(Tcell.query, metadata = predictions)
Idents(Tcell.query) <- "predicted.id"
dittoDimPlot(Tcell.query, "naive_cluster", reduction.use = "umap")

p <- dittoDimPlot(Tcell.query, c("predicted.id"), reduction.use = "umap" , size = 1.5, main = "GSE130157" , legend.show = F)

p + theme(
  panel.border = element_blank(),  # 박스 숨기기
  axis.line = element_blank(),     # 축 라인 숨기기
  axis.text = element_blank(),     # 축 텍스트 숨기기
  axis.ticks = element_blank(),    # 축 눈금 숨기기
  axis.title = element_blank(),    # 축 제목 숨기기
  plot.title = element_text(hjust = 0.5)  # 타이틀 가운데 정렬
)

Tcell.query_meta <- as.data.frame(Tcell.query$predicted.id)
Tcell.query_meta[Tcell.query_meta[ , 1] == "0" ,1] <- "naive_c0"
Tcell.query_meta[Tcell.query_meta[ , 1] == "1" ,1] <- "naive_c1"
Tcell.query_meta[Tcell.query_meta[ , 1] == "2" ,1] <- "naive_c2"
Tcell.query_meta[Tcell.query_meta[ , 1] == "3" ,1] <- "naive_c3"

Tcell.query$naive_cluster <- Tcell.query_meta$`Tcell.query$predicted.id`

##### Divide into before/after treatment patient #####

Before_Tcell_query <- subset(Tcell.query , subset = Time.Point =="C1")
After_Tcell_query <- subset(Tcell.query , subset = Time.Point !="C1")

After_C5_Tcell_query <- subset(Tcell.query , subset = Time.Point =="C5")
After_C13_15_Tcell_query <- subset(Tcell.query , subset = Time.Point %in% c("C13", "C15"))

patient_C5_table <- as.data.frame(cbind(as.matrix(table(After_C5_Tcell_query$Patient.ID , After_C5_Tcell_query$predicted.id)), 'C5'))
patient_C13_C15_table <- as.data.frame(cbind(as.matrix(table(After_C13_15_Tcell_query$Patient.ID , After_C13_15_Tcell_query$predicted.id)), 'C13,15'))

##### Cell composition #####

Before_Tcell_query_metadata <- Before_Tcell_query[[]]
After_Tcell_query_metadata <- After_Tcell_query[[]]


# 1. Drug Response
# Before
Before_responder <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$Responder=="Responder" ,]
Before_non_responder <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$Responder=="Non.Responder" ,]

Before_response <- rbind(cbind(as.matrix(table(Before_non_responder$predicted.id)), c(0:3), 'Non-responder'),
                     cbind(as.matrix(table(Before_responder$predicted.id)), c(0:3), 'Responder'))

Before_response <- as.data.frame(Before_response)
colnames(Before_response) <- c('Proportion',  'cluster', 'Response')
Before_response$Proportion <- as.numeric(Before_response$Proportion)
Before_response$cluster <- factor(Before_response$cluster, levels=c(0:3))

# plot Stacked + percent bargraph
ggplot(Before_response, aes(fill=cluster, y=Proportion, x=Response)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = dittoColors()[1:9])+theme(axis.text=element_text(size=14),
                                                                                                    axis.title=element_text(size=16,face="bold"))


# After
After_responder <- After_Tcell_query_metadata[After_Tcell_query_metadata$Responder=="Responder" ,]
After_non_responder <- After_Tcell_query_metadata[After_Tcell_query_metadata$Responder=="Non.Responder" ,]

After_response <- rbind(cbind(as.matrix(table(After_non_responder$predicted.id)), c(0:3), 'Non-responder'),
                         cbind(as.matrix(table(After_responder$predicted.id)), c(0:3), 'Responder'))

After_response <- as.data.frame(After_response)
colnames(After_response) <- c('Proportion',  'cluster', 'Response')
After_response$Proportion <- as.numeric(After_response$Proportion)
After_response$cluster <- factor(After_response$cluster, levels=c(0:3))

# plot Stacked + percent bargraph
ggplot(After_response, aes(fill=cluster, y=Proportion, x=Response)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = dittoColors()[1:9])+theme(axis.text=element_text(size=14),
                                                                                                    axis.title=element_text(size=16,face="bold"))



# 2. C1/Naive T proportion per patient

patient_before <- Tcell.query_metadata[Tcell.query_metadata$Time.Point=="C1" ,]
patient_after  <- Tcell.query_metadata[Tcell.query_metadata$Time.Point!="C1" ,]


patient_before_table <- as.data.frame(cbind(as.matrix(table(patient_before$Patient.ID , patient_before$predicted.id)), 'Before'))
patient_after_table  <- as.data.frame(cbind(as.matrix(table(patient_after$Patient.ID , patient_after$predicted.id)), 'After'))

setwd('/home/sbm/kyungmin/HYIBB/scRNA/Data/GSE130157/Result/R/')

write.table(patient_before_table , "patient_before_table.tsv", sep = '\t' , col.names = TRUE, row.names = TRUE)
write.table(patient_after_table , "patient_after_table.tsv", sep = '\t' , col.names = TRUE, row.names = TRUE)

# 1. Drug Response (Th1 , Treg CD4)

# subset only CD4 cells


GI_seurat_QC$tosubset_CD4 <- grepl("CD4", GI_seurat_QC$DICE)
GI_seurat_QC_CD4 <- subset(GI_seurat_QC, subset = tosubset_CD4=="TRUE")

# Add predicted id
Tcell.query <- GI_seurat_QC_CD4 
Tcell.anchors <- FindTransferAnchors(reference = Naive_ortholog_final_mouse, query = GI_naive_CD4_T,
                                     dims = 1:30, reference.reduction = "pca" )

predictions <- TransferData(anchorset = Tcell.anchors, refdata = Naive_ortholog_final_mouse$RNA_snn_res.0.2,
                            dims = 1:30)

Tcell.query <- AddMetaData(Tcell.query, metadata = predictions)
Idents(Tcell.query) <- "predicted.id"





  
  

#### DEG analysis ####
DEG_respondse  <- FindMarkers(Tcell.query, ident.1 = 'Non.Responder', ident.2 = 'Responder', group.by = 'Responder', subset.ident = "1",
                             test.use = "wilcox" , min.pct = 0.1, logfc.threshold = 0)

DEG_naive_all  <-FindAllMarkers(Tcell.query, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, assay="RNA")


VlnPlot(Tcell.query , c("IL7R" , "ADGRE5" , "PIM1" , "KLF3", "RIPOR2" , "ETS1") , group.by = "Responder" , pt.size =0 , ncol =3)

#write.csv(DEG_naive_all,'DEG_naive_all.csv',  quote = F)



##### Generate contingency table for predicted clusters #####

testresult <- data.frame()
odd_ratio_result  <- data.frame()

for(i in 1:length(unique(Idents(Tcell.query)))-1){
  
  A <- nrow(Before_Tcell_query_metadata[Before_Tcell_query_metadata$Responder=="Responder" & Before_Tcell_query_metadata$predicted.id==i,])
  B <- nrow(Before_Tcell_query_metadata[Before_Tcell_query_metadata$Responder=="Responder" & Before_Tcell_query_metadata$predicted.id!=i,])
  C <- nrow(Before_Tcell_query_metadata[Before_Tcell_query_metadata$Responder=="Non.Responder" & Before_Tcell_query_metadata$predicted.id==i,])
  D <- nrow(Before_Tcell_query_metadata[Before_Tcell_query_metadata$Responder=="Non.Responder" & Before_Tcell_query_metadata$predicted.id!=i,])
  
  
  contingency_table_tmp    <- matrix(c(A, B, C, D), nrow = 2)
  colnames(contingency_table_tmp) <- c("Responder", "Non.Responder")
  rownames(contingency_table_tmp) <- c("In cluster", "Not in cluster")
  if(A/(A+B)<C/(C+D)){
    fisher_result <- fisher.test(contingency_table_tmp, alternative = "less")
  }
  else{
    fisher_result <- fisher.test(contingency_table_tmp, alternative = "greater")
  }
  p_value <- fisher_result$p.value
  log10pvalue <- -log10(p_value)
  testresult[i+1,1] <- log10pvalue
  
  odd_ratio <- (A/B)/(C/D)
  odd_ratio_result[i+1,1] <- odd_ratio
}

testresult$cluster <- paste0("naive_c",0:3)
colnames(testresult)[1] <- "-log10p-value"

# Plot the barplot
bp <- barplot(testresult$`-log10p-value`, names.arg = testresult$cluster, xlab = "naive T cell cluster",
              ylab = "-log10p-value", main = "Significance score", ylim = c(0, max(testresult$`-log10p-value`)))
abline(h = 1.3 , col = "red")



# 1. Proportion of all T cells
# Before
Before_Treg       <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD4+, memory TREG" ,]
Before_naive_CD4  <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD4+, naive" ,      ]
Before_naive_Treg <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD4+, naive TREG" , ]
Before_TFH        <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD4+, TFH" ,        ]
Before_Th1        <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD4+, Th1" ,        ]
Before_Th1_17     <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD4+, Th1_17" ,     ]
Before_Th17       <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD4+, Th17" ,       ]
Before_Th2        <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD4+, Th2" ,        ]
Before_naive_CD8  <- Before_Tcell_query_metadata[Before_Tcell_query_metadata$DICE=="T cells, CD8+, naive" ,      ]



Before_T_cells <- rbind(cbind(as.matrix(table(Before_non_responder$DICE)), c(0:6), 'Non-responder'),
                        cbind(as.matrix(table(Before_responder$DICE))    , c(0:8),    'Responder'))

Before_T_cells <- as.data.frame(Before_T_cells)
colnames(Before_T_cells) <- c('Proportion',  'cluster', 'Response')
Before_T_cells$Proportion <- as.numeric(Before_T_cells$Proportion)
Before_T_cells$cluster <- factor(Before_T_cells$cluster, levels=c(0:8))

# plot Stacked + percent bargraph
ggplot(Before_T_cells, aes(fill=cluster, y=Proportion, x=Response)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual(values = dittoColors()[1:9])+theme(axis.text=element_text(size=14),
                                                                                                    axis.title=element_text(size=16,face="bold"))


# Test
tmp <- table(Idents(After_Tcell_query), After_Tcell_query$Responder)
tmp <- as.data.frame(tmp)
tmp$Var1 <- as.character(tmp$Var1)
colnames(tmp) <- c('Cluster',  'Response', 'Freq')

ggplot(tmp, aes(x = Response, y = Freq, fill = Cluster)) +
  theme_bw(base_size = 10) +
  geom_col(position = "fill", width = 0.7) +
  xlab("Response State") +
  ylab("Proportion") +
  scale_fill_manual(values = dittoColors()[1:9]) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold")
  )

  