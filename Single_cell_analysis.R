
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
adj.matrix <- Read10X("data/update/soupX_pbmc10k_filt")


########load and process single cell T cell data set of GSE99254 (log2 transformed)########
load("gse99254.RData")
patient_information<-read.csv("patient_information.csv")
patient_information_process<-patient_information[!is.na(patient_information$Lymphocyte.count..10e9.L.),]
dim(patient_information_process)
gse99254$geneSymbol[1:10]
rownames(gse99254)<-gse99254$geneSymbol
names(table(gse99254$geneSymbol))[table(gse99254$geneSymbol)>1]

#extract gene name
gene_name<-as.data.frame(table(gse99254$geneSymbol))
gene_name[gene_name$Freq>1,]
gse99254$geneID[is.na(gse99254$geneSymbol)]
gse99254<-gse99254[!is.na(gse99254$geneSymbol),]

#extract tissue cell subsets
gse99254<-gse99254[,c(-1,-2)]
gse99254_tissue<-gse99254[,grepl(c("TTC|TTH|TTR|TTY"),colnames(gse99254))]
dim(gse99254_tissue)
unique(colnames(gse99254_tissue))
colnames(gse99254_tissue)<-gsub("\\.","_",colnames(gse99254_tissue))

#explore the characteristics of gse99254_tissue data
srat <- CreateSeuratObject(gse99254_tissue,project = "single_cell")
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(srat), 10)
plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)

#unsupervised clustering by SNN
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
DimHeatmap(srat, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(srat, reduction = "pca")
srat <- FindNeighbors(srat, dims = 1:10)
srat <- FindClusters(srat, resolution = 0.8)
srat <- RunUMAP(srat, dims = 1:10, verbose = F)
table(srat@meta.data$seurat_clusters)
DimPlot(srat,label.size = 4,repel = T,label = T)
FeaturePlot(srat, features = c("IL-4"))
"CD8" %in% rownames(srat)
FeaturePlot(srat,"IL7R") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral"))) + ggtitle("CD8B: CD8 T cells")
all.markers <- FindAllMarkers(srat, only.pos = T, min.pct = 0.5, logfc.threshold = 0.5)
top3_markers <- as.data.frame(all.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC))

# processing patient information
colnames(gse99254_tissue)[1:10]
meta_data_gse99254<-data.frame(Cell=colnames(gse99254_tissue))
head(meta_data_gse99254)

patient_id<-c()
cell_type<-c()
for(i in 1:nrow(meta_data_gse99254)){
  cellid<-meta_data_gse99254$Cell[[i]]
  cell_id<-strsplit(cellid,"_")
  cell<-cell_id[[1]][[1]]
  pid<-cell_id[[1]][[2]]
  cell_t<-substr(cell, 1, 3)
  cell_type<-c(cell_type,cell_t)
  for (patientid in as.character(patient_information_process$Patient_ID)){
    if (grepl(patientid,cellid)){
      patient_id<-c(patient_id,patientid)
    }
  }
}

colnames(patient_information_process)[[6]]<-"Lymphocyte_count"
patient_information_process$Patient_ID<-gsub("P","",patient_information_process$Patient.ID)
patient_information_process$Patient_ID[[7]]<-"0616P"
meta_data_gse99254$Patient_ID<-patient_id
meta_data_gse99254$cell_type<-cell_type
patient_information_process<-merge(patient_information_process,meta_data_gse99254,by="Patient_ID")
View(patient_information_process)
patient_information_process$Cell[patient_information_process$Patient_ID=="0616P"]



################# analysis based on the raw count data set of GSE99254######################

install.packages("remotes")
library(remotes)
remotes::install_github("carmonalab/scGate")
remotes::install_github("carmonalab/ProjecTILs")

install.packages("renv")
install.packages("patchwork")
library(renv)
renv::restore()
library(patchwork)
library(ggplot2)
library(reshape2)
library(ProjecTILs)

#read raw count data
gse99254_count<-read.csv("GSE99254_NSCLC.TCell.S12346.count.txt",sep="\t")
gse99254_count<-gse99254_count[!is.na(gse99254_count$symbol),]
rownames(gse99254_count)<-gse99254_count$symbol
gse99254_count_tissue<-gse99254_count[,grepl(c("TTC|TTH|TTR|TTY"),colnames(gse99254_count))]
dim(gse99254_count_tissue)
meta_data_gse99254_count<-data.frame(Cell=colnames(gse99254_count_tissue))
colnames(gse99254_count_tissue)<-gsub("\\.","_",colnames(gse99254_count_tissue))
sum(gsub("\\.","_",colnames(gse99254_count_tissue)) %in% patient_information_process$Cell)
gse99254_count_tissue<-gse99254_count_tissue[,colnames(gse99254_count_tissue) %in% patient_information_process$Cell]
dim(gse99254_count_tissue)

# keep only genes with the average count>1
gse99254_count_tissue_12415_gene<-gse99254_count_tissue[rownames(gse99254_count_tissue) %in% gene_name$Var1,]
dim(gse99254_count_tissue_12415_gene)

#construct seurat object
gse99254_count_tissue_12415_gene.seurate <- CreateSeuratObject(gse99254_count_tissue_12415_gene, project = "GSE_raw_count", meta.data = patient_information_process)
rm(gse99254_count_tissue_12415_gene)

#log2 normalization
gse99254_count_tissue_12415_gene.seurate <- NormalizeData(gse99254_count_tissue_12415_gene.seurate, normalization.method = "LogNormalize", scale.factor = 10000)

#Project local data to reference atlas using ProjectTILs package

ncores = 4
ref <- readRDS("Z:/xiaojun/xiaojun_0330_2022/ky744_single_cell/Bassez_Tcell_data/ref_TILAtlas_mouse_v1.rds")
class(ref)
dim(ref)
colnames(ref)
rownames(ref)
table(ref$functional.cluster)

# split data by clinical stage
data.split <- SplitObject(gse99254_count_tissue_12415_gene.seurate, split.by = "Stage")
query.projected <- make.projection(data.split, ref, filter.cells = F, ncores = ncores)
query.projected <- lapply(query.projected, function(x) {
  cellstate.predict(ref = ref, query = x, reduction = "umap", ndim = 2)
})
query.projected.merged <- suppressMessages(Reduce(ProjecTILs:::merge.Seurat.embeddings,
                                                            query.projected))
genes <- c("Foxp3", "Cd4", "Cd8a", "Tcf7", "Ccr7", "Sell", "Il7r", "Gzmb", "Gzmk",
           "Pdcd1", "Lag3", "Havcr2", "Tox", "Entpd1", "Cxcr3", "Cxcr5", "Ifng", "Cd200",
           "Xcl1", "Itgae")

pdf("gene_expression_radar.pdf",width=15,height=15)
plot.states.radar(ref, query = query.projected.merged, min.cells = 30, genes4radar = genes)
dev.off()

query.sub.bystage <- SplitObject(query.projected.merged, split.by = "Stage")

plot.states.radar(query = list(I = query.sub.bystage$I, III=query.sub.bystage$III,IV = query.sub.bystage$IV))
plot.states.radar(ref, query = list(I = query.sub.bystage$I, III=query.sub.bystage$III,IV = query.sub.bystage$IV))

a <- plot.projection(ref, query = query.sub.bystage$I, linesize = 0.5, pointsize = 0.2)+ ggtitle("Stage I")
b <- plot.projection(ref, query = query.sub.bystage$III, linesize = 0.5, pointsize = 0.2)+ ggtitle("Stage III")
c <- plot.projection(ref, query = query.sub.bystage$IV, linesize = 0.5, pointsize = 0.2)+ ggtitle("Stage IV")

pdf("uamp_by_stage.pdf",width=16,height=12)
a | b | c
dev.off()

cell_functional_cluster<-as.data.frame(query.projected.merged$functional.cluster)
cell_functional_cluster$Cell<-rownames(cell_functional_cluster)
colnames(cell_functional_cluster)[[1]]<-"functional_cluster"
save(cell_functional_cluster,file="cell_functional_cluster.RData")

which.types <- table(query.projected.merged$functional.cluster) > 20  # disregard very small populations
stateColors_func <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc",
                      "#FF0000", "#87f6a5", "#e812dd")
states_all <- levels(ref$functional.cluster)
names(stateColors_func) <- states_all
cols_use <- stateColors_func[names(which.types)][which.types]
query.projected.merged$functional.cluster <- factor(query.projected.merged$functional.cluster,levels = states_all)

query.list <- SplitObject(query.projected.merged, split.by = "Stage")
d <- plot.statepred.composition(ref, query.list$I, metric = "percent") + ylim(0,45)+ ggtitle("Stage I")
e <- plot.statepred.composition(ref, query.list$III, metric = "percent") + ylim(0,45)+ ggtitle("Stage III")
f <- plot.statepred.composition(ref, query.list$IV, metric = "percent") + ylim(0,45)+ggtitle("Stage IV")

pdf("percentage_by_stage.pdf",width=12,height=3)
d | e |f
dev.off()

# fold change stage IV vs. stage I
norm.I <- table(query.list[["I"]]$functional.cluster)/sum(table(query.list[["I"]]$functional.cluster))
norm.III <- table(query.list[["III"]]$functional.cluster)/sum(table(query.list[["III"]]$functional.cluster))
#norm.I <- table(query.list[["I"]]$functional_cluster)/sum(table(query.list[["I"]]$functional_cluster))
norm.IV <- table(query.list[["IV"]]$functional.cluster)/sum(table(query.list[["IV"]]$functional.cluster))

foldchange <- norm.IV[which.types]/norm.I[which.types]
foldchange <- sort(foldchange, decreasing = T)

tb.m <- melt(foldchange)
colnames(tb.m) <- c("Cell_type", "Fold_change")
pll <- list()

pdf("fold_change.pdf",width=6,height=6)
ggplot(tb.m, aes(x = Cell_type, y = Fold_change, fill = Cell_type)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = cols_use) + geom_hline(yintercept = 1) + scale_y_continuous(trans = "log2") +
  theme(axis.text.x = element_blank(), legend.position = "left") + ggtitle("Stage IV vs. Stage I")+theme_light()
dev.off()



dim(cell_functional_cluster)
colnames(gse99254_count_tissue_12415_gene.seurate) %in% cell_functional_cluster$Cell
colnames(cell_functional_cluster)[[2]]<-"cell_name"
functional_cluster<-cell_functional_cluster$functional_cluster
names(functional_cluster)<-cell_functional_cluster$cell_name

gse99254_count_tissue_12415_gene.seurate <- AddMetaData(
  object = gse99254_count_tissue_12415_gene.seurate,
  metadata = functional_cluster,
  col.name = 'functional_cluster'
)

gse99254_count_tissue_12415_gene.seurate <- FindVariableFeatures(gse99254_count_tissue_12415_gene.seurate, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gse99254_count_tissue_12415_gene.seurate)
gse99254_count_tissue_12415_gene.seurate<- ScaleData(gse99254_count_tissue_12415_gene.seurate, features =all.genes)
gse99254_count_tissue_12415_gene.seurate <- RunPCA(gse99254_count_tissue_12415_gene.seurate, features = VariableFeatures(object = gse99254_count_tissue_12415_gene.seurate))
ElbowPlot(gse99254_count_tissue_12415_gene.seurate)
gse99254_count_tissue_12415_gene.seurate <- FindNeighbors(gse99254_count_tissue_12415_gene.seurate, dims = 1:20)
gse99254_count_tissue_12415_gene.seurate <- FindClusters(gse99254_count_tissue_12415_gene.seurate, resolution = 0.4)
gse99254_count_tissue_12415_gene.seurate <- RunUMAP(gse99254_count_tissue_12415_gene.seurate, dims = 1:20)
gse99254_count_tissue_12415_gene.seurate <- RunTSNE(gse99254_count_tissue_12415_gene.seurate, dims = 1:20)
DimPlot(gse99254_count_tissue_12415_gene.seurate, reduction = "umap",group.by="functional_cluster")


######################## identifying differential expression genes by DAseq package####################
devtools::install_github("KlugerLab/DAseq")
install.packages("githubinstall")
library(githubinstall)

gse99254_pca<-gse99254_count_tissue_12415_gene.seurate@reductions$pca@cell.embeddings[,1:20]
gse99254_umap<-gse99254_count_tissue_12415_gene.seurate@reductions$umap@cell.embeddings
gse99254_tsne<-gse99254_count_tissue_12415_gene.seurate@reductions$tsne@cell.embeddings

X.label.gse99254<-paste(gse99254_count_tissue_12415_gene.seurate$Patient_ID,gse99254_count_tissue_12415_gene.seurate$stage_binary,sep="_")
labels.I<-X.label.gse99254[which(gse99254_count_tissue_12415_gene.seurate$stage_binary=="I")]
labels.III_IV<-X.label.gse99254[which(gse99254_count_tissue_12415_gene.seurate$stage_binary=="III_IV")]

# stage I vs. stage III_IV
da_cells_gse99254 <- getDAcells(
  X = gse99254_pca,
  cell.labels = X.label.gse99254,
  labels.1 = labels.I,
  labels.2 = labels.III_IV,
  k.vector = seq(50, 500, 50),
  plot.embedding = gse99254_umap
)

da_cells_gse99254$pred.plot

da_cells_gse99254_update <- updateDAcells(
  X = da_cells_gse99254, pred.thres = c(-0.9,0.9),
  plot.embedding = gse99254_umap
)

da_cells_gse99254_update$da.cells.plot
table(gse99254_count_tissue_12415_gene.seurate$functional_cluster[da_cells_gse99254_update$da.up])
table(gse99254_count_tissue_12415_gene.seurate$functional_cluster[da_cells_gse99254_update$da.down])

da_regions_gse99254 <- getDAregion(
  X = gse99254_pca,
  da.cells = da_cells_gse99254_update,
  cell.labels = X.label.gse99254,
  labels.1 = labels.I,
  labels.2 = labels.III_IV,
  resolution = 0.01,
  plot.embedding = gse99254_umap
)

da_regions_gse99254$da.region.plot


##########Analysis based on centered log2 transformed data set of GSE99254#########

gse99254_tissue_stageI_IV<-gse99254_tissue[,colnames(gse99254_tissue) %in% patient_information_process$Cell[patient_information_process$Stage %in% c("I","IV")]]
dim(gse99254_tissue_stageI_IV)
table(patient_information_process$Stage)
patient_information_process_stageI_IV<-patient_information_process[colnames(gse99254_tissue_stageI_IV),]

library(FactoMineR)
library(factoextra)
patient_information_process$stage_binary<-"I"
patient_information_process$stage_binary[patient_information_process$Stage %in% c("III","IV")]<-"III_IV"

patient_information_process_stageI_IV$Cell==colnames(gse99254_tissue_stageI_IV)
gse99254_tissue.pca <- PCA(t(gse99254_tissue_stageI_IV), ncp =10, graph = FALSE)
dim(gse99254_tissue.pca$ind$coord)
gse99254_tissue.pca.all.cell <- PCA(t(gse99254_tissue), ncp =10, graph = FALSE)
dim(gse99254_tissue.pca.all.cell$ind$coord)

fviz_pca_ind(gse99254_tissue.pca,
             axes = c(1, 2),
             geom.ind = "point", # show points only (nbut not "text")
             pointshape = 21,
             pointsize = 3,
             col.ind = "black",
             fill.ind=patient_information_process_stageI_IV$cell_type,
             palette="npj",
             legend.title = "Groups"
)+ggpubr::fill_palette("npj")


library(umap)
gse99254.umap <- umap(t(gse99254_tissue_stageI_IV))
gse99254.umap$layout
class(gse99254.umap$layout)
gse99254.2d.umap<-as.matrix(gse99254.umap$layout)
colnames(gse99254.2d.umap)<-c("tSNE_1","tSNE_2")

gse99254.umap.all.cells <- umap(t(gse99254_tissue))

cell_type_table<-read.csv("cell_type_table.csv",stringsAsFactors = F)

da_cells <- getDAcells(
  X = gse99254_tissue.pca$ind$coord,
  cell.labels = rownames(gse99254_tissue.pca$ind$coord),
  labels.1 = patient_information_process_stageI_IV$Cell[patient_information_process_stageI_IV$Stage=="IV"],
  labels.2 = patient_information_process_stageI_IV$Cell[patient_information_process_stageI_IV$Stage=="I"],
  k.vector = seq(50, 500, 50),
  plot.embedding = gse99254.2d.umap
)
da_cells$da.cells.plot

da_cells.all <- getDAcells(
  X = gse99254_tissue.pca.all.cell$ind$coord,
  cell.labels = rownames(gse99254_tissue.pca.all.cell$ind$coord),
  labels.1 = patient_information_process$Cell[patient_information_process$stage_binary=="I"],
  labels.2 = patient_information_process$Cell[patient_information_process$stage_binary=="III_IV"],
  k.vector = seq(50, 500, 50),
  plot.embedding = tsne.embeds
)
da_cells.all$da.cells.plot

da_cells_select <- updateDAcells(
  X = da_cells, pred.thres = c(-0.6,0.6),
  plot.embedding = gse99254.2d.umap
)

da_cells_select.all <- updateDAcells(
  X = da_cells.all, pred.thres = c(-0.7,0.7),
  plot.embedding = tsne.embeds
)
da_cells_select.all$da.cells.plot

da_cells_select_0_7 <- updateDAcells(
  X = da_cells, pred.thres = c(-0.8,0.8),
  plot.embedding = gse99254.2d.umap
)
da_cells_select_0_7$da.cells.plot

da_regions <- getDAregion(
  X = gse99254_tissue.pca$ind$coord,
  da.cells = da_cells_select,
  cell.labels = rownames(gse99254_tissue.pca$ind$coord),
  labels.1 = patient_information_process_stageI_IV$Cell[patient_information_process_stageI_IV$Stage=="IV"],
  labels.2 = patient_information_process_stageI_IV$Cell[patient_information_process_stageI_IV$Stage=="I"],
  resolution = 0.08,
  size = 0.8,
  plot.embedding = gse99254.2d.umap
)

da_regions.all <- getDAregion(
  X = gse99254_tissue.pca.all.cell$ind$coord,
  da.cells = da_cells_select.all,
  cell.labels = rownames(gse99254_tissue.pca.all.cell$ind$coord),
  labels.1 = patient_information_process$Cell[patient_information_process$stage_binary=="I"],
  labels.2 = patient_information_process$Cell[patient_information_process$stage_binary=="III_IV"],
  resolution = 0.08,
  size = 0.8,
  plot.embedding =tsne.embeds
)

da_regions.all$da.region.plot
da_regions$da.region.plot

python2use <-"C:/Users/jxf43/Anaconda3/python.exe"
STG_markers <- STGmarkerFinder(
  X = gse99254_tissue_stageI_IV,
  da.regions = da_regions,
  lambda = 1.5, n.runs = 5, return.model = T,
  python.use = python2use, GPU = GPU
)

save(gse99254_tissue_stageI_IV,da_regions,file="DA_gene_search.RData")

gse99254.2d.umap.dataframe<-as.data.frame(gse99254.2d.umap)
gse99254.2d.umap.dataframe$Cell<-rownames(gse99254.2d.umap.dataframe)
gse99254.2d.umap.dataframe<-merge(gse99254.2d.umap.dataframe,patient_information_process_stageI_IV,by="Cell")

library(dplyr)
gse99254.2d.umap.dataframe %>%
  ggplot(aes(x = tSNE_1,
             y = tSNE_2,
             color = cell_type,
             shape = Stage)) +
  geom_point() +
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")

################### original cell clusters of GSE99254########################
cell_type_table<-read.csv("cell_type_table.csv",stringsAsFactors = F)
cell_type_table$Cell<-gsub("-","_",cell_type_table$Cell)
rownames(patient_information_process)<-patient_information_process$Cell
patient_information_process<-merge(patient_information_process,cell_type_table[,c("Cell","cluster")],by="Cell")

gse99254_seurat <- CreateSeuratObject(counts = gse99254_tissue, project = "gse99254", min.cells = 3, min.features = 200,meta.data = patient_information_process)

plot(GetAssayData(object = gse99254_seurat, slot = "counts")[1,],GetAssayData(object = gse99254_seurat, slot = "scale.data")[1,])
gse99254_seurat <- FindVariableFeatures(gse99254_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(gse99254_seurat)
gse99254_seurat<- ScaleData(gse99254_seurat, features =all.genes)
gse99254_seurat <- RunPCA(gse99254_seurat, features = VariableFeatures(object = gse99254_seurat))
ElbowPlot(gse99254_seurat)
gse99254_seurat <- FindNeighbors(gse99254_seurat, dims = 1:10)
gse99254_seurat <- FindClusters(gse99254_seurat, resolution = 0.7)
gse99254_seurat <- RunUMAP(gse99254_seurat, dims = 1:10)
gse99254_seurat <- RunTSNE(gse99254_seurat, dims = 1:10)
DimPlot(gse99254_seurat, reduction = "tsne",group.by="cluster")
average_expression_12_cluster<-AverageExpression(gse99254_seurat,slot="scale.data")

cell_class<-as.data.frame(gse99254_seurat$seurat_clusters)
colnames(cell_class)[[1]]<-"seurat_cluster"
cell_class$Cell<-rownames(cell_class)
cell_class$Cell %in% patient_information_process$Cell
patient_information_process_12cluster<-merge(patient_information_process,cell_class,by="Cell")
table(patient_information_process_12cluster$Stage,patient_information_process_12cluster$cluster)
table(patient_information_process_12cluster$Stage)
DimPlot(gse99254_seurat, reduction="tsne",label = TRUE, split.by = "stage_binary", ncol = 2)
DimPlot(gse99254_seurat, reduction = "umap",label = TRUE, split.by = "stage_binary", ncol = 2)
cluster1.markers <- FindMarkers(gse99254_seurat, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

gse99254_seurat.markers <- FindAllMarkers(gse99254_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
marker<-as.data.frame(gse99254_seurat.markers %>%
                        group_by(cluster) %>%
                        slice_min(n = 10, order_by =  p_val_adj))
marker[marker$cluster==8,]

marker[marker$gene %in% c("CCR3","CCR4","CCR8","CXCR4","STAT6"),]


######cell type assigning by SCINA package and add cell type to patient_information meta data#####
install.packages('SCINA')
library('SCINA')

log2_expression<-GetAssayData(object = gse99254_seurat, slot = "counts")
T_signature<-read.csv("Tcell_13type.csv",stringsAsFactors = F)
signatures=preprocess.signatures("Tcell_13type.csv")
results = SCINA(log2_expression, signatures, max_iter = 50, convergence_n = 10,
                convergence_rate = 0.999, sensitivity_cutoff = 0.9, rm_overlap=F, allow_unknown=F, log_file='SCINA.log')
table(results$cell_labels)
View(results$probabilities)
cell_class<-as.data.frame(results$cell_labels)
colnames(cell_class)[[1]]<-"seurat_cluster"
cell_class$Cell<-rownames(cell_class)
cell_class$Cell %in% patient_information_process$Cell
patient_information_process$test_class<-results$cell_labels
table(patient_information_process$Stage,patient_information_process$test_class)
table(patient_information_process$Stage)
