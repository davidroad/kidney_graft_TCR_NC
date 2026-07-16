# Canonical manuscript-scale integration and TCR analysis program.
#
# This is the original monolithic entry point. For a portable rebuild, first
# run the scripts in raw_to_intermediate/ and source_data_export/ with the
# documented manifest. Set SCVDJ_RAW_ROOT and SCVDJ_WORK_ROOT to local roots;
# no server-specific paths are required by this public copy.

library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(sctransform)
#library(monocle)
library(ComplexHeatmap)
library(circlize)
library(future)
library(Platypus)
library(tidyverse)
library(DoubletFinder)
# Use the HDF5 support provided by the active R environment. The original
# server-specific dyn.load calls are intentionally omitted from the public
# copy because they are not portable across R installations.
library(hdf5r)
library(scRepertoire)

scvdj_raw_root <- Sys.getenv("SCVDJ_RAW_ROOT", unset = "data/raw")
scvdj_work_root <- Sys.getenv("SCVDJ_WORK_ROOT", unset = ".")
###setting functions
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
    tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))

    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    tcr <- tcr[!duplicated(tcr$barcode), ]

    # Only keep the barcode and clonotype columns. 
    # We'll get additional clonotype info from the clonotype table.
    tcr <- tcr[,c("barcode", "raw_clonotype_id")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

    # Clonotype-centric info.
    clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))

    # Slap the AA sequences onto our original table by clonotype_id.
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
    names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"

    # Reorder so barcodes are first column and set them as rownames.
    tcr <- tcr[, c(2,1,3)]
    rownames(tcr) <- tcr[,1]
    tcr[,1] <- NULL
    colnames(tcr) <- paste(type, colnames(tcr), sep="_")
    # Add to the Seurat object's metadata.
    clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
    return(clono_seurat)
}
find_doublet <- function(samp,sample_name){
tmp_object  <- NormalizeData(samp)
tmp_object  <- FindVariableFeatures(tmp_object,selection.method = "vst", nfeatures = 2000)
tmp_object  <- ScaleData(tmp_object)
tmp_object  <- RunPCA(tmp_object)
tmp_object  <- FindNeighbors(tmp_object, dims = 1:20)
tmp_object  <- FindClusters(tmp_object, resolution = 0.5)
tmp_object  <- RunUMAP(tmp_object, dims = 1:10)
#DimPlot(tmp_object)
#----#v1-------------
sweep.res.list <- paramSweep(tmp_object, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
print("Testing modelHomotypic...")

homotypic.prop <- modelHomotypic(tmp_object@meta.data$seurat_clusters)          ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
##infer doublet rate from empirical data https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x
print("modelHomotypic OK")
doublet_rate <- nrow(tmp_object@meta.data)/500*0.4/100
#nExp_poi <- round(0.055*nrow(MEL3_ECM@meta.data)) 
#nExp_poi <- round(0.01*nrow(MEL3_BM@meta.data)) 
#nExp_poi <- round(0.008*nrow(RCC1_ECM@meta.data)) 
#nExp_poi <- round(0.003*nrow(RCC2_P@meta.data))
#nExp_poi <- round(0.02*nrow(tmp_object@meta.data))
nExp_poi <- round(doublet_rate*nrow(tmp_object@meta.data))
print("nExp_poi OK")

##500 0.4%
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

print("nExp_poi.adj OK")

print("Testing pK_optimal extraction...")
print("Checking bcmvn structure...")
print(str(bcmvn))
print(class(bcmvn))
print("BCmetric column:")
print(class(bcmvn$BCmetric))
print(bcmvn$BCmetric)

bc_values <- as.numeric(bcmvn$BCmetric)
pk_values <- as.character(bcmvn$pK)

# 手动找最大值的位置
max_value <- max(bc_values)
for(i in 1:length(bc_values)) {
    if(bc_values[i] == max_value) {
        max_idx <- i
        break
    }
}
# 在调用doubletFinder之前，检查所有参数
pK_optimal <- as.numeric(pk_values[max_idx])

print("=== Checking doubletFinder parameters ===")
print(paste("pK_optimal:", pK_optimal))
print(paste("nExp_poi:", nExp_poi))
print(paste("Class of pK_optimal:", class(pK_optimal)))
print(paste("Class of nExp_poi:", class(nExp_poi)))
print(paste("Is pK_optimal finite:", is.finite(pK_optimal)))
print(paste("Is nExp_poi finite:", is.finite(nExp_poi)))

print("test3")
tryCatch({
  tmp_object <- doubletFinder(
    tmp_object,
    PCs = 1:10,
    pN = 0.25,
    pK = as.numeric(pK_optimal),
    nExp = nExp_poi,
    sct = FALSE
  )
  print("test4 - SUCCESS!")
}, error = function(e) {
  print("=== DETAILED ERROR INFO ===")
  print(paste("Error message:", e$message))
  print(paste("Error class:", class(e)))
  print("Full error:"); print(e)
  print("=== TRACEBACK ==="); traceback()
  stop(e)   # 建议第一趟失败直接停止，避免后续再错
})
print("test4")

# 抓取 pANN 列名（字符向量）
pann_col <- grep("^pANN", colnames(tmp_object@meta.data), value = TRUE)
if (length(pann_col) < 1) stop("No pANN column found after first pass.")
pann_col <- pann_col[1]

# ---- Pass 2：复用 pANN ----
tmp_object <- doubletFinder(
  tmp_object,
  PCs = 1:10,
  pN = 0.25,
  pK = as.numeric(pK_optimal),
  nExp = nExp_poi.adj,
  reuse.pANN = pann_col,  # 这里必须是“列名字符串”
  sct = FALSE
)
doubletfinder<-tmp_object@meta.data
colnames(doubletfinder)[ncol(doubletfinder)]<-"doubletfinder"
head(doubletfinder)
write.table(doubletfinder,paste0(sample_name,"illumina_doubletfinder_",doublet_rate,".txt"),sep="\t",quote = F)
return(tmp_object)
}
options(future.globals.maxSize= 8912896000)


###load PBMC data
output_directory_1 <- file.path(scvdj_work_root, "merge_09042022_PBMC_R2_only")
setwd(output_directory_1)
load("pbmc_merge_9_5_2022_data.Rdata")
ls()
##reload the find_doublet function here to avoid error
###adding VDJ data
PBMC_0712_data <- add_clonotype(input_VDJ_directory_pbmc_3, PBMC_0712_data, "t")
PBMC_0825_data <- add_clonotype(input_VDJ_directory_pbmc_2, PBMC_0825_data, "t")
PBMC_0708_data <- add_clonotype(input_VDJ_directory_pbmc, PBMC_0708_data, "t")
PBMC_0701_data <- add_clonotype(input_VDJ_directory_pbmc_1, PBMC_0701_data, "t")
###doublet finding
PBMC1obj <- find_doublet(PBMC_0825_data,"PBMC_0825")
PBMC2obj <- find_doublet(PBMC_0712_data,"PBMC_0712")
PBMC3obj <- find_doublet(PBMC_0701_data,"PBMC_0701")
PBMC4obj <- find_doublet(PBMC_0708_data,"PBMC_0708")
###filtering doublets
PBMC1obj_filter <- subset(PBMC1obj,cells  = colnames(PBMC1obj)[PBMC1obj@meta.data[,ncol(PBMC1obj@meta.data)]=="Singlet"])
PBMC2obj_filter <- subset(PBMC2obj,cells  = colnames(PBMC2obj)[PBMC2obj@meta.data[,ncol(PBMC2obj@meta.data)]=="Singlet"])
PBMC3obj_filter <- subset(PBMC3obj,cells  = colnames(PBMC3obj)[PBMC3obj@meta.data[,ncol(PBMC3obj@meta.data)]=="Singlet"])
PBMC4obj_filter <- subset(PBMC4obj,cells  = colnames(PBMC4obj)[PBMC4obj@meta.data[,ncol(PBMC4obj@meta.data)]=="Singlet"])
###load Graft data
output_directory <- file.path(scvdj_work_root, "merge_09042022_Graft_R2_only")
setwd(output_directory)
load("graft_merge_9_5_2022_data.Rdata")
input_directory <- file.path(scvdj_raw_root, "multiplex_Graft_0708_2022_VDJ/per_sample_outs/Graft_0708_2022/count")
input_VDJ_directory <- file.path(scvdj_raw_root, "multiplex_Graft_0708_2022_VDJ/per_sample_outs/Graft_0708_2022/vdj_t")
#0701
input_directory_1 <- file.path(scvdj_raw_root, "multiplex_Graft_0701_2022_VDJ/per_sample_outs/Graft_0701_2022/count")
input_VDJ_directory_1 <- file.path(scvdj_raw_root, "multiplex_Graft_0701_2022_VDJ/per_sample_outs/Graft_0701_2022/vdj_t")
#0825
input_directory_2 <- file.path(scvdj_work_root, "08252021_Graft_out_R2only/08252021_Graft_R2only/outs/per_sample_outs/08252021_Graft_R2only/count/sample_filtered_feature_bc_matrix")
input_VDJ_directory_2 <- file.path(scvdj_work_root, "08252021_Graft_out_R2only/08252021_Graft_R2only/outs/per_sample_outs/08252021_Graft_R2only/vdj_t")
#0712 misslabeled as 0721
input_directory_3 <- file.path(scvdj_work_root, "07212021_Graft_out_R2only/07212021_Graft_R2only/outs/per_sample_outs/07212021_Graft_R2only/count/sample_filtered_feature_bc_matrix")
input_VDJ_directory_3 <- file.path(scvdj_work_root, "07212021_Graft_out_R2only/07212021_Graft_R2only/outs/per_sample_outs/07212021_Graft_R2only/vdj_t")
##0825
graft_0825_data <- add_clonotype(input_VDJ_directory_2, graft_0825_data, "t")
##0712
graft_0712_data <- add_clonotype(input_VDJ_directory_3, graft_0712_data, "t")
##0708
graft_0708_data <- add_clonotype(input_VDJ_directory, graft_0708_data, "t")
##0701
graft_0701_data <- add_clonotype(input_VDJ_directory_1, graft_0701_data, "t")
##reload the find_doublet function here to avoid error
graft1obj <- find_doublet(graft_0825_data,"graft_0825")
graft2obj <- find_doublet(graft_0712_data,"graft_0712")
graft3obj <- find_doublet(graft_0701_data,"graft_0701")
graft4obj <- find_doublet(graft_0708_data,"graft_0708")

graft1obj_filter <- subset(graft1obj,cells  = colnames(graft1obj)[graft1obj@meta.data[,ncol(graft1obj@meta.data)]=="Singlet"])
graft2obj_filter <- subset(graft2obj,cells  = colnames(graft2obj)[graft2obj@meta.data[,ncol(graft2obj@meta.data)]=="Singlet"])
graft3obj_filter <- subset(graft3obj,cells  = colnames(graft3obj)[graft3obj@meta.data[,ncol(graft3obj@meta.data)]=="Singlet"])
graft4obj_filter <- subset(graft4obj,cells  = colnames(graft4obj)[graft4obj@meta.data[,ncol(graft4obj@meta.data)]=="Singlet"])

#####merge without batch correction
object.list <- list(PBMC1obj, PBMC2obj, PBMC3obj, PBMC4obj, graft1obj, graft2obj, graft3obj, graft4obj)
merged.all <- merge(x = object.list[[1]],
                    y = object.list[2:length(object.list)], 
                    add.cell.ids = c("PBMC_0825","PBMC_0712","PBMC_0701","PBMC_0708","graft_0825","graft_0712","graft_0701","graft_0708"),project = "merged_all_eight")
###merge and batch correction using rpca
features <- SelectIntegrationFeatures(object.list = object.list)
anchors <- FindIntegrationAnchors(object.list = object.list, 
                                  anchor.features = features,
                                  reduction = "rpca") 
merged.all.batchcorrected <- IntegrateData(anchorset = anchors)

#merged_cca <- NormalizeData(merged_cca)
DefaultAssay(merged.all.batchcorrected) <- "integrated"
merged.all.batchcorrected <- ScaleData(merged.all.batchcorrected, verbose = FALSE)
merged.all.batchcorrected <- RunPCA(merged.all.batchcorrected, npcs = 50, verbose = FALSE)
merged.all.batchcorrected <- RunUMAP(merged.all.batchcorrected, dims = 1:30, reduction = "pca")
sample_names <- c("PBMC_0825","PBMC_0712","PBMC_0701","PBMC_0708",
                  "graft_0825","graft_0712","graft_0701","graft_0708")

stopifnot(length(sample_names) == 8)  # 保险起见

## 从 cell 名里取出最后一个下划线后的数字（例如 "...-1_8" 抽取到 8）
cell_ids  <- colnames(merged.all.batchcorrected)  # 或 rownames(merged.all.batchcorrected@meta.data)
suffix_id <- sub(".*_(\\d+)$", "\\1", cell_ids)
suffix_id_num <- suppressWarnings(as.integer(suffix_id))

## 建立映射并写入 meta.data
id2name <- setNames(sample_names, 1:8)
merged.all.batchcorrected$batch <- unname(id2name[as.character(suffix_id_num)])

## 快速查看分布
table(merged.all.batchcorrected$batch)

DimPlot(merged.all.batchcorrected, group.by = "batch", label = TRUE)
merged.all.batchcorrected <- FindNeighbors(object = merged.all.batchcorrected, dims = 1:20, verbose = FALSE)
merged.all.batchcorrected <- FindClusters(object = merged.all.batchcorrected, verbose = FALSE, resolution = 0.5 )

merged.all.batchcorrected$celltype.sample.indent <- paste(merged.all.batchcorrected$integrated_snn_res.0.5, merged.all.batchcorrected$batch, sep = "_")
DimPlot(merged.all.batchcorrected, group.by = "integrated_snn_res.0.5", label = TRUE)


###merged.all
rn <- rownames(merged.all@meta.data)

# 抽取到 "graft_0708" / "PBMC_0825" 这类前缀
batch_lab <- sub("^([^_]+_[^_]+)_.*$", "\\1", rn)

# 可选：检查不匹配的细胞名（防止格式不一致）
bad <- !grepl("^[^_]+_[^_]+_", rn)
if (any(bad)) {
  warning("warning",
          paste(head(rn[bad]), collapse = ", "))
  batch_lab[bad] <- NA_character_
}

# 写入 meta
merged.all$batch <- batch_lab

# 快速查看分布
table(merged.all$batch)
DefaultAssay(merged.all) <- "RNA"
merged.all[['SCT']] <- NULL
merged_combined_regressed <- SCTransform(object = merged.all,verbose = FALSE)
#merged_cca <- NormalizeData(merged_cca)
merged_combined_regressed <- FindVariableFeatures(merged_combined_regressed, selection.method = "vst", nfeatures = 2000)#2000
merged_combined_regressed <- RunPCA(object = merged_combined_regressed, ndims.print = 1:20,verbose = FALSE)
merged_combined_regressed <- FindNeighbors(object = merged_combined_regressed, dims = 1:20, verbose = FALSE)
merged_combined_regressed <- FindClusters(object = merged_combined_regressed, verbose = FALSE, resolution = 0.5 )
merged_combined_regressed <- RunUMAP(object = merged_combined_regressed, dims = 1:20,reduction.name = "UMAP")
DimPlot(merged_combined_regressed, group.by="batch", label = TRUE, reduction = "UMAP")
merged_combined_regressed$celltype.sample.indent <- paste(merged_combined_regressed$SCT_snn_res.0.5, merged_combined_regressed$batch, sep = "_")


###help me to save the merged.all and merged.all.batchcorrected object to Rdata file together.
setwd(file.path(scvdj_work_root, "merge_4_with_vdj"))

###standard_workflow
pdf("merged_batchcorrected.pdf",10,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = TRUE,group.by="batch"))
dev.off()
pdf("merged_all_regressed.pdf",10,8)
print(DimPlot(merged_combined_regressed, reduction = "UMAP", label = TRUE,group.by="batch"))
dev.off()

pdf("merged_batchcorrected_no_label.pdf",8,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = FALSE,group.by="batch"))
dev.off()
pdf("merged_all_regressed_no_label.pdf",8,8)
print(DimPlot(merged_combined_regressed, reduction = "UMAP", label = FALSE,group.by="batch"))
dev.off()

Idents(merged.all.batchcorrected) <- merged.all.batchcorrected$integrated_snn_res.0.5
pdf("merged_batchcorrected_by_cluster.pdf",12,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = TRUE))
dev.off()
pdf("merged_batchcorrected_by_cluster_no_label.pdf",10,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = FALSE))
dev.off()
pdf("merged_batchcorrected_by_cluster_split.pdf",16,8)
DimPlot(object = merged.all.batchcorrected, split.by = "batch",label = TRUE,reduction = "umap", pt.size = 0.3,ncol=4) + ggtitle(label = "UMAP")
dev.off()
pdf("merged_batchcorrected_by_cluster_split_no_label.pdf",16,8)
DimPlot(object = merged.all.batchcorrected, split.by = "batch",label = FALSE,reduction = "umap", pt.size = 0.3,ncol=4) + ggtitle(label = "UMAP")
dev.off()
Idents(merged_combined_regressed) <- merged_combined_regressed$SCT_snn_res.0.5
pdf("merged_all_regressed_by_cluster.pdf",10,8)
print(DimPlot(merged_combined_regressed, reduction = "UMAP", label = TRUE))
dev.off()

Idents(merged_combined_regressed) <- merged_combined_regressed$SCT_snn_res.0.5
pdf("merged_all_regressed_by_cluster_no_label.pdf",10,8)
print(DimPlot(merged_combined_regressed, reduction = "UMAP", label = FALSE))
dev.off()

pdf("Regress_0.5_cluster_umap_split.pdf",16,8)
DimPlot(object = merged_combined_regressed, split.by = "batch",label = TRUE,reduction = "UMAP", pt.size = 0.3,ncol=4) + ggtitle(label = "UMAP")
dev.off()

##6/7/2022 routin analysis
gene_list_0 <- c("CXCR6","SELL")
gene_list_1 <- c("IFNG", "TBX21", "MKI67", "KLRG1", "GZMA", "GZMB", "FASL", "ZEB2", "ID2", "IRF4", "TOX", "IKZF2", "LDHA", "TCF7", "FOXO1", "ID3")
gene_list_2 <- c("SELL", "IL7R", "BACH2", "PDCD1", "HAVCR2", "PRF1", "STAT3", "HIF1A", "IL1B", "NR4A1", "NR4A2", "NR4A3", "IFNGR1", "CXCR1")
gene_list_3 <- c("CR6", "CD28", "ICOS", "TNFSF4", "ID2", "S100A4", "S100A10", "ID3", "BCL2", "CD200", "TOX2", "CD69", "TNFSF8", "CXCR5")
gene_list_4 <- c("CD27", "CD40L", "KLF2", "CTLA4", "CD38", "TIGIT", "LAG3", "CD160", "CD44", "CD244", "BATF", "BATF3", "PKM", "EZH2", "BRD4")
gene_list_5 <- c("APEX1","DNMT1","DNMT3A","EOMES","GATA3","RORA","FOXP1","FOXO1","MYC","CDKN1A","FOXP3","CX3CR1","TRAC","TRBC2","TRDC","TRDV2")
gene_list_6 <- c("TRGV9","CD3D","CD3E","CD3G","CD4","CD8A","CD8B","NCR1","KLRC1","KLRB1C","CD56","CD7","CD16","CD1C","CD200R1","FCER1A","FCGR3A")
gene_list_7 <- c("SIGLEC6","CSF1R","LY6H","ITGAM","ADGRE1","ITGAX","LYZ","CD68","CD86","CD80","NKG7","GNLY","CLEC4C","CXCR2","CCR2","CD15","CD32","CCR3","CD45")
gene_list_8 <- c("EMR1","CD125","CD40","HLA- DR","CF2RA","ARG1","NOS2","CXCL22","CXCL3","CD274","CD14","MARCO","RETNLB","TGFB1","STAT1")
gene_list_9 <- c("STAT3","STAT5A","MRC1","EGR2","FN1","VEGFA","IL1B","NFKB1","HIF1A","TOX","NR4A1","NR4A22","NR4A3","IFNGR1,CD19","CD79A")
gene_list_10 <- c("VCAM11","CXCR44","TLR4","CCL18","CCL3","CCL22","IL6","IL10","MMP9","CSF2RA","IL3RA","APOE","ICOS","TNFRSF18","CXCR3","MS4A1","NCAM1","CD163","CD33","CD34","CCR5","CXCR5")
gene_list_11 <- c(
  "CCL1","CCL2","CCL3","CCL3L1","CCL4","CCL4L1","CCL5","CCL7","CCL8","CCL11",
  "CCL13","CCL14","CCL15","CCL16","CCL17","CCL18","CCL19","CCL20","CCL21","CCL22",
  "CCL23","CCL24","CCL25","CCL26","CCL27","CCL28","CXCL1","CXCL2"
)

gene_list_12 <- c(
  "CXCL3","CXCL4","CXCL5","CXCL6","CXCL7","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12",
  "CXCL13","CXCL14","CXCL16","CXCL17","XCL1","XCL2","CX3CL1","CCR1","CCR2","CCR3",
  "CCR4","CCR5","CCR6","CCR7","CCR8","CCR9","CCR10","CXCR1"
)

gene_list_13 <- c(
  "CXCR2","CXCR3","CXCR4","CXCR5","CXCR6","CXCR7","XCR1","CX3CR1","ACKR1","ACKR2",
  "ACKR4","ACKR5","KIR2DL1","KIR2DL2","KIR2DL3","KIR2DL4","KIR2DL5A","KIR2DL5B","KIR2DS1","KIR2DS2",
  "KIR2DS3","KIR2DS4","KIR2DS5","KIR3DL1","KIR3DL2","KIR3DL3","KIR3DS1","KLRC1"
)

gene_list_14 <- c(
  "KLRC2","KLRC3","KLRC4","KLRK1","KLRD1","KLRB1","KLRF1","KLRG1","NCR1","NCR2",
  "NCR3","SLAMF1","SLAMF2","SLAMF3","SLAMF4","SLAMF5","SLAMF6","SLAMF7","SLAMF8","SLAMF9",
  "FCGR3A","LILRB1","LILRB2","BCL11B","PIK3IP1","PIK3AP1"
)

#################
Idents(object = merged.all.batchcorrected)<- merged.all.batchcorrected$integrated_snn_res.0.5
merged.all.batchcorrected.markers <- FindAllMarkers(object = merged.all.batchcorrected, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

merged.all.batchcorrected.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
write.table(test,"cluster_top10_specific_genes.xls",quote=F,row.names=F,sep="\t")
merged.all.batchcorrected.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
write.table(test,"cluster_top50_specific_genes.xls",quote=F,row.names=F,sep="\t")

merged.all.batchcorrected.markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
write.table(test,"cluster_top200_specific_genes.xls",quote=F,row.names=F,sep="\t")

merged.all.batchcorrected.markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
write.table(test,"cluster_top500_specific_genes.xls",quote=F,row.names=F,sep="\t")

merged.all.batchcorrected.markers %>% group_by(cluster) %>% top_n(n = 1000, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
write.table(test,"cluster_top1000_specific_genes.xls",quote=F,row.names=F,sep="\t")
#######

##get DEGs from each cluster##
Idents(object = merged.all.batchcorrected)<- merged.all.batchcorrected$integrated_snn_res.0.5
cluster_list <- names(table(merged.all.batchcorrected$integrated_snn_res.0.5)[table(merged.all.batchcorrected$integrated_snn_res.0.5)>1])
merged.all.batchcorrected$celltype.sample.indent <- paste(Idents(object = merged.all.batchcorrected), merged.all.batchcorrected$batch, sep = "_")
Idents(object = merged.all.batchcorrected) <- "celltype.sample.indent"

##generate DEGs in each clusters##
lab <- rep(NA_character_, ncol(merged.all.batchcorrected))

lab[grepl("^graft", merged.all.batchcorrected$batch, ignore.case = TRUE)] <- "Graft"
lab[grepl("^pbmc",  merged.all.batchcorrected$batch, ignore.case = TRUE)] <- "PBMC"

## 快速检查
table(merged.all.batchcorrected$label)
merged.all.batchcorrected$label <- factor(lab, levels = c("Graft","PBMC"))

Idents(object = merged.all.batchcorrected) <- "label"
pdf("merged_batchcorrected_umap_split_graft_pbmc.pdf",16,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = TRUE,group.by="label",split.by = "label"))
dev.off()

pdf("merged_batchcorrected_umap_graft_pbmc.pdf",9,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = TRUE,group.by="label"))
dev.off()

pdf("merged_batchcorrected_umap_split_graft_pbmc_no_label.pdf",16,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = FALSE,group.by="label",split.by = "label"))
dev.off()

pdf("merged_batchcorrected_umap_graft_pbmc_no_label.pdf",9,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = FALSE,group.by="label"))
dev.off()

merged.all.batchcorrected$celltype.label.indent <- paste(merged.all.batchcorrected$integrated_snn_res.0.5, merged.all.batchcorrected$label, sep = "_")
Idents(object = merged.all.batchcorrected)<- merged.all.batchcorrected$celltype.label.indent

pdf("merged_batchcorrected_umap_celltype.graft_pbmc.pdf",16,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = TRUE,group.by="integrated_snn_res.0.5",split.by = "label"))
dev.off()

pdf("merged_batchcorrected_umap_celltype.graft_pbmc_no_label.pdf",16,8)
print(DimPlot(merged.all.batchcorrected, reduction = "umap", label = FALSE,group.by="integrated_snn_res.0.5",split.by = "label"))
dev.off()

###DEGs in each cluster between Graft and PBMC
for (i in names(table(merged.all.batchcorrected$integrated_snn_res.0.5)) ){
print(i)
if (( length(try(WhichCells(object = merged.all.batchcorrected, idents = paste0(i,"_Graft")),TRUE))>=3) & (length( try(WhichCells(object = merged.all.batchcorrected, idents = paste0(i,"_PBMC")),TRUE))>=3) ){
cluster.markers <- FindMarkers(object = merged.all.batchcorrected, ident.1 = paste0(i,"_Graft"), ident.2 = paste0(i,"_PBMC") , only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster.markers, paste0("DEGs_in_each_cluster/cluster.",i,".DEGs_Graft_vs_PMBC.xls"),sep="\t",row.names=T,col.names=T,quote=F)
}
}

##generate top expressed 3000 gene in each case and control##
Idents(object = merged.all.batchcorrected)<- merged.all.batchcorrected$integrated_snn_res.0.5
for (i in cluster_list ){
cluster <- WhichCells(object = merged.all.batchcorrected, idents = i)
expression_cluster <- GetAssayData(object = merged.all.batchcorrected, slot = "scale.data")[,cluster]
means_cluster <-  rowMeans(expression_cluster)[order(rowMeans(expression_cluster),decreasing =T)]
top3000 <- head(means_cluster,3000)
output <- data.frame(top3000)
rownames(output) <- names(top3000)
colnames(output) <- "average_scaled_expression"
write.table(output, paste0("top3000_scaled_expressed_genes/cluster.",i,".top_3000_expressed_genes.xls"),sep="\t",row.names=T,col.names=T,quote=F)
}

##cell_prop
cell_prop <- prop.table(x = table(Idents(object = merged.all.batchcorrected),merged.all.batchcorrected$batch),margin =2)
cell_prop <- round(cell_prop,3)
write.table(cell_prop,"cell_proportion.xls",quote=F,sep="\t")

#####heatmap#####
##plot top 3 DEG in each clusters###
merged.all.batchcorrected.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
DefaultAssay(merged.all.batchcorrected) <- "integrated"
heatmap_matrix <- GetAssayData(object = merged.all.batchcorrected, slot = "counts")
all_markers_need <- unique(test[,7])
pdf("Heatmap_All_scale_data_v2.pdf",20,12)
DoHeatmap(subset(merged.all.batchcorrected, downsample = 200),features = all_markers_need, size = 3,label=F)
dev.off()

DefaultAssay(merged.all.batchcorrected) <- "RNA"
pdf("feature_plot_gene_list0.pdf",16,length(gene_list_0)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list1.pdf",16,length(gene_list_1)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list2.pdf",16,length(gene_list_2)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list3.pdf",16,length(gene_list_3)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list4.pdf",16,length(gene_list_4)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list5.pdf",16,length(gene_list_5)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list6.pdf",16,length(gene_list_6)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list7.pdf",16,length(gene_list_7)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list8.pdf",16,length(gene_list_8)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list9.pdf",16,length(gene_list_9)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list10.pdf",16,length(gene_list_10)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list11.pdf",16,length(gene_list_11)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list12.pdf",16,length(gene_list_12)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list13.pdf",16,length(gene_list_13)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list14.pdf",16,length(gene_list_14)*4)
FeaturePlot(object = merged.all.batchcorrected, features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()


DefaultAssay(merged.all.batchcorrected) <- "RNA"
###split by batch
pdf("feature_plot_gene_list0_split.pdf",64,length(gene_list_0)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list1_split.pdf",64,length(gene_list_1)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list2_split.pdf",64,length(gene_list_2)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list3_split.pdf",64,length(gene_list_3)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list4_split.pdf",64,length(gene_list_4)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list5_split.pdf",64,length(gene_list_5)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list6_split.pdf",64,length(gene_list_6)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list7_split.pdf",64,length(gene_list_7)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list8_split.pdf",64,length(gene_list_8)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list9_split.pdf",64,length(gene_list_9)*8)
FeaturePlot(object = merged.all.batchcorrected, split.by="batch",features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list10_split.pdf",64,length(gene_list_10)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list11_split.pdf",64,length(gene_list_11)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list12_split.pdf",64,length(gene_list_12)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list13_split.pdf",64,length(gene_list_13)*8)
FeaturePlot(object = merged.all.batchcorrected, split.by="batch",features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("feature_plot_gene_list14_split.pdf",64,length(gene_list_14)*8)
FeaturePlot(object = merged.all.batchcorrected,split.by="batch", features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()

DefaultAssay(merged.all.batchcorrected) <- "RNA"
###split by label
pdf("feature_plot_gene_list0_split_by_tissue.pdf",8,length(gene_list_0)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list1_split_by_tissue.pdf",8,length(gene_list_1)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list2_split_by_tissue.pdf",8,length(gene_list_2)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list3_split_by_tissue.pdf",8,length(gene_list_3)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list4_split_by_tissue.pdf",8,length(gene_list_4)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list5_split_by_tissue.pdf",8,length(gene_list_5)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list6_split_by_tissue.pdf",8,length(gene_list_6)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list7_split_by_tissue.pdf",8,length(gene_list_7)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list8_split_by_tissue.pdf",8,length(gene_list_8)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list9_split_by_tissue.pdf",8,length(gene_list_9)*4)
FeaturePlot(object = merged.all.batchcorrected, split.by="label",features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list10_split_by_tissue.pdf",8,length(gene_list_10)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list11_split_by_tissue.pdf",8,length(gene_list_11)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list12_split_by_tissue.pdf",8,length(gene_list_12)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list13_split_by_tissue.pdf",8,length(gene_list_13)*4)
FeaturePlot(object = merged.all.batchcorrected, split.by="label",features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("feature_plot_gene_list14_split_by_tissue.pdf",8,length(gene_list_14)*4)
FeaturePlot(object = merged.all.batchcorrected,split.by="label", features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
# DefaultAssay(merged.all.batchcorrected) <- "RNA"


##update 9/22/2025
cluster2celltype <- c(
  "0"  = "Classical monocyte",
  "1"  = "Naive/Memory-like CD4+",
  "2"  = "KLRG1+ effector CD8+",
  "3"  = "Naïve B cell",
  "4"  = "NK cell",
  "5"  = "CXCL13+ effector CD8+",
  "6"  = "Resting monocyte",
  "7"  = "Naive/Memory-like CD4+",
  "8"  = "OASL+ effector CD8+",
  "9"  = "FOXP3+ Treg",
  "10" = "γδ T cell",
  "11" = "Effector-like CD4+",
  "12" = "GC B cell",
  "13" = "XCL1+ NK cell",
  "14" = "HLA-DR+ myeloid cell",
  "15" = "Plasma cell",
  "16" = "C1Q+ macrophage",
  "17" = "Neutrophil",
  "18" = "Cycling cell",
  "19" = "Naive B cell",
  "20" = "Memory B cell",
  "21" = "Plasmablast",
  "22" = "Transitional B cell",
  "23" = "cDC",
  "24" = "pDC",
  "25" = "Erythroid cell",
  "26" = "Mast cell/Basophil",
  "27" = "Plasmablast"
)

# 把 Idents 设置成你的 cluster 列（比如 integrated_snn_res.0.5）
Idents(merged.all.batchcorrected) <- merged.all.batchcorrected$integrated_snn_res.0.5

# 用 RenameIdents 改名
merged.all.batchcorrected <- RenameIdents(
  merged.all.batchcorrected,
  cluster2celltype
)

# 如果需要，把新的 Idents 存一份到 meta.data
merged.all.batchcorrected$cell_type_annotated <- Idents(merged.all.batchcorrected)

# 检查结果
table(Idents(merged.all.batchcorrected))
genes_use <- c(
  "LYZ","CD14","TCF7","LEF1","CD3D","CD8A","TCL1A","MS4A1","GNLY","NKG7",
  "LAG3","GZMB","FOXP3","CTLA4","TRDV1","TRDC","TSHZ2","RALGPS2","XCL1","XCL2",
  "CDKN1C","CX3CR1","JCHAIN","MZB1","C1QA","C1QB","CSF3R","NAMPT","TYMS","STMN1",
  "AFF3","IGHM","IGHD","CD1C","FCER1A","IL3RA","CLEC4C","HBB","ALAS2","CPA3","HDC"
)

# 确保 Idents 是 cell_type_annotated
Idents(merged.all.batchcorrected) <- merged.all.batchcorrected$cell_type_annotated

# 画 bubble plot (交换坐标轴)
p <- DotPlot(
  merged.all.batchcorrected,
  features = genes_use,
  assay   = "RNA",
  scale   = FALSE,            
  cols = c("lightgrey", "red"),
  dot.scale = 6
) + 
  coord_flip() +
  RotatedAxis() +
  theme(axis.text.y = element_text(size = 8))

# 保存 PDF
pdf("celltype_marker_bubbleplot_swapped.pdf", width = 12, height = 16)
print(p)
dev.off()

# 设置 Idents
Idents(merged.all.batchcorrected) <- merged.all.batchcorrected$cell_type_annotated

library(ggplot2)

p <- DotPlot(
  merged.all.batchcorrected,
  features = genes_use,
  assay    = "RNA",
  scale    = TRUE,                     
  cols     = c("lightgrey", "red"),
  dot.scale = 6
)

# 把 z-score 映射到 0-1 区间，避免负值
p$data$avg.exp.scaled <- scales::rescale(p$data$avg.exp.scaled, to = c(0, 1))

# 修改 legend 名称：颜色映射 & 点大小映射
p <- p +
  coord_flip() +
  RotatedAxis() +
  theme(axis.text.y = element_text(size = 8)) +
  scale_color_gradient(
    name = "Scaled Average Expression",   # 改颜色 legend 标题
    low = "lightgrey", high = "red"
  ) +
  guides(size = guide_legend(title = "Percent Expressed"))  # 改点大小 legend 标题

# 保存 PDF
pdf("celltype_marker_bubbleplot_rescaled.pdf", width = 12, height = 14)
print(p)
dev.off()

# 确保 Idents 是新的人群注释
Idents(merged.all.batchcorrected) <- merged.all.batchcorrected$cell_type_annotated

# DimPlot（标签自动排斥，避免重叠；shuffle+raster 改善密点渲染）
p_dim <- DimPlot(
  merged.all.batchcorrected,
  reduction = "umap",
  group.by  = "cell_type_annotated",
  label     = TRUE,
  repel     = TRUE,
  label.size = 3,
  shuffle    = TRUE,
  raster     = TRUE
) + ggplot2::theme(legend.position = "right")

# 保存 PDF
ggplot2::ggsave(
  filename = "DimPlot_cell_type_annotated.pdf",
  plot = p_dim,
  width = 12, height = 8, units = "in"
)
# 确保 Idents 是新的人群注释
Idents(merged.all.batchcorrected) <- merged.all.batchcorrected$cell_type_annotated

# DimPlot（标签自动排斥，避免重叠；shuffle+raster 改善密点渲染）
p_dim <- DimPlot(
  merged.all.batchcorrected,
  reduction = "umap",
  group.by  = "cell_type_annotated",
  label     = FALSE,
  repel     = TRUE,
  label.size = 3,
  shuffle    = TRUE,
  raster     = TRUE
) + ggplot2::theme(legend.position = "right")

# 保存 PDF
ggplot2::ggsave(
  filename = "DimPlot_cell_type_annotated_no_label.pdf",
  plot = p_dim,
  width = 12, height = 8, units = "in"
)

save(merged_combined_regressed, merged.all.batchcorrected, file = "merged_4_with_vdj_8_24_2025.RData")

##stratify to CD4 and CD8 together to conduct same analyses and clone/subtype analyses 10/8/2025
Idents(merged.all.batchcorrected) <- "cell_type_annotated"
CD4_CD8_object <- subset(merged.all.batchcorrected,idents =c("Naive/Memory-like CD4+","Effector-like CD4+", "FOXP3+ Treg","KLRG1+ effector CD8+","CXCL13+ effector CD8+","OASL+ effector CD8+"))

##redo batch effect correction 
#####merge without batch correction
cells_all <- colnames(CD4_CD8_object)
sample_id  <- sub(".*_(\\d+)$", "\\1", cells_all)      # 提取末尾样本号 1..8
cells_base <- sub("_(\\d+)$", "", cells_all)           # 去掉 _数字 后缀，得到原始对象里的cell名

map_df <- data.frame(
  cell_merged = cells_all,
  sample_id   = sample_id,
  cell_orig   = cells_base,
  stringsAsFactors = FALSE
)

# 2) 原始对象列表（顺序必须与样本号一致）
object.list <- list(PBMC1obj, PBMC2obj, PBMC3obj, PBMC4obj,
                    graft1obj, graft2obj, graft3obj, graft4obj)
target_names <- c("PBMC1obj_T","PBMC2obj_T","PBMC3obj_T","PBMC4obj_T",
                  "graft1obj_T","graft2obj_T","graft3obj_T","graft4obj_T")

# 3) 对每个样本：筛CD4/CD8的cells -> 去后缀 -> 在对应原始对象中subset
for (i in seq_along(object.list)) {
  keep_in_sample <- map_df$cell_orig[map_df$sample_id == as.character(i)]
  # 与原始对象取交集，避免潜在不匹配
  keep_in_sample <- intersect(keep_in_sample, colnames(object.list[[i]]))
  assign(target_names[i], subset(object.list[[i]], cells = keep_in_sample))
}

# 4) 简单核对一下每个样本匹配了多少个细胞（可选）
matched_counts <- sapply(seq_along(object.list), function(i) {
  sum(colnames(object.list[[i]]) %in% map_df$cell_orig[map_df$sample_id == as.character(i)])
})
names(matched_counts) <- target_names
matched_counts

object_T.list <- list(PBMC1obj_T, PBMC2obj_T, PBMC3obj_T, PBMC4obj_T,
                    graft1obj_T, graft2obj_T, graft3obj_T, graft4obj_T)
merged_T.all <- merge(x = object_T.list[[1]],
                    y = object_T.list[2:length(object.list)], 
                    add.cell.ids = c("PBMC_0825","PBMC_0712","PBMC_0701","PBMC_0708","graft_0825","graft_0712","graft_0701","graft_0708"),project = "merged_all_eight")



###merge and batch correction using rpca
features <- SelectIntegrationFeatures(object.list = object_T.list)
anchors <- FindIntegrationAnchors(object.list = object_T.list, 
                                  anchor.features = features,
                                  reduction = "rpca") 
merged_all_T_batchcorrected <- IntegrateData(anchorset = anchors)

#merged_cca <- NormalizeData(merged_cca)
DefaultAssay(merged_all_T_batchcorrected) <- "integrated"
merged_all_T_batchcorrected<- ScaleData(merged_all_T_batchcorrected, verbose = FALSE)
merged_all_T_batchcorrected<- RunPCA(merged_all_T_batchcorrected, npcs = 50, verbose = FALSE)
merged_all_T_batchcorrected <- RunUMAP(merged_all_T_batchcorrected, dims = 1:30, reduction = "pca")
sample_names <- c("PBMC_0825","PBMC_0712","PBMC_0701","PBMC_0708","graft_0825","graft_0712","graft_0701","graft_0708")


## 从 cell 名里取出最后一个下划线后的数字（例如 "...-1_8" 抽取到 8）
cell_ids  <- colnames(merged_all_T_batchcorrected) 
suffix_id <- sub(".*_(\\d+)$", "\\1", cell_ids)
suffix_id_num <- suppressWarnings(as.integer(suffix_id))

## 建立映射并写入 meta.data
id2name <- setNames(sample_names, 1:8)
merged_all_T_batchcorrected$batch <- unname(id2name[as.character(suffix_id_num)])

## 快速查看分布
table(merged_all_T_batchcorrected$batch)

DimPlot(merged_all_T_batchcorrected, group.by = "batch", label = TRUE)
merged_all_T_batchcorrected <- FindNeighbors(object = merged_all_T_batchcorrected, dims = 1:20, verbose = FALSE)
merged_all_T_batchcorrected <- FindClusters(object = merged_all_T_batchcorrected, verbose = FALSE, resolution = 0.5 )

merged_all_T_batchcorrected$celltype.sample.indent <- paste(merged_all_T_batchcorrected$integrated_snn_res.0.5, merged.all.batchcorrected$batch, sep = "_")
#DimPlot(merged_all_T_batchcorrected, group.by = "integrated_snn_res.0.5", label = TRUE)
merged_all_T_batchcorrected$cell_type_annotated_previous <- CD4_CD8_object$cell_type_annotated[match(rownames(CD4_CD8_object@meta.data),rownames(merged_all_T_batchcorrected@meta.data))]

pdf("merged_batchcorrected_T_cells.pdf",10,8)
print(DimPlot(merged_all_T_batchcorrected, reduction = "umap", label = TRUE,group.by="batch"))
dev.off()

pdf("merged_batchcorrected_T_cells_no_label.pdf",8,8)
print(DimPlot(merged_all_T_batchcorrected, reduction = "umap", label = FALSE,group.by="batch"))
dev.off()

Idents(merged_all_T_batchcorrected) <- merged_all_T_batchcorrected$integrated_snn_res.0.5
pdf("merged_batchcorrected_T_cells_by_cluster.pdf",12,8)
print(DimPlot(merged_all_T_batchcorrected,reduction = "umap", label = TRUE))
dev.off()
pdf("merged_batchcorrected_T_cells_by_cluster_no_label.pdf",10,8)
print(DimPlot(merged_all_T_batchcorrected, reduction = "umap", label = FALSE))
dev.off()
pdf("merged_batchcorrected_T_cells_by_cluster_split.pdf",16,8)
DimPlot(object = merged_all_T_batchcorrected, split.by = "batch",label = TRUE,reduction = "umap", pt.size = 0.3,ncol=4) + ggtitle(label = "UMAP")
dev.off()
pdf("merged_batchcorrected_T_cells_by_cluster_split_no_label.pdf",16,8)
DimPlot(object = merged_all_T_batchcorrected, split.by = "batch",label = FALSE,reduction = "umap", pt.size = 0.3,ncol=4) + ggtitle(label = "UMAP")
dev.off()

Idents(merged_all_T_batchcorrected) <- merged_all_T_batchcorrected$cell_type_annotated_previous
pdf("merged_batchcorrected_T_cells_by_previous_annotation_no_label.pdf",10,8)
print(DimPlot(merged_all_T_batchcorrected, reduction = "umap", label = FALSE))
dev.off()
pdf("merged_batchcorrected_T_cells_by_cluster_by_previous_annotation_split_no_label.pdf",16,8)
DimPlot(object = merged_all_T_batchcorrected, split.by = "batch",label = FALSE,reduction = "umap", pt.size = 0.3,ncol=4) + ggtitle(label = "UMAP")
dev.off()

Idents(merged_all_T_batchcorrected) <- merged_all_T_batchcorrected$integrated_snn_res.0.5

###DEGs in each cluster between Graft and PBMC
for (i in names(table(merged_all_T_batchcorrected$integrated_snn_res.0.5)) ){
print(i)
if (( length(try(WhichCells(object = merged_all_T_batchcorrected, idents = paste0(i,"_Graft")),TRUE))>=3) & (length( try(WhichCells(object = merged_all_T_batchcorrected, idents = paste0(i,"_PBMC")),TRUE))>=3) ){
cluster.markers <- FindMarkers(object = merged_all_T_batchcorrected, ident.1 = paste0(i,"_Graft"), ident.2 = paste0(i,"_PBMC") , only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(cluster.markers, paste0("T_cells_DEGs_in_each_cluster/cluster.",i,".DEGs_Graft_vs_PMBC.xls"),sep="\t",row.names=T,col.names=T,quote=F)
}
}

##generate top expressed 3000 gene in each case and control##
Idents(object = merged_all_T_batchcorrected)<- merged_all_T_batchcorrected$integrated_snn_res.0.5
cluster_list <- names(table(merged_all_T_batchcorrected$integrated_snn_res.0.5)[table(merged_all_T_batchcorrected$integrated_snn_res.0.5)>1])

for (i in cluster_list ){
cluster <- WhichCells(object = merged_all_T_batchcorrected, idents = i)
expression_cluster <- GetAssayData(object = merged_all_T_batchcorrected, slot = "scale.data")[,cluster]
means_cluster <-  rowMeans(expression_cluster)[order(rowMeans(expression_cluster),decreasing =T)]
top3000 <- head(means_cluster,3000)
output <- data.frame(top3000)
rownames(output) <- names(top3000)
colnames(output) <- "average_scaled_expression"
write.table(output, paste0("T_cells_top3000_scaled_expressed_genes/cluster.",i,".T_cells_top_3000_expressed_genes.xls"),sep="\t",row.names=T,col.names=T,quote=F)
}

##cell_prop
cell_prop <- prop.table(x = table(Idents(object = merged_all_T_batchcorrected),merged_all_T_batchcorrected$batch),margin =2)
cell_prop <- round(cell_prop,3)
write.table(cell_prop,"T_cells_cell_proportion.xls",quote=F,sep="\t")
write.table(table(merged_all_T_batchcorrected$batch),"T_cells_cell_number.xls",quote=F,sep="\t")


#####heatmap#####
##plot top 3 DEG in each clusters###
merged_all_T_batchcorrected.markers <- FindAllMarkers(object = merged_all_T_batchcorrected, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
merged_all_T_batchcorrected.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
DefaultAssay(merged_all_T_batchcorrected) <- "integrated"
heatmap_matrix <- GetAssayData(object = merged_all_T_batchcorrected, slot = "counts")
all_markers_need <- unique(test[,7])
pdf("T_cells_Heatmap_All_scale_data_v2.pdf",20,12)
DoHeatmap(subset(merged_all_T_batchcorrected, downsample = 200),features = all_markers_need, size = 3,label=F)
dev.off()

DefaultAssay(merged_all_T_batchcorrected) <- "RNA"
pdf("T_cells_feature_plot_gene_list0.pdf",16,length(gene_list_0)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list1.pdf",16,length(gene_list_1)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list2.pdf",16,length(gene_list_2)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list3.pdf",16,length(gene_list_3)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list4.pdf",16,length(gene_list_4)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list5.pdf",16,length(gene_list_5)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list6.pdf",16,length(gene_list_6)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list7.pdf",16,length(gene_list_7)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list8.pdf",16,length(gene_list_8)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list9.pdf",16,length(gene_list_9)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list10.pdf",16,length(gene_list_10)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list11.pdf",16,length(gene_list_11)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list12.pdf",16,length(gene_list_12)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list13.pdf",16,length(gene_list_13)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list14.pdf",16,length(gene_list_14)*4)
FeaturePlot(object = merged_all_T_batchcorrected, features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()


DefaultAssay(merged_all_T_batchcorrected) <- "RNA"
###split by batch
pdf("T_cells_feature_plot_gene_list0_split.pdf",64,length(gene_list_0)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list1_split.pdf",64,length(gene_list_1)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list2_split.pdf",64,length(gene_list_2)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list3_split.pdf",64,length(gene_list_3)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list4_split.pdf",64,length(gene_list_4)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list5_split.pdf",64,length(gene_list_5)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list6_split.pdf",64,length(gene_list_6)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list7_split.pdf",64,length(gene_list_7)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list8_split.pdf",64,length(gene_list_8)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list9_split.pdf",64,length(gene_list_9)*8)
FeaturePlot(object = merged_all_T_batchcorrected, split.by="batch",features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list10_split.pdf",64,length(gene_list_10)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list11_split.pdf",64,length(gene_list_11)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list12_split.pdf",64,length(gene_list_12)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list13_split.pdf",64,length(gene_list_13)*8)
FeaturePlot(object = merged_all_T_batchcorrected, split.by="batch",features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list14_split.pdf",64,length(gene_list_14)*8)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="batch", features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()

DefaultAssay(merged_all_T_batchcorrected) <- "RNA"
merged_all_T_batchcorrected$label <- CD4_CD8_object$label
###split by label
pdf("T_cells_feature_plot_gene_list0_split_by_tissue.pdf",8,length(gene_list_0)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list1_split_by_tissue.pdf",8,length(gene_list_1)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list2_split_by_tissue.pdf",8,length(gene_list_2)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list3_split_by_tissue.pdf",8,length(gene_list_3)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list4_split_by_tissue.pdf",8,length(gene_list_4)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list5_split_by_tissue.pdf",8,length(gene_list_5)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list6_split_by_tissue.pdf",8,length(gene_list_6)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list7_split_by_tissue.pdf",8,length(gene_list_7)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list8_split_by_tissue.pdf",8,length(gene_list_8)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list9_split_by_tissue.pdf",8,length(gene_list_9)*4)
FeaturePlot(object = merged_all_T_batchcorrected, split.by="label",features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list10_split_by_tissue.pdf",8,length(gene_list_10)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list11_split_by_tissue.pdf",8,length(gene_list_11)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list12_split_by_tissue.pdf",8,length(gene_list_12)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list13_split_by_tissue.pdf",8,length(gene_list_13)*4)
FeaturePlot(object = merged_all_T_batchcorrected, split.by="label",features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list14_split_by_tissue.pdf",8,length(gene_list_14)*4)
FeaturePlot(object = merged_all_T_batchcorrected,split.by="label", features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()

save(merged_all_T_batchcorrected, file = "merged_4_with_vdj_T_cells_10_8_2025.RData")

###plot top 3 CDR3 sequencee in graft
merged_all_T_batchcorrected$clone_by_label <- paste0(merged_all_T_batchcorrected$t_cdr3s_aa,"_",merged_all_T_batchcorrected$label)
tmp <- table(merged_all_T_batchcorrected$t_cdr3s_aa, merged_all_T_batchcorrected$label)
t_idx <- {}
##get shared clone in either Graft and PBMC
for (t in 1:nrow(tmp)){
##with at least three clone
	if (tmp[t,1]>3 && tmp[t,2] >3){
		t_idx <- c(t_idx,t)
	}
}
tmp = tmp[t_idx,]
idx <- order(tmp[,1],decreasing=T)
tmp <- tmp[idx,]
system("mkdir -p shared_clone")

for (i in rownames(tmp)){
print(i)
Idents(merged_all_T_batchcorrected) <- "clone_by_label"
if (( length(try(WhichCells(object = merged_all_T_batchcorrected, idents = paste0(i,"_Graft")),TRUE))>=1) & (length( try(WhichCells(object = merged_all_T_batchcorrected, idents = paste0(i,"_PBMC")),TRUE))>=1) ){
cluster.markers <- FindMarkers(object = merged_all_T_batchcorrected,min.cells.group=1, ident.1 = paste0(i,"_Graft"), ident.2 = paste0(i,"_PBMC") , only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
colnames(cluster.markers)[2] <- "avg_logFC"
pdf(paste0("shared_clone/",i,".DEGs_volcano_plot.pdf"))
print(GEX_volcano(cluster.markers,input.type= "findmarkers",condition.1="Graft", condition.2="PBMC",n.label.up=10,n.label.down=10))
dev.off()

write.table(cluster.markers, paste0("shared_clone/cluster.",i,".DEGs_Graft_vs_PMBC.xls"),sep="\t",row.names=T,col.names=T,quote=F)
g1_Graft <- WhichCells(merged_all_T_batchcorrected, idents = paste0(i,"_Graft"))
g1_PBMC <- WhichCells(merged_all_T_batchcorrected, idents = paste0(i,"_PBMC"))
pdf(paste0("shared_clone/",i,".DEGs_graft_vs_PBMC.pdf"),6,5)
print(DimPlot(merged_all_T_batchcorrected, label=T, group.by="integrated_snn_res.0.5", cells.highlight= list(g1_Graft, g1_PBMC), cols.highlight = c("darkblue", "red"), cols= "grey"))
dev.off()
}
}

###update 10/27/2025 remove and rename cell clusters
# Cluster annotation record from the 10/27/2025 update:
# 0: Memory-like CD4+
# 1: CXCL13+CXCR6+ effector CD8+
# 2: CXCR6+ effector CD8+
# 3: CCR2+ CD4+
# 4: Naive/Memory-like CD8+
# 5: Naive CD4+
# 6: Treg
# 7: CX3CR1+ effector CD8+
# 9: ZNF683+TCF7hi CD8+
# 10: CXCL13+ CD4+
# 13: ISG-high CD8+
# Clusters 8, 11, and 12 are excluded from downstream analysis.
# 1. 把分群列设成当前身份
Idents(merged_all_T_batchcorrected) <- "integrated_snn_res.0.5"

# 2. 定义要丢掉的 cluster
clusters_to_drop <- c("8", "11", "12")

# 3. 计算要保留哪些 cluster
clusters_keep <- setdiff(levels(Idents(merged_all_T_batchcorrected)), clusters_to_drop)

# 4. 找到这些 cluster 里所有的细胞
cells_keep <- WhichCells(merged_all_T_batchcorrected, idents = clusters_keep)

# 5. 子集
merged_all_T_batchcorrected_filter <- subset(
  merged_all_T_batchcorrected,
  cells = cells_keep
)

# 6. 现在对 filtered 对象重命名 cluster
merged_all_T_batchcorrected_filter <- RenameIdents(
  merged_all_T_batchcorrected_filter,
  `0`  = "Memory-like CD4+",
  `1`  = "CXCL13+CXCR6+ effector CD8+",
  `2`  = "CXCR6+ effector CD8+",
  `3`  = "CCR2+ CD4+",
  `4`  = "Naive/Memory-like CD8+",
  `5`  = "Naive CD4+",
  `6`  = "Treg",
  `7`  = "CX3CR1+ effector CD8+",
  `9`  = "ZNF683+TCF7hi CD8+",
  `10` = "CXCL13+ CD4+",
  `13` = "ISG-high CD8+"
)
merged_all_T_batchcorrected_filter$cell_type_annotated_granularity <- Idents(merged_all_T_batchcorrected_filter)

save(merged_all_T_batchcorrected_filter, file = "merged_4_with_vdj_T_cells_filter_10_27_2025.RData")

pdf("merged_batchcorrected_T_cells_by_final_annotation_no_label.pdf",10,8)
print(DimPlot(merged_all_T_batchcorrected_filter, reduction = "umap", label = FALSE))
dev.off()
pdf("merged_batchcorrected_T_cells_by_cluster_by_final_annotation_split_no_label.pdf",16,8)
DimPlot(object = merged_all_T_batchcorrected_filter, split.by = "batch",label = FALSE,reduction = "umap", pt.size = 0.3,ncol=4) + ggtitle(label = "UMAP")
dev.off()

##cell_prop
cell_prop <- prop.table(x = table(Idents(object = merged_all_T_batchcorrected_filter),merged_all_T_batchcorrected_filter$batch),margin =2)
cell_prop <- round(cell_prop,3)
write.table(cell_prop,"T_cells_cell_proportion_final.xls",quote=F,sep="\t")
write.table(table(merged_all_T_batchcorrected_filter$batch),"T_cells_cell_number_final.xls",quote=F,sep="\t")


#####heatmap#####
##plot top 3 DEG in each clusters###
merged_all_T_batchcorrected_filter.markers <- FindAllMarkers(object = merged_all_T_batchcorrected_filter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
merged_all_T_batchcorrected_filter.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
DefaultAssay(merged_all_T_batchcorrected_filter) <- "integrated"
heatmap_matrix <- GetAssayData(object = merged_all_T_batchcorrected_filter, slot = "counts")
all_markers_need <- unique(test[,7])
pdf("T_cells_Heatmap_All_scale_data_final.pdf",20,12)
DoHeatmap(subset(merged_all_T_batchcorrected_filter, downsample = 200),features = all_markers_need, size = 3,label=F)
dev.off()

DefaultAssay(merged_all_T_batchcorrected_filter) <- "RNA"
pdf("T_cells_feature_plot_gene_list0_final.pdf",16,length(gene_list_0)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list1_final.pdf",16,length(gene_list_1)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list2_final.pdf",16,length(gene_list_2)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list3_final.pdf",16,length(gene_list_3)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list4_final.pdf",16,length(gene_list_4)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list5_final.pdf",16,length(gene_list_5)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list6_final.pdf",16,length(gene_list_6)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list7_final.pdf",16,length(gene_list_7)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list8_final.pdf",16,length(gene_list_8)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list9_final.pdf",16,length(gene_list_9)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list10_final.pdf",16,length(gene_list_10)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list11_final.pdf",16,length(gene_list_11)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list12_final.pdf",16,length(gene_list_12)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list13_final.pdf",16,length(gene_list_13)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list14_final.pdf",16,length(gene_list_14)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()


DefaultAssay(merged_all_T_batchcorrected_filter) <- "RNA"
###split by batch
pdf("T_cells_feature_plot_gene_list0_split_final.pdf",64,length(gene_list_0)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list1_split_final.pdf",64,length(gene_list_1)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list2_split_final.pdf",64,length(gene_list_2)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list3_split_final.pdf",64,length(gene_list_3)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list4_split_final.pdf",64,length(gene_list_4)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list5_split_final.pdf",64,length(gene_list_5)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list6_split_final.pdf",64,length(gene_list_6)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list7_split_final.pdf",64,length(gene_list_7)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list8_split_final.pdf",64,length(gene_list_8)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list9_split_final.pdf",64,length(gene_list_9)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter, split.by="batch",features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list10_split_final.pdf",64,length(gene_list_10)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list11_split_final.pdf",64,length(gene_list_11)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list12_split_final.pdf",64,length(gene_list_12)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list13_split_final.pdf",64,length(gene_list_13)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter, split.by="batch",features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list14_split_final.pdf",64,length(gene_list_14)*8)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="batch", features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()

DefaultAssay(merged_all_T_batchcorrected_filter) <- "RNA"
###split by label
pdf("T_cells_feature_plot_gene_list0_split_by_tissue_final.pdf",8,length(gene_list_0)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list1_split_by_tissue_final.pdf",8,length(gene_list_1)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list2_split_by_tissue_final.pdf",8,length(gene_list_2)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list3_split_by_tissue_final.pdf",8,length(gene_list_3)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list4_split_by_tissue_final.pdf",8,length(gene_list_4)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list5_split_by_tissue_final.pdf",8,length(gene_list_5)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list6_split_by_tissue_final.pdf",8,length(gene_list_6)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list7_split_by_tissue_final.pdf",8,length(gene_list_7)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list8_split_by_tissue_final.pdf",8,length(gene_list_8)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list9_split_by_tissue_final.pdf",8,length(gene_list_9)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, split.by="label",features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list10_split_by_tissue_final.pdf",8,length(gene_list_10)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list11_split_by_tissue_final.pdf",8,length(gene_list_11)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list12_split_by_tissue_final.pdf",8,length(gene_list_12)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list13_split_by_tissue_final.pdf",8,length(gene_list_13)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter, split.by="label",features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("T_cells_feature_plot_gene_list14_split_by_tissue_final.pdf",8,length(gene_list_14)*4)
FeaturePlot(object = merged_all_T_batchcorrected_filter,split.by="label", features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()

##with cell type annotation
#all_clone$new_cluster_name <- all_clone$SCT_snn_res.0.6
#Idents(all_clone) <- "new_cluster_name"
#all_clone<- RenameIdents(all_clone, `0` = "CD8", `3` = "CD8", `5` = "CD8", `6` = "CD8", `7` = "CD8",`10` = "CD8", `1`= "CD4", `2`= "CD4",`4`= "CD4")
#all_clone$new_cluster_name <- Idents(all_clone)
#cluster_list <- names(table(all_clone$new_cluster_name)[table(all_clone$new_cluster_name)>1])
#all_clone$celltype.batch <- paste(Idents(object = all_clone), all_clone$batch, sep = "_")
#Idents(object = all_clone) <- "celltype.batch"

setwd(output_directory)

##T-cell clone in the cluster
#clone_in_cluster <- table(!is.na(merged_all_T_batchcorrected_filter$t_clonotype_id),merged_all_T_batchcorrected_filter$SCT_snn_res.0.6)
#write.table(clone_in_cluster,file="clone_in_clusters.txt",quote=F,sep="\t")

#clone_in_cluster_detail <- table(merged_all_T_batchcorrected_filter$t_clonotype_id,merged_all_T_batchcorrected_filter$SCT_snn_res.0.6)
#idx <- order( as.numeric(as.vector(gsub("clonotype","",rownames(clone_in_cluster_detail)))))
#clone_in_cluster_detail <- clone_in_cluster_detail[idx,]
#write.table(clone_in_cluster_detail,file="clone_in_clusters_detail.txt",quote=F,sep="\t")

##test whether t_clonotype_id are generalizable across samples， the answers are NO!!!
clono100_obj <- subset(
  merged_all_T_batchcorrected_filter,
  subset = t_clonotype_id == "clonotype100"
)
##https://ncborcherding.github.io/vignettes/vignette.html
# Historical single-sample input examples are omitted from the public copy;
# use the current sample manifest in raw_to_intermediate/ instead.
#c("PBMC_0825","PBMC_0712","PBMC_0701","PBMC_0708","graft_0825","graft_0712","graft_0701","graft_0708") in these order to merge TCRs
PBMC_0825_vdj <- read.csv(paste0(input_VDJ_directory_pbmc_2,"/filtered_contig_annotations.csv"))
PBMC_0712_vdj <- read.csv(paste0(input_VDJ_directory_pbmc_3,"/filtered_contig_annotations.csv"))
PBMC_0701_vdj <- read.csv(paste0(input_VDJ_directory_pbmc_1,"/filtered_contig_annotations.csv"))
PBMC_0708_vdj <- read.csv(paste0(input_VDJ_directory_pbmc,"/filtered_contig_annotations.csv"))
graft_0825_vdj <- read.csv(paste0(input_VDJ_directory_2,"/filtered_contig_annotations.csv"))
graft_0712_vdj <- read.csv(paste0(input_VDJ_directory_3,"/filtered_contig_annotations.csv"))
graft_0701_vdj <- read.csv(paste0(input_VDJ_directory_1,"/filtered_contig_annotations.csv"))
graft_0708_vdj <- read.csv(paste0(input_VDJ_directory,"/filtered_contig_annotations.csv"))

contig_list <- list(PBMC_0825_vdj, PBMC_0712_vdj,PBMC_0701_vdj,PBMC_0708_vdj,graft_0825_vdj,graft_0712_vdj,graft_0701_vdj,graft_0708_vdj)
##T-AB - T cells, alpha-beta TCR
combined <- combineTCR(contig_list, 
                samples = c("PBMC_0825","PBMC_0712","PBMC_0701","PBMC_0708","graft_0825","graft_0712","graft_0701","graft_0708"), 
                ID = c("1","2","3","4","5","6","7","8"))
for (i in seq_along(combined)) {
  id <- combined[[i]]$ID[1]  # 从每个样本里取对应 ID
  combined[[i]]$barcode <- gsub("^.*_", "", combined[[i]]$barcode)  # 去掉前缀直到最后一个 "_"
  combined[[i]]$barcode <- paste0(combined[[i]]$barcode, "_", id)   # 在后面加 "_ID"
}
##clonotype analysis
##Visualizing Contigs
##Clonotype Grouping  use "gene+nt" method, or "gene+aa"
##https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/algorithms/clonotyping
# 克隆扩增 vs 稀有克隆（分组到不同"size bin"）
pdf("unique_clone_proportion.pdf")
clonalQuant(
    input.data = combined,
    cloneCall  = "gene+nt",  # same meaning as before
    chain      = "both",     # default in v2; keep explicit for clarity
    scale      = TRUE        # TRUE -> relative %
)
dev.off()
combined <- addVariable(combined, 
                            variable.name = "patient", 
                            variables = rep(c("0825","0712","0701","0708"),2))

table(combined[[1]]$patient)


# 1b. 也拿到表格 (原 quantContig exportTable=TRUE, 
#     v2 的 clonalQuant() 仍然支持返回 data.frame)
quantContig_output <- clonalQuant(
    input.data   = combined,
    cloneCall    = "gene+nt",
    chain        = "both",
    scale        = TRUE
    # 在 v2, clonalQuant() 直接返回一个 data.frame 可再用。
    # 不再用 exportTable=TRUE 这种老参数；现在默认 return 是可抓取对象。:contentReference[oaicite:2]{index=2}
)

# 1c. clone uniqueness / expansion distribution
# abundanceContig -> clonalAbundance
pdf("clone_uniqueness.pdf")
clonalAbundance(
    input.data = combined,
    cloneCall  = "gene+nt",
    chain      = "both",
    scale      = TRUE      # same idea: show proportions not raw counts
)
dev.off()

# 1d. CDR3 length distribution
# lengthContig(s) -> clonalLength
pdf("length_distribution_CDR3_sequences_aa_both_chain.pdf")
clonalLength(
    input.data = combined,
    cloneCall  = "aa",     # amino acid sequence length
    chain      = "both"    # both TRA/TRB chains
)
dev.off()


### 2. Compare clonotypes between samples / alluvial flows ###

# compareClonotypes -> clonalCompare
# graph="alluvial" still supported, but they renamed "alluvialClonotypes" (per-sample flows) 
# to "alluvialClones". Here you were using compareClonotypes(..., graph="alluvial"),
# so use clonalCompare(..., graph="alluvial").
# samples = c("PBMC_0825","PBMC_0712","PBMC_0701","PBMC_0708","graft_0825","graft_0712","graft_0701","graft_0708"), 
# ID = c("1","2","3","4","5","6","7","8"))
                
pdf("clone_transfer_0825.pdf")
clonalCompare(
    input.data = combined,
    cloneCall  = "aa",
    samples    = c("PBMC_0825_1", "graft_0825_5"),
    top.clones    = 10,               # top N clones
    graph      = "alluvial"        # produce alluvial flow
)
dev.off()

pdf("clone_transfer_0712.pdf")
clonalCompare(
    input.data = combined,
    cloneCall  = "aa",
    samples    = c("PBMC_0712_2", "graft_0712_6"),
    top.clones    = 10,               # top N clones
    graph      = "alluvial"        # produce alluvial flow
)
dev.off()

pdf("clone_transfer_0701.pdf")
clonalCompare(
    input.data = combined,
    cloneCall  = "aa",
    samples    = c("PBMC_0701_3", "graft_0701_7"),
    top.clones    = 10,               # top N clones
    graph      = "alluvial"        # produce alluvial flow
)
dev.off()

pdf("clone_transfer_0708.pdf")
clonalCompare(
    input.data = combined,
    cloneCall  = "aa",
    samples    = c("PBMC_0708_4", "graft_0708_8"),
    top.clones    = 10,               # top N clones
    graph      = "alluvial"        # produce alluvial flow
)
dev.off()

### 3. Clonal Space Homeostasis / Proportion / Size Distribution ###

# clonalHomeostasis: name is unchanged in v2 (still clonalHomeostasis) :contentReference[oaicite:3]{index=3}
pdf("unique_clone_proportion_advanced.pdf")
clonalHomeostasis(
    input.data = combined,
    cloneCall  = "gene+nt",
    chain      = "both",
    cloneSize  = c(
      Rare          = 1e-04,
      Small         = 0.001,
      Medium        = 0.01,
      Large         = 0.1,
      Hyperexpanded = 1
    )
)
dev.off()

# clonalProportion: unchanged name in v2
pdf("unique_clone_number_advanced.pdf")
clonalProportion(
    input.data   = combined,
    cloneCall    = "gene+nt",                 # 你之前用的定义
    chain        = "both",
    clonalSplit  = c(10, 100, 1000, 10000, 30000, 1e5),
    group.by     = NULL,                      # 可选: 比如 "sample" 或 "batch"
    order.by     = NULL,                      # 可选: 手动排序
    exportTable  = FALSE                      # 如果想拿矩阵，把它设 TRUE
)
dev.off()

# clonesizeDistribution() got renamed in changelog:
# clonotypeSizeDistribution -> clonalSizeDistribution. 
# Your code called clonesizeDistribution(), which is basically the same idea.
# We'll use clonalSizeDistribution() with method="ward.D2".
pdf("clone_distribution.pdf",10,8)
clonalSizeDistribution(
    input.data = combined,
    cloneCall  = "gene+nt",
    chain      = "both",
    method     = "ward.D2"
)
dev.off()


### 4. Overlap Analysis (Morisita / Jaccard / Overlap / Raw) ###

# clonalOverlap(): name unchanged in v2 and still supports the same methods. :contentReference[oaicite:4]{index=4}
pdf("cloneoverlap_all_moristia_nt.pdf",10,10)
clonalOverlap(combined,
              cloneCall = "gene+nt",
              chain     = "both",
              method    = "morisita")
dev.off()

pdf("cloneoverlap_all_jaccard_nt.pdf",10,10)
clonalOverlap(combined,
              cloneCall = "gene+nt",
              chain     = "both",
              method    = "jaccard")
dev.off()

pdf("cloneoverlap_all_raw_nt.pdf",10,10)
clonalOverlap(combined,
              cloneCall = "gene+nt",
              chain     = "both",
              method    = "raw")
dev.off()

pdf("cloneoverlap_all_overlap_nt.pdf",10,10)
clonalOverlap(combined,
              cloneCall = "gene+nt",
              chain     = "both",
              method    = "overlap")
dev.off()

pdf("cloneoverlap_all_moristia_aa.pdf",10,10)
clonalOverlap(combined,
              cloneCall = "aa",
              chain     = "both",
              method    = "morisita")
dev.off()

pdf("cloneoverlap_all_jaccard_aa.pdf",10,10)
clonalOverlap(combined,
              cloneCall = "aa",
              chain     = "both",
              method    = "jaccard")
dev.off()

pdf("cloneoverlap_all_raw_aa.pdf",10,10)
clonalOverlap(combined,
              cloneCall = "aa",
              chain     = "both",
              method    = "raw")
dev.off()

pdf("cloneoverlap_all_overlap_aa.pdf",10,10)
clonalOverlap(combined,
              cloneCall = "aa",
              chain     = "both",
              method    = "overlap")
dev.off()

### 5. Scatter between samples ###

# scatterClonotype(s) -> clonalScatter
pdf("comparison_between_samples_0825.pdf")
clonalScatter(
    input.data = combined,
    cloneCall  = "gene+nt",
    x.axis     = "PBMC_0825_1",
    y.axis     = "graft_0825_5",
    dot.size   = "total",        # size aesthetic
    graph      = "proportion"    # same meaning as before
)
dev.off()

pdf("comparison_between_samples_0712.pdf")
clonalScatter(
    input.data = combined,
    cloneCall  = "gene+nt",
    x.axis     = "PBMC_0712_2",
    y.axis     = "graft_0712_6",
    dot.size   = "total",        # size aesthetic
    graph      = "proportion"    # same meaning as before
)
dev.off()

pdf("comparison_between_samples_0701.pdf")
clonalScatter(
    input.data = combined,
    cloneCall  = "gene+nt",
    x.axis     = "PBMC_0701_3",
    y.axis     = "graft_0701_7",
    dot.size   = "total",        # size aesthetic
    graph      = "proportion"    # same meaning as before
)
dev.off()

pdf("comparison_between_samples_0708.pdf")
clonalScatter(
    input.data = combined,
    cloneCall  = "gene+nt",
    x.axis     = "PBMC_0708_4",
    y.axis     = "graft_0708_8",
    dot.size   = "total",        # size aesthetic
    graph      = "proportion"    # same meaning as before
)
dev.off()


### 6. 把克隆信息投到 Seurat 对象 ###

# combineExpression() is still called combineExpression() in v2 and same idea:
# attach clone info to Seurat meta.data (cloneType, etc.). :contentReference[oaicite:5]{index=5}
# 重要：你仍然需要先把 Seurat 的细胞名改成 "Sample_ID-CellBarcode" 这种和 combineTCR 输出一致的前缀，否则不会对上。


merged_all_T_batchcorrected_filter <- combineExpression(
    input.data = combined,
    sc         = merged_all_T_batchcorrected_filter,
    cloneCall  = "gene+nt",
    chain      = "both",   # 你的对象是TCR α/β，所以还是 both
    proportion = FALSE,
    cloneSize  = c(
        Single        = 1,
        Small         = 5,
        Medium        = 20,
        Large         = 100,
        Hyperexpanded = 500
    )
)

colorblind_vector <- colorRampPalette(rev(c(
  "#0D0887FF", "#47039FFF", "#7301A8FF", "#9C179EFF", "#BD3786FF",
  "#D8576BFF", "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF"
)))

pdf("Clone_expansion_in_all_Cell.pdf",16,8)
DimPlot(
    merged_all_T_batchcorrected_filter,
    group.by = "cloneSize",   # <- 改这里
    split.by = "batch"
) +
    scale_color_manual(
        values   = colorblind_vector(5),
        na.value = "grey"
    ) +
    theme(plot.title = element_blank())
dev.off()


### 7. “shared clones” 分析对象 all_clone_test ###

all_clone_test <- all_clone
Idents(all_clone_test) <- all_clone_test$SCT_snn_res.0.8

NEW_NAME_clone <- Cells(all_clone_test)
NEW_NAME_clone <- gsub("^Graft", "Graft_07082022", NEW_NAME_clone)
NEW_NAME_clone <- gsub("^PBMC",  "PBMC_07082022",  NEW_NAME_clone)
all_clone_test <- RenameCells(
    object    = all_clone_test,
    new.names = NEW_NAME_clone
)

all_clone_test <- combineExpression(
    input.data  = combined,
    sc          = all_clone_test,
    cloneCall   = "gene+nt",
    chain       = "both",
    proportion  = FALSE,
    cloneTypes  = c(
        Single        = 1,
        Small         = 5,
        Medium        = 20,
        Large         = 100,
        Hyperexpanded = 500
    )
)

pdf("Clone_expansion_in_all_shared_clone.pdf",16,8)
DimPlot(
    all_clone_test,
    group.by = "cloneType",
    split.by = "batch"
) +
    scale_color_manual(values = colorblind_vector(5), na.value="grey") +
    theme(plot.title = element_blank())
dev.off()


### 8. Overlay, Network, Occupancy, Alluvial, Circos ###

# clonalOverlay(): name unchanged in changelog (still clonalOverlay)
pdf("clonalOverlay_all_shared_clone.pdf")
clonalOverlay(
    sc             = all_clone_test,
    reduction      = "UMAP",
    freq.cutpoint  = 10,
    bins           = 10,
    facet          = "batch"
) + guides(color = "none")
dev.off()

pdf("clonalOverlay_all_Cell.pdf",12,6)
clonalOverlay(
    sc.data     = merged_all_T_batchcorrected_filter,
    reduction   = "umap",
    cut.category= "clonalFrequency",    # 或者你的 meta.data 列名，例如 "cloneSize" 或 "cloneType"
    cutpoint    = 10,
    bins        = 10,
    facet.by    = "batch"
)
dev.off()

# clonalNetwork(): name unchanged
clonalNetwork(
    sc         = merged_all_T_batchcorrected_filter,
    reduction  = "umap",
    identity   = "SCT_snn_res.0.6",
    filter.clones   = NULL,
    filter.identity = NULL,
    cloneCall  = "aa",
    chain      = "both"
)

DimPlot(merged_all_T_batchcorrected_filter,label =TRUE)
clonalNetwork(
     sc = merged_all_T_batchcorrected_filter,
     reduction = "umap",
     chain = "both",
     cloneCall = "aa",
     group.by = "cell_type_annotated_granularity",
     filter.clones = NULL,
     filter.identity = "Naive CD4+"
)

clonalNetwork(
    sc         = merged_all_T_batchcorrected_filter,
    reduction  = "umap",
    identity   = "cell_type_annotated_granularity",
    filter.clones   = NULL,
    filter.identity = NULL,
    cloneCall  = "aa",
    chain      = "both"
)

clonalNetwork(
    sc         = all_clone_test,
    reduction  = "umap",
    identity   = "SCT_snn_res.0.8",
    filter.clones   = NULL,
    filter.identity = NULL,
    cloneCall  = "clone_by_batch",
    chain      = "both"
)

library(dplyr)
library(purrr)

# 假设 combined 是 combineTCR(...) 的结果
# 比如 combined[[1]] 是 PBMC, combined[[2]] 是 Graft, etc.

get_top20_for_one_sample <- function(df_one_sample) {
    # df_one_sample: one element of combined (one sample's contig table)

    # 我们先 check 它的 sample 名。很多版本里有列叫 "sample" 或 "Samples" 或 "orig.ident"
    # 如果没有，就用你在 combineTCR 里给的名字（即这个 df 的 attr）
    sample_name <- unique(df_one_sample$sample)
    if (length(sample_name) == 0 || all(is.na(sample_name))) {
        # fallback: try "Samples" or "ID"
        if ("Samples" %in% colnames(df_one_sample)) {
            sample_name <- unique(df_one_sample$Samples)
        } else if ("ID" %in% colnames(df_one_sample)) {
            sample_name <- unique(df_one_sample$ID)
        } else {
            # last resort: name not in column, we just mark sample_name = NA and fill later
            sample_name <- NA
        }
    }
    # 如果 sample_name 是向量长度>1，取第一个
    sample_name <- sample_name[1]

    # 这里定义“克隆是谁”
    # 我们用 CTstrict 当克隆ID主键（严格定义的克隆），
    # 因为这是最细的定义，一对一映射到 CTaa / CTnt / CTgene
    # 然后算这个克隆在这个样本里出现多少细胞(barcodes)
    top20_df <- df_one_sample %>%
        filter(!is.na(CTstrict)) %>%
        group_by(CTstrict, CTaa, CTnt, CTgene) %>% 
        summarise(
            clone_size = n(),   # 这个克隆有多少细胞
            .groups = "drop"
        ) %>%
        arrange(desc(clone_size)) %>%
        slice_head(n = 20) %>%
        mutate(batch = sample_name) %>%
        # 统一列顺序
        select(
            batch,
            CTstrict,
            clone_size,
            CTaa,
            CTnt,
            CTgene
        )

    return(top20_df)
}

# 对 combined 里每个样本跑一遍上面的函数，然后绑在一起
top20_all_samples <- map_df(combined, get_top20_for_one_sample)

write.csv(top20_all_samples,
          file = "top20_clones_combined.csv",
          row.names = FALSE)


# 1. UMAP with labels from Seurat
p_dim <- DimPlot(
    merged_all_T_batchcorrected_filter,
    reduction = "umap",
    group.by  = "cell_type_annotated_granularity",
    label     = TRUE,
    repel     = TRUE
) + theme(plot.title = element_blank())

# 2. clonalNetwork (你的参数原样)
p_net <- clonalNetwork(
    sc              = merged_all_T_batchcorrected_filter,
    reduction       = "umap",
    chain           = "both",
    cloneCall       = "aa",
    group.by        = "cell_type_annotated_granularity",
    filter.clones   = NULL,
    filter.identity = "Naive CD4+"
)

# 3. 并排拼图
combined_plot <- p_dim | p_net

pdf("DimPlot_plus_clonalNetwork_NAIVE_CD4.pdf", width = 14, height = 6)
print(combined_plot)
dev.off()



# 1. 高亮指定氨基酸序列的克隆




scRep_example <- highlightClones(
  merged_all_T_batchcorrected_filter,
  cloneCall = "aa",
  sequence  = c("CALSENRDDKIIF_CASSQGSSDTQYF")
)
pdf("highlight_0825_CALSENRDDKIIF_CASSQGSSDTQYF.pdf")
DimPlot(
  scRep_example,
  group.by = "highlight"
) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(
    plot.title = element_blank(),
    legend.position = "bottom"
  )
dev.off()

library(dplyr)


# 统计每个 cell type 中目标克隆的数量
clone_table <- merged_all_T_batchcorrected_filter@meta.data %>%
  filter(
    CTaa == "CALSENRDDKIIF_CASSQGSSDTQYF",
    patient == "0825"
  ) %>%
  count(cell_type_annotated_granularity, name = "cell_count") %>%
  arrange(desc(cell_count))

# 打印看一下结果
print(clone_table)

# 保存到当前工作目录
write.csv(clone_table, "clone_CALSENRDDKIIF_CASSQGSSDTQYF_0825_counts.csv", row.names = FALSE)
p_alluvial <- alluvialClones(
    merged_all_T_batchcorrected_filter,
    cloneCall = "cTaa",   # 如果 "aa" 报错，改成 "CTaa"（你真实的AA克隆列名）
    y.axes    = c("batch", "cell_type_annotated_granularity", "label"), 
    color     = c("CALSENRDDKIIF_CASSQGSSDTQYF")  # 要高亮的克隆
)
pdf("riverplot_CALSENRDDKIIF_CASSQGSSDTQYF.pdf",15,9)
print(p_alluvial)
dev.off()


###another example for 0712 CAEMYSGGGADGLTF_CASSLFSGNNEQFF ##
scRep_example2 <- highlightClones(
  merged_all_T_batchcorrected_filter,
  cloneCall = "aa",
  sequence  = c("CAEMYSGGGADGLTF_CASSLFSGNNEQFF")
)
pdf("highlight_0712_CAEMYSGGGADGLTF_CASSLFSGNNEQFF.pdf")
DimPlot(
  scRep_example2,
  group.by = "highlight"
) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(
    plot.title = element_blank(),
    legend.position = "bottom"
  )
dev.off()

library(dplyr)


# 统计每个 cell type 中目标克隆的数量
clone_table2 <- merged_all_T_batchcorrected_filter@meta.data %>%
  filter(
    CTaa == "CAEMYSGGGADGLTF_CASSLFSGNNEQFF",
    patient == "0712"
  ) %>%
  count(cell_type_annotated_granularity, name = "cell_count") %>%
  arrange(desc(cell_count))

# 打印看一下结果
print(clone_table2)

# 保存到当前工作目录
write.csv(clone_table2, "clone_CAEMYSGGGADGLTF_CASSLFSGNNEQFF_0712_counts.csv", row.names = FALSE)
p_alluvial2 <- alluvialClones(
    merged_all_T_batchcorrected_filter,
    cloneCall = "cTaa",   # 如果 "aa" 报错，改成 "CTaa"（你真实的AA克隆列名）
    y.axes    = c("batch", "cell_type_annotated_granularity", "label"), 
    color     = c("CAEMYSGGGADGLTF_CASSLFSGNNEQFF")  # 要高亮的克隆
)
pdf("riverplot_0712_CAEMYSGGGADGLTF_CASSLFSGNNEQFF.pdf",15,9)
print(p_alluvial2)
dev.off()

save(combined, merged_all_T_batchcorrected_filter,file="merged_4_with_vdj_T_cells_filter_10_28_2025.RData")
# occupiedscRepertoire() -> clonalOccupy()
pdf("occupiedscRepertoire_all.pdf",12,10)
clonalOccupy(
    sc      = merged_all_T_batchcorrected_filter,
    x.axis  = "cell_type_annotated_granularity"
)
dev.off()

pdf("occupiedscRepertoire_shared_clone.pdf")
clonalOccupy(
    sc      = all_clone_test,
    x.axis  = "cell_type_annotated_granularity"
)
dev.off()

# alluvialClonotypes() -> alluvialClones()
pdf("riverplot_clone_info_all_clone_by_cluster.pdf")
alluvialClones(
    sc        = all_clone_test,
    cloneCall = "gene",
    # y.axes renamed? In v2 docs this stays y.axes or y.axes/by? 
    # 最新说明里 alluvialClones() retains y.axes-style arg to pick columns for the strata.:contentReference[oaicite:6]{index=6}
    y.axes    = c("batch", "seurat_clusters"),
    color     = "seurat_clusters"
)
dev.off()

###plot shared clones between each patients' graft and PBMC
merged_all_T_batchcorrected_filter$CTaa



# Circos: getCirclize() still exists and returns matrix for circlize::chordDiagram()
library(circlize)
library(scales)

circles <- getCirclize(
    sc       = all_clone_test,
    group.by = "ident"
)

grid.cols <- scales::hue_pal()(length(unique(all_clone_test@active.ident)))
names(grid.cols) <- levels(all_clone_test@active.ident)

circlize::chordDiagram(
    circles,
    self.link = 1,
    grid.col  = grid.cols
)

pdf("circlize_clone_info_all_clone_by_cluster.pdf")
circlize::chordDiagram(
    circles,
    self.link = 1,
    grid.col  = grid.cols
)
dev.off()

circlize::chordDiagram(
    circles,
    self.link        = 1,
    grid.col         = grid.cols,
    directional      = 1,
    direction.type   = "arrows",
    link.arr.type    = "big.arrow"
)


### 9. clusterTCR, StartracDiversity, clonotypeBias ###

# clusterTCR(): still clusterTCR() in v2, with similar args.
sub_combined <- clusterTCR(
    input.data = combined[[2]],
    chain      = "TRA",
    sequence   = "aa",
    threshold  = 0.85,
    group.by   = NULL
)

pdf("StartracDiversity_shared_clone.pdf")
StartracDiversity(
    sc      = all_clone_test,
    type    = "seurat_clusters",
    sample  = "batch",
    by      = "overall"
)
dev.off()
star_tab <- StartracDiversity(
    sc   = merged_all_T_batchcorrected_filter,
    type = "batch"
)
pdf("StartracDiversity_all.pdf")
StartracDiversity(
    sc      = merged_all_T_batchcorrected_filter,
    type    = "seurat_clusters",
    sample  = "batch",
    by      = "overall"
)
dev.off()
merged_all_T_batchcorrected_filter$patient <- merged_all_T_batchcorrected_filter$batch
tmp <- gsub("PBMC_","",merged_all_T_batchcorrected_filter$patient)
tmp <- gsub("graft_","",tmp)
merged_all_T_batchcorrected_filter$patient <- tmp
Idents(merged_all_T_batchcorrected_filter) <- merged_all_T_batchcorrected_filter$patient

pdf("clone_sparsity_by_sample.pdf")
clonalRarefaction(merged_all_T_batchcorrected_filter,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)
dev.off()
pdf("clone_sparsity_by_sample_by_tissue.pdf")				  
clonalRarefaction(combined,
                  plot.type = 1,
                  hill.numbers = 0,
                  n.boots = 2)
dev.off()
pdf("clone_bias_shared_clone.pdf")
clonotypeBias(
    sc          = all_clone_test,
    cloneCall   = "aa",
    group.by    = "batch",
    n.boots     = 20,
    min.expand  = 10
)
dev.off()

pdf("clone_bias_all.pdf")
clonalBias(
    sc        = merged_all_T_batchcorrected_filter,
    cloneCall = "aa",
    group.by  = "batch",
    n.boots   = 20,
    min.expand= 2
)
dev.off()

### 10. expression2List -> cluster-level analyses ###

combined2 <- expression2List(
    sc       = merged_all_T_batchcorrected_filter,
    split.by = "SCT_snn_res.0.6"
)
# clonalDiversity(), clonalHomeostasis(), clonalProportion(), clonalOverlap()
# still keep those names in v2 (unchanged).:contentReference[oaicite:7]{index=7}

pdf("clonalDiversity_all1.pdf")
clonalDiversity(
    input.data = combined2,
    split.by   = "ident",
    cloneCall  = "nt",
    chain      = "both"
)
dev.off()

pdf("clonalDiversity_all2.pdf",10,6)
clonalDiversity(
    input.data   = merged_all_T_batchcorrected_filter,
    cloneCall    = "nt",
    chain        = "both",
    group.by     = "cell_type_annotated_granularity",
    exportTable  = FALSE
)
dev.off()

pdf("clonalDiversity_all_by_patient.pdf",10,6)
clonalDiversity(
    input.data   = merged_all_T_batchcorrected_filter,
    cloneCall    = "nt",
    chain        = "both",
    group.by     = "patient",
    exportTable  = FALSE
)
dev.off()

pdf("clonalDiversity_all_by_patient_tissue.pdf",10,6)
clonalDiversity(
    input.data   = merged_all_T_batchcorrected_filter,
    cloneCall    = "nt",
    chain        = "both",
    group.by     = "batch",
    exportTable  = FALSE
)
dev.off()

pdf("clonalHomeostasis_all.pdf")
clonalHomeostasis(
    input.data = combined2,
    cloneCall  = "nt",
    chain      = "both"
)
dev.off()

pdf("clonalProportion_all.pdf")
clonalProportion(
    input.data = combined2,
    cloneCall  = "nt",
    chain      = "both"
)
dev.off()

pdf("clonalOverlap_all.pdf")
clonalOverlap(
    input.data = combined2,
    cloneCall  = "aa",
    chain      = "both",
    method     = "overlap"
)
dev.off()


combined3 <- expression2List(
    sc       = all_clone_test,
    split.by = "SCT_snn_res.0.8"
)

pdf("clonalDiversity_shared_clone1.pdf")
clonalDiversity(
    input.data = combined3,
    split.by   = "ident",
    cloneCall  = "nt",
    chain      = "both"
)
dev.off()

pdf("clonalDiversity_shared_clone2.pdf")
clonalDiversity(
    input.data = all_clone_test,
    split.by   = "ident",
    cloneCall  = "nt",
    chain      = "both"
)
dev.off()

pdf("clonalHomeostasis_shared_clone.pdf")
clonalHomeostasis(
    input.data = combined3,
    cloneCall  = "nt",
    chain      = "both"
)
dev.off()

pdf("clonalProportion_shared_clone.pdf")
clonalProportion(
    input.data = combined3,
    cloneCall  = "nt",
    chain      = "both"
)
dev.off()

pdf("clonalOverlap_shared_clone.pdf")
clonalOverlap(
    input.data = combined3,
    cloneCall  = "aa",
    chain      = "both",
    method     = "overlap"
)
dev.off()

###update shared clone plot
###plot top 3 CDR3 sequencee in graft

tmp <- table(merged_all_T_batchcorrected_filter$t_cdr3s_aa, merged_all_T_batchcorrected_filter$label)
t_idx <- {}
##get shared clone in either Graft and PBMC
for (t in 1:nrow(tmp)){
##with at least three clone
	if (tmp[t,1]>3 && tmp[t,2] >3){
		t_idx <- c(t_idx,t)
	}
}
tmp = tmp[t_idx,]
idx <- order(tmp[,2],decreasing=T)
tmp <- tmp[idx,]
system("mkdir -p shared_clone_update")

for (i in rownames(tmp)){
print(i)
Idents(merged_all_T_batchcorrected_filter) <- "clone_by_label"
if (( length(try(WhichCells(object = merged_all_T_batchcorrected_filter, idents = paste0(i,"_Graft")),TRUE))>=1) & (length( try(WhichCells(object = merged_all_T_batchcorrected_filter, idents = paste0(i,"_PBMC")),TRUE))>=1) ){
cluster.markers <- FindMarkers(object = merged_all_T_batchcorrected_filter,min.cells.group=1, ident.1 = paste0(i,"_Graft"), ident.2 = paste0(i,"_PBMC") , only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0)
colnames(cluster.markers)[2] <- "avg_logFC"
pdf(paste0("shared_clone_update/",i,".DEGs_volcano_plot.pdf"))
print(GEX_volcano(cluster.markers,input.type= "findmarkers",condition.1="Graft", condition.2="PBMC",n.label.up=10,n.label.down=10))
dev.off()

write.table(cluster.markers, paste0("shared_clone_update/cluster.",i,".DEGs_Graft_vs_PMBC.xls"),sep="\t",row.names=T,col.names=T,quote=F)
g1_Graft <- WhichCells(merged_all_T_batchcorrected_filter, idents = paste0(i,"_Graft"))
patient_tab <- table(merged_all_T_batchcorrected_filter$patient[g1_Graft])
if (length(patient_tab) == 1) {
  g1_Graft <- g1_Graft
  top_patient <- names(which.max(patient_tab))
} else {
  top_patient <- names(which.max(patient_tab))
  message("⚠ 多个 patient，保留数量最多的：", top_patient)
  g1_Graft <- g1_Graft[merged_all_T_batchcorrected_filter$patient[g1_Graft] == top_patient]
}
g1_PBMC <- WhichCells(merged_all_T_batchcorrected_filter, idents = paste0(i,"_PBMC"))
patient_tab <- table(merged_all_T_batchcorrected_filter$patient[g1_PBMC])
if (length(patient_tab) == 1) {
  g1_PBMC <- g1_PBMC
  top_patient <- names(which.max(patient_tab))
} else {
  top_patient <- names(which.max(patient_tab))
  message("⚠ 多个 patient，保留数量最多的：", top_patient)
  g1_PBMC <- g1_PBMC[merged_all_T_batchcorrected_filter$patient[g1_PBMC] == top_patient]
}
total_cells = sum(length(g1_PBMC)+length(g1_Graft))
pdf(paste0("shared_clone_update/",top_patient,"_",i,".DEGs_graft_vs_PBMC.pdf"),6,5)
print(DimPlot(merged_all_T_batchcorrected_filter, label=T, label.size = 2,group.by="cell_type_annotated_granularity", cells.highlight= list("Graft" = g1_Graft, "PBMC" = g1_PBMC), cols.highlight = c("darkblue", "red"), cols= "grey")+ 
  ggtitle(paste0("patient: ",top_patient," clone expension: ",total_cells," cells")) +       # ← 把标题改成变量 top_patient 的内容
  theme(
    plot.title = element_text(hjust = 0.5)  # 居中显示标题，可选
  ))
dev.off()
}
}

###parse all cells with shared clones from whole objects 10/28/2025
unique_clone_graft <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$label == "Graft")])
unique_clone_PBMC <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$label == "PBMC")])
shared_clone_list <- intersect(unique_clone_graft,unique_clone_PBMC)
shared_clone_list <-  shared_clone_list[!is.na(shared_clone_list)]
Idents(object = merged_all_T_batchcorrected_filter)<- merged_all_T_batchcorrected_filter$label
merged_all_T_batchcorrected_filter$clone_by_label <- paste(merged_all_T_batchcorrected_filter$t_cdr3s_aa,Idents(object = merged_all_T_batchcorrected_filter), sep = "_")
##
unique_clone_graft <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$batch == "graft_0825")])
unique_clone_PBMC <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$batch == "PBMC_0825")])
shared_clone_list <- intersect(unique_clone_graft,unique_clone_PBMC)
shared_clone_list <-  shared_clone_list[!is.na(shared_clone_list)]
cells_shared_0825 <- rownames(
  merged_all_T_batchcorrected_filter@meta.data[
    merged_all_T_batchcorrected_filter$t_cdr3s_aa %in% shared_clone_list &
    merged_all_T_batchcorrected_filter$patient == "0825",
  ]
)
unique_clone_graft <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$batch == "graft_0701")])
unique_clone_PBMC <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$batch == "PBMC_0701")])
shared_clone_list <- intersect(unique_clone_graft,unique_clone_PBMC)
shared_clone_list <-  shared_clone_list[!is.na(shared_clone_list)]
cells_shared_0701 <- rownames(
  merged_all_T_batchcorrected_filter@meta.data[
    merged_all_T_batchcorrected_filter$t_cdr3s_aa %in% shared_clone_list &
    merged_all_T_batchcorrected_filter$patient == "0701",
  ]
)
unique_clone_graft <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$batch == "graft_0708")])
unique_clone_PBMC <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$batch == "PBMC_0708")])
shared_clone_list <- intersect(unique_clone_graft,unique_clone_PBMC)
shared_clone_list <-  shared_clone_list[!is.na(shared_clone_list)]
cells_shared_0708 <- rownames(
  merged_all_T_batchcorrected_filter@meta.data[
    merged_all_T_batchcorrected_filter$t_cdr3s_aa %in% shared_clone_list &
    merged_all_T_batchcorrected_filter$patient == "0708",
  ]
)
unique_clone_graft <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$batch == "graft_0712")])
unique_clone_PBMC <- unique(merged_all_T_batchcorrected_filter$t_cdr3s_aa[which(merged_all_T_batchcorrected_filter$batch == "PBMC_0712")])
shared_clone_list <- intersect(unique_clone_graft,unique_clone_PBMC)
shared_clone_list <-  shared_clone_list[!is.na(shared_clone_list)]
cells_shared_0712 <- rownames(
  merged_all_T_batchcorrected_filter@meta.data[
    merged_all_T_batchcorrected_filter$t_cdr3s_aa %in% shared_clone_list &
    merged_all_T_batchcorrected_filter$patient == "0712",
  ]
)
cells_shared_all <- c(cells_shared_0712,cells_shared_0708,cells_shared_0701,cells_shared_0825)

##get shared clones all and conduct the DEG together
merged_all_T_batchcorrected_filter$shared_clone_list <- "non_shared_clone"
shared_clone_list_graft_cells <- {}
shared_clone_list_pbmc_cells <- {}
shared_clone_list_graft_cells <- intersect(cells_shared_all,  rownames(merged_all_T_batchcorrected_filter@meta.data[which(merged_all_T_batchcorrected_filter$label == "Graft"),] ))
shared_clone_list_pbmc_cells <- intersect(cells_shared_all, rownames(merged_all_T_batchcorrected_filter@meta.data[which(merged_all_T_batchcorrected_filter$label == "PBMC"),] ))

merged_all_T_batchcorrected_filter$shared_clone_list[match(shared_clone_list_graft_cells,rownames(merged_all_T_batchcorrected_filter@meta.data))] <- "shared_clone_Graft"
merged_all_T_batchcorrected_filter$shared_clone_list[match(shared_clone_list_pbmc_cells,rownames(merged_all_T_batchcorrected_filter@meta.data))] <- "shared_clone_PBMC"
Idents(merged_all_T_batchcorrected_filter) <- "shared_clone_list"
markers <- FindMarkers(object = merged_all_T_batchcorrected_filter, ident.1 = "shared_clone_Graft", ident.2 = "shared_clone_PBMC" , only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0)
write.table(markers, "./shared_clone_update/shared_clone_DEGs_Graft_vs_PMBC.xls",sep="\t",row.names=T,col.names=T,quote=F)
colnames(markers)[2] <- "avg_logFC"
pdf("./shared_clone_update/comparison_between_all_shared_clones_graft_vs_PBMC.pdf")
GEX_volcano(markers,input.type= "findmarkers",condition.1="Graft", condition.2="PBMC",n.label.up=10,n.label.down=10)
dev.off()
Idents(merged_all_T_batchcorrected_filter) <- "shared_clone_list"
g1_Graft <- WhichCells(merged_all_T_batchcorrected_filter, idents = c("shared_clone_Graft"))
g1_PBMC <- WhichCells(merged_all_T_batchcorrected_filter, idents = c( "shared_clone_PBMC"))
pdf("./shared_clone_update/pseudo_bulk_clone_vs.pdf",8,8)
DimPlot(merged_all_T_batchcorrected_filter, label=T, label.size =2, group.by="cell_type_annotated_granularity", cells.highlight= list(g1_Graft, g1_PBMC), cols.highlight = c("darkblue", "red"), cols= "grey")
dev.off()

###
shared_clone_in_batch_detail <- table(merged_all_T_batchcorrected_filter$t_cdr3s_aa[match(cells_shared_all,rownames(merged_all_T_batchcorrected_filter@meta.data))],merged_all_T_batchcorrected_filter$batch[match(cells_shared_all,rownames(merged_all_T_batchcorrected_filter@meta.data))])
idx_tmp <- order(shared_clone_in_batch_detail[,1],decreasing =T)
shared_clone_in_batch_detail <- shared_clone_in_batch_detail[idx_tmp,]
write.table(shared_clone_in_batch_detail,file="./shared_clone_update/clone_in_batch_detail.txt",quote=F,sep="\t")
all_clone_test_filter <- subset( merged_all_T_batchcorrected_filter, cells = cells_shared_all)

pdf("all_clone_test_filter_by_source.pdf",9,6)
DimPlot(all_clone_test_filter,pt.size=1)
dev.off()
pdf("all_clone_test_filter_by_cell_type.pdf",9,6)
DimPlot(all_clone_test_filter,pt.size=1,group.by = "cell_type_annotated_granularity")
dev.off()
# 生成一个数据框，只包含 cell_type_annotated_granularity 列
celltype_table <- data.frame(
  cell_type_annotated_granularity = all_clone_test_filter$cell_type_annotated_granularity
)

# 保存成 CSV 文件

celltype_summary <- as.data.frame(table(all_clone_test_filter$cell_type_annotated_granularity))
write.csv(celltype_summary, "cell_type_summary.csv", row.names = FALSE)
save(all_clone_test_filter,"all_clone_filter_10_28.Rdata")

##monocle3
library("Seurat")
library("Moncole3")
##monocle3
#input embedding from seurat
setwd("D:\\Dropbox\\collaboration\\for dawei\\manuscript\\transplant\\result_update_8_24_2025\\merge_4_T_cells_with_vdj\\")
load("all_clone_filter_10_28.Rdata")

library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(tibble)

## 1. 准备 cds
data <- GetAssayData(all_clone_test_filter, assay = "RNA", slot = "count")  # 或 layer="counts" 如果你是Seurat v5
cell_metadata <- all_clone_test_filter@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# 确保行名/列名对应：colnames(data) 应该等于 rownames(cell_metadata)
stopifnot(identical(colnames(data), rownames(cell_metadata)))

cds <- new_cell_data_set(
  data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 现在 cds 里还没有我们想要的 UMAP；我们马上塞进去

## 2. 拿到 Seurat 的 UMAP
seurat_umap <- Embeddings(all_clone_test_filter, reduction = "umap")
# 这是一个矩阵，行名是细胞名，列通常是 UMAP_1, UMAP_2

## 3. 对齐细胞顺序
cds_cells <- colnames(cds)

# 检查交集（调试用，不是必须，但能确认没有大问题）
common_cells <- intersect(cds_cells, rownames(seurat_umap))
message("length(cds_cells) = ", length(cds_cells))
message("length(rownames(seurat_umap)) = ", nrow(seurat_umap))
message("length(common_cells) = ", length(common_cells))

# 如果 length(common_cells) != length(cds_cells)，
# 说明有的细胞名不匹。那下面的对齐会产生 NA 行，这会告诉我们是不是命名问题。

# 按照 monocle3 的细胞顺序重排 Seurat 的 UMAP
seurat_umap_aligned <- seurat_umap[cds_cells, , drop = FALSE]

# 看看有没有 NA（大量 NA 说明名字没对齐）
print(colSums(is.na(seurat_umap_aligned)))

## 4. 把这个 UMAP 塞进 monocle3 的 reducedDims
# 注意：是 "UMAP" (大写)，不是 "umap"
cds@int_colData$reducedDims$UMAP <- seurat_umap_aligned

## 5. 画图，用 monocle3 但坐标就是 Seurat 的
p <- plot_cells(
  cds,
  reduction_method = "UMAP",  # 这里会用我们刚刚放进去的 cds@int_colData$reducedDims$UMAP
  color_cells_by = "cell_type_annotated_granularity",  # 这个列必须存在于 cds@colData
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
) +
  theme(
    plot.title = element_blank()
  )

print(p)


##monocle3 cluster
cds <- cluster_cells(cds,resolution=0.00105)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

cds$clusters <- cds$cell_type_annotated_granularity
pdf("monocle3_by_clusters.pdf",8,6)
plot_cells(cds,
           color_cells_by = "clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


##CD4+, CD8+ separate 
##monocle3 cluster
cds <- cluster_cells(cds,resolution=0.005)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_split_CD8.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

gene_list_curves = c("TCF7", "LEF1", "KLF2", "IL7R", "CISH", "JUN", "KLRD1", "CXCR6", "IFNG", "GZMB", "PRF1", "ID2")
##plot gene curves
pdf("genes_pseudotime_curves_split_CD8.pdf", 8, 35)
plot_genes_in_pseudotime(cds[gene_list_curves,],
                         color_cells_by = "pseudotime",
                         min_expr = 0.5)
dev.off()


cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_split_CD4_update.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

pdf("genes_pseudotime_curves_split_CD4.pdf", 8, 35)
plot_genes_in_pseudotime(cds[gene_list_curves,],
                         color_cells_by = "pseudotime",
                         min_expr = 0.5)
dev.off()


cds$clusters <- cds$cell_type_annotated_granularity
pdf("monocle3_by_clusters_split_CD4_CD8.pdf",8,6)
plot_cells(cds,
           color_cells_by = "clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

###split the individuals
all_clone_test_filter
# 用 SplitObject 一次性拆分
lst <- SplitObject(all_clone_test_filter, split.by = "patient")

# 取出四个对象（若某个分组不存在会报错，便于及时发现）
seurat_0701 <- lst[["0701"]]
seurat_0708 <- lst[["0708"]]
seurat_0712 <- lst[["0712"]]
seurat_0825 <- lst[["0825"]]

## 1. 准备 cds
data <- GetAssayData(seurat_0701, assay = "RNA", slot = "count")  # 或 layer="counts" 如果你是Seurat v5
cell_metadata <- seurat_0701@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# 确保行名/列名对应：colnames(data) 应该等于 rownames(cell_metadata)
stopifnot(identical(colnames(data), rownames(cell_metadata)))

cds <- new_cell_data_set(
  data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 现在 cds 里还没有我们想要的 UMAP；我们马上塞进去

## 2. 拿到 Seurat 的 UMAP
seurat_umap <- Embeddings(seurat_0701, reduction = "umap")
# 这是一个矩阵，行名是细胞名，列通常是 UMAP_1, UMAP_2

## 3. 对齐细胞顺序
cds_cells <- colnames(cds)

# 检查交集（调试用，不是必须，但能确认没有大问题）
common_cells <- intersect(cds_cells, rownames(seurat_umap))
message("length(cds_cells) = ", length(cds_cells))
message("length(rownames(seurat_umap)) = ", nrow(seurat_umap))
message("length(common_cells) = ", length(common_cells))

# 如果 length(common_cells) != length(cds_cells)，
# 说明有的细胞名不匹。那下面的对齐会产生 NA 行，这会告诉我们是不是命名问题。

# 按照 monocle3 的细胞顺序重排 Seurat 的 UMAP
seurat_umap_aligned <- seurat_umap[cds_cells, , drop = FALSE]

# 看看有没有 NA（大量 NA 说明名字没对齐）
print(colSums(is.na(seurat_umap_aligned)))

## 4. 把这个 UMAP 塞进 monocle3 的 reducedDims
# 注意：是 "UMAP" (大写)，不是 "umap"
cds@int_colData$reducedDims$UMAP <- seurat_umap_aligned

## 5. 画图，用 monocle3 但坐标就是 Seurat 的
p <- plot_cells(
  cds,
  reduction_method = "UMAP",  # 这里会用我们刚刚放进去的 cds@int_colData$reducedDims$UMAP
  color_cells_by = "cell_type_annotated_granularity",  # 这个列必须存在于 cds@colData
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
) +
  theme(
    plot.title = element_blank()
  )

print(p)

#CD4+ trajectory source should be naive CD4+
##monocle3 cluster
cds <- cluster_cells(cds,resolution=0.09)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_0701_CD4+.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


#CD8+ trajectory source should be tcf7+ and naive/memory-like CD8+
##monocle3 cluster
cds <- cluster_cells(cds,resolution=0.09)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_0701_CD8+.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

cds$clusters <- cds$cell_type_annotated_granularity

pdf("monocle3_0701_by_clusters.pdf",8,6)
plot_cells(cds,
           color_cells_by = "clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


pdf("monocle3_0701_by_batches.pdf",8,6)
plot_cells(cds,
           color_cells_by = "batch",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

##0708
## 1. 准备 cds
data <- GetAssayData(seurat_0708, assay = "RNA", slot = "count")  # 或 layer="counts" 如果你是Seurat v5
cell_metadata <- seurat_0708@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# 确保行名/列名对应：colnames(data) 应该等于 rownames(cell_metadata)
stopifnot(identical(colnames(data), rownames(cell_metadata)))

cds <- new_cell_data_set(
  data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 现在 cds 里还没有我们想要的 UMAP；我们马上塞进去

## 2. 拿到 Seurat 的 UMAP
seurat_umap <- Embeddings(seurat_0708, reduction = "umap")
# 这是一个矩阵，行名是细胞名，列通常是 UMAP_1, UMAP_2

## 3. 对齐细胞顺序
cds_cells <- colnames(cds)

# 检查交集（调试用，不是必须，但能确认没有大问题）
common_cells <- intersect(cds_cells, rownames(seurat_umap))
message("length(cds_cells) = ", length(cds_cells))
message("length(rownames(seurat_umap)) = ", nrow(seurat_umap))
message("length(common_cells) = ", length(common_cells))

# 如果 length(common_cells) != length(cds_cells)，
# 说明有的细胞名不匹。那下面的对齐会产生 NA 行，这会告诉我们是不是命名问题。

# 按照 monocle3 的细胞顺序重排 Seurat 的 UMAP
seurat_umap_aligned <- seurat_umap[cds_cells, , drop = FALSE]

# 看看有没有 NA（大量 NA 说明名字没对齐）
print(colSums(is.na(seurat_umap_aligned)))

## 4. 把这个 UMAP 塞进 monocle3 的 reducedDims
# 注意：是 "UMAP" (大写)，不是 "umap"
cds@int_colData$reducedDims$UMAP <- seurat_umap_aligned

## 5. 画图，用 monocle3 但坐标就是 Seurat 的
p <- plot_cells(
  cds,
  reduction_method = "UMAP",  # 这里会用我们刚刚放进去的 cds@int_colData$reducedDims$UMAP
  color_cells_by = "cell_type_annotated_granularity",  # 这个列必须存在于 cds@colData
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
) +
  theme(
    plot.title = element_blank()
  )

print(p)

#CD4+ trajectory source should be naive CD4+
##monocle3 cluster
cds <- cluster_cells(cds,resolution=1)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_0708_CD4+.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


#CD8+ trajectory source should be tcf7+ and naive/memory-like CD8+
##monocle3 cluster
cds <- cluster_cells(cds,resolution=1)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_0708_CD8+.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

cds$clusters <- cds$cell_type_annotated_granularity

pdf("monocle3_0708_by_clusters.pdf",8,6)
plot_cells(cds,
           color_cells_by = "clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


pdf("monocle3_0708_by_batches.pdf",8,6)
plot_cells(cds,
           color_cells_by = "batch",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


##0712
## 1. 准备 cds
data <- GetAssayData(seurat_0712, assay = "RNA", slot = "count")  # 或 layer="counts" 如果你是Seurat v5
cell_metadata <- seurat_0712@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# 确保行名/列名对应：colnames(data) 应该等于 rownames(cell_metadata)
stopifnot(identical(colnames(data), rownames(cell_metadata)))

cds <- new_cell_data_set(
  data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 现在 cds 里还没有我们想要的 UMAP；我们马上塞进去

## 2. 拿到 Seurat 的 UMAP
seurat_umap <- Embeddings(seurat_0712, reduction = "umap")
# 这是一个矩阵，行名是细胞名，列通常是 UMAP_1, UMAP_2

## 3. 对齐细胞顺序
cds_cells <- colnames(cds)

# 检查交集（调试用，不是必须，但能确认没有大问题）
common_cells <- intersect(cds_cells, rownames(seurat_umap))
message("length(cds_cells) = ", length(cds_cells))
message("length(rownames(seurat_umap)) = ", nrow(seurat_umap))
message("length(common_cells) = ", length(common_cells))

# 如果 length(common_cells) != length(cds_cells)，
# 说明有的细胞名不匹。那下面的对齐会产生 NA 行，这会告诉我们是不是命名问题。

# 按照 monocle3 的细胞顺序重排 Seurat 的 UMAP
seurat_umap_aligned <- seurat_umap[cds_cells, , drop = FALSE]

# 看看有没有 NA（大量 NA 说明名字没对齐）
print(colSums(is.na(seurat_umap_aligned)))

## 4. 把这个 UMAP 塞进 monocle3 的 reducedDims
# 注意：是 "UMAP" (大写)，不是 "umap"
cds@int_colData$reducedDims$UMAP <- seurat_umap_aligned

## 5. 画图，用 monocle3 但坐标就是 Seurat 的
p <- plot_cells(
  cds,
  reduction_method = "UMAP",  # 这里会用我们刚刚放进去的 cds@int_colData$reducedDims$UMAP
  color_cells_by = "cell_type_annotated_granularity",  # 这个列必须存在于 cds@colData
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
) +
  theme(
    plot.title = element_blank()
  )

print(p)

#CD4+ trajectory source should be naive CD4+
##monocle3 cluster
cds <- cluster_cells(cds,resolution=0.0065)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_0712_CD4+.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


#CD8+ trajectory source should be tcf7+ and naive/memory-like CD8+
##monocle3 cluster
cds <- cluster_cells(cds,resolution=0.0065)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_0712_CD8+.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

cds$clusters <- cds$cell_type_annotated_granularity

pdf("monocle3_0712_by_clusters.pdf",8,6)
plot_cells(cds,
           color_cells_by = "clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


pdf("monocle3_0712_by_batches.pdf",8,6)
plot_cells(cds,
           color_cells_by = "batch",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


##0825
seurat_0825
# 定义各组对应的标签
cd4_clusters <- c("Memory-like CD4+", "CCR2+ CD4+", "CXCL13+ CD4+")
cd8_clusters <- c("CXCL13+CXCR6+ effector CD8+", 
                  "CXCR6+ effector CD8+", 
                  "Naive/Memory-like CD8+", 
                  "CX3CR1+ effector CD8+", 
                  "ZNF683+TCF7hi CD8+", 
                  "ISG-high CD8+")

# 子集化生成新的 Seurat 对象
seurat_0825_CD4 <- subset(seurat_0825, subset = cell_type_annotated_granularity %in% cd4_clusters)
seurat_0825_CD8 <- subset(seurat_0825, subset = cell_type_annotated_granularity %in% cd8_clusters)

## 1. 准备 cds seurat_0825_CD4
data <- GetAssayData(seurat_0825_CD4, assay = "RNA", slot = "count")  # 或 layer="counts" 如果你是Seurat v5
cell_metadata <- seurat_0825_CD4@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# 确保行名/列名对应：colnames(data) 应该等于 rownames(cell_metadata)
stopifnot(identical(colnames(data), rownames(cell_metadata)))

cds <- new_cell_data_set(
  data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 现在 cds 里还没有我们想要的 UMAP；我们马上塞进去

## 2. 拿到 Seurat 的 UMAP
seurat_umap <- Embeddings(seurat_0825_CD4, reduction = "umap")
# 这是一个矩阵，行名是细胞名，列通常是 UMAP_1, UMAP_2

## 3. 对齐细胞顺序
cds_cells <- colnames(cds)

# 检查交集（调试用，不是必须，但能确认没有大问题）
common_cells <- intersect(cds_cells, rownames(seurat_umap))
message("length(cds_cells) = ", length(cds_cells))
message("length(rownames(seurat_umap)) = ", nrow(seurat_umap))
message("length(common_cells) = ", length(common_cells))

# 如果 length(common_cells) != length(cds_cells)，
# 说明有的细胞名不匹。那下面的对齐会产生 NA 行，这会告诉我们是不是命名问题。

# 按照 monocle3 的细胞顺序重排 Seurat 的 UMAP
seurat_umap_aligned <- seurat_umap[cds_cells, , drop = FALSE]

# 看看有没有 NA（大量 NA 说明名字没对齐）
print(colSums(is.na(seurat_umap_aligned)))

## 4. 把这个 UMAP 塞进 monocle3 的 reducedDims
# 注意：是 "UMAP" (大写)，不是 "umap"
cds@int_colData$reducedDims$UMAP <- seurat_umap_aligned

## 5. 画图，用 monocle3 但坐标就是 Seurat 的
p <- plot_cells(
  cds,
  reduction_method = "UMAP",  # 这里会用我们刚刚放进去的 cds@int_colData$reducedDims$UMAP
  color_cells_by = "cell_type_annotated_granularity",  # 这个列必须存在于 cds@colData
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
) +
  theme(
    plot.title = element_blank()
  )

print(p)

#CD4+ trajectory source should be naive CD4+
##monocle3 cluster
cds <- cluster_cells(cds,resolution=0.01)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_0825_CD4+.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

####
## 1. 准备 cds seurat_0825_CD8
data <- GetAssayData(seurat_0825_CD8, assay = "RNA", slot = "count")  # 或 layer="counts" 如果你是Seurat v5
cell_metadata <- seurat_0825_CD8@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)

# 确保行名/列名对应：colnames(data) 应该等于 rownames(cell_metadata)
stopifnot(identical(colnames(data), rownames(cell_metadata)))

cds <- new_cell_data_set(
  data,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")
# 现在 cds 里还没有我们想要的 UMAP；我们马上塞进去

## 2. 拿到 Seurat 的 UMAP
seurat_umap <- Embeddings(seurat_0825_CD8, reduction = "umap")
# 这是一个矩阵，行名是细胞名，列通常是 UMAP_1, UMAP_2

## 3. 对齐细胞顺序
cds_cells <- colnames(cds)

# 检查交集（调试用，不是必须，但能确认没有大问题）
common_cells <- intersect(cds_cells, rownames(seurat_umap))
message("length(cds_cells) = ", length(cds_cells))
message("length(rownames(seurat_umap)) = ", nrow(seurat_umap))
message("length(common_cells) = ", length(common_cells))

# 如果 length(common_cells) != length(cds_cells)，
# 说明有的细胞名不匹。那下面的对齐会产生 NA 行，这会告诉我们是不是命名问题。

# 按照 monocle3 的细胞顺序重排 Seurat 的 UMAP
seurat_umap_aligned <- seurat_umap[cds_cells, , drop = FALSE]

# 看看有没有 NA（大量 NA 说明名字没对齐）
print(colSums(is.na(seurat_umap_aligned)))

## 4. 把这个 UMAP 塞进 monocle3 的 reducedDims
# 注意：是 "UMAP" (大写)，不是 "umap"
cds@int_colData$reducedDims$UMAP <- seurat_umap_aligned

## 5. 画图，用 monocle3 但坐标就是 Seurat 的
p <- plot_cells(
  cds,
  reduction_method = "UMAP",  # 这里会用我们刚刚放进去的 cds@int_colData$reducedDims$UMAP
  color_cells_by = "cell_type_annotated_granularity",  # 这个列必须存在于 cds@colData
  show_trajectory_graph = FALSE,
  label_cell_groups = FALSE
) +
  theme(
    plot.title = element_blank()
  )

print(p)

#CD8+ trajectory source should be tcf7+ and naive/memory-like CD8+
##monocle3 cluster
cds <- cluster_cells(cds,resolution=0.002)
plot_cells(cds,show_trajectory_graph = FALSE)

cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "cluster", group_label_size = 
             8, cell_size = 1)
plot_cells(cds, label_branch_points = FALSE, label_groups_by_cluster = FALSE, label_leaves=FALSE,cell_size = 1)
cds <- order_cells(cds)

pdf("monocle3_by_pseudotime_0825_CD8+.pdf",8,6)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size =0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

pdf("seurat_0825_DimPlot.pdf",8,6)
DimPlot(seurat_0825,group.by = "cell_type_annotated_granularity")
dev.off()


cds$clusters <- cds$cell_type_annotated_granularity

pdf("monocle3_0825_by_clusters.pdf",8,6)
plot_cells(cds,
           color_cells_by = "clusters",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()


pdf("monocle3_0825_by_batches.pdf",8,6)
plot_cells(cds,
           color_cells_by = "batch",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.9,
           graph_label_size=1.5,label_roots=F)+theme(axis.line.x  = element_blank(),axis.line.y  = element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),axis.title  = element_blank())
dev.off()

###
all_clone_test_filter.markers <- FindAllMarkers(object = all_clone_test_filter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all_clone_test_filter.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC) -> test
ttt <- data.frame(test)
##reorder the result##
ttt <- ttt[order(ttt[,6],ttt[,2],decreasing =T),]
test <- ttt[order(ttt[,6]),]
DefaultAssay(all_clone_test_filter) <- "integrated"
heatmap_matrix <- GetAssayData(object = all_clone_test_filter, slot = "counts")
all_markers_need <- unique(test[,7])
pdf("shared_clone_T_cells_Heatmap_All_scale_data_final.pdf",20,12)
DoHeatmap(subset(all_clone_test_filter, downsample = 200),features = all_markers_need, size = 3,label=F)
dev.off()

DefaultAssay(all_clone_test_filter) <- "RNA"
pdf("shared_clone_T_cells_feature_plot_gene_list0_final.pdf",16,length(gene_list_0)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list1_final.pdf",16,length(gene_list_1)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list2_final.pdf",16,length(gene_list_2)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list3_final.pdf",16,length(gene_list_3)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list4_final.pdf",16,length(gene_list_4)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list5_final.pdf",16,length(gene_list_5)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list6_final.pdf",16,length(gene_list_6)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list7_final.pdf",16,length(gene_list_7)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list8_final.pdf",16,length(gene_list_8)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list9_final.pdf",16,length(gene_list_9)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list10_final.pdf",16,length(gene_list_10)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list11_final.pdf",16,length(gene_list_11)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list12_final.pdf",16,length(gene_list_12)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list13_final.pdf",16,length(gene_list_13)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list14_final.pdf",16,length(gene_list_14)*4)
FeaturePlot(object = all_clone_test_filter, features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()


DefaultAssay(all_clone_test_filter) <- "RNA"
###split by batch
pdf("shared_clone_T_cells_feature_plot_gene_list0_split_final.pdf",64,length(gene_list_0)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list1_split_final.pdf",64,length(gene_list_1)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list2_split_final.pdf",64,length(gene_list_2)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list3_split_final.pdf",64,length(gene_list_3)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list4_split_final.pdf",64,length(gene_list_4)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list5_split_final.pdf",64,length(gene_list_5)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list6_split_final.pdf",64,length(gene_list_6)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list7_split_final.pdf",64,length(gene_list_7)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list8_split_final.pdf",64,length(gene_list_8)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list9_split_final.pdf",64,length(gene_list_9)*8)
FeaturePlot(object = all_clone_test_filter, split.by="batch",features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list10_split_final.pdf",64,length(gene_list_10)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list11_split_final.pdf",64,length(gene_list_11)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list12_split_final.pdf",64,length(gene_list_12)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list13_split_final.pdf",64,length(gene_list_13)*8)
FeaturePlot(object = all_clone_test_filter, split.by="batch",features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list14_split_final.pdf",64,length(gene_list_14)*8)
FeaturePlot(object = all_clone_test_filter,split.by="batch", features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.5,order=TRUE)
dev.off()

DefaultAssay(all_clone_test_filter) <- "RNA"
###split by label
pdf("shared_clone_T_cells_feature_plot_gene_list0_split_by_tissue_final.pdf",8,length(gene_list_0)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_0,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list1_split_by_tissue_final.pdf",8,length(gene_list_1)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_1,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list2_split_by_tissue_final.pdf",8,length(gene_list_2)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_2,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list3_split_by_tissue_final.pdf",8,length(gene_list_3)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_3,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list4_split_by_tissue_final.pdf",8,length(gene_list_4)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_4,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list5_split_by_tissue_final.pdf",8,length(gene_list_5)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_5,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list6_split_by_tissue_final.pdf",8,length(gene_list_6)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_6,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list7_split_by_tissue_final.pdf",8,length(gene_list_7)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_7,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list8_split_by_tissue_final.pdf",8,length(gene_list_8)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_8,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list9_split_by_tissue_final.pdf",8,length(gene_list_9)*4)
FeaturePlot(object = all_clone_test_filter, split.by="label",features = gene_list_9,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list10_split_by_tissue_final.pdf",8,length(gene_list_10)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_10,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list11_split_by_tissue_final.pdf",8,length(gene_list_11)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_11,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list12_split_by_tissue_final.pdf",8,length(gene_list_12)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_12,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list13_split_by_tissue_final.pdf",8,length(gene_list_13)*4)
FeaturePlot(object = all_clone_test_filter, split.by="label",features = gene_list_13,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()
pdf("shared_clone_T_cells_feature_plot_gene_list14_split_by_tissue_final.pdf",8,length(gene_list_14)*4)
FeaturePlot(object = all_clone_test_filter,split.by="label", features = gene_list_14,cols = c("lightgrey", "firebrick3"),ncol=2,label = FALSE,pt.size=0.2,order=TRUE)
dev.off()

##Run CellChat on merged.all.batchcorrected
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

## ========= 0. 建 CellChat =========
cellchat <- createCellChat(
  object  = merged.all.batchcorrected,
  group.by = "cell_type_annotated"
)

CellChatDB <- CellChatDB.human   # 如果是鼠，用 CellChatDB.mouse
cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))

## 尝试找 CXCL16–CXCR6 这对
db.cxcl16 <- subsetDB(cellchat@DB, search = "CXCL16")
pair.cxcl16.cxcr6 <- db.cxcl16$interaction_name_2[
  grep("CXCR6", db.cxcl16$interaction_name_2, ignore.case = TRUE)
][1]
# 如果没找到，就手动写，比如：
# pair.cxcl16.cxcr6 <- "CXCL16 - CXCR6"

## 尝试找它所属的 pathway
signaling.name <- db.cxcl16$pathway_name[1]
if (length(signaling.name) == 0) signaling.name <- NA

## ========= 1. 总体 overview：circle (count) =========
pdf("01_overview_circle_count.pdf", width = 15, height = 15)
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)
dev.off()

## ========= 2. 总体 overview：circle (weight) =========
pdf("02_overview_circle_weight.pdf", width = 15, height = 15)
netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights"
)
dev.off()
library(ggplot2)



## 0. 先拿到所有有名字的通路
comm.all <- subsetCommunication(cellchat)
all.pw <- sort(unique(comm.all$pathway_name))
all.pw <- all.pw[!is.na(all.pw) & all.pw != ""]


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht.out.all <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern   = "outgoing",
  signaling = all.pw
)
ht.in.all <- netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern   = "incoming",
  signaling = all.pw
)

pdf("03_signalingRole_heatmap.pdf", width = 12, height = 7)
print(ht.out.all + ht.in.all)
dev.off()


## ========= 4. 发送者/接受者角色分析（看谁说得多、收得多） =========
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

pdf("04_signaling_role_scatter.pdf", width = 10, height = 10)
netAnalysis_signalingRole_scatter(cellchat)  # 不再加 type 参数
dev.off()

## ========5.circle plot ==
pdf("05_aggregate_circle.pdf", width = 8, height = 7)
netVisual_aggregate(
  cellchat,
  signaling = all.pw,
  layout   = "circle"
)
dev.off()

## 提取全部显著通信（所有通路）
comm.all <- subsetCommunication(cellchat)

## 查看表结构（可选）
print(head(comm.all))

## 导出到 CSV
write.csv(
  comm.all,
  file = "All_significant_signaling_results.csv",
  row.names = FALSE
)

message("✅ 已导出所有显著的 signaling 通路到 All_significant_signaling_results.csv")

## ============ 0. 所有显著通信 ============
comm.all <- subsetCommunication(cellchat)
comm.all <- comm.all[!is.na(comm.all$pathway_name) & comm.all$pathway_name != "", ]

all.pw <- sort(unique(comm.all$pathway_name))

## ============ 1. 分组：cDC / pDC 单独，其它进大类 ============
all.clusters <- levels(cellchat@idents)

group.cellType <- setNames(rep("OTHERS", length(all.clusters)), all.clusters)

# 单独的
group.cellType[all.clusters %in% "cDC"] <- "cDC"
group.cellType[all.clusters %in% "pDC"] <- "pDC"

# Myeloid 其它
group.cellType[all.clusters %in% c(
  "Classical monocyte",
  "FCGR3A+ monocyte",
  "HLA-DR+ myeloid cell",
  "MRC1+ macrophage",
  "Neutrophil",
  "Mast cell/Basophil"
)] <- "Myeloid"

# T 系
group.cellType[all.clusters %in% c(
  "Naive/Memory-like CD4+",
  "Effector-like CD4+",
  "CXCL13+ effector CD8+",
  "OASL+ effector CD8+",
  "FOXP3+ Treg",
  "γδ T cell"
)] <- "T"

# B 系
group.cellType[all.clusters %in% c(
  "Naive B cell",
  "GC B cell",
  "Memory B cell",
  "Transitional B cell",
  "Plasma cell",
  "Plasmablast"
)] <- "B"

# NK
group.cellType[all.clusters %in% c("NK cell", "XCL1+ NK cell")] <- "NK"

# Erythroid
group.cellType[all.clusters %in% "Erythroid cell"] <- "Erythroid"

# Cycling
group.cellType[all.clusters %in% "MKI67+ cycling cell"] <- "Cycling"

print(table(group.cellType))

# 计算中心性
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# 生成热图对象
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

# 保存成 pdf
pdf("00_signalingRole_heatmaps.pdf", width = 10, height = 8)
print(ht1 + ht2)
dev.off()

## ============ 2. 对每个 pathway 循环 ============
for (pw in all.pw) {
  message(">>> processing pathway: ", pw)
  safe_pw <- gsub("[^A-Za-z0-9_]+", "_", pw)

  # 这个通路的显著通信
  pw.comm <- subsetCommunication(cellchat, signaling = pw)

  ## 2.1 Hierarchy
  tryCatch({
    pdf(paste0(safe_pw, "_01_hierarchy.pdf"), width = 15, height = 15)
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      vertex.receiver = seq(1, 4)
    )
    dev.off()
  }, error = function(e) {
    message("hierarchy failed for ", pw, ": ", e$message)
    if (dev.cur() > 1) dev.off()
  })

  ## 2.2 Circle
  tryCatch({
    pdf(paste0(safe_pw, "_02_circle.pdf"), width = 15, height = 15)
    par(mfrow = c(1, 1))
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      layout = "circle"
    )
    dev.off()
  }, error = function(e) {
    message("circle failed for ", pw, ": ", e$message)
    if (dev.cur() > 1) dev.off()
  })

  ## 2.3 Chord
  tryCatch({
    pdf(paste0(safe_pw, "_03_chord.pdf"), width = 15, height = 15)
    par(mfrow = c(1, 1))
    netVisual_aggregate(
      cellchat,
      signaling = pw,
      layout = "chord"
    )
    dev.off()
  }, error = function(e) {
    message("chord failed for ", pw, ": ", e$message)
    if (dev.cur() > 1) dev.off()
  })

  ## 2.4 Heatmap —— 你说能跑的版本，这里加保护
  tryCatch({
    pdf(paste0(safe_pw, "_04_heatmap.pdf"), width = 15, height = 15)
    par(mfrow = c(1, 1))
    print(netVisual_heatmap(
      cellchat,
      signaling = pw,
      color.heatmap = "Reds"
    ))
    dev.off()
  }, error = function(e) {
    message("heatmap failed for ", pw, ": ", e$message)
    if (dev.cur() > 1) dev.off()
  })

  ## 2.5 分组 chord（cDC / pDC 单独）
  tryCatch({
    pdf(paste0(safe_pw, "_05_chord_grouped.pdf"), width = 15, height = 15)
    netVisual_chord_cell(
      cellchat,
      signaling = pw,
      group = group.cellType,
      title.name = paste0(pw, " signaling network")
    )
    dev.off()
  }, error = function(e) {
    message("grouped chord failed for ", pw, ": ", e$message)
    if (dev.cur() > 1) dev.off()
  })

  ## 2.6 contribution —— 先官方，失败就自己算
  contrib.file <- paste0(safe_pw, "_06_contribution.pdf")
  contrib.data <- NULL
  official.ok <- TRUE

  tryCatch({
    pdf(contrib.file, width = 7, height = 4)
    contrib.data <- print(netAnalysis_contribution(cellchat, signaling = pw))
    dev.off()
  }, error = function(e) {
    official.ok <<- FALSE
    message("netAnalysis_contribution failed for ", pw, ": ", e$message)
    if (dev.cur() > 1) dev.off()
  })

  if (!official.ok) {
    if (!is.null(pw.comm) && nrow(pw.comm) > 0) {
      contrib.data <- aggregate(
        prob ~ interaction_name_2,
        data = pw.comm,
        FUN = sum
      )
      contrib.data <- contrib.data[order(-contrib.data$prob), ]

      write.csv(
        contrib.data,
        file = paste0(safe_pw, "_contribution_table.csv"),
        row.names = FALSE
      )

      pdf(contrib.file, width = 7, height = 4)
      barplot(
        contrib.data$prob,
        names.arg = contrib.data$interaction_name_2,
        las = 2,
        cex.names = 0.6,
        main = paste0(pw, " contribution (sum prob)"),
        ylab = "sum(prob)"
      )
      dev.off()
    }
  } else {
    if (!is.null(contrib.data)) {
      try(
        write.csv(
          contrib.data,
          file = paste0(safe_pw, "_contribution_table.csv"),
          row.names = FALSE
        ),
        silent = TRUE
      )
    }
  }

  ## 2.7 single LR（如果能提取到）
  pairLR.pw <- tryCatch(
    extractEnrichedLR(cellchat, signaling = pw, geneLR.return = FALSE),
    error = function(e) NULL
  )

  if (!is.null(pairLR.pw) && nrow(pairLR.pw) > 0) {
    LR.show <- pairLR.pw[1, ]

    # hierarchy
    tryCatch({
      pdf(paste0(safe_pw, "_07_singleLR_hierarchy.pdf"), width = 15, height = 15)
      netVisual_individual(
        cellchat,
        signaling = pw,
        pairLR.use = LR.show,
        vertex.receiver = seq(1, 4)
      )
      dev.off()
    }, error = function(e) {
      message("singleLR hierarchy failed for ", pw, ": ", e$message)
      if (dev.cur() > 1) dev.off()
    })

    # circle
    tryCatch({
      pdf(paste0(safe_pw, "_08_singleLR_circle.pdf"), width = 15, height = 15)
      netVisual_individual(
        cellchat,
        signaling = pw,
        pairLR.use = LR.show,
        layout = "circle"
      )
      dev.off()
    }, error = function(e) {
      message("singleLR circle failed for ", pw, ": ", e$message)
      if (dev.cur() > 1) dev.off()
    })
  } else {
    message("(", pw, ") no enriched LR, skip singleLR plots.")
  }

  ## 2.8 导出这个 pathway 的 communication
  if (!is.null(pw.comm) && nrow(pw.comm) > 0) {
    write.csv(
      pw.comm,
      file = paste0(safe_pw, "_communication.csv"),
      row.names = FALSE
    )
  }
	## 10. Bubble plot —— 通路下所有显著通信 (source × target)
	if (!is.null(pw.comm) && nrow(pw.comm) > 0) {
	  pdf(paste0(safe_pw, "_09_bubble.pdf"), width = 12, height = 12)
	  p <- ggplot(
		pw.comm,
		aes(x = target, y = source, size = prob, fill = prob)
	  ) +
		geom_point(shape = 21, color = "black") +
		scale_size(range = c(2, 10)) +
		labs(
		  title = paste0(pw, " pathway communications"),
		  x = "Target cell",
		  y = "Source cell"
		) +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1))
	  print(p)
	  dev.off()
	}
	## 11. incoming/outgoing 
	## ========= signaling role heatmaps for this pathway =========
	tryCatch({
	  ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = pw)
	  ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = pw)
	  
	  pdf(paste0(safe_pw, "_10_signalingRole_heatmap.pdf"), width = 14, height = 9)
	  print(ht1 + ht2)
	  dev.off()
	}, error = function(e) {
	  message("signalingRole heatmap failed for ", pw, ": ", e$message)
	  if (dev.cur() > 1) dev.off()
	})

}

## 1. 拿所有显著的通信
comm.all <- subsetCommunication(cellchat)

## 2. 只保留有通路名的
comm.all <- comm.all[!is.na(comm.all$pathway_name) & comm.all$pathway_name != "", ]

## 3. prob 强制转数值（关键！）
comm.all$prob <- as.numeric(comm.all$prob)

## 4. 按通路做一个 list，后面好算
pw.list <- split(comm.all, comm.all$pathway_name)

## 5. 对每个通路做统计
pw.summary.list <- lapply(names(pw.list), function(pw) {
  df <- pw.list[[pw]]

  # 通路里的所有基因 = 所有 ligand + 所有 receptor 去重
  genes <- unique(c(df$ligand, df$receptor))

  # 这条通路的概率向量
  pr <- df$prob

  data.frame(
    pathway_name = pw,
    mean_prob = mean(pr, na.rm = TRUE),
    sum_prob  = sum(pr,  na.rm = TRUE),
    n_pairs   = length(pr),
    genes     = paste(genes, collapse = ", "),
    stringsAsFactors = FALSE
  )
})

## 6. 绑回一个大表
pw.summary.full <- do.call(rbind, pw.summary.list)

## 7. 导出
write.csv(
  pw.summary.full,
  file = "All_significant_pathways_summary.csv",
  row.names = FALSE
)

saveRDS(cellchat, file = "cellchat_merged_all_batchcorrected.rds")

###volcano plot update 2/8/2026
library(ggplot2)
library(ggrepel)
library(dplyr)

# 1. 读取并修正数据对齐
data <- read.table("DEG.txt", sep = "\t", skip = 1, header = FALSE)
colnames(data) <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")

# 2. 设定阈值与基础处理
lfc_thresh <- 0.5
p_thresh <- 0.05
data <- data %>%
  mutate(
    logP = -log10(p_val_adj + 1e-300),
    group = case_when(
      p_val_adj > p_thresh ~ "NS",
      avg_log2FC >= lfc_thresh ~ "Up",
      avg_log2FC <= -lfc_thresh ~ "Down",
      TRUE ~ "NS"
    )
  )

# 3. 筛选目标基因
target_genes <- c("KLF2", "JUN", "CXCR6", "IFNG", "GZMB", "CISH", "TCF7")
highlight_df <- data %>% filter(gene %in% target_genes)

# 4. 绘图
ggplot(data, aes(x = avg_log2FC, y = logP)) +
  # 背景点
  geom_point(data = filter(data, group == "NS"), color = "grey85", alpha = 0.4, size = 1) +
  geom_point(data = filter(data, group == "Up"), color = "#E41A1C", alpha = 0.4, size = 1) +
  geom_point(data = filter(data, group == "Down"), color = "#377EB8", alpha = 0.4, size = 1) +
  
  # 高亮大点
  geom_point(data = highlight_df, aes(fill = group), shape = 21, size = 3, color = "black", stroke = 0.8) +
  
  # 极短指引线标注 (nudge_x = 0.25)
  geom_text_repel(data = highlight_df,
                  aes(label = gene),
                  size = 5, fontface = "italic",
                  nudge_x = ifelse(highlight_df$avg_log2FC > 0, 0.25, -0.25), 
                  direction = "y", 
                  hjust = ifelse(highlight_df$avg_log2FC > 0, 0, 1),
                  segment.size = 0.4,
                  segment.color = "black",
                  min.segment.length = 0,
                  max.overlaps = Inf) +
  
  # 辅助线
  geom_vline(xintercept = c(-lfc_thresh, 0, lfc_thresh), 
             linetype = c("dashed", "solid", "dashed"), color = c("grey70", "black", "grey70"), linewidth = 0.4) +
  geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "grey70", linewidth = 0.4) +
  
  # 设置 X 轴对称范围：-3 到 3
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  
  # 样式美化
  scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
  theme_classic() +
  labs(title = "Graft vs PBMC (shared-TCR clones)", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P value") +
  theme(legend.position = "none", 
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("Volcano_Plot.pdf", width = 7, height = 6)


###
library(Seurat)

  setwd(file.path(scvdj_work_root, "merge_4_with_vdj"))

  load("merged_4_with_vdj_T_cells_filter_10_28_2025.RData")

  cd8_clusters <- c(
    "CXCL13+CXCR6+ effector CD8+",
    "CXCR6+ effector CD8+",
    "Naive/Memory-like CD8+",
    "CX3CR1+ effector CD8+",
    "ZNF683+TCF7hi CD8+",
    "ISG-high CD8+"
  )

  merged_all_T_batchcorrected_filter_CD8 <- subset(
    merged_all_T_batchcorrected_filter,
    subset = cell_type_annotated_granularity %in% cd8_clusters
  )

  table(merged_all_T_batchcorrected_filter_CD8$cell_type_annotated_granularity)
  table(merged_all_T_batchcorrected_filter_CD8$batch)

  save(
    merged_all_T_batchcorrected_filter_CD8,
    file = "merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData"
  )
