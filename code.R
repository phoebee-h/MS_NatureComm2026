## R version 4.1.0
library(Seurat)  #v4.1.1
library(SingleR)  #v1.8.1
library(SeuratWrappers) #v0.3.0
library(clusterProfiler) #v4.2.2
library(CellChat) #v1.1.3

h5file <- Read10X_h5("filtered_feature_bc_matrix.h5")
seurat_obj <- CreateSeuratObject(h5file, min.cells=1, min.features=1)

##---------- [Quality Control] ----------##
# Subset "each sample" and do the filtering QC steps
seurat_obj <- subset(seurat_obj, idents = "certain_sample")

# Mitochondrial gene percentage count
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
nFeature_Q3 <- as.numeric(quantile(seurat_obj$nFeature_RNA, 0.75))
nFeature_Q1 <- as.numeric(quantile(seurat_obj$nFeature_RNA, 0.25))
IQR <- nFeature_Q3 - nFeature_Q1

# Remove cells by filtering threshold
upthred <- round(nFeature_Q3 + 1.5 * IQR)
seurat_obj <- seurat_obj[, which(seurat_obj$nFeature_RNA < upthred)]
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 20)

# Merge all samples after quality control
seurat_obj_merge <- merge(obj1, y = c(obj2, obj3, ...), 
	add.cell.ids = c("obj1", "obj2", "obj3", ...), 
	project = "irAEs_SJS")
	
##---------- [Seurat pipeline] ----------##
# Seurat pipeline, Normalization
seurat_obj_merge <- NormalizeData(seurat_obj_merge) %>% FindVariableFeatures() 
notzerogenes <- which(Matrix::rowSums(seurat_obj_merge[["RNA"]]@counts)!=0)
	
vargene <- VariableFeatures(seurat_obj_merge)
seurat_obj_mnn <- seurat_obj_merge %>%
	subset(features = notzerogenes) %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
vargene <- seurat_obj_mnn[["RNA"]]@var.features

## Batch effect correction
seurat_obj_mnn <- RunFastMNN(object.list = SplitObject(seurat_obj_mnn, split.by = "batch")) 
seurat_obj_mnn <- ScaleData(seurat_obj_mnn) %>%
    RunPCA(features=vargene) %>%
    RunTSNE(reduction="mnn", dims = 1:20, perplexity = 30, seed.use = 100) %>%
    RunUMAP(reduction="mnn", dims = 1:20, n.components = 2, seed.use = 100) %>%
    FindNeighbors(reduction="mnn", dims = 1:20) %>%
    FindClusters(resolution = 0.4)


##---------- [SingleR] ----------##
obj_SCE <- as.SingleCellExperiment(seurat_obj_mnn)
db <- loadRData("celldexdb.BE.Rdata")
obj_SCE.main <- SingleR(test = obj_SCE, ref = db, assay.type.test=1, labels = db$label.main)
seurat_obj_mnn[["BE.main"]] <- obj_SCE.main$labels


##---------- [DEGs/ GSEA] ----------##
# group1 and group2 are defined by different condition
data.GMT_H <- clusterProfiler::read.gmt("h.all.v2022.1.Hs.symbols.gmt")
DEGlist <- FindMarkers(seurat_obj_mnn, only.pos = FALSE, 
		       min.pct = 0.25, logfc.threshold = 0, 
		       ident.1 = "group1",
		       ident.2 = "group2")

genelist <- DEGlist$avg_log2FC
names(genelist) <- DEGlist$symbol
genelist <- sort(genelist, decreasing=TRUE) 
gsea.gmt <- GSEA(genelist, 
		 maxGSSize = 500, 
		 minGSSize = 10,
		 nPerm = 1000,
		 TERM2GENE = data.GMT, 
		 pvalueCutoff = 1, 
		 by="fgsea", seed = TRUE)
				 
				 
##---------- [CellChat] ----------##
# See the Cellchat official website for more details
# We used the older version of CellChat v1.1.3 back then, some functions might be changed.
# https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html
