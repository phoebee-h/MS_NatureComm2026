library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sctransform)

## dataset processing
skin1 <- Load10X_Spatial(data.dir = "path of skin1", slice = "skin1")

skin1 <- AddMetaData(object = skin1, metadata = "skin1", col.name = "sampleID")
skin1 <- AddMetaData(object = skin1, metadata = "group_name", col.name = "Group")

skin1 <- RenameCells(object = skin1, add.cell.id = "skin1")

skin1 <- SCTransform(skin1, assay = "Spatial", vst.flavor = "v2")

## working with multiple slices in Seurat
skin_merge <- merge(x = skin1, y = c(skin2, skin3, ...))

skin_merge <- SCTransform(skin_merge, assay = "Spatial", vst.flavor = "v2")

skin_merge <- RunPCA(skin_merge, assay = "SCT", verbose = FALSE)
skin_merge <- FindNeighbors(skin_merge, reduction = "pca", dims = 1:30)
skin_merge <- FindClusters(skin_merge, verbose = FALSE)
skin_merge <- RunUMAP(skin_merge, reduction = "pca", dims = 1:30)

## integrated layer 
# https://satijalab.org/seurat/articles/integration_rpca
skin_merge <- IntegrateLayers(object = skin_merge, method = RPCAIntegration, 
                          normalization.method = "SCT")

skin_merge <- FindNeighbors(skin_merge, dims = 1:30, reduction = "integrated.dr")
skin_merge <- FindClusters(skin_merge, resolution = 0.4)
skin_merge <- RunUMAP(skin_merge, dims = 1:30, reduction = "integrated.dr")

skin_merge <- PrepSCTFindMarkers(skin_merge, assay = "SCT")
Allmarkers <- FindAllMarkers(skin_merge, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = TRUE)


##---------- [DEGs/ GSEA] ----------##
DEG <- FindMarkers(skin_merge, min.pct = 0.1, logfc.threshold = 0,
                   ident.1 = "group1",
                   ident.2 = "group2")


library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)

genelist <- DEG[order(DEG$avg_log2FC, decreasing = T),]
genelist$gene <- rownames(genelist)

genelist_list <- as.numeric(genelist$avg_log2FC)
names(genelist_list) <- genelist$gene

# MSigDB (c2cp)
C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2") %>%
  dplyr::select(gs_name, gene_symbol)

region1_gsea <- GSEA(genelist_list, TERM2GENE = C2_t2g, pvalueCutoff = 1)


##---------- [Spatial plots/ UMAP plots/ dot plots] ----------##
SpatialDimPlot(skin_merge)
SpatialDimPlot(skin_merge, cells.highlight = WhichCells(skin_merge, idents = "region"), facet.highlight = TRUE)

SpatialFeaturePlot(skin_merge, features = "genes")

DimPlot(skin_merge, reduction = "umap")

DotPlot(skin_merge,
        features = c("genes"),
        cols = "RdYlBu") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_y_discrete(limits = rev)


##---------- [RCTD (spacexr)] ----------##
# See the official spacexr website for more details:
# https://raw.githack.com/dmcable/spacexr/master/vignettes/visium_full_regions.html
# https://satijalab.org/seurat/articles/spatial_vignette
# RCTD in full mode was used for the analysis
#
# The normalized RCTD weight data were then visualized using scatterpie. 
# https://cran.r-project.org/web/packages/scatterpie/vignettes/scatterpie.html

