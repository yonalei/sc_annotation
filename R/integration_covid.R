setwd("/gpfs/gibbs/pi/csei/users/yl2485/")
getwd()

devtools::install_github(repo = "satijalab/seurat", ref = "develop")
library(Seurat)

covid_flu <- readRDS(file = "data/COVID_Flu_CITE_Seurat_Object.rds")
acute_covid <- readRDS("/gpfs/gibbs/pi/csei/Projects/CFLU/metadata/AcuteCovid_GSE161918_AllBatches_SeuratObj.rds")
acute_covid = UpdateSeuratObject(object = acute_covid)

library(ggplot2)
library(cowplot)

## standard flow
acute_covid <- ScaleData(acute_covid, verbose = FALSE)
acute_covid <- RunPCA(acute_covid, npcs = 30, verbose = FALSE)
acute_covid <- RunUMAP(acute_covid, reduction = "pca", dims = 1:30)
# p3 <- DimPlot(acute_covid, reduction = "umap", group.by = "WCTcoursecelltype", raster = FALSE)
p4 <- DimPlot(dataset.list[["AcuteCovid2"]], reduction = "umap", group.by = "WCTcoursecelltype", label = TRUE,
               repel = TRUE, raster = FALSE) + NoLegend()
p5 <- DimPlot(dataset.list[["AcuteCovid2"]], reduction = "umap", group.by = "predicted.id", label = TRUE,
              repel = TRUE, raster = FALSE) + NoLegend()

plot_grid(p4, p5)

p1 <- DimPlot(acute_covid, reduction = "umap", group.by = "mergedcelltype", raster = FALSE)
p2 <- DimPlot(acute_covid, reduction = "umap", group.by = "mergedcelltype", label = TRUE,
              repel = TRUE, raster = FALSE) + NoLegend()
plot_grid(p1, p2)

# covid.query <- list(covid_flu)
# covid.reference <- list(acute_covid)
# covid.anchors <- FindTransferAnchors(reference = covid.reference, query = covid.query, 
#                                         dims = 1:30)
# predictions <- TransferData(anchorset = covid.anchors, refdata = covid.integrated$celltype, 
#                             dims = 1:30)
# covid.query <- AddMetaData(covid.query, metadata = predictions)
# 


# covid_flu <- FindVariableFeatures(covid_flu, selection.method = "vst", 
#                                       nfeatures = 2000, verbose = FALSE)
# 

dataset.list<-list()
dataset.list[["AcuteCovid1"]]<-acute_covid #reference
dataset.list[["AcuteCovid2"]]<-acute_covid #query
#dataset.anchors <- FindIntegrationAnchors(object.list = dataset.list, dims = 1:30)


options(future.globals.maxSize = 1000 * 1024^2)
covid.anchors <- FindTransferAnchors(reference = dataset.list[["AcuteCovid1"]], query = dataset.list[["AcuteCovid2"]], 
                                        dims = 1:30)


predictions <- TransferData(anchorset = covid.anchors, refdata = dataset.list[["AcuteCovid1"]]$WCTcoursecelltype, 
                            dims = 1:30)
# Warning messages:
# 1: ggrepel: 10 unlabeled data points (too many overlaps). Consider increasing max.overlaps 
# 2: UNRELIABLE VALUE: One of the ‘future.apply’ iterations (‘future_lapply-1’) unexpectedly generated random numbers without declaring so. 
# There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'future.seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced via the L'Ecuyer-CMRG method. To disable this check, use 'future.seed = NULL', or set option 'future.rng.onMisuse' to "ignore". 

dataset.list[["AcuteCovid2"]] <- AddMetaData(dataset.list[["AcuteCovid2"]], metadata = predictions)

dataset.list[["AcuteCovid2"]]$prediction.match <- dataset.list[["AcuteCovid2"]]$predicted.id == dataset.list[["AcuteCovid2"]]$WCTcoursecelltype


table(dataset.list[["AcuteCovid2"]]$prediction.match)
 # output
# FALSE   TRUE 
# 47307 364595 


