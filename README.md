# LncPairs
LncPairs: lncRNA-mRNA pairs identification at single-cell level

uncompress LncPairs.tar.gz and move LncPairs directory to R package library (e.g. C:\Users\Administrator\AppData\Local\R\win-library\4.2)

### Usage
```{r}
library(Seurat)
libary(LncPairs)

data('lncRNA')
data('mRNA')

pbmc_small <- NormalizeData(pbmc_small, normalization.method = 'LogNormalize', scale.factor = 10000)
pbmc_small <- FindVariableFeatures(pbmc_small, selection.method = 'vst', nfeatures = 2000)
pbmc_small <- ScaleData(pbmc_small, features = VariableFeatures(pbmc_small),do.center = T)
pbmc_small <- RunPCA(pbmc_small, features = VariableFeatures(pbmc_small))

pct <- pbmc_small[['pca']]@stdev / sum(pbmc_small[['pca']]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
print(paste0('The best of PC: ', pcs))
pbmc_small <- FindNeighbors(pbmc_small, dims = 1:pcs, reduction = 'pca')
pbmc_small <- RunUMAP(pbmc_small, dims = 1:pcs, reduction = 'pca')
pbmc_small <- FindClusters(pbmc_small, resolution = 1, algorithm = 1)


LncPairs(
    seurat_obj = pbmc_small,
    groupBy = 'seurat_clusters',
    sampleName = 'pbmc',
    top_genes = 50
)
```
