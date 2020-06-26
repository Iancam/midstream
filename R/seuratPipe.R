library(purrr)
library(Seurat)
library(patchwork)
library(zeallot)
library(magrittr)

# source("./plotting.r")
# source("./fileManip.r")

tenX2PCA <- function(
        path,
        cache_dir = "pcaCache",
        percent.mt.thresh = 5,
        minNFeature=200,
        maxNFeature=2500,
        numTop =10
) {
    print(path)

    name <- getFname(path)
    c(plots, save_plot) %<-% plotsFor(name)
    # Initialize the Seurat object with the raw (non-normalized data).
    if(path %>% endsWith("rds")) dataset <- readRDS(path)
    else dataset <- CreateSeuratObject(get10x(path))
    print("percent.mt")
    percent.mt <- PercentageFeatureSet(dataset, pattern = "^[Mm][tT]-")
    AddMetaData(object=dataset, metadata = percent.mt, col.name = "percent.mt")
    dataset[["percent.mt"]] <- percent.mt
    save_plot(function() VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),"beforeThresh.mito.vln")
    plot1 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(function() plot1 + plot2, "mito.scatter")

    if (!is.null(minNFeature)) dataset <- subset(dataset, subset = nFeature_RNA > minNFeature)
    if (!is.null(maxNFeature)) dataset <- subset(dataset, subset = nFeature_RNA < maxNFeature)
    if (!is.null(percent.mt))  dataset <- subset(dataset, subset = percent.mt < percent.mt.thresh)
    save_plot(function() VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),"afterThresh.mito.vln")
    print("Normalize")
    dataset <- NormalizeData(dataset)
    
    print("variable features")
    dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 2000)
    top <- head(VariableFeatures(dataset), numTop)
    plot1 <- VariableFeaturePlot(dataset)
    plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE)
    save_plot(function() plot2, "expressionVariance")
    
    print("scaling data")
    all.genes <- rownames(dataset)
    dataset <- ScaleData(dataset)
    save_plot(function() VizDimLoadings(dataset, dims = 1:2, reduction = "pca"),
     "dim.loadings")
    
    print("PCA")
    dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset))
    save_plot(function() DimPlot(dataset, reduction = "pca"), "dim.pca")
    
    # save_plot(function() DimHeatmap(dataset, dims = 1:15, cells = 500, balanced = TRUE),
    #  "dim.heatmap")

    dir.create(cache_dir, showWarnings= FALSE)
    op_path <- paste(cache_dir, paste0(name, ".rds"), sep = "/")
    saveRDS(dataset, file = op_path)
    list(plots, dataset, name, op_path)
}

getdimsInfo = function(dir_info) getDims(dir_info = dir_info)
getDimsCached = function(dir) getDims(dir_path = dir)

getDims <- function(dir_info = NULL, dir_path = NULL) {
    if (!is.null(dir_path)) {
        name = getFname(dir_path)
        print(dir_path)
        dataset <- readRDS(dir_path)
    } else {
        c(dataset, name, dir_path) %<-% dir_info
    }
    c(plots, save_plot) %<-% plotsFor(name)
    elbow <- ElbowPlot(dataset)
    print(elbow)
    save_plot(function() elbow, "elbow")()

    dims <- as.numeric(readline(prompt = "what shall be your dimensions? "))
    c(dataset, name, dims)
}

PCA2Cluster <- function(
    info,
    output_dir = "seuratOutput",
    clusterResolution= 0.5,
    reductionTypes= c("umap"),
    min.pct = 0.25,
    logfc.threshold = 0.25
) {
    print("toCluster")
    print(reductionTypes)
    c(dataset, name, dims) %<-% info
    c(plots, save_plot) %<-% plotsFor(name)
    markerFN <- paste(output_dir, paste0(name, ".markers.csv"), sep = "/")
    rdsFN <- paste(output_dir, paste0(name, ".rds"), sep = "/")
    dataset <- FindNeighbors(dataset, dims = 1:dims)
    dataset <- FindClusters(dataset, resolution = clusterResolution)
    
    reductionMapping = c("umap" = RunUMAP, "tsne" = RunTSNE)
    dataset <- reduce(reductionTypes, function(dataset,type) {
      dataset <- reductionMapping[[type]](dataset, dims = 1:dims)
      save_plot(function() DimPlot(dataset, reduction = type), type)
      dataset
    }, .init=dataset)
    # if("tsne" %in% reductionTypes) {
    #     dataset <- reductionMapping[["tsne"]](dataset, dims = 1:dims)
    #     save_plot(function() DimPlot(dataset, reduction = "tsne"), "tsne") 
    # }
    # if("umap" %in% reductionTypes) {
    #     dataset <- reductionMapping[["umap"]](dataset, dims = 1:dims)
    #     save_plot(function() DimPlot(dataset, reduction = "umap"), "umap") 
    # }
    print("markers")
    dataset.markers <- FindAllMarkers(
        dataset,
        only.pos = TRUE,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold
    )
    print(markerFN)
    list(plots, dataset.markers, markerFN, dataset, rdsFN, name)
}