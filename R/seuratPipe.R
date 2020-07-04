#' @import Seurat
QC = function(
    dataset,
    percent.mt.thresh = 5,
    minNFeature=200,
    maxNFeature=2500,
    numTop =10,
    nfeatures = 2000,
    plotSaver = NULL
) {
    name = getNameFromSeurat(dataset)
    plotStuff = plotSaver(name)
    plots = plotStuff[[1]]
    save_plot = plotStuff[[2]]
    analysisFile = plotStuff[[3]]
    # Initialize the Seurat object with the raw (non-normalized data).
    
    print("percent.mt")
    percent.mt <- PercentageFeatureSet(dataset, pattern = "^[Mm][tT]-")
    AddMetaData(object=dataset, metadata = percent.mt, col.name = "percent.mt")
    dataset[["percent.mt"]] <- percent.mt

    save_plot(function() VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),"before.mito.vln")
    plot1 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(function() plot1, "mito-rnaCount.scatter")
    save_plot(function() plot2, "rnaFeature-rnaCount.scatter")

    sink(analysisFile)
    print(Sys.time())
    print(dataset)
    if (!is.null(minNFeature)) {
        dataset <- reportChanges(
            dataset,
            function() subset(dataset, subset = nFeature_RNA > minNFeature),
            "subsetMinNFeature",
            c(minNFeature)
        )
    }
    if (!is.null(maxNFeature)) {
        dataset <- reportChanges(
            dataset,
            function() subset(dataset, subset = nFeature_RNA < maxNFeature),
            "subsetMaxNFeature",
            c(maxNFeature)
        )
    }
    if (!is.null(percent.mt)) {
        dataset <- reportChanges(
            dataset,
            function() subset(dataset, subset = percent.mt < percent.mt.thresh),
            "subsetPercent.mt",
            c(percent.mt.thresh)
        )
    }
    sink()
    save_plot(function() VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),"afterThresh.vln")
    dataset
}

toPCA <- function(
        dataset,
        cache_dir = "pcaCache",
        plotSaver = NULL,
        percent.mt.thresh = 5,
        minNFeature=200,
        maxNFeature=2500,
        numTop =10,
        nfeatures = 2000,
        plotHeatMap = F
) {
    name <- getNameFromSeurat(dataset)
    plotStuff = plotSaver(name)
    plots = plotStuff[[1]]
    save_plot = plotStuff[[2]]
    analysisFile = plotStuff[[3]]
    
    print("Normalize")
    dataset <- NormalizeData(dataset, verbose = F)
    
    print("variable features")
    sink(analysisFile, append = T)
    dataset <- reportChanges(dataset,
        function() FindVariableFeatures(dataset, selection.method = "vst", nfeatures = nfeatures),
        "findVariableFeatures",
        c(selection.method = "vst", nfeatures = nfeatures)
    )
    sink()
    top <- head(VariableFeatures(dataset), numTop)
    plot2 <- LabelPoints(plot = VariableFeaturePlot(dataset), points = top, repel = TRUE)
    save_plot(function() plot2, "expressionVariance")
    
    
    print("scaling data")
    all.genes <- rownames(dataset)
    dataset <- ScaleData(dataset)
    save_plot(function() VizDimLoadings(dataset, dims = 1:2, reduction = "pca"),
     "dim.loadings")
    
    print("PCA")
    dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset, nfeatures = nfeatures))
    save_plot(function() DimPlot(dataset, reduction = "pca"), "dim.pca")
    
    if(plotHeatMap) save_plot(function() DimHeatmap(dataset, dims = 1:15, cells = 500, balanced = TRUE),
     "dim.heatmap")

    dir.create(cache_dir, showWarnings= FALSE)
    op_path <- paste(cache_dir, paste0(name, ".rds"), sep = "/")
    saveRDS(dataset, file = op_path)
    list(plots, dataset, name, op_path)
}

getDims <- function(dir_info = NULL, dir_path = NULL, plotSaver = NULL, dims=20, getManually= F) {
    print("getting dims")
    if (!is.null(dir_path)) {
        name = getFname(dir_path)
        print(dir_path)
        dataset <- readRDS(dir_path)
    } else {
        dataset = dir_info[[1]]
        name = dir_info[[2]]
        dir_path = dir_info[[3]]
    }
    plotStuff = plotSaver(name)
    plots = plotStuff[[1]]
    save_plot =plotStuff[[2]]
    analysisFile = plotStuff[[3]]
    elbow <- ElbowPlot(dataset)
    # print(elbow)
    save_plot(function() elbow, "elbow")
    sink(analysisFile, append = T)
    if(getManually) dims <- as.numeric(readline(prompt = "what shall be your dimensions? "))
    cat("dims: ")
    print(dims)
    sink()
    c(dataset, name, dims)
}

toCluster <- function(
    info,
    plotSaver,
    output_dir = "seuratOutput",
    clusterResolution= 0.5,
    reductionTypes= c("umap"),
    min.pct = 0.25,
    logfc.threshold = 0.25
) {
    dataset = info[[1]]
    name = info[[2]]
    dims = info[[3]]
    plotStuff = plotSaver(name)
    plots = plotStuff[[1]]
    save_plot = plotStuff[[2]]
    analysisFile = plotStuff[[3]]
    print("toCluster")

    markerFN <- paste(output_dir, paste0(name, ".markers.csv"), sep = "/")
    rdsFN <- paste(output_dir, paste0(name, ".rds"), sep = "/")
    dataset <- FindNeighbors(dataset, dims = 1:dims)
    dataset <- FindClusters(dataset, resolution = clusterResolution)
    
    sink(analysisFile, append = T)
    print("added reductions: ")
    print(reductionTypes)
    print(table(Idents(dataset)))
    sink()

    reductionMapping = c("umap" = RunUMAP, "tsne" = RunTSNE)
    lapply(reductionTypes, function(type) {
        dataset <<- reductionMapping[[type]](dataset, dims = 1:dims)
        save_plot(function() DimPlot(dataset, reduction = type, label = T), type)
    })
    print(dataset)
    
    print("markers")
    dataset.markers <- FindAllMarkers(
        dataset,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold
    )
    print(markerFN)
    
    list(plots, dataset.markers, markerFN, dataset, rdsFN, name)
}

tenX2Combined <- function(
    input_dir = "input",
    comb_delim = "_",
    sub_delim="-",
    min.cells= 3,
    min.features= 200,
    force = F) {
    if(file.exists(paste0(input_dir, "/comb.rds"))) {
        print("comb exists, skipping")
        return()
    }
    dirs = list.dirs(input_dir, recursive = FALSE)
    print(input_dir)
    print(dirs)
    print("toCombined")
    named_data <- lapply(dirs, function(dir_path) {
        dir_name <- sapply(stringr::str_split(dir_path, "/"), tail, 1)
        data <- get10x(dir_path)
        id <- gsub(comb_delim, sub_delim, dir_name)
        colnames(data) <- paste0(id, "_", colnames(data))
        data
    })
    data <- CreateSeuratObject(do.call(cbind, named_data),
        min.cells =min.cells,
        min.features =min.features,
        project = paste0(getExperimentName(paste0(input_dir, "/comb.rds")), "/comb")
    )
    saveRDS(data, file = paste0(input_dir, "/comb.rds"))
}