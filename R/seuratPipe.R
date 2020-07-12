#' @import Seurat
#' @import dplyr
#' @importFrom utils head tail
#' 
QC = function(
    props,
    percent.mt.thresh = 5,
    minNFeature=200,
    maxNFeature=2500,
    numTop =10
) {
    dataset = props$dataset
    report = props$report
    name = getNameFromSeurat(dataset)
    # Initialize the Seurat object with the raw (non-normalized data).
    report("startTime", Sys.time())
    print("percent.mt")
    percent.mt <- PercentageFeatureSet(dataset, pattern = "^[Mm][tT]-")
    AddMetaData(object=dataset, metadata = percent.mt, col.name = "percent.mt")
    dataset[["percent.mt"]] <- percent.mt
    report("plot::before.mito.vln", g(VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)))
    report("plot::mito.rnaCount.scatter", g(FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "percent.mt")))
    report("plot::rnaFeature.rnaCount.scatter", g(FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")))

    if (!is.null(minNFeature)) {
        dataset = reportChanges(
            dataset,
            function() subset(dataset, subset = nFeature_RNA > minNFeature),
            "subsetMinNFeature",
            c(minNFeature),
            report
        )
    }
    if (!is.null(maxNFeature)) {
        dataset <- reportChanges(
            dataset,
            function() subset(dataset, subset = nFeature_RNA < maxNFeature),
            "subsetMaxNFeature",
            c(maxNFeature),
            report
        )
    }
    if (!is.null(percent.mt)) {
        dataset <- reportChanges(
            dataset,
            function() subset(dataset, subset = percent.mt < percent.mt.thresh),
            "subsetPercent.mt",
            c(percent.mt.thresh),
            report
        )
    }
    report("plot::afterThresh.vln", g(VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)))
    list("dataset" = dataset,  "report" = report)
}

toPCA <- function(
        input,
        numTop =10,
        nfeatures = 2000,
        plotHeatMap = F
) {
    dataset = input$dataset
    report  = input$report
    name <- getNameFromSeurat(dataset)
    
    print("Normalize")
    dataset <- NormalizeData(dataset, verbose = F)
    
    print("variable features")
    dataset <- reportChanges(
        dataset,
        function() FindVariableFeatures(dataset, selection.method = "vst", nfeatures = nfeatures),
        "findVariableFeatures",
        c(selection.method = "vst", nfeatures = nfeatures),
        report
    )
    top <- head(VariableFeatures(dataset), numTop)
    plot2 <- LabelPoints(plot = VariableFeaturePlot(dataset), points = top, repel = TRUE)
    report("plot::expressionVariance", g(plot2))
    
    print("scaling data")
    all.genes <- rownames(dataset)
    dataset <- ScaleData(dataset)
    
    print("PCA")
    dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset, nfeatures = nfeatures))
    report("plot::dim.pca", g(DimPlot(dataset, reduction = "pca")))
    report("plot::dim.loadings", g(VizDimLoadings(dataset, dims = 1:2, reduction = "pca")))
    if(plotHeatMap) report("dim.heatmap", g(DimHeatmap(dataset, dims = 1:15, cells = 500, balanced = TRUE)))

    list("dataset" = dataset,  "report" = report)
}

getDims <- function(
    input = NULL,
    dir_path = NULL,
    dims=20,
    getManually= F
) {
    print("getting dims")
    dataset = input$dataset
    report = input$report
    elbow <- ElbowPlot(dataset)
    report("plot::elbow", g(elbow))
    if(getManually) {
        print(elbow)
        dims <- as.numeric(readline(prompt = "what shall be your dimensions? "))
    }
    report("dims", dims)
    list("dataset" = dataset,  "report" = report, "dims" = dims)
}

toCluster <- function(
    input,
    output_dir = "seuratOutput",
    clusterResolution= 0.5,
    reductionTypes= c("umap"),
    min.pct = 0.25,
    logfc.threshold = 0.25
) {
    dataset = input$dataset
    dims = input$dims
    report = input$report
    name = getNameFromSeurat(dataset)

    print("toCluster")
    dataset <- FindNeighbors(dataset, dims = 1:dims)
    dataset <- FindClusters(dataset, resolution = clusterResolution)
    
    report("clusterCount", table(Idents(dataset)))


    reductionMapping = c("umap" = RunUMAP, "tsne" = RunTSNE)
    lapply(reductionTypes, function(type) {
        dataset <<- reductionMapping[[type]](dataset, dims = 1:dims)
        p = DimPlot(dataset, reduction = type, label = T)
        p$layers[[2]]$aes_params$size = 10
        report(paste0("plot::",type), p)
    })
    
    print("markers")
    dataset.markers <- FindAllMarkers(
        dataset,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold
    )
    
    list(
        "dataset" = dataset,
        "report" = report,
        "dataset.markers" = dataset.markers
    )
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