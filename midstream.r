library(dplyr)
library(Seurat)
library(patchwork)
library(parallel)
library(zeallot)
library(magrittr)
library(ClusterMap)


#' @import zeallot
num_cores <- detectCores()

##### FILE MANIP

getFname = function (path) {
    fName <- sapply(stringr::str_split(path, "/"), tail, 1)
    sapply(stringr::str_split(fName, "\\."), head, 1)
}

makeOutputDir = function(path, output_dir = "seuratOutput"){
    output_path = paste(output_dir, getFname(path), sep = "/")
    dir.create(output_path, showWarnings= FALSE)
    output_path
}

getInputPaths <- function(input_dir = "input") {
    rdss = list.files(input_dir, full.names = TRUE) %>% Filter(function (x) endsWith(x, ".rds"), .)
    c(list.dirs(input_dir, full.names = TRUE, recursive = FALSE), rdss)
}

#' @param input_dir the directory of the experiment's 10x output
#' @param output_dir the directory to output to, default seuratOutput
#' updates the output directory to receive
#' 
prepInput <- function(input_dir = "input", output_dir = "seuratOutput") {
    input_paths <- getInputPaths(input_dir)
    output_paths <- lapply(input_paths, makeOutputDir)
    list("input_paths"=input_paths,"output_paths"=output_paths)
}



saveDataset <- function(clustered) {
    c(dataset.markers, markerFN, dataset, rdsFN, name) %<-% clustered

    write.csv(dataset.markers, markerFN)
    saveRDS(dataset, file = rdsFN)
}

get10x <- function(dir_path, useFiltered = T) {
    if (useFiltered)
        geneName <- "filtered_gene_bc_matrices"
    else
        geneName <- "raw_gene_bc_matrices"

    Read10X(
        data.dir = paste(
            dir_path,
            "outs",
            geneName,
            "mm10",
            sep = "/"
            )
        )
}


##### PLOTTING 

mySavePlot <- function(plotFx, fname, dir) {
    print(fname)
    curwd <- getwd()
    if (!is.null(dir)) {
        setwd(dir)
    }
    plotted <- plotFx()
    png(paste0(fname, ".png"))
        print(plotted)
    dev.off()
    pdf(paste0(fname, ".pdf"), onefile=FALSE)
        print(plotted)
    dev.off()
    setwd(curwd)
}

plotsFor <- function(datasetName, output_dir = "seuratOutput") {
    plots = list()
    save = function(plot, plotType) {
        print(plotType)
        fx = function() mySavePlot(plot,
            paste(datasetName, plotType, sep = "."),
            paste(output_dir, datasetName, sep = "/"))
        
        plots <<- append(plots, fx)
        fx
    }
    list(function() plots, save)
}

dumpPrintPlots = function(args) {
    lapply(args[[1]](), function(fn) fn())
    args[2:length(args)]
}

##### SEURAT 

toPCA <- function(
        path,
        cache_dir = "pcaCache",
        percent.mt = 5,
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
    save_plot(function() VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3),"mito.vln")
    plot1 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(dataset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    save_plot(function() plot1+plot2, "mito.scatter")
    print("Normalize")
    if (!is.null(minNFeature)) {
        if(!is.null(maxNFeature)) dataset <- subset(dataset, subset = nFeature_RNA > minNFeature & nFeature_RNA < maxNFeature & percent.mt < percent.mt)
        else dataset <- subset(dataset, subset = nFeature_RNA > minNFeature & percent.mt < percent.mt)   
    } else if(!is.null(maxNFeature)) dataset <- subset(dataset, subset = nFeature_RNA < maxNFeature & percent.mt < percent.mt)
    else dataset <- subset(dataset, subset = percent.mt < percent.mt)
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

getDims <- function(dir_info) {
    c(dataset, name, dir_path) %<-% dir_info
    name = getFname(dir_path)

    c(plots, save_plot) %<-% plotsFor(name)
    dataset <- readRDS(dir_path)

    elbow <- ElbowPlot(dataset)
    print(elbow)
    save_plot(function() elbow, "elbow")()

    dims <- as.numeric(readline(prompt = "what shall be your dimensions? "))
    c(dataset, name, dims)
}

toCluster <- function(
    info,
    output_dir = "seuratOutput",
    clusterResolution= 0.5,
    reductionTypes= c("umap"),
    min.pct = 0.25,
    logfc.threshold = 0.25) {
    
    print("toCluster")
    c(dataset, name, dims) %<-% info
    c(plots, save_plot) %<-% plotsFor(name)
    markerFN <- paste(output_dir, paste0(name, ".markers.csv"), sep = "/")
    rdsFN <- paste(output_dir, paste0(name, ".rds"), sep = "/")
    dataset <- FindNeighbors(dataset, dims = 1:dims)
    dataset <- FindClusters(dataset, resolution = clusterResolution)
    
    reductionMapping = c("umap"=RunUMAP, "tsne"=RunTSNE)
    lapply(reductionTypes, function(name) datset <- reductionMapping[name](dataset, dims = 1:dims))
    save_plot(function() DimPlot(dataset, reduction = reductionType), reductionType)
    # print(DimPlot(dataset, reduction = "umap"))
    print("markers")
    dataset.markers <- FindAllMarkers(
        dataset,
        only.pos = TRUE,
        min.pct =min.pct,
        logfc.threshold =logfc.threshold
    )
    print(markerFN)
    list(plots, dataset.markers, markerFN, dataset, rdsFN, name)
}

###### CLUSTERMAP

tenX2Combined <- function(
    input_dir = "input",
    comb_delim = "_",
    targets = NULL,
    sub_delim="-",
    min.cells= 3,
    min.features= 200) {

    dirs = list.dirs(input_dir, recursive = FALSE)
    print("toCombined")
    named_data <- lapply(dirs, function(dir_path) {
        dir_name <- sapply(stringr::str_split(dir_path, "/"), tail, 1)
        data <- getnumTopx(dir_path)
        id <- gsub(comb_delim, sub_delim, dir_name)
        colnames(data) <- paste0(id, "_", colnames(data))
        data
    })
    data <- CreateSeuratObject(do.call(cbind, named_data),
        min.cells =min.cells,
        min.features =min.features
    )
    saveRDS(data, file = paste0(input_dir, "/comb.rds"))
}

seurat2ClusterMap <- function(
        seurat_dir = "seuratOutput",
        output_dir = "clusterMapOutput",
        output_name = NULL,
        comb_delim = "_",
        sub_delim="-") 
    {
        dir.create(output_dir,showWarnings= FALSE)
        curr_dir = getwd()
        setwd(seurat_dir)
        files = list.files(full.names = TRUE)
        markers <- Filter(function(x) endsWith(x, ".markers.csv") && !startsWith(basename(x), "comb"), files)
        comb <- Filter(function(x) endsWith(x, ".rds") && startsWith(basename(x), "comb"), files)
        rdsFiles <- Filter(function(x) endsWith(x, ".rds") && !startsWith(basename(x), "comb"), files)
        obj.names = lapply(rdsFiles, function(path) gsub(comb_delim, sub_delim, getFname(path)))
        single_obj_list <- lapply(rdsFiles,readRDS)
        single_obj_list <- setNames(single_obj_list, obj.names)
        markers = setNames(markers, obj.names)
        comb = readRDS(comb)  
        
        results <- cluster_map(
            markers,
            edge_cutoff = 0.1,
            output = paste(curr_dir, output_dir, basename(curr_dir), sep="/"),
            single_obj_list = single_obj_list,
            comb_obj = comb,
            comb_delim = comb_delim
        )
        dir.create(paste(curr_dir, output_dir,"figures",sep="/"))
        figures <- Filter(function(x) endsWith(x, "png") || endsWith(x, "pdf"), list.files())
        print(figures)
        print(output_dir)
        lapply(figures, function (file) file.rename(file, paste(curr_dir, output_dir,"figures", file, sep="/")))
        setwd(curr_dir)
}

tenX2Seurat = function (
    experimentPath,
    targets = NULL,
    output_dir = "seuratOutput",
    input_dir = "input",
    useFiltered = TRUE,
    combine = TRUE,

    percent.mt = 5,
    minNFeature=200,
    maxNFeature=2500,
    numTop =10,

    clusterResolution= 0.5,
    reductionTypes= c("umap"),
    min.pct = 0.25,
    logfc.threshold = 0.25
) {
    setwd(experimentPath)

    if(combine) {
        tenX2Combined(targets= targets)
    }
    vc <- prepInput()
    paths <- vc$input_paths
    if(!is.null(targets)){
        paths <- lapply(targets, function(t) paste(input_dir,t, sep="/") )
    }

    mclapply(paths, toPCA, mc.cores = num_cores) %>%
        lapply(dumpPrintPlots) %>%
        lapply(getDims) %>%
        mclapply(toCluster,
            output = output_dir, 
            clusterResolution = clusterResolution ,
            reductionTypes= reductionTypes,
            min.pct = min.pct,
            logfc.threshold = logfc.threshold,
            mc.cores = num_cores) %>%
        lapply(dumpPrintPlots) %>%
        lapply(saveDataset)
}


# , targets= c("comb.rds"), combine=F)
tenX2Seurat(experimentPath = "~/experiments/IM", reductionTypes = c("umap", "tsne"), maxNFeature=NULL)
seurat2ClusterMap()