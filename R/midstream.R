#' @import parallel
#' @import magrittr
#' @import Seurat
num_cores <- parallel::detectCores()
#'
#' @param experimentPath path to a dir containing input_dir (defaults to input), will contain output directories
#' @param output_dir = "seuratOutput":
#' @param input_dir = "input": directory name in experimentPath containing 10x output, should look like "outs/filtered_gene_bc_matrices/mm10"
#' @param cache_dir = "pcaCache": directory name to store intermediate data from PCA step
#' @param combine = TRUE: set to false if you have already created a comb.rds input file. Useful if rerunning the analysis
#' @param targets = NULL: list of individual 10x filepaths on which to run the analysis. Useful if you're just testing the set-up
#' @param useFiltered = TRUE: if True, selects filtered_gene_bc_matrices, else selects raw_gene_bc_matrices from the 10x folders
#' @param percent.mt = 5: exclude cells with more than percent.mt mitochondrial dna. 
#' @param minNFeature=200; exclude cells with less than minNFeature features 
#' @param maxNFeature=2500, exclude cells with more than maxNFeature features 
#' @param numTop =10: number of genes to label as contributing to variability
#' @param clusterResolution= 0.5:
#' @param reductionTypes= c("umap"): currently supports umap, and tsne
#' @param min.pct = 0.25: 
#' @param logfc.threshold = 0.25
tenX2Seurat = function (
    experimentPath,
    targets = NULL,
    output_dir = "seuratOutput",
    input_dir = "input",
    useFiltered = TRUE,
    combine = TRUE,
    cache_dir = "pcaCache",
    ignore_cache = F,
    getDimsManually = F,
    suppressIO = F,
    defaultDims = 20,
    percent.mt = 5,
    minNFeature=200,
    maxNFeature=2500,
    numTop =10,
    nfeatures = 2000,

    clusterResolution= 0.5,
    reductionTypes= c("umap"),
    min.pct = 0.25,
    logfc.threshold = 0.25
) {
    sink()
    
    input_dir  <- paste(experimentPath, input_dir, sep = "/")
    if(!file.exists(input_dir)) {
        warning("input dir missing from path: ", experimentPath, "\n",
        "dir: ", input_dir, "\nskipping")

        return()
    }
    output_dir <- paste(experimentPath, output_dir, sep = "/")
    cache_dir  <- paste(experimentPath, cache_dir, sep = "/")
    if(combine) {
        tenX2Combined(input_dir = input_dir)
    }
    vc <- prepInput(experimentPath)
    
    plotSaver = Logger(output_dir)
    if(suppressIO) dumpPrintPlots = function(args) {args[2:length(args)]}
    cache_paths = c()
    if (dir.exists(cache_dir) && ! ignore_cache) cache_paths = list.files(cache_dir, full.names = T)
    cached_names = lapply(cache_paths, getFname)
    input_paths = Filter(function(path) not(getFname(path) %in% cached_names), vc$input_paths)
    print(input_paths)
    if (!is.null(targets)) {
        input_paths <- Filter(function(x) getFname(x) %in% targets, input_paths)
        cache_paths <- Filter(function(x) getFname(x) %in% targets, cache_paths)
    }
    
    print("cache_paths")
    print(cache_paths)
    print("input_paths")
    print(input_paths)
    
    from_cache = lapply(
            cache_paths, 
            function(path) getDims(
                dir_path = path,
                plotSaver=plotSaver,
                getManually = getDimsManually,
                dims = defaultDims
                )
        ) %>%
        future_lapply(toCluster,
            plotSaver = plotSaver,
            output_dir = output_dir,
            clusterResolution = clusterResolution,
            reductionTypes = reductionTypes,
            min.pct = min.pct,
            logfc.threshold = logfc.threshold
            # mc.cores = num_cores
        ) %>%
        lapply(dumpPrintPlots) %>%
        lapply(saveDataset)

    from_input =
    future_lapply(input_paths, loadFile) %>%
    future_lapply(toPCA,
        cache_dir = cache_dir,
        plotSaver = plotSaver,
        percent.mt.thresh = percent.mt,
        minNFeature = minNFeature,
        maxNFeature = maxNFeature,
        numTop = numTop,
        nfeatures = nfeatures
        ) %>%
        lapply(dumpPrintPlots) %>%
        lapply(function(info) getDims(
            dir_info = info,
            plotSaver = plotSaver,
            getManually = getDimsManually,
            dims = defaultDims)
        ) %>%
        future_lapply(toCluster,
            plotSaver = plotSaver,
            output_dir = output_dir,
            clusterResolution = clusterResolution,
            reductionTypes = reductionTypes,
            min.pct = min.pct,
            logfc.threshold = logfc.threshold
            ) %>%
        lapply(dumpPrintPlots) %>%
        lapply(saveDataset)
    
    datasets = c(from_cache, from_input)
    lapply(datasets, function(clustered) c(
         "dataset"= clustered[[3]] ,
         "filename"= clustered[[4]] ,
         "output"= output_dir     
        )
    )
}