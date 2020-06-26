#' @import zeallot
#' @import parallel
#' @import magrittr
# source("./seuratPipe.r")
# source("./clusterMapPipe.r")
# source("./plotting")
# source("./fileManip.r")
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
    vc <- prepInput(experimentPath)
    paths <- vc$input_paths
    if (!is.null(targets)) {
        paths <- lapply(targets, function(t) paste(input_dir, t, sep = "/"))
    }
    if (dir.exists(paste0("./", cache_dir))) {
        list.files(cache_dir, full.names = T) %>%
        lapply(getDimsCached) %>%
        mclapply(PCA2Cluster,
            output = output_dir, 
            clusterResolution = clusterResolution ,
            reductionTypes= reductionTypes,
            min.pct = min.pct,
            logfc.threshold = logfc.threshold,
            mc.cores = num_cores) %>%
        lapply(dumpPrintPlots) %>%
        lapply(saveDataset)
    }
    else
    mclapply(paths, tenX2PCA,
    percent.mt.thresh = percent.mt,
    minNFeature = minNFeature,
    maxNFeature = maxNFeature,
    numTop = numTop,
    mc.cores = num_cores
    ) %>%
        lapply(dumpPrintPlots) %>%
        lapply(getDimsInfo) %>%
        mclapply(toCluster,
            output = output_dir,
            clusterResolution = clusterResolution,
            reductionTypes = reductionTypes,
            min.pct = min.pct,
            logfc.threshold = logfc.threshold,
            mc.cores = num_cores) %>%
        lapply(dumpPrintPlots) %>%
        lapply(saveDataset)
}