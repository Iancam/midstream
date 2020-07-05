#' @import magrittr
#' @import Seurat
#' @import future.apply
#' @import drake

num_cores <- parallel::detectCores()

#' @param experimentPath path to a dir containing input_dir (defaults to input), will contain output directories
#' @param output_dir = "seuratOutput":
#' @param input_dir = "input": directory name in experimentPath containing 10x output, should look like "outs/filtered_gene_bc_matrices/mm10"
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
#' @export
tenX2Seurat = function (
    experimentPath,
    targets = NULL,
    output_dir = "seuratOutput",
    input_dir = "input",
    useFiltered = TRUE,
    combine = TRUE,
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
    if(combine) {
        tenX2Combined(input_dir = input_dir)
    }
    vc <- prepInput(experimentPath)
    plotSaver = Logger(output_dir)
    if (suppressIO) dumpPrintPlots = function(args) args[2:length(args)]

    if (!is.null(targets)) {
        vc$input_paths <- Filter(function(x) getFname(x) %in% targets, vc$input_paths)
    }

    num_cores <- future::availableCores()
    future::plan(future::multiprocess, workers = num_cores)

    dumpFixer = function(fn) {
        function(...) fn(...) %>% dumpPrintPlots()
    }

    fixedFuncs = lapply(list(
        QC,
        toPCA,
        getDims,
        toCluster
     ), dumpFixer)

    QC = fixedFuncs[[1]]
    toPCA = fixedFuncs[[2]]
    getDims = fixedFuncs[[3]]
    toCluster = fixedFuncs[[4]]

    pipeLine = drake_plan(
            input = target(
                loadFile(file_in(assayPath)),
                transform = map(assayPath = !!vc$input_paths)
            ),
            qc = target(QC(input, plotSaver = plotSaver,
                    percent.mt.thresh = percent.mt,
                    minNFeature = minNFeature,
                    maxNFeature = maxNFeature
                    ), transform = map(input)),
            pca = target(toPCA(qc, 
                        plotSaver = plotSaver,
                        numTop = numTop,
                        nfeatures = nfeatures
                        ), transform = map(qc)),
            withDims = target(getDims(pca, 
                        plotSaver = plotSaver,
                        getManually = getDimsManually,
                        dims = defaultDims
                        ), transform = map(pca)),
            clustered = target(toCluster(withDims, 
                        plotSaver = plotSaver,
                        output_dir = output_dir,
                        clusterResolution = clusterResolution,
                        reductionTypes = reductionTypes,
                        min.pct = min.pct,
                        logfc.threshold = logfc.threshold
                    ), transform = map(withDims)),
            fin = target(saveDataset(clustered), transform = map(clustered))
            
        )
    

    vis_drake_graph(pipeLine)
    make(pipeLine)
}