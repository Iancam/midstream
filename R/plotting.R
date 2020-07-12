##### PLOTTING 
#' @import Seurat
reportChanges = function(oldDataset, fx, functionName, args, report){
    newDataset = fx()
    datasetReports = lapply(list(oldDataset, newDataset), function(dataset) {
        c("features" = nrow(dataset),
        "cells" = ncol(dataset),
        "variable features" = length(VariableFeatures(dataset)))
    })
    after = datasetReports[[2]]
    before = datasetReports[[1]]
    diff = after - before

    changes = list()
    changes$args = args
    changes$comparison = cbind(before, after)
    changes$difference = diff
    report(functionName, changes)
    newDataset
}

g = function(plot){
    plot
}

reportFactory = function() {
    report = list()
    function(key=NULL, value=NULL) {
        if (!is.null(key) && !is.null(value)) report[[key]] <<- value
        return(report)
    }
}

mySavePlot <- function(plotted, fname, dir) {
    path = paste0(dir, "/", fname)
    p_ng = paste0(path, ".png")
    p_df = paste0(path, ".pdf")
    try(R.devices::suppressGraphics(ggsave(p_ng, plot = plotted)))
    print(paste0(path, ".pdf"))
    try(R.devices::suppressGraphics(ggsave(p_df, plot = plotted)))
}

LoggerFor <- function(datasetName, output_dir = "seuratOutput") {
    dir = paste(output_dir, datasetName, sep = "/")

    plots = list()
    save = function(plot, plotType) {
        print(plotType)
        fx = function() mySavePlot(plot,
            paste(datasetName, plotType, sep = "."),
            dir)
        
        plots <<- append(plots, fx)
        fx
    }
    list(function() plots, save, paste0(dir, "/run.txt"))
}


Logger = function(output_dir) {
    function(datasetName) LoggerFor(datasetName, output_dir)
}

dumpPrintPlots = function(args) {
    lapply(args[[1]](), function(fn) fn())
    args[2:length(args)]
}

