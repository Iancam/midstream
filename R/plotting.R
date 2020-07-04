##### PLOTTING 
#' @import Seurat
reportChanges = function(oldDataset, fx, functionName, args){
    newDataset = fx()
    datasetReports = lapply(list(oldDataset, newDataset), function(dataset) {
        c("features" = nrow(dataset),
        "cells" = ncol(dataset),
        "variable features" = length(VariableFeatures(dataset)))
    })
    after = datasetReports[[2]]
    before = datasetReports[[1]]
    diff = after - before
    cat("\nafter applying "); cat(functionName); cat(" with args: (\n"); print(args); cat(")\n")
    print(cbind(before, after))
    print("difference: ")
    print(diff)
    newDataset
}

mySavePlot <- function(plotFx, fname, dir) {
    path = paste0(dir, "/", fname)
    plotted <- plotFx()
    print(paste0(path, ".png"))
    try(ggsave(paste0(path, ".png"), plot = plotted))
    print(paste0(path, ".pdf"))
    try(ggsave(paste0(path, ".pdf"), plot = plotted))
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

