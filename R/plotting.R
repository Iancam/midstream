##### PLOTTING 

mySavePlot <- function(plotFx, fname, dir) {
    print(fname)
    plotted <- plotFx()
    png(paste0(dir, "/", fname, ".png"))
        print(plotted)
    dev.off()
    pdf(paste0(dir, "/", fname, ".pdf"), onefile=FALSE)
        print(plotted)
    dev.off()
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
