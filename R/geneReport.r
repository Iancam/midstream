#' @import Seurat
#' @import ggplot2
splitToLen = function(list,len) {
    split(list, ceiling(seq_along(list)/len))
}

#' @export
#' @param geneInfo first element is the name, second is the featureList. eg: list("macID", c("Cd68", "H2-Eb1"))
#' @param filename the filename of the dataset, required
#' @param experimentPath the path to the experiment containing the file
geneReport = function(
    geneInfo,
    filename,
    experimentPath,
    output_dir = "geneReportOutput",
    dataset = NULL,
    plots = c("dot")
) {
  if(is.null(dataset)) dataset = readRDS(filename)
  name = geneInfo[[1]]
  features = geneInfo[[2]]
  dir = paste(experimentPath, output_dir, getFname(filename), sep = "/")
  dir.create(dir, recursive= T, showWarnings = F)

  plotMap = list(
    "ridge"   = list(9, function(features) RidgePlot(dataset, features = features, ncol = 3)),
    "vln"     = list(9, function(features) VlnPlot(dataset, features = features)),
    "feature" = list(4, function(features) FeaturePlot(dataset,features = features)),
    "dot"     = list(20, function(features) DotPlot(dataset, features = features))
  )

  lapply(plots, function(plotType){
    foundPlot = plotMap[[plotType]]
    if (is.null(foundPlot)) {
        warning("no plot of name ", plotType, "\n", "try one of: ", names(plotMap))
        return()
    }
    length = foundPlot[1]
    plotter = foundPlot[2]
    plotted = lapply(splitToLen(features, length), plotter)
    lapply(seq_along(plotted), function(plotIndex) {
      lapply(c("pdf", "png"), function(fType){
        fname = paste(name, plotType, plotIndex, fType, sep = ".")
        ggsave(
          paste(dir, fname, sep = "/"),
          plot = plotted[[plotIndex]]
        )
      })
    })
  })    
}