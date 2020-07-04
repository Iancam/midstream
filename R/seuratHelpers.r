#' @import Seurat

getNameFromSeurat = function(dataset) {
    unique(dataset$datasetName)
    # stringr::str_split(levels(dataset), "/")[[1]][2]
}
getExperimentFromSeurat = function(dataset) {
    unique(dataset$experimentName)
    # stringr::str_split(levels(dataset), "/")[[1]][1]
}