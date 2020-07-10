#' @import Seurat

getNameFromSeurat = function(dataset) {
    unique(dataset$datasetName)
    # stringr::str_split(levels(dataset), "/")[[1]][2]
}
getExperimentFromSeurat = function(dataset) {
    unique(dataset$experimentName)
    # stringr::str_split(levels(dataset), "/")[[1]][1]
}

splitClusters = function(markers, path, name, inp_p_val = 0.05) {
    dir = paste(path, "clusters", name, sep = "/")
    dir.create(dir, recursive = T)
    clusters = unique(markers$clust)
    lapply(clusters, function(clusterName) {
        subcluster = subset(subset(markers, cluster == clusterName), p_val_adj <= inp_p_val)

        write.csv(subcluster, paste0(dir, "/", clusterName, ".csv"))})
}