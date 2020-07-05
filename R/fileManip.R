#' @import stringr
#' @import dplyr
#' @import Seurat
#' @importFrom utils head tail write.csv
getFname = function (path) {
    fName <- sapply(stringr::str_split(path, "/"), tail, 1)
    sapply(stringr::str_split(fName, "\\."), head, 1)
}

getExperimentName = function(inputFile) {
    basename(normalizePath(paste0(inputFile, "/../../")))
}

makeOutputDir = function(path, experimentDir, output_dir = "seuratOutput" ){
    output_path = paste(experimentDir, output_dir, getFname(path), sep = "/")
    dir.create(output_path, recursive = T, showWarnings = F)
    output_path
}

getInputPaths <- function(experimentDir, input_dir = "input") {
    dir = paste(experimentDir,input_dir, sep="/")
    rdss = list.files(dir, full.names = TRUE) %>% Filter(function (x) endsWith(x, ".rds"), .)
    c(list.dirs(dir, full.names = TRUE, recursive = FALSE), rdss)
}

#' @param input_dir the directory of the experiment's 10x output
#' @param output_dir the directory to output to, default seuratOutput
#' updates the output directory to receive
#' 
prepInput <- function(experimentDir, input_dir = "input", output_dir = "seuratOutput") {
    input_paths <- getInputPaths(experimentDir)
    output_paths <- lapply(input_paths, makeOutputDir, experimentDir)
    list("input_paths"=input_paths,"output_paths"=output_paths)
}



saveDataset <- function(clustered) {
    dataset = clustered[[1]]
    dataset.markers = clustered[[2]]
    markerFN = clustered[[3]]
    rdsFN = clustered[[4]]
    write.csv(dataset.markers, markerFN)
    dir = paste(dirname(rdsFN), "clusters", getFname(rdsFN), sep = "/")
    dir.create(dir)
    clusters = unique(dataset.markers$clust)
    lapply(clusters, function(cluster) {
        subcluster = subset(dataset.markers, cluster == cluster)
        write.csv(subcluster, paste0(dir, "/", cluster, ".csv"))})
    saveRDS(dataset, file = rdsFN)
    clustered
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

loadFile = function(path) {
  name <- getFname(path)
  experimentName = getExperimentName(path)
  if(path %>% endsWith("rds")) dataset <- readRDS(path)
  else dataset <- CreateSeuratObject(get10x(path))
  dataset$experimentName = experimentName
  dataset$datasetName = name
  dataset$orig.ident <- paste(experimentName, name, sep = "/")
  dataset
}
