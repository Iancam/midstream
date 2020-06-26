#' @import stringr
getFname = function (path) {
    fName <- sapply(stringr::str_split(path, "/"), tail, 1)
    sapply(stringr::str_split(fName, "\\."), head, 1)
}

makeOutputDir = function(path, experimentDir, output_dir = "seuratOutput" ){
    output_path = paste(experimentDir, output_dir, getFname(path), sep = "/")
    dir.create(output_path, recursive = T, showWarnings = F)
    output_path
}

getInputPaths <- function(input_dir = "input") {
    rdss = list.files(input_dir, full.names = TRUE) %>% Filter(function (x) endsWith(x, ".rds"), .)
    c(list.dirs(input_dir, full.names = TRUE, recursive = FALSE), rdss)
}

#' @param input_dir the directory of the experiment's 10x output
#' @param output_dir the directory to output to, default seuratOutput
#' updates the output directory to receive
#' 
prepInput <- function(experimentDir, input_dir = "input", output_dir = "seuratOutput") {
    input_paths <- getInputPaths(input_dir)
    output_paths <- lapply(input_paths, makeOutputDir, experimentDir)
    list("input_paths"=input_paths,"output_paths"=output_paths)
}



saveDataset <- function(clustered) {
    c(dataset.markers, markerFN, dataset, rdsFN, name) %<-% clustered

    write.csv(dataset.markers, markerFN)
    saveRDS(dataset, file = rdsFN)
}

get10x <- function(dir_path, useFiltered = T) {
    if (useFiltered)
        geneMatrices <- "filtered_gene_bc_matrices"
    else
        geneMatrices <- "raw_gene_bc_matrices"

    Read10X(
        data.dir = paste(
            dir_path,
            "outs",
            geneMatrices,
            "mm10",
            sep = "/"
            )
        )
}