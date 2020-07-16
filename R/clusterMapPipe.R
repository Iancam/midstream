
#' @import stringr
#' @import Seurat
#' @import ClusterMap
#' @importFrom stats setNames
#' @export
tenX2Combined <- function(
    input_dir = "input",
    comb_delim = "_",
    targets = NULL,
    sub_delim="-",
    min.cells= 3,
    min.features= 200) {

    dirs = list.dirs(input_dir, recursive = FALSE)
    print("toCombined")
    named_data <- lapply(dirs, function(dir_path) {
        dir_name <- sapply(stringr::str_split(dir_path, "/"), tail, 1)
        data <- get10x(dir_path)
        id <- gsub(comb_delim, sub_delim, dir_name)
        colnames(data) <- paste0(id, "_", colnames(data))
        data
    })
    data <- CreateSeuratObject(do.call(cbind, named_data),
        min.cells =min.cells,
        min.features =min.features
    )
    saveRDS(data, file = paste0(input_dir, "/comb.rds"))
}
#' @export

seurat2ClusterMap <- function(
        experimentPath,
        seurat_dir = "seuratOutput",
        output_dir = "clusterMapOutput",
        output_name = NULL,
        comb_delim = "_",
        sub_delim="-") 
    {
        seurat_dir = paste(experimentPath, seurat_dir, sep = "/")
        output_dir = paste(experimentPath, output_dir, sep = "/")
        dir.create(output_dir,showWarnings= FALSE)
        files = list.files(seurat_dir, full.names = TRUE)
        markers <- Filter(function(x) endsWith(x, ".markers.csv") && !startsWith(basename(x), "comb"), files)
        comb <- Filter(function(x) endsWith(x, ".rds") && startsWith(basename(x), "comb") && !endsWith(x, ".report.rds"), files)
        rdsFiles <- Filter(function(x) endsWith(x, ".rds") &&
            !startsWith(basename(x), "comb") &&
            !endsWith(x, ".report.rds"), files)
        obj.names = lapply(rdsFiles, function(path) gsub(comb_delim, sub_delim, getFname(path)))
        single_obj_list <- lapply(rdsFiles, readRDS)
        single_obj_list <- setNames(single_obj_list, obj.names)
        markers = setNames(markers, obj.names)
        comb = readRDS(comb)  
        
        results <- cluster_map(
            markers,
            edge_cutoff = 0.1,
            output = paste(output_dir, basename(experimentPath), sep="/"),
            single_obj_list = single_obj_list,
            comb_obj = comb,
            comb_delim = comb_delim
        )
        dir.create(paste(output_dir, "figures", sep = "/"))
        figures <- Filter(function(x) endsWith(x, "png") || endsWith(x, "pdf"), list.files())
        print(figures)
        print(output_dir)
        lapply(figures, function (file) file.rename(file, paste(output_dir, "figures", file, sep="/")))
}