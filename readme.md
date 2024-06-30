# Midstream

([Seurat](https://satijalab.org/seurat/) for biologists)

```r
devtools::install_github("iancam/midstream")
```

## What

A "one-click" solution for biologists doing single cell RNA seq clustering analyses.

## Using midstream

- create a folder for your experiment. For the example, we'll call it experimentFolder. The name of the folder will be used in the project, so name it thoughtfully.

- Insie of the experiment folder, create a folder called input.

- Copy the relevant 10X directories into that folder.

- Rename **those** directories to reflect their experimental names, ie. control, day 1, etc.

- create an analysis.r file in the experiment folder
  Your analysis might look something like this.

```r
library(midstream)
# run the standard seurat pipeline on each of the files in input
# and on the combination of all the files together
tenX2Seurat("./")

# collect your clusters together and compare them
seurat2ClusterMap("./")

# create reports for specific genes
geneReport("./", name = "interestingGenes", geneList = c("caspr", "alg1", "pepp2")))

```

- To run your analysis, open terminal and cd to your experiment directory, then run `rscript analysis.r`

## Why

- someone with little coding experience can take days to perform an analysis on a single dataset
- results aren't easily reproducible or shareable

## But Seurat is Great!

It really is. But if you're just following the vignettes, you might end up doing a lot of repetitive work without realizing it.

A simple Seurat run through (with plots) looks like this:

```r
inputPath = paste(dir_path,"outs",geneName,"mm10",sep = "/")
dataset = Read10X(data.dir = inputPath)
percent.mt = PercentageFeatureSet(dataset, pattern = "^[Mm][Tt]-")
dataset[["percent.mt"]] = percent.mt
VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, combine = T)
dataset = subset(dataset, subset = nFeature_RNA > minNFeature)
dataset = subset(dataset, subset = nFeature_RNA < maxNFeature)
dataset = subset(dataset, subset = percent.mt < percent.mt.thresh)
VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, combine = T)
dataset = SCTransform(dataset)
dataset = RunPCA(dataset)
ElbowPlot(dataset)
if(plotHeatmap) DimHeatmap(dataset, dims = 1:20, cells = 500, balanced = TRUE))
dataset = RunUMAP(dataset, dims = 1:20)
dataset = FindNeighbors(dataset, dims = 1:20)
dataset = FindClusters(dataset)
UMAPPlot(dataset)
markers = FindAllMarkers(
        dataset,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold
    )
clusterCount = table(Idents(dataset))
```

Imagine typing all of that into the console every time you want to run an analysis. Also, you still have to manage saving the plots and data. This gets overwhelming when you've got:

- multiple assays,
- multiple filtering parameters that you'd like to try
- people asking for reports every two weeks.

You'd spend a lot of time choosing names and locations for those files, and it gets hard to tell which analysis produced which results.
# Future Directions

- command line tool
- include functional analyses, IE GSEA, GAGE
- include cell trajectory analyses, IE SCVelo
