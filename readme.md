# Midstream

Midstream is a wrapper designed to simplify single-cell RNA sequencing analysis, minimizing manual work and making results more reproducible and shareable.

## Overview

"Midstream" provides a streamlined, single-click solution for biologists performing single-cell RNA sequencing (scRNA-seq) clustering analyses. Built on Seurat, a widely used R package, Midstream aims to automate repetitive processes, allowing users with minimal coding experience to obtain reproducible, shareable insights with ease.

## Why Midstream?

Using Seurat for scRNA-seq analysis can involve manually running extensive code and managing multiple parameters and output files, which becomes time-consuming and difficult to track across experiments. "Midstream" alleviates these issues by:

- **Reducing Repetitive Work**: Automates common scRNA-seq workflows, minimizing user input.
- **Improving Reproducibility**: Standardizes processes, making it easy to replicate and share analyses.
- **Simplifying Analysis Management**: Organizes output files consistently to avoid confusion, especially when working with multiple datasets and parameters.

## Installation

```r
devtools::install_github("iancam/midstream")
```

## Using midstream

- create a folder for your experiment. For the example, we'll call it experimentFolder. The name of the folder will be used in the project, so name it thoughtfully.

- Inside of the experiment folder, create a folder called input.

- Copy the relevant 10X directories into that folder.

- Rename **those** directories to reflect their experimental names, ie. control, day 1, etc.

- create an analysis.r file in the experiment folder.
  Your analysis might look something like this.

```
├── experiment_name
│   ├── input
│   │   ├── control
│   |   ├── day_1
│   |   ├── day_5
│   |   ├── day_30
|   |–– analysis.r
```

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

## Comparison to Seurat Workflow.

Seurat greatly simplifies single cell RNA sequencing already, but if you're just following the vignettes, you might end up doing a lot of repetitive work without realizing it.

A simple Seurat run through (with plots) looks like this:

```r
# Read in the data
inputPath = paste(dir_path,"outs",geneName,"mm10",sep = "/")
dataset = Read10X(data.dir = inputPath)

# visualize and filter on percent mt threshold
percent.mt = PercentageFeatureSet(dataset, pattern = "^[Mm][Tt]-")
dataset[["percent.mt"]] = percent.mt
VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, combine = T)
dataset = subset(dataset, subset = nFeature_RNA > minNFeature)
dataset = subset(dataset, subset = nFeature_RNA < maxNFeature)
dataset = subset(dataset, subset = percent.mt < percent.mt.thresh)
VlnPlot(dataset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, combine = T)

# run and analyze pca
dataset = SCTransform(dataset)
dataset = RunPCA(dataset)
ElbowPlot(dataset)
if(plotHeatmap) DimHeatmap(dataset, dims = 1:20, cells = 500, balanced = TRUE))

# Cluster using umap space
dataset = RunUMAP(dataset, dims = 1:20)
dataset = FindNeighbors(dataset, dims = 1:20)
dataset = FindClusters(dataset)
UMAPPlot(dataset)

# find markers
markers = FindAllMarkers(
        dataset,
        min.pct = min.pct,
        logfc.threshold = logfc.threshold
    )
clusterCount = table(Idents(dataset))
```

That's at least 18 steps per assay, not to mention managing the plot and summary statistic output. Imagine typing all of that into the console every time you want to run an analysis. This gets overwhelming when you've got:

- multiple assays (our lab had 4-5 per experiment),
- multiple filtering parameters that you'd like to try
- people asking for reports every two weeks.

You'd spend a lot of time choosing names and locations for those files, and it gets hard to tell which analysis produced which results.

midstream takes what could be a day or even multi-day analysis and makes it take a few minutes.

# Associated Paper

you can find an associated analysis [here](https://www.researchgate.net/publication/360846890_Single_Cell_RNA_Sequencing_and_Binary_Hierarchical_Clustering_Defines_Lung_Interstitial_Macrophage_Heterogeneity_in_Response_to_Hypoxia)

# Future Directions

- command line tool
- include functional analyses, IE GSEA, GAGE
- include cell trajectory analyses, IE SCVelo
