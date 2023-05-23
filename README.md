# Leader introduction
Cel**L** s**E**letion vi**A** a**D**aptiv**E** **R**eweighting (LEADER) is a method to subset cells aiming to capture rare cell types with the goal of creating more balanced datasets. Leader requires a processed seurat object with gene expression data and an optional marker file specifying known marker genes for known cell types.

Leader generates subsampled datasets by clustering all cells and selecting the desired number of cells from these clusters. In the marker aware version, each cluster is putatively assigned to a cell type using the average expression of marker genes for each expected cell type provided by the user. The desired number of cells is then selected from the putatively assigned cell types:

<img width="751" alt="image" src="https://github.com/camlab-bioml/leader/assets/45369908/4337c2ef-ba82-418b-8f96-c776dd1c674d">


## Installing leader
```
BiocManager::install("leader")
```

## Running leader
### Data processing
The input seurat oject must be log normalized and the most variable features must be identified and scaled prior to cell selection with LEADER as follows:

```
library(Seurat)
devtools::load_all()

# Load example data
data(seu, package = "leader")

# Normalize expression and find variable features
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 2000)

# Scale gene expression
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)
```

Note that the seurat object must also have unique cell ID's as the rownames of the dataframe in the `seu@meta.data` slot. These ID's will then be used to output a subset of cells that can be used as a more balanced dataset.

### Cell selection

Next cells are selected using the LEADER package. This can be in a cell type marker aware way which requires a list of positive and optional negative markers for each cell type expected in the data.

### Cluster based selection

Running LEADER without markers simply requires setting `mode = "cluster_based"`, and determining the required number of cells to select as follows:

```
adaptive_reweighting(seu, n_requested_cells = 40, mode = "cluster_based")
```


### Marker informed selection

The marker informed cell selection requires the creation of a `yaml` file containing marker genes for the expected cell types. These could be positive markers (e.g. CD3 for T-cells) or negative markers (e.g. B220 for T-cells).

The marker yaml file must follow the following structure:

    cell_types:
      CD8 T cells:
        positive:
          - CD3
          - CD8
        negative:
          - B220
          - CD4
      
      CD4 T cells:
        positive:
          - CD3
          - CD4
        negative:
          - CD8
          - B220
          
      B cells:
        positive:
          - B220
          - IgM
      
      ...

The first row must contain the word `cell_types`. Each indented element is an expected cell type with required positive labels and optional negative labels.

To run LEADER using the marker informed selection simply read in a marker file using the `get_markers()` function and add this object to the `adaptive_reweighting` call along with the `mode = "marker_based"` parameter.

```
markers <- get_markers(system.file("extdata", "markers.yml", package = "leader"))
adaptive_reweighting(seu, markers, n_requested_cells = 40, mode = "marker_based")
```

The output of the `adaptive_reweighting` call can then be saved and the cell_id's can be used as more balanced training data.
