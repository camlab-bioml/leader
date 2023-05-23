#' Read in user specified marker file.
#'
#' @param yaml_path Path to a yaml file. The file extension must be yml
#'
#' @return A list of positive and negative markers for each cell type
#'
#' @importFrom yaml read_yaml
#' @importFrom dplyr bind_rows
#' @export
get_markers <- function(yaml_path){
  if(!file.exists(yaml_path)){
    stop("Error: it looks like the marker file does not exist in the directory specified.")
  }
  # TODO: try reading in file and throw error if e.g. duplicate keys
  markers <- read_yaml(yaml_path)

  if(!("cell_types" %in% names(markers))){
    stop("The first key of the yaml file must be cell_types.")
  }

  markers <- markers$cell_types
  num_cell_types <- length(markers)

  if(num_cell_types < 1){
    stop("No markers were found in the provided marker file")
  }

  # Find cell types for which positive and negative have been misspelled
  missing_markers <- lapply(seq_along(markers), function(x){
    if(!('positive' %in% names(markers[[x]])) & !('negative' %in% names(markers[[x]]))){
      names(markers[x])
    }
  }) |> unlist()

  if(!is.null(missing_markers)){
    stop(paste0("The following cell types: ", missing_markers, " have neither positive nor negative markers specified in the input marker file. Did you spell positive and negative correctly and all lower case?" ))
  }

  # Check that no cell type has the same marker for positive and negative for the
  # same cell type
  pos_neg_markers <- lapply(seq_along(markers), function(x){
    if('positive' %in% names(markers[[x]]) & 'negative' %in% names(markers[[x]])){
      pos_neg_markers <- markers[[x]]$positive %in% markers[[x]]$negative
      if(any(pos_neg_markers)){
        tibble(cell_type = names(markers[x]),
               markers = paste(markers[[x]]$positive[pos_neg_markers],
                               collapse = ", ")
        )
      }
    }
  }) |>
    bind_rows()

  if(nrow(pos_neg_markers) > 0){
    stop(
      paste0(
        "The following markers: ",
        paste(pos_neg_markers$markers, collapse = ", "),
        " were found to be positive and negative marker(s) for the following cell type(s): ",
        paste(pos_neg_markers$cell_type, collapse = ", "),
        ". A marker can only be a positive or negative marker for a single cell type, but never both.")
      )
  }

  # Calculate the number of positive and negative markers for the each cell type
  marker_counts <- lapply(seq_along(markers), function(x){
    tibble(
      cell_type = names(markers[x]),
      pos_num = length(markers[[x]]$positive),
      neg_num = length(markers[[x]]$negative))
  }) |>
    bind_rows() |>
    mutate(total = pos_num + neg_num)

  # Check if any cell types have no markers
  missing_all_markers <- filter(marker_counts, total == 0) |>
    pull(cell_type)

  if(length(missing_all_markers) > 0){
    stop(paste0("The following cell types have no markers specified: ",
               paste(missing_all_markers, collapse = ", "), "."))
  }

  # Check which cell types only have positive/negative markers
  no_pos_markers <- filter(marker_counts, pos_num == 0) |>
    pull(cell_type)
  no_neg_markers <- filter(marker_counts, neg_num == 0) |>
    pull(cell_type)

  if(length(no_pos_markers) > 0){
    message(paste0("No positive cell types were specified in the input marker file for the following cell types: ",
                   paste(no_pos_markers, collapse = ", "), "."))
  }
  if(length(no_neg_markers) > 0){
    message(paste0("No negative cell types were specified in the input marker file for the following cell types: ",
                   paste(no_neg_markers, collapse = ", "), "."))
  }

  # Get a list of unique markers
  unique_markers <- unlist(markers) |>
    unique()
  message(paste0(length(unique_markers), " unique markers for ", num_cell_types,
                 " cell type(s) were found."))
  markers
}


#' Select cells using adaptive reweighting.
#'
#' @param seu A processed \code{Seurat} object (must be normalized, the most variable features should have already been found using \code{FindVariableFeatures} and should be scaled using \code{ScaleData})
#' @param markers A \code{yaml} file with the set of markers to be used. The name of the marker must match that in \code{seu}, else an error is thrown.
#' @param n_requested_cells Number of cells to sample from dataset. Cannot be bigger than total number of cells in `seu`.
#' @param mode Either \code{marker_based} or \code{cluster_based}. When marker_based is selected, cells are sampled from aggregated clusters that have been putatively assigned to a cell type based on the average gene expression. When cluster_based is used cells are sampled from the individual clusters.
#' @param n_of_pcs Defines the number of PCA dimensions to use for clustering. Defaults to 30.
#' @param knn Defines the number of k for the k-nearest neighbor algorithm during clustering. Defaults to 20.
#' @param res Defines the clustering resulution. Defaults to 0.8.
#'
#' @return None
#'
#' @export
adaptive_reweighting <- function(seu, markers, n_requested_cells,
                                 mode = "marker_based", n_of_pcs = 30,
                                 knn = 20, res = 0.8){
  # Warnings/errors if number of requested cells is too distinct
  if(n_requested_cells >= dim(seu)[2]){
    error("You are trying to select as many or more cells as are present in the dataset.")
  }
  else if(dim(seu)[2] * 0.7 < n_requested_cells){
    warning("You are trying to select more than 70% of the cells in the dataset.")
  }

  cluster_output <- seurat_cluster(seu, n_of_pcs = n_of_pcs, knn = knn, res = res)

  if(mode == "cluster_based"){
    sel_cells <- sample_cells(cluster_output, n_requested_cells = n_requested_cells)
  }else{
    ar <- reweighting(seu, markers, cluster_output)
    sel_cells <- sample_cells(ar, n_requested_cells = n_requested_cells)
  }

  sel_cells
}

#' Creates PCA embedding if missing, clusters cells and returns a cluster - cell id mapping.
#' @returns A dataframe with two columns: cell_id and cluster_id for each cell.
#' @importFrom Seurat RunPCA FindNeighbors FindClusters VariableFeatures Idents
#' @importFrom tibble rownames_to_column
seurat_cluster <- function(seu, n_of_pcs, knn, res){
  # Check if the number of PC's requested is higher than the number of features/cells
  n_features <- length(seu@assays[[1]]@var.features)
  n_cells <- nrow(seu@meta.data)
  
  # If PCA has not been computed go ahead now
  if(!('pca' %in% names(seu@reductions))){
    # use min(n_features, n_cells) as npcs if requested higher
    if(n_of_pcs > n_features | n_of_pcs > n_cells){
      n_of_pcs <- min(n_features, n_cells) - 1
      warning(paste("The number of selected principal components is larger than the number of features or cells in the dataset. Computing", n_of_pcs, "total PCs."))
    }
    seu <- RunPCA(seu, npcs = n_of_pcs, features = VariableFeatures(object = seu))
  }
  # If the number of PCs requested is higher than that calculated in existing PCA
  # calculate pca embedding again
  else if(ncol(seu@reductions$pca@feature.loadings) < n_of_pcs){
    # use min(n_features, n_cells) as npcs if requested higher
    if(n_of_pcs > n_features | n_of_pcs > n_cells){
      n_of_pcs <- min(n_features, n_cells) - 1
      warning(paste("The number of selected principal components is larger than the number of features or cells in the dataset. Computing", n_of_pcs, "total PCs."))
    }
    seu <- RunPCA(seu, npcs = n_of_pcs, features = VariableFeatures(object = seu))
  }

  seu <- FindNeighbors(seu, dims = 1:n_of_pcs, k.param = knn)
  seu <- FindClusters(seu, resolution = res)

  # Get clusters
  clusters <- Idents(seu) |>
    as.data.frame() |>
    rownames_to_column('cell_id')
  colnames(clusters)[2] <- "cluster_id"

  clusters
}


#' Samples the requested number of cells from identified clusters
#' @param cluster_output cell id - cluster/putative cell type mapping dataframe
#' @param n_requested_cells Number of cells to select
#'
#' @importFrom splitstackshape stratified
#' @importFrom dplyr pull
#' @importFrom tibble as_tibble
#' @importFrom dplyr slice_head tally
sample_cells <- function(cluster_output, n_requested_cells){
  n_cells <- 0

  while (n_requested_cells > n_cells){
    # Calculate the maximum number of cells to select per cluster
    n_clusters <- unique(cluster_output$cluster_id) |>
      length()

    if(exists("selected_cells")){
      no_to_select <- ceiling((n_requested_cells - nrow(selected_cells)) / n_clusters)
    }else{
      no_to_select <- ceiling(n_requested_cells / n_clusters)
    }

    # Sample cells from clusters
    new_cells <- stratified(cluster_output, "cluster_id", no_to_select)

    n_cells <- n_cells + nrow(new_cells)

    # Calculate how many cells too many were selected
    oversampled <- n_cells - n_requested_cells

    while(oversampled > 0){
      rem_from_cluster <- new_cells |>
        group_by(cluster_id) |>
        tally() |>
        filter(n == max(n)) |>
        slice_head(n = 1) |>
        pull(cluster_id)

      rem_cell_id <- new_cells |>
        filter(cluster_id == rem_from_cluster) |>
        slice_head(n = 1) |>
        pull(cell_id)

      new_cells <- filter(new_cells, cell_id != rem_cell_id)
      oversampled <- oversampled - 1
    }

    if(exists("selected_cells")){
      selected_cells <- bind_rows(selected_cells, new_cells)
    }else{
      selected_cells <- new_cells
    }

    cluster_output <- filter(cluster_output, !(cell_id %in% selected_cells$cell_id))
  }

  as_tibble(selected_cells)
}

#' Takes into account average marker gene expression to identify overlapping clusters
#' @param seu a seurat object
#' @param markers an object with cell type markers read in by get_markers()
#' @param cluster_output
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join select rename
reweighting <- function(seu, markers, cluster_output){
  expression <- seu@assays$RNA@counts |>
    as.matrix() |>
    t() |>
    as.data.frame() |>
    rownames_to_column('cell_id')

  unique_markers <- unlist(markers) |>
    unique()
  expression <- expression[, c('cell_id', unique_markers)] |>
    left_join(cluster_output, by = "cell_id") |>
    select(-cell_id)

  cluster_assignments <- cluster_enrich(markers, expression)

  # Join cell id's with cell type label
  assignments <- left_join(cluster_output, cluster_assignments, by = "cluster_id") |>
    select(-cluster_id) |>
    rename(cluster_id = predicted_cell_type)

  assignments
}


#' calculates the average expression of marker genes in each cluster
#' @param markers an object with cell type markers read in by get_markers()
#' @param expression dataframe with expression values for all cells and marker genes. rows are cells, columns are markers. The first column should be `cell_id`.
#'
#' @importFrom dplyr group_by summarize mutate filter select
#' @importFrom tidyr pivot_longer
#' @importFrom tibble tibble
cluster_enrich <- function(markers, expression){
  # Calculate enrichments for each cell type/cluster
  enrichments <- lapply(1:length(markers), function(x){

    # Positive enrichment
    p <- expression[, c(markers[[x]]$positive, 'cluster_id')] |>
      pivot_longer(-cluster_id, names_to = "marker", values_to = "expression") |>
      group_by(cluster_id) |>
      summarize(mean = mean(expression))

    if(!is.null(markers[[x]]$negative)){
      # If this cell types has negative markers - calculate negative enrichment
      n <- expression[, c(markers[[x]]$negative, 'cluster_id')] |>
        pivot_longer(-cluster_id, names_to = "marker", values_to = "expression") |>
        group_by(cluster_id) |>
        summarize(mean = mean(expression))

      diff <- p$mean - n$mean
    }else{
      diff <- p$mean
    }

    assignment <- tibble(cluster_id = p$cluster_id,
                         avrg_expression = diff,
                         predicted_cell_type = names(markers[x]))
    assignment
  }) |> bind_rows()

  # Find cell type with max value for each cluster
  cluster_assignments <- enrichments |>
    group_by(cluster_id) |>
    mutate(max_score = avrg_expression == max(avrg_expression)) |>
    filter(max_score == TRUE) |>
    select(-max_score, -avrg_expression)

  cluster_assignments
}
