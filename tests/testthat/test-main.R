test_that("get_markers() reads in a yaml marker file", {
  expect_error(get_markers("marker_files/wrong-key.yml"))
  expect_error(get_markers("marker_files/no-celltypes.yml"))
  expect_error(get_markers("marker_files/misspelled-pos-neg.yml"))
  expect_error(get_markers("marker_files/no-markers-per-celltype.yml"))
  expect_error(get_markers("non-existent-marker-file.yml"))
  expect_error(get_markers("marker_files/same-neg-pos-marker.yml"))
  expect_message(get_markers("marker_files/only-positive.yml"))

  expect_equal(get_markers("marker_files/correct-markers.yml"),
               list("CD4+ T cell" = list("positive" = c("CD3", "CD4"),
                                         "negative" = c("CD8"))))
})


test_that("seurat_cluster", {
  data(seu, package = "leader")
  
  # Using too many PCs
  expect_warning(seurat_cluster(seu, 2101, knn = 20, res = 0.8))
  
  # check clustering output is correct
  cluster_out <- seurat_cluster(seu, 30, knn = 20, res = 1.2)
  
  expect_equal(dim(cluster_out), c(100, 2))
  expect_equal(length(unique(cluster_out$cluster_id)), 10)
})


test_that("sample_cells", {
  data(seu, package = "leader")
  
  # check clustering output is correct
  cluster_out <- seurat_cluster(seu, 30, knn = 20, res = 1.2)
  sampled_cells <- sample_cells(cluster_out, 20)
  
  expect_equal(dim(sampled_cells), c(20, 2))
  expect_equal(length(unique(cluster_out$cluster_id)),
               length(unique(sampled_cells$cluster_id)))
})


test_that("reweighting", {
  data(seu, package = "leader")
  markers <- get_markers(system.file("extdata", "markers.yml", package = "leader")) 
  cluster_out <- seurat_cluster(seu, 30, knn = 20, res = 1.2)
  
  reweighted_assignment <- reweighting(seu, markers, cluster_out)
  
  expect_equal(dim(reweighted_assignment), dim(cluster_out))
  expect_true(all(reweighted_assignment$cluster_id %in% names(markers)))
})


test_that("cluster_enrichment", {
  data(seu, package = "leader")
  markers <- get_markers(system.file("extdata", "markers.yml", package = "leader")) 
  cluster_out <- seurat_cluster(seu, 30, knn = 20, res = 1.2)
  
  expression <- seu@assays$RNA@counts |>
    as.matrix() |>
    t() |>
    as.data.frame() |>
    rownames_to_column('cell_id')
  
  unique_markers <- unlist(markers) |>
    unique()
  expression <- expression[, c('cell_id', unique_markers)] |>
    left_join(cluster_out, by = "cell_id") |>
    select(-cell_id)
  
  cluster_assignments <- cluster_enrich(markers, expression)
  
  expect_false(any(duplicated(cluster_assignments$cluster_id)))
  expect_true(all(cluster_assignments$predicted_cell_type %in% names(markers)))
  
  # there should only be one assigned cell type per cluster
  cell_type_per_cluster <- group_by(cluster_assignments, 
                                    predicted_cell_type, cluster_id) |> 
    tally() |> 
    filter(n > 1)
  expect_true(nrow(cell_type_per_cluster) == 0)
})


test_that("adaptive_reweighting", {
  # seu <- read in seurat object
  markers <- get_markers(system.file("extdata", "markers.yml", package = "leader")) #"extdata/markers.yml")

  # Selecting more cells than are present in dataset
  expect_error(adaptive_reweighting(seu, markers, 101))
  
  # selecting more than 70% of cells
  expect_warning(adaptive_reweighting(seu, markers, 71))
  
  ar_res <- adaptive_reweighting(seu, markers, 20)
  expect_equal(dim(ar_res), c(20, 2))
})




# Todo
# read in seurat object properly for seurat_clustering function
# computing too large a percentage, use a standard svd instead