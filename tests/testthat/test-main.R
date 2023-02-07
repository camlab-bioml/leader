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


test_that("seurat_clustering", {
  seu <- readRDS("data/example-seurat.rds")
  seu2 <- FindVariableFeatures(seu, selection.method = 'vst', nfeatures = 20)

  expect_warning(seurat_cluster(seu2, n_of_pcs = 20, knn = 10, res = 0.8))
})
