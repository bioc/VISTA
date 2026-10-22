# Deprecation warnings are throttled for users (see .vista_deprecate_arg), which
# would make them invisible to the second and later assertions in a run. Tests
# must see every occurrence.
options(vista.deprecation_frequency = "always")

make_small_vista <- function() {
  example_vista()
}

# PCA and MDS axes come from an eigendecomposition, and the sign of an
# eigenvector is arbitrary: prcomp() and cmdscale() may return v or -v for the
# same input, and which one you get depends on the LAPACK/BLAS build and on the
# order of the input. A runner image update flipped both MDS axes and turned
# the suite red without anything in VISTA changing; later, shuffling the genes
# passed to get_pca_plot() flipped both PCA axes on Linux CI only. Magnitudes
# are identical, so canonicalise each axis to a fixed orientation before
# comparing -- a real change still moves the numbers.
canon_axis <- function(v) {
  if (!is.numeric(v) || !length(v)) return(v)
  i <- which.max(abs(v))                       # largest coordinate is unambiguous
  if (length(i) == 1L && !is.na(v[i]) && v[i] < 0) -v else v
}

make_benchmark_inputs <- function(n_genes = 150) {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  target_groups <- c("control", "treatment1")
  sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, , drop = FALSE]
  counts_subset <- count_data[seq_len(n_genes), c("gene_id", sample_subset$sample_names)]

  list(
    counts = counts_subset,
    sample_info = sample_subset,
    column_geneid = "gene_id",
    group_column = "cond_long",
    group_numerator = "treatment1",
    group_denominator = "control",
    min_counts = 5,
    min_replicates = 1,
    log2fc_cutoff = 1,
    pval_cutoff = 0.05,
    p_value_type = "padj"
  )
}
