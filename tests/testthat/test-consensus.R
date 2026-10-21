consensus_object <- function(mode) {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  counts <- count_data[seq_len(600), c("gene_id", sample_metadata$sample_names)]
  # One extreme count gives DESeq2 an untested gene (pvalue NA).
  sample <- sample_metadata$sample_names[1]
  gene <- which(rowSums(counts[, -1]) > 500)[1]
  counts[[sample]][gene] <- counts[[sample]][gene] * 50
  create_vista(
    counts = counts, sample_info = sample_metadata, column_geneid = "gene_id",
    group_column = "cond_long", group_numerator = "treatment1", group_denominator = "control",
    method = "both", consensus_mode = mode, result_source = "consensus",
    min_counts = 5, log2fc_cutoff = 1, pval_cutoff = 0.05
  )
}

no_evidence <- function(p) replace(p, is.na(p), 1)

test_that("intersection consensus: padj is the less-significant backend's for every gene", {
  cons <- comparisons(consensus_object("intersection"), source = "consensus")[[1]]
  expect_true(any(cons$support %in% c("deseq2_only", "edger_only")))
  expect_equal(cons$padj, pmax(no_evidence(cons$padj_deseq2), no_evidence(cons$padj_edger)))
  # So filtering the table on its own columns cannot return single-backend genes.
  passing <- cons$padj <= 0.05
  expect_true(all(cons$support[passing] %in% c("both", "discordant", "none")))
})

test_that("union consensus: single-backend calls carry the calling backend's values", {
  cons <- comparisons(consensus_object("union"), source = "consensus")[[1]]
  d_only <- cons$support == "deseq2_only"
  e_only <- cons$support == "edger_only"
  expect_true(any(d_only | e_only))
  expect_equal(cons$padj[d_only], no_evidence(cons$padj_deseq2[d_only]))
  expect_equal(cons$padj[e_only], no_evidence(cons$padj_edger[e_only]))
  expect_equal(cons$log2fc[d_only], cons$log2fc_deseq2[d_only])
  expect_true(all(cons$regulation[d_only | e_only] %in% c("Up", "Down")))
})

test_that("consensus keeps each backend's own p-values, missing ones included", {
  v <- consensus_object("intersection")
  cons <- comparisons(v, source = "consensus")[[1]]
  d <- comparisons(v, source = "deseq2")[[1]]
  e <- comparisons(v, source = "edger")[[1]]
  expect_true(anyNA(d$pvalue))
  expect_identical(cons$pvalue_deseq2, d$pvalue[match(cons$gene_id, d$gene_id)])
  expect_identical(cons$pvalue_edger, e$pvalue[match(cons$gene_id, e$gene_id)])
})

test_that("provenance records the VISTA, R and engine versions", {
  # A freshly built object; example_vista() returns a stored one.
  prov <- S4Vectors::metadata(consensus_object("intersection"))$provenance
  expect_identical(prov$vista_version, as.character(utils::packageVersion("VISTA")))
  expect_identical(prov$r_version, R.version.string)
  expect_named(prov$engine_versions, c("DESeq2", "edgeR", "limma"))
})

test_that("example_vista() default object carries full provenance", {
  prov <- S4Vectors::metadata(example_vista())$provenance
  expect_true(all(c("constructor", "created_at", "vista_version", "r_version", "engine_versions") %in% names(prov)))
  expect_named(prov$engine_versions, c("DESeq2", "edgeR", "limma"))
})
