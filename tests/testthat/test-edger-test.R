edger_inputs <- function() {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  list(
    counts = count_data[seq_len(400), c("gene_id", sample_metadata$sample_names)],
    sample_info = sample_metadata,
    column_geneid = "gene_id",
    group_column = "dex",
    group_numerator = "trt",
    group_denominator = "untrt",
    min_counts = 5
  )
}

build <- function(...) do.call(create_vista, c(edger_inputs(), list(...)))

test_that("edgeR runs the QL F-test on request and records the test used", {
  v_lrt <- suppressMessages(build(method = "edger", edger_test = "LRT"))
  v_ql <- build(method = "edger", edger_test = "QL")
  expect_identical(cutoffs(v_lrt)$edger_test, "LRT")
  expect_identical(cutoffs(v_ql)$edger_test, "QL")
  lrt <- comparisons(v_lrt)[[1]]
  ql <- comparisons(v_ql)[[1]]
  expect_true("stat" %in% names(ql))
  expect_false(isTRUE(all.equal(lrt$pvalue, ql$pvalue)))
  # QL treats the dispersion as estimated, so it is the more conservative test.
  expect_lte(sum(ql$regulation != "Other"), sum(lrt$regulation != "Other"))
})

# The vista_default_change messages a call emits, letting every other message
# through untouched.
notices <- function(expr) {
  seen <- character()
  withCallingHandlers(
    expr,
    vista_default_change = function(m) {
      seen <<- c(seen, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  seen
}

test_that("only calls that run edgeR on the default are told the default is changing", {
  expect_length(notices(build(method = "edger")), 1L)
  expect_length(notices(build(method = "both")), 1L)
  expect_length(notices(build(method = "edger", edger_test = NULL)), 1L)
  expect_length(notices(build(method = "edger", edger_test = "LRT")), 0L)
  expect_length(notices(build(method = "edger", edger_test = "QL")), 0L)
  expect_length(notices(build(method = "limma")), 0L)
  expect_length(notices(build(method = "deseq2")), 0L)
  expect_null(cutoffs(suppressMessages(build(method = "limma")))$edger_test)
  expect_length(notices(example_vista(method = "edger")), 0L)
})

test_that("the edgeR entry points outside create_vista() give the notice once", {
  inputs <- edger_inputs()
  direct <- inputs[c("counts", "sample_info", "column_geneid", "group_column", "group_numerator", "group_denominator", "min_counts")]
  expect_length(notices(do.call(run_edger_analysis, direct)), 1L)
  expect_length(notices(do.call(run_edger_analysis, c(direct, list(edger_test = "QL")))), 0L)
  expect_length(notices(do.call(benchmark_vista_equivalence, c(inputs, list(methods = "edger")))), 1L)
  expect_length(notices(do.call(benchmark_vista_equivalence, c(inputs, list(methods = "limma")))), 0L)
  expect_length(notices(do.call(validate_vista_deep, c(inputs, list(methods = "edger", error = FALSE)))), 1L)
})

test_that("negative values are refused and unmistakable non-counts draw a warning", {
  inputs <- edger_inputs()
  raw <- as.matrix(inputs$counts[, -1])
  tpm <- inputs; tpm$counts[, -1] <- sweep(raw + 0.5, 2, colSums(raw + 0.5), "/") * 1e6
  logged <- inputs; logged$counts[, -1] <- log2(raw + 1.5)
  negative <- inputs; negative$counts[1, 2] <- -3
  for (m in c("deseq2", "edger", "limma")) {
    expect_warning(try(suppressMessages(do.call(create_vista, c(tpm, list(method = m, edger_test = "LRT")))), silent = TRUE),
                   class = "vista_not_counts_warning")
    expect_warning(try(suppressMessages(do.call(create_vista, c(logged, list(method = m, edger_test = "LRT")))), silent = TRUE),
                   class = "vista_not_counts_warning")
    expect_error(do.call(create_vista, c(negative, list(method = m, edger_test = "LRT"))), "non-negative", info = m)
  }
})

test_that("estimated counts near one million reads and small fractional panels are accepted quietly", {
  inputs <- edger_inputs()
  raw <- as.matrix(inputs$counts[, -1])
  # Estimated counts from libraries of about a million reads: totals near, but
  # not exactly, one million.
  depth <- inputs
  depth$counts[, -1] <- sweep(raw + 0.3, 2, colSums(raw + 0.3), "/") * (1e6 + c(-6000, 2500, 7000, -1200, 4000, -3000, 500, 6500))
  expect_no_warning(suppressMessages(do.call(create_vista, c(depth, list(method = "edger", edger_test = "QL")))),
                    class = "vista_not_counts_warning")
  # A small panel of lowly expressed genes with fractional estimated counts.
  panel <- inputs
  low <- utils::head(which(rowMeans(raw) > 2 & rowMeans(raw) < 15), 40)
  expect_gte(length(low), 10)
  panel$counts <- panel$counts[low, ]
  panel$counts[, -1] <- as.matrix(panel$counts[, -1]) * 1.07
  panel$min_counts <- 1
  expect_no_warning(try(suppressMessages(do.call(create_vista, c(panel, list(method = "limma")))), silent = TRUE),
                    class = "vista_not_counts_warning")
})

test_that("every engine keeps the same genes from non-integer counts", {
  inputs <- edger_inputs()
  inputs$counts[, -1] <- as.matrix(inputs$counts[, -1]) * 1.13
  # A gene just below the threshold must be dropped by DESeq2 too, not rounded in.
  inputs$counts[1, -1] <- 9.6
  inputs$min_counts <- 10
  v_cons <- suppressMessages(do.call(create_vista, c(inputs, list(method = "both", edger_test = "QL"))))
  v_edger <- suppressMessages(do.call(create_vista, c(inputs, list(method = "both", edger_test = "QL", result_source = "edger"))))
  expect_identical(sort(rownames(v_cons)), sort(rownames(v_edger)))
  # counts(v) holds DESeq2's rounded counts whichever source is active.
  expect_identical(counts(v_cons)[rownames(v_edger), ], counts(v_edger))
  expect_true(all(counts(v_edger) == round(counts(v_edger))))
  expect_false(inputs$counts$gene_id[1] %in% rownames(v_cons))
  expect_false(anyNA(comparisons(v_cons, source = "consensus")[[1]]$pvalue_edger))
})

test_that("the consensus uses the requested edgeR test", {
  v <- build(method = "both", edger_test = "QL")
  ql_alone <- comparisons(build(method = "edger", edger_test = "QL"))[[1]]
  e <- comparisons(v, source = "edger")[[1]]
  expect_equal(e$pvalue, ql_alone$pvalue[match(e$gene_id, ql_alone$gene_id)])
})

test_that("the equivalence check reproduces the QL test from an independent reference", {
  report <- do.call(benchmark_vista_equivalence, c(edger_inputs(), list(methods = "edger", edger_test = "QL")))
  expect_true(report$valid)
  expect_identical(report$comparison_summary$max_abs_pvalue_diff, 0)
})
