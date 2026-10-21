test_that("benchmark_vista_equivalence matches standalone DESeq2 and edgeR", {
  inputs <- make_benchmark_inputs()

  report <- do.call(
    benchmark_vista_equivalence,
    c(inputs, list(methods = c("deseq2", "edger"), tolerance = 1e-8, return_plots = FALSE))
  )

  expect_true(report$valid)
  expect_true(all(report$comparison_summary$deg_sets_identical))
  expect_true(all(report$comparison_summary$regulation_identical))
  expect_true(all(report$comparison_summary$norm_counts_identical))
  expect_true(all(report$comparison_summary$log2fc_within_tolerance))
  expect_true(all(report$comparison_summary$pvalue_within_tolerance))
  expect_true(all(report$comparison_summary$padj_within_tolerance))
  expect_true(all(report$visual_summary$pass))
})

test_that("benchmark_vista_equivalence matches standalone limma", {
  skip_if_not_installed("limma")
  inputs <- make_benchmark_inputs()

  report <- do.call(
    benchmark_vista_equivalence,
    c(inputs, list(methods = "limma", tolerance = 1e-8, return_plots = FALSE))
  )

  expect_true(report$valid)
  expect_true(all(report$comparison_summary$deg_sets_identical))
  expect_true(all(report$comparison_summary$regulation_identical))
  expect_true(all(report$comparison_summary$norm_counts_identical))
  expect_true(all(report$comparison_summary$log2fc_within_tolerance))
  expect_true(all(report$comparison_summary$pvalue_within_tolerance))
  expect_true(all(report$comparison_summary$padj_within_tolerance))
  expect_true(all(report$visual_summary$pass))
})

test_that("validate_vista_deep returns a passing report on benchmark inputs", {
  inputs <- make_benchmark_inputs()

  report <- do.call(
    validate_vista_deep,
    c(inputs, list(methods = "deseq2", tolerance = 1e-8, error = FALSE))
  )

  expect_true(report$valid)
})

# The reference side must not share code with what it checks, or a defect in
# that code is reproduced on both sides and passes. These tests enforce the
# separation and show that the benchmark catches defects injected into the
# helpers VISTA itself uses.

reference_functions <- c(
  ".run_reference_backend", ".reference_sample_ids", ".reference_inputs",
  ".reference_contrast", ".reference_regulation", ".reference_fit_deseq2",
  ".reference_fit_edger", ".reference_fit_limma",
  ".ma_payload_from_reference", ".volcano_payload_from_reference",
  ".deg_count_payload_from_reference", ".pca_payload_from_reference"
)

# Every string constant in a function body, so a call by name
# (do.call("helper"), match.fun("helper"), get("helper")) is seen too.
string_literals <- function(expr) {
  if (is.character(expr)) return(expr)
  if (is.call(expr) || is.expression(expr) || is.pairlist(expr) || is.list(expr)) {
    return(unlist(lapply(as.list(expr), string_literals), use.names = FALSE))
  }
  character()
}

# Names reached as VISTA::name or VISTA:::name, which findGlobals() reports
# only as a call to `::` or `:::`.
qualified_names <- function(expr, pkg = "VISTA") {
  if (is.call(expr) &&
      (identical(expr[[1]], as.name("::")) || identical(expr[[1]], as.name(":::"))) &&
      identical(as.character(expr[[2]]), pkg)) {
    return(as.character(expr[[3]]))
  }
  if (is.call(expr) || is.expression(expr) || is.pairlist(expr) || is.list(expr)) {
    return(unlist(lapply(as.list(expr), qualified_names, pkg = pkg), use.names = FALSE))
  }
  character()
}

test_that("reference side calls no VISTA code", {
  skip_if_not_installed("codetools")
  ns <- asNamespace("VISTA")
  vista_functions <- Filter(function(nm) is.function(get(nm, envir = ns)), ls(ns, all.names = TRUE))
  forbidden <- setdiff(vista_functions, reference_functions)
  for (fn in reference_functions) {
    f <- get(fn, envir = ns)
    code <- c(body(f), as.list(formals(f)))
    # Globals used as calls or as values (lapply(x, helper), h <- helper),
    # names reached through VISTA:: / VISTA:::, and names given as strings.
    used <- c(
      codetools::findGlobals(f, merge = TRUE),
      qualified_names(code),
      string_literals(code)
    )
    expect_identical(intersect(used, forbidden), character(0), info = fn)
  }
})

benchmark_with_defect <- function(method, ..., inputs = make_benchmark_inputs()) {
  report <- do.call(
    benchmark_vista_equivalence,
    c(inputs, list(methods = method, ...))
  )
  report$comparison_summary
}

test_that("a sign error in VISTA's contrast builder is caught", {
  original <- VISTA:::.edgeR_group_contrast
  local_mocked_bindings(.edgeR_group_contrast = function(...) -original(...))
  cs <- benchmark_with_defect("edger")
  expect_false(cs$log2fc_within_tolerance)
  expect_false(cs$all_checks_pass)
})

test_that("swapped Up/Down calls in VISTA are caught", {
  original <- VISTA:::.categorize_deg_results
  local_mocked_bindings(.categorize_deg_results = function(...) {
    out <- original(...)
    out$regulation <- c(Up = "Down", Down = "Up", Other = "Other")[out$regulation]
    out
  })
  cs <- benchmark_with_defect("deseq2")
  expect_false(cs$regulation_identical)
  expect_false(cs$all_checks_pass)
})

test_that("values moved under the wrong gene labels are caught", {
  original <- VISTA:::.align_de_to_counts
  local_mocked_bindings(.align_de_to_counts = function(df, ref_rn, ...) {
    out <- original(df, ref_rn, ...)
    value_cols <- setdiff(names(out), "gene_id")
    out[value_cols] <- out[rev(seq_len(nrow(out))), value_cols, drop = FALSE]
    out
  })
  cs <- benchmark_with_defect("edger")
  expect_false(cs$log2fc_within_tolerance)
  expect_false(cs$all_checks_pass)
})

test_that("samples matched to the wrong metadata rows are caught", {
  original <- VISTA:::.normalize_sample_info
  # Normalisation runs more than once on the way in, so the injected defect
  # must survive repetition: a descending sort does, a reversal undoes itself.
  local_mocked_bindings(.normalize_sample_info = function(...) {
    out <- original(...)
    out$sample_names <- sort(out$sample_names, decreasing = TRUE)
    out
  })
  cs <- benchmark_with_defect("limma")
  expect_false(cs$all_checks_pass)
})

test_that("a covariate dropped from VISTA's design is caught", {
  original <- VISTA:::.resolve_design_formula
  local_mocked_bindings(.resolve_design_formula = function(sample_info, group_column, covariates = NULL, ...) {
    original(sample_info, group_column, covariates = NULL, ...)
  })
  cs <- benchmark_with_defect("edger", covariates = "cell")
  expect_false(cs$all_checks_pass)
})

test_that("benchmark agrees with covariates and a contrast between two non-reference groups", {
  inputs <- make_benchmark_inputs(n_genes = 400)
  si <- inputs$sample_info
  si$grp <- si$cond_long
  si$grp[si$cond_long == "treatment1"][1:2] <- "treatment2"
  inputs$sample_info <- si
  inputs$group_column <- "grp"
  inputs$group_numerator <- c("treatment1", "treatment2")
  inputs$group_denominator <- c("control", "treatment1")

  report <- do.call(
    benchmark_vista_equivalence,
    c(inputs, list(methods = c("deseq2", "edger", "limma"), covariates = "cell"))
  )
  expect_true(report$valid)
  expect_true(all(report$comparison_summary$all_columns_within_tolerance))
})

test_that("DE tables keep a p-value the engine did not compute as NA", {
  de <- data.frame(
    gene_id = c("g1", "g2", "g3"),
    log2FoldChange = c(2, -1, 0.5),
    pvalue = c(0.001, NA, 0.4),
    padj = c(0.01, NA, 0.6)
  )
  out <- VISTA:::.tidy_de_results(de, rowname_col = "gene_id")
  expect_identical(out$pvalue, c(0.001, NA, 0.4))
  expect_identical(out$padj, c(0.01, NA, 0.6))
})

# One sample with an extreme count makes DESeq2 flag the gene as an outlier
# (Cook's distance) and report no p-value for it.
inputs_with_outlier <- function() {
  inputs <- make_benchmark_inputs(n_genes = 300)
  sample <- inputs$sample_info$sample_names[1]
  gene <- which(rowSums(inputs$counts[, -1]) > 500)[1]
  inputs$counts[[sample]][gene] <- inputs$counts[[sample]][gene] * 50
  inputs
}

test_that("an untested gene's missing p-value is kept, and replacing it is caught", {
  inputs <- inputs_with_outlier()
  report <- do.call(benchmark_vista_equivalence, c(inputs, list(methods = "deseq2")))
  expect_true(any(is.na(comparisons(report$methods$deseq2$vista)[[1]]$pvalue)))
  expect_true(report$valid)

  # The 1.1.5 behaviour: a missing p-value stored as 1.
  original <- VISTA:::.tidy_de_results
  local_mocked_bindings(.tidy_de_results = function(...) {
    out <- original(...)
    if ("pvalue" %in% names(out)) out$pvalue[is.na(out$pvalue)] <- 1
    out
  })
  cs <- benchmark_with_defect("deseq2", inputs = inputs)
  expect_false(cs$pvalue_within_tolerance)
  expect_identical(cs$max_abs_pvalue_diff, Inf)
})

test_that("a design with no DEGs at all passes", {
  inputs <- modifyList(make_benchmark_inputs(), list(log2fc_cutoff = 10))
  report <- do.call(benchmark_vista_equivalence, c(inputs, list(methods = c("deseq2", "edger", "limma"))))
  expect_true(report$valid)
})

test_that("data-dependent formula terms, an explicit intercept and a Run-only sample sheet pass", {
  for (f in c("~ factor(cell) + cond_long", "~ scale(avgLength) + cond_long", "~ 1 + cell + cond_long")) {
    inputs <- modifyList(make_benchmark_inputs(), list(design_formula = f))
    report <- do.call(benchmark_vista_equivalence, c(inputs, list(methods = c("edger", "limma"))))
    expect_true(report$valid, info = f)
  }

  inputs <- make_benchmark_inputs()
  inputs$sample_info$sample_names <- NULL
  report <- do.call(benchmark_vista_equivalence, c(inputs, list(methods = "edger")))
  expect_true(report$valid)
})

test_that("a comparison dropped from DEG-count plot data is caught", {
  original <- VISTA:::.collect_deg_count_data
  local_mocked_bindings(.collect_deg_count_data = function(...) {
    out <- original(...)
    out[out$sample_comparisons != out$sample_comparisons[1], , drop = FALSE]
  })
  report <- do.call(benchmark_vista_equivalence, c(make_benchmark_inputs(), list(methods = "edger")))
  expect_false(report$valid)
  expect_false(report$visual_summary$pass[report$visual_summary$visual == "deg_count"])
})

test_that("the MA check reads what get_ma_plot() draws", {
  original <- VISTA:::.align_de_to_counts
  local_mocked_bindings(.align_de_to_counts = function(df, ref_rn, ...) {
    out <- original(df, ref_rn, ...)
    value_cols <- setdiff(names(out), "gene_id")
    out[value_cols] <- out[rev(seq_len(nrow(out))), value_cols, drop = FALSE]
    out
  })
  report <- do.call(benchmark_vista_equivalence, c(make_benchmark_inputs(), list(methods = "edger")))
  expect_false(report$visual_summary$pass[report$visual_summary$visual == "ma"])
})
