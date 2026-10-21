#' Benchmark VISTA against standalone differential-expression backends
#'
#' Run VISTA and independent DESeq2, edgeR, and limma pipelines with the same
#' filtering, design, contrast definitions, and DEG thresholds, then compare
#' the retained genes, normalized matrices, every numeric column the engine
#' reports, DEG calls, and critical plot inputs.
#'
#' @details
#' The reference pipelines are written independently of the code they check.
#' They call only base R, \pkg{stats} and the engines' own functions -- never
#' VISTA's internal helpers -- and re-implement VISTA's documented contract
#' from scratch:
#'
#' * **Samples** are matched to count columns by a `sample_names` column, else
#'   by non-default rownames, else by the first of `sample_name`, `sample`,
#'   `sample_id`, `Run`, `run`, `Sample` and `SampleName` whose values name
#'   count columns (the rule documented for `sample_info` in [create_vista()]).
#' * **Filtering** treats missing counts as zero and keeps genes whose total
#'   reaches `min_counts` and that reach `min_counts` in at least
#'   `min_replicates` samples. For DESeq2, non-integer counts (estimated counts
#'   from tximport or RSEM) are then rounded, as DESeq2 requires.
#' * **Tests** are DESeq2's Wald test via `results()`, edgeR's `glmLRT()` or,
#'   with `edger_test = "QL"`, `glmQLFit()` and `glmQLFTest()`, and limma's
#'   `eBayes()` after `voom()`.
#' * **Design** is `~ covariates + group` for DESeq2 and the same terms without
#'   an intercept for edgeR and limma-voom, or `design_formula` when supplied.
#' * **Contrasts** for edgeR and limma are the difference between the design
#'   rows of two samples that differ only in group, so they do not depend on
#'   how design columns are named or which level is the reference.
#' * **DEG calls** are `Up` when log2FC >= `log2fc_cutoff` and the selected
#'   p-value <= `pval_cutoff`, `Down` symmetrically, and `Other` otherwise,
#'   including when the engine reports no p-value.
#'
#' VISTA's DE tables and normalized counts are read through its exported
#' accessors ([comparisons()], [norm_counts()]). The MA, PCA and DEG-count
#' checks read the data [get_ma_plot()], [get_pca_plot()] and the DEG-count
#' plots draw; the volcano check rebuilds its inputs from the stored table.
#' Rows are matched by gene identifier, not position, and the set of
#' comparisons must agree as well. A
#' defect in sample matching, design construction, contrast signs, DEG calling
#' or row alignment therefore shows up as a mismatch instead of being
#' reproduced on both sides. The two sides share only the engines and the
#' contract above: the benchmark checks that VISTA runs DESeq2, edgeR and limma
#' faithfully, not that those packages are correct.
#'
#' @param counts Raw counts (matrix/data.frame) with a gene-id column and sample columns.
#' @param sample_info Data frame with sample metadata.
#' @param column_geneid Column name in `counts` that contains gene identifiers.
#' @param group_column Column in `sample_info` used to group samples.
#' @param group_numerator Character vector of numerator groups for pairwise comparisons.
#' @param group_denominator Character vector of denominator groups.
#' @param methods Character vector of backends to benchmark. Any subset of
#'   `c("deseq2", "edger", "limma")`.
#' @param min_counts Minimum count a gene must reach in a sample to count as
#'   detected there; also the minimum row total.
#' @param min_replicates Minimum number of samples, across the whole
#'   experiment, in which a gene must reach `min_counts`.
#' @param log2fc_cutoff Absolute log2 fold-change threshold for DEG calling.
#' @param pval_cutoff P-value (or adjusted p-value) threshold.
#' @param p_value_type Either `"padj"` or `"pvalue"`.
#' @param covariates Optional character vector of additional sample_info columns.
#' @param design_formula Optional model formula (or formula string) including `group_column`.
#' @param tolerance Numeric tolerance used for floating-point comparisons.
#' @param return_plots Logical; if `TRUE`, return paired VISTA/reference plots
#'   for MA, volcano, DEG count, and PCA views.
#' @param edger_test The edgeR test, `"LRT"` (the default until VISTA 1.4.0) or
#'   `"QL"`; see [create_vista()]. VISTA and the reference both use it.
#'
#' @return A list with fields `valid`, `comparison_summary`, `visual_summary`,
#'   and `methods`. `comparison_summary` has one row per method and comparison
#'   recording whether the retained gene set, normalized counts, every numeric
#'   column the engine reports (`columns_compared`), regulation calls and DEG
#'   sets agree, with the largest absolute difference for the main columns. A
#'   `*_within_tolerance` field is `NA` when the engine does not report that
#'   column (edgeR reports `logCPM`, not `baseMean`). Each element of `methods`
#'   contains the VISTA object, the reference results, structural validation
#'   output, self-consistency checks, and optional plot objects.
#' @examples
#' v <- example_vista()
#' si <- as.data.frame(sample_info(v))
#' data("count_data", package = "VISTA")
#' count_subset <- count_data[seq_len(500), c("gene_id", si$sample_names), drop = FALSE]
#'
#' bm <- benchmark_vista_equivalence(
#'   counts = count_subset,
#'   sample_info = si,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   methods = "limma",
#'   min_counts = 5,
#'   min_replicates = 1
#' )
#'
#' bm$comparison_summary
#'
#' \donttest{
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' target_groups <- c("control", "treatment1")
#' sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
#' count_subset <- count_data[seq_len(150), c("gene_id", sample_subset$sample_names)]
#'
#' bm <- benchmark_vista_equivalence(
#'   counts = count_subset,
#'   sample_info = sample_subset,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   methods = c("deseq2", "edger"),
#'   min_counts = 5,
#'   min_replicates = 1
#' )
#'
#' bm$comparison_summary
#' bm$visual_summary
#' }
#' @export
benchmark_vista_equivalence <- function(counts,
                                        sample_info,
                                        column_geneid,
                                        group_column,
                                        group_numerator,
                                        group_denominator,
                                        methods = c("deseq2", "edger", "limma"),
                                        min_counts = 10,
                                        min_replicates = 1,
                                        log2fc_cutoff = 1,
                                        pval_cutoff = 0.05,
                                        p_value_type = "padj",
                                        covariates = NULL,
                                        design_formula = NULL,
                                        tolerance = 1e-8,
                                        return_plots = FALSE,
                                        edger_test = c("LRT", "QL")) {
  methods <- match.arg(methods, c("deseq2", "edger", "limma"), several.ok = TRUE)
  p_value_type <- match.arg(p_value_type, c("padj", "pvalue"))
  if ((missing(edger_test) || is.null(edger_test)) && "edger" %in% methods) {
    .vista_default_change_notice("benchmark_vista_equivalence", "edger_test")
  }
  edger_test <- match.arg(edger_test)

  method_results <- stats::setNames(vector("list", length(methods)), methods)
  comparison_summary <- vector("list", length(methods))
  visual_summary <- vector("list", length(methods))

  for (i in seq_along(methods)) {
    method <- methods[[i]]

    # Both sides receive the inputs exactly as the user supplied them; neither
    # sees the other's preprocessing.
    vista_obj <- create_vista(
      counts = counts,
      sample_info = sample_info,
      column_geneid = column_geneid,
      group_column = group_column,
      group_numerator = group_numerator,
      group_denominator = group_denominator,
      method = method,
      min_counts = min_counts,
      min_replicates = min_replicates,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type,
      covariates = covariates,
      design_formula = design_formula,
      validate = TRUE,
      edger_test = edger_test
    )

    reference <- .run_reference_backend(
      method = method,
      counts = counts,
      sample_info = sample_info,
      column_geneid = column_geneid,
      group_column = group_column,
      group_numerator = group_numerator,
      group_denominator = group_denominator,
      covariates = covariates,
      design_formula = design_formula,
      min_counts = min_counts,
      min_replicates = min_replicates,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      p_value_type = p_value_type,
      edger_test = edger_test
    )

    structural_validation <- validate_vista(vista_obj, level = "full", error = FALSE)
    self_consistency <- .validate_vista_self_consistency(
      x = vista_obj,
      p_value_type = p_value_type,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff
    )

    method_check <- .compare_vista_and_reference(
      vista_obj = vista_obj,
      reference = reference,
      method = method,
      p_value_type = p_value_type,
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff,
      tolerance = tolerance,
      return_plots = return_plots
    )

    method_results[[method]] <- c(
      list(
        vista = vista_obj,
        reference = reference,
        structural_validation = structural_validation,
        self_consistency = self_consistency
      ),
      method_check
    )

    comparison_summary[[i]] <- method_check$comparison_summary |>
      dplyr::mutate(
        structural_valid = structural_validation$valid,
        self_consistency_valid = length(self_consistency$issues) == 0,
        all_checks_pass = .data$all_checks_pass &
          .data$structural_valid &
          .data$self_consistency_valid
      )

    visual_summary[[i]] <- method_check$visual_summary |>
      dplyr::mutate(
        structural_valid = structural_validation$valid,
        self_consistency_valid = length(self_consistency$issues) == 0,
        pass = .data$pass & .data$structural_valid & .data$self_consistency_valid
      )
  }

  comparison_summary <- dplyr::bind_rows(comparison_summary)
  visual_summary <- dplyr::bind_rows(visual_summary)
  valid <- all(comparison_summary$all_checks_pass) && all(visual_summary$pass)

  list(
    valid = valid,
    comparison_summary = comparison_summary,
    visual_summary = visual_summary,
    methods = method_results
  )
}

#' Deep validation of VISTA differential-expression fidelity
#'
#' This validator combines `validate_vista()` with backend-to-backend numerical
#' equivalence checks against standalone DESeq2, edgeR, and limma runs.
#'
#' @inheritParams benchmark_vista_equivalence
#' @param error Logical; if `TRUE`, abort when any discrepancy is detected.
#'
#' @return Invisibly returns the full benchmark report.
#' @examples
#' v <- example_vista()
#' si <- as.data.frame(sample_info(v))
#' data("count_data", package = "VISTA")
#' count_subset <- count_data[seq_len(500), c("gene_id", si$sample_names), drop = FALSE]
#'
#' report <- validate_vista_deep(
#'   counts = count_subset,
#'   sample_info = si,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   methods = "limma",
#'   min_counts = 5,
#'   min_replicates = 1,
#'   error = FALSE
#' )
#'
#' report$valid
#'
#' \donttest{
#' data("count_data", package = "VISTA")
#' data("sample_metadata", package = "VISTA")
#'
#' target_groups <- c("control", "treatment1")
#' sample_subset <- sample_metadata[sample_metadata$cond_long %in% target_groups, ]
#' count_subset <- count_data[seq_len(150), c("gene_id", sample_subset$sample_names)]
#'
#' validate_vista_deep(
#'   counts = count_subset,
#'   sample_info = sample_subset,
#'   column_geneid = "gene_id",
#'   group_column = "cond_long",
#'   group_numerator = "treatment1",
#'   group_denominator = "control",
#'   methods = c("deseq2", "edger"),
#'   min_counts = 5,
#'   min_replicates = 1
#' )
#' }
#' @export
validate_vista_deep <- function(counts,
                                sample_info,
                                column_geneid,
                                group_column,
                                group_numerator,
                                group_denominator,
                                methods = c("deseq2", "edger", "limma"),
                                min_counts = 10,
                                min_replicates = 1,
                                log2fc_cutoff = 1,
                                pval_cutoff = 0.05,
                                p_value_type = "padj",
                                covariates = NULL,
                                design_formula = NULL,
                                tolerance = 1e-8,
                                return_plots = FALSE,
                                error = TRUE,
                                edger_test = c("LRT", "QL")) {
  methods <- match.arg(methods, c("deseq2", "edger", "limma"), several.ok = TRUE)
  if ((missing(edger_test) || is.null(edger_test)) && "edger" %in% methods) {
    .vista_default_change_notice("validate_vista_deep", "edger_test")
  }
  edger_test <- match.arg(edger_test)
  report <- benchmark_vista_equivalence(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    group_numerator = group_numerator,
    group_denominator = group_denominator,
    methods = methods,
    min_counts = min_counts,
    min_replicates = min_replicates,
    log2fc_cutoff = log2fc_cutoff,
    pval_cutoff = pval_cutoff,
    p_value_type = p_value_type,
    covariates = covariates,
    design_formula = design_formula,
    tolerance = tolerance,
    return_plots = return_plots,
    edger_test = edger_test
  )

  if (!report$valid) {
    issues <- c(
      .collect_benchmark_issues(report$comparison_summary, type = "comparison"),
      .collect_benchmark_issues(report$visual_summary, type = "visual"),
      unlist(lapply(report$methods, function(x) x$self_consistency$issues), use.names = FALSE),
      unlist(lapply(report$methods, function(x) x$structural_validation$issues), use.names = FALSE)
    )
    issues <- unique(stats::na.omit(issues))
    bullets <- stats::setNames(as.character(issues), rep("x", length(issues)))

    if (isTRUE(error)) {
      cli::cli_abort(c("Deep VISTA validation failed.", bullets))
    } else {
      cli::cli_warn(c("Deep VISTA validation reported issues.", bullets))
    }
  }

  invisible(report)
}

#' @keywords internal
.collect_benchmark_issues <- function(df, type = c("comparison", "visual")) {
  type <- match.arg(type)
  if (is.null(df) || !nrow(df)) {
    return(character())
  }

  if (type == "comparison") {
    failed <- df[!df$all_checks_pass, , drop = FALSE]
    if (!nrow(failed)) {
      return(character())
    }
    apply(failed, 1, function(row) {
      paste0(
        row[["method"]], " / ", row[["comparison"]], ": ",
        "comparison_sets_identical=", row[["comparison_sets_identical"]],
        ", gene_sets_identical=", row[["gene_sets_identical"]],
        ", regulation_identical=", row[["regulation_identical"]],
        ", deg_sets_identical=", row[["deg_sets_identical"]],
        ", norm_counts_identical=", row[["norm_counts_identical"]],
        ", all_columns_within_tolerance=", row[["all_columns_within_tolerance"]],
        " (", row[["columns_compared"]], ")",
        ", log2fc_within_tolerance=", row[["log2fc_within_tolerance"]],
        ", pvalue_within_tolerance=", row[["pvalue_within_tolerance"]],
        ", padj_within_tolerance=", row[["padj_within_tolerance"]]
      )
    })
  } else {
    failed <- df[!df$pass, , drop = FALSE]
    if (!nrow(failed)) {
      return(character())
    }
    apply(failed, 1, function(row) {
      comp <- if (!is.na(row[["comparison"]]) && nzchar(row[["comparison"]])) {
        paste0(" / ", row[["comparison"]])
      } else {
        ""
      }
      paste0(
        row[["method"]], comp, " / ", row[["visual"]], ": ",
        row[["detail"]]
      )
    })
  }
}

# Independent reference pipelines --------------------------------------------
#
# Nothing in this section may call VISTA's analysis helpers: sample matching,
# design construction, contrasts, DEG calling, or table standardisation and
# alignment. A reference that shared that code would reproduce any defect in
# it, and the comparison would pass. Only base R, stats and the engines' own
# functions are used; the documented contract is re-implemented here from
# scratch (see the Details of ?benchmark_vista_equivalence).

#' @keywords internal
.run_reference_backend <- function(method,
                                   counts,
                                   sample_info,
                                   column_geneid,
                                   group_column,
                                   group_numerator,
                                   group_denominator,
                                   covariates = NULL,
                                   design_formula = NULL,
                                   min_counts = 10,
                                   min_replicates = 1,
                                   log2fc_cutoff = 1,
                                   pval_cutoff = 0.05,
                                   p_value_type = "padj",
                                   edger_test = "LRT") {
  if (method %in% c("edger", "limma") && !requireNamespace("edgeR", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg edgeR} is required for {method} benchmarking.")
  }
  if (identical(method, "limma") && !requireNamespace("limma", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg limma} is required for limma benchmarking.")
  }

  inputs <- .reference_inputs(
    counts = counts,
    sample_info = sample_info,
    column_geneid = column_geneid,
    group_column = group_column,
    covariates = covariates,
    design_formula = design_formula,
    min_counts = min_counts,
    min_replicates = min_replicates,
    intercept = identical(method, "deseq2"),
    round_counts = identical(method, "deseq2")
  )
  fit <- switch(
    method,
    deseq2 = .reference_fit_deseq2(inputs),
    edger = .reference_fit_edger(inputs, edger_test = edger_test),
    limma = .reference_fit_limma(inputs)
  )

  comparisons <- list()
  for (i in seq_along(group_numerator)) {
    num <- as.character(group_numerator[[i]])
    den <- as.character(group_denominator[[i]])
    if (identical(num, den)) {
      next
    }
    tbl <- fit$test(num, den)
    tbl$regulation <- .reference_regulation(
      log2fc = tbl$log2fc,
      p = tbl[[p_value_type]],
      log2fc_cutoff = log2fc_cutoff,
      pval_cutoff = pval_cutoff
    )
    comparisons[[paste(num, den, sep = "_VS_")]] <- tbl
  }

  deg_summary <- lapply(comparisons, function(tbl) {
    data.frame(
      regulation = c("Up", "Down"),
      n = c(sum(tbl$regulation == "Up"), sum(tbl$regulation == "Down")),
      stringsAsFactors = FALSE
    )
  })

  list(
    norm_counts = fit$norm_counts,
    sample_info = inputs$col_data,
    group_column = group_column,
    comparisons = comparisons,
    deg_summary = deg_summary
  )
}

#' @keywords internal
.reference_sample_ids <- function(sample_info, count_columns) {
  if ("sample_names" %in% names(sample_info)) {
    return(as.character(sample_info$sample_names))
  }
  rn <- rownames(sample_info)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(sample_info))))) {
    return(as.character(rn))
  }
  # Otherwise the first of the documented ID columns (see ?create_vista,
  # `sample_info`) whose values name count columns.
  documented <- c("sample_name", "sample", "sample_id", "Run", "run", "Sample", "SampleName")
  for (nm in intersect(documented, names(sample_info))) {
    v <- as.character(sample_info[[nm]])
    if (!anyNA(v) && all(nzchar(v)) && !anyDuplicated(v) && all(v %in% count_columns)) {
      return(v)
    }
  }
  cli::cli_abort(c(
    "Cannot match {.arg sample_info} rows to count columns.",
    "i" = "Add a {.field sample_names} column naming each sample's count column."
  ))
}

#' @keywords internal
.reference_inputs <- function(counts,
                              sample_info,
                              column_geneid,
                              group_column,
                              covariates = NULL,
                              design_formula = NULL,
                              min_counts = 10,
                              min_replicates = 1,
                              intercept = TRUE,
                              round_counts = FALSE) {
  counts <- as.data.frame(counts, stringsAsFactors = FALSE)
  sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
  ids <- .reference_sample_ids(sample_info, setdiff(names(counts), column_geneid))
  if (anyNA(ids) || anyDuplicated(ids) || !all(ids %in% names(counts))) {
    cli::cli_abort("Sample identifiers must be unique and name columns of {.arg counts}.")
  }
  genes <- as.character(counts[[column_geneid]])
  if (anyNA(genes) || anyDuplicated(genes)) {
    cli::cli_abort("Gene identifiers in {.field {column_geneid}} must be unique and non-missing.")
  }

  mat <- as.matrix(counts[, ids, drop = FALSE])
  storage.mode(mat) <- "double"
  mat[is.na(mat)] <- 0
  dimnames(mat) <- list(genes, ids)
  keep <- rowSums(mat) >= min_counts & rowSums(mat >= min_counts) >= min_replicates
  mat <- mat[keep, , drop = FALSE]
  if (round_counts) {
    mat <- round(mat)
  }

  col_data <- sample_info
  col_data$sample_names <- NULL
  rownames(col_data) <- ids
  col_data[[group_column]] <- factor(col_data[[group_column]])

  if (is.null(design_formula)) {
    terms <- unique(c(setdiff(covariates, group_column), group_column))
    design_formula <- stats::reformulate(terms, intercept = intercept)
  } else {
    if (is.character(design_formula)) {
      design_formula <- stats::as.formula(design_formula)
    }
    if (!intercept && attr(stats::terms(design_formula), "intercept") == 1L) {
      design_formula <- stats::update(design_formula, ~ . - 1)
    }
  }

  vars <- all.vars(design_formula)
  if (!group_column %in% vars) {
    cli::cli_abort("{.arg design_formula} must include {.field {group_column}}.")
  }
  if (!all(vars %in% names(col_data))) {
    cli::cli_abort("Design variables missing from {.arg sample_info}: {.val {setdiff(vars, names(col_data))}}")
  }
  for (vv in vars) {
    if (is.character(col_data[[vv]])) {
      col_data[[vv]] <- factor(col_data[[vv]])
    }
  }

  list(
    counts = mat,
    col_data = col_data,
    design_formula = design_formula,
    group_column = group_column
  )
}

#' @keywords internal
.reference_contrast <- function(inputs, design, numerator, denominator) {
  # The design built with every sample in the numerator group, minus the
  # design with every sample in the denominator group, all other columns as
  # observed: each row of the difference is the numerator-minus-denominator
  # contrast, whatever the coding or the reference level. Covariate terms that
  # depend on the data (factor(), scale(), cut(), I(x / sd(x))) see the same
  # values as in the fit, so they evaluate exactly as they did there.
  col_data <- inputs$col_data
  group <- inputs$group_column
  lv <- levels(col_data[[group]])
  if (!all(c(numerator, denominator) %in% lv)) {
    cli::cli_abort("Groups {.val {setdiff(c(numerator, denominator), lv)}} are not levels of {.field {group}}.")
  }
  all_in <- function(level) {
    d <- col_data
    d[[group]] <- factor(rep(level, nrow(d)), levels = lv)
    stats::model.matrix(inputs$design_formula, data = d)
  }
  rows <- all_in(numerator) - all_in(denominator)
  if (!identical(colnames(rows), colnames(design))) {
    cli::cli_abort("Could not rebuild the design columns for the {numerator} vs {denominator} contrast.")
  }
  # Additive designs give the same contrast for every sample; a difference
  # between rows would mean the group effect depends on the covariates.
  if (any(abs(sweep(rows, 2, rows[1L, ])) > 1e-12)) {
    cli::cli_abort("The {numerator} vs {denominator} contrast depends on covariates; only additive designs are supported.")
  }
  rows[1L, ]
}

#' @keywords internal
.reference_regulation <- function(log2fc, p, log2fc_cutoff, pval_cutoff) {
  sig <- !is.na(log2fc) & !is.na(p) & p <= pval_cutoff
  ifelse(
    sig & log2fc >= log2fc_cutoff, "Up",
    ifelse(sig & log2fc <= -log2fc_cutoff, "Down", "Other")
  )
}

#' @keywords internal
.reference_fit_deseq2 <- function(inputs) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = inputs$counts,
    colData = inputs$col_data,
    design = inputs$design_formula
  )
  dds <- DESeq2::DESeq(dds, quiet = TRUE)

  test <- function(num, den) {
    res <- DESeq2::results(dds, contrast = c(inputs$group_column, num, den))
    data.frame(
      gene_id = rownames(res),
      baseMean = res$baseMean,
      log2fc = res$log2FoldChange,
      lfcSE = res$lfcSE,
      stat = res$stat,
      pvalue = res$pvalue,
      padj = res$padj,
      row.names = rownames(res),
      stringsAsFactors = FALSE
    )
  }

  list(norm_counts = DESeq2::counts(dds, normalized = TRUE), test = test)
}

#' @keywords internal
.reference_fit_edger <- function(inputs, edger_test = "LRT") {
  design <- stats::model.matrix(inputs$design_formula, data = inputs$col_data)
  dge <- edgeR::DGEList(counts = inputs$counts, group = inputs$col_data[[inputs$group_column]])
  dge <- edgeR::calcNormFactors(dge)
  if (identical(edger_test, "QL")) {
    # edgeR v4 quasi-likelihood pipeline, as the edgeR user's guide runs it.
    fit <- edgeR::glmQLFit(dge, design)
    run_test <- function(contrast) edgeR::glmQLFTest(fit, contrast = contrast)
    stat_col <- "F"
  } else {
    dge <- edgeR::estimateDisp(dge, design)
    fit <- edgeR::glmFit(dge, design)
    run_test <- function(contrast) edgeR::glmLRT(fit, contrast = contrast)
    stat_col <- "LR"
  }

  test <- function(num, den) {
    contrast <- .reference_contrast(inputs, design, num, den)
    tab <- edgeR::topTags(run_test(contrast), n = Inf, sort.by = "none")$table
    data.frame(
      gene_id = rownames(tab),
      log2fc = tab$logFC,
      logCPM = tab$logCPM,
      stat = tab[[stat_col]],
      pvalue = tab$PValue,
      padj = tab$FDR,
      row.names = rownames(tab),
      stringsAsFactors = FALSE
    )
  }

  list(norm_counts = edgeR::cpm(dge), test = test)
}

#' @keywords internal
.reference_fit_limma <- function(inputs) {
  design <- stats::model.matrix(inputs$design_formula, data = inputs$col_data)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = inputs$counts))
  fit <- limma::lmFit(limma::voom(dge, design = design, plot = FALSE), design)

  test <- function(num, den) {
    contrast <- matrix(
      .reference_contrast(inputs, design, num, den),
      ncol = 1,
      dimnames = list(colnames(design), "contrast")
    )
    fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrasts = contrast))
    tab <- limma::topTable(fit2, coef = 1, number = Inf, sort.by = "none")
    data.frame(
      gene_id = rownames(tab),
      baseMean = tab$AveExpr,
      log2fc = tab$logFC,
      stat = tab$t,
      pvalue = tab$P.Value,
      padj = tab$adj.P.Val,
      B = tab$B,
      # VISTA's limma tables keep limma's own names beside pvalue/padj.
      P.Value = tab$P.Value,
      adj.P.Val = tab$adj.P.Val,
      row.names = rownames(tab),
      stringsAsFactors = FALSE
    )
  }

  list(norm_counts = edgeR::cpm(dge), test = test)
}

# Comparison ------------------------------------------------------------------
#
# VISTA is read through its exported accessors, and every match between the
# two sides is by identifier. Positional comparison would let a label/value
# misalignment on VISTA's side pass whenever both tables share an order.

#' @keywords internal
.compare_vista_and_reference <- function(vista_obj,
                                         reference,
                                         method,
                                         p_value_type = "padj",
                                         log2fc_cutoff = 1,
                                         pval_cutoff = 0.05,
                                         tolerance = 1e-8,
                                         return_plots = FALSE) {
  vista_norm <- as.matrix(norm_counts(vista_obj))
  ref_norm <- as.matrix(reference$norm_counts)
  gene_sets_identical <- nrow(vista_norm) == nrow(ref_norm) &&
    setequal(rownames(vista_norm), rownames(ref_norm))
  norm_check <- .compare_numeric_matrix(vista_norm, ref_norm, tolerance = tolerance)

  vista_comps <- comparisons(vista_obj)
  comp_names <- names(reference$comparisons)

  comparison_rows <- vector("list", length(comp_names))
  visual_rows <- vector("list", length(comp_names) + 2L)
  plot_bundle <- if (isTRUE(return_plots)) list(ma = list(), volcano = list()) else NULL

  for (i in seq_along(comp_names)) {
    comp_name <- comp_names[[i]]
    has_comp <- comp_name %in% names(vista_comps)

    comparison_rows[[i]] <- .compare_de_tables(
      vista_tbl = if (has_comp) vista_comps[[comp_name]] else NULL,
      ref_tbl = reference$comparisons[[comp_name]],
      method = method,
      comparison = comp_name,
      norm_check = norm_check,
      gene_sets_identical = gene_sets_identical,
      tolerance = tolerance
    )$summary

    if (!has_comp) {
      missing_row <- tibble::tibble(
        method = method,
        comparison = comp_name,
        visual = c("ma", "volcano"),
        pass = FALSE,
        detail = "comparison missing from VISTA object"
      )
      visual_rows[[i]] <- missing_row[1, ]
      visual_rows[[length(comp_names) + i]] <- missing_row[2, ]
      next
    }

    visual_rows[[i]] <- .compare_visual_payload(
      method = method,
      comparison = comp_name,
      visual = "ma",
      vista_df = .ma_payload_from_vista(vista_obj, comp_name),
      ref_df = .ma_payload_from_reference(reference, comp_name),
      key = "gene_id",
      tolerance = tolerance
    )

    visual_rows[[length(comp_names) + i]] <- .compare_visual_payload(
      method = method,
      comparison = comp_name,
      visual = "volcano",
      vista_df = .volcano_payload_from_vista(
        vista_obj,
        comp_name,
        p_value_type = p_value_type,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff
      ),
      ref_df = .volcano_payload_from_reference(
        reference,
        comp_name,
        p_value_type = p_value_type,
        log2fc_cutoff = log2fc_cutoff,
        pval_cutoff = pval_cutoff
      ),
      key = "gene_id",
      tolerance = tolerance
    )

    if (isTRUE(return_plots)) {
      plot_bundle$ma[[comp_name]] <- list(
        vista = get_ma_plot(vista_obj, sample_comparison = comp_name),
        reference = .plot_reference_ma(.ma_payload_from_reference(reference, comp_name))
      )
      plot_bundle$volcano[[comp_name]] <- list(
        vista = .plot_reference_volcano(
          .volcano_payload_from_vista(
            vista_obj,
            comp_name,
            p_value_type = p_value_type,
            log2fc_cutoff = log2fc_cutoff,
            pval_cutoff = pval_cutoff
          )
        ),
        reference = .plot_reference_volcano(
          .volcano_payload_from_reference(
            reference,
            comp_name,
            p_value_type = p_value_type,
            log2fc_cutoff = log2fc_cutoff,
            pval_cutoff = pval_cutoff
          )
        )
      )
    }
  }

  visual_rows[[length(comp_names) * 2 + 1L]] <- .compare_visual_payload(
    method = method,
    comparison = NA_character_,
    visual = "deg_count",
    vista_df = .deg_count_payload_from_vista(vista_obj),
    ref_df = .deg_count_payload_from_reference(reference),
    key = c("sample_comparisons", "regulation"),
    tolerance = tolerance
  )
  visual_rows[[length(comp_names) * 2 + 2L]] <- .compare_visual_payload(
    method = method,
    comparison = NA_character_,
    visual = "pca",
    vista_df = .pca_payload_from_vista(vista_obj),
    ref_df = .pca_payload_from_reference(reference, vista_obj),
    key = "sample",
    tolerance = tolerance,
    allow_sign_flip = TRUE
  )

  if (isTRUE(return_plots)) {
    plot_bundle$deg_count <- list(
      vista = get_deg_count_barplot(vista_obj, label = TRUE),
      reference = .plot_reference_deg_count(.deg_count_payload_from_reference(reference))
    )
    plot_bundle$pca <- list(
      vista = get_pca_plot(vista_obj),
      reference = .plot_reference_pca(.pca_payload_from_reference(reference, vista_obj), vista_obj)
    )
  }

  # A comparison VISTA holds that the reference never asked for is a mismatch
  # too, not something to skip.
  comparison_summary <- dplyr::bind_rows(comparison_rows) |>
    dplyr::mutate(
      comparison_sets_identical = setequal(names(vista_comps), comp_names),
      all_checks_pass = .data$all_checks_pass & .data$comparison_sets_identical
    )

  list(
    comparison_summary = comparison_summary,
    visual_summary = dplyr::bind_rows(visual_rows),
    plots = plot_bundle
  )
}

#' @keywords internal
.compare_de_tables <- function(vista_tbl,
                               ref_tbl,
                               method,
                               comparison,
                               norm_check,
                               gene_sets_identical,
                               tolerance = 1e-8) {
  # Every numeric column the engine reports is compared, not just the four
  # that feed the summary fields.
  value_cols <- setdiff(names(ref_tbl), c("gene_id", "regulation"))

  vista_tbl <- if (is.null(vista_tbl)) {
    data.frame(gene_id = character(), stringsAsFactors = FALSE)
  } else {
    as.data.frame(vista_tbl, stringsAsFactors = FALSE)
  }
  vista_ids <- if ("gene_id" %in% names(vista_tbl)) {
    as.character(vista_tbl$gene_id)
  } else {
    rownames(vista_tbl)
  }
  idx <- match(ref_tbl$gene_id, vista_ids)
  rows_matched <- !anyNA(idx) && nrow(vista_tbl) == nrow(ref_tbl)

  checks <- lapply(value_cols, function(col) {
    if (!rows_matched || !col %in% names(vista_tbl)) {
      return(list(identical = FALSE, max_abs_diff = Inf))
    }
    .compare_numeric_vectors(vista_tbl[[col]][idx], ref_tbl[[col]], tolerance = tolerance)
  })
  names(checks) <- value_cols
  # A numeric column VISTA reports that the reference does not cannot be
  # checked, so it counts against the comparison rather than being skipped.
  vista_numeric <- names(vista_tbl)[vapply(vista_tbl, is.numeric, logical(1))]
  unchecked <- setdiff(vista_numeric, value_cols)
  for (col in unchecked) {
    checks[[col]] <- list(identical = FALSE, max_abs_diff = Inf)
  }
  within <- function(col) if (col %in% value_cols) checks[[col]]$identical else NA
  max_diff <- function(col) if (col %in% value_cols) checks[[col]]$max_abs_diff else NA_real_
  all_columns <- all(vapply(checks, `[[`, logical(1), "identical"))

  vista_reg <- if ("regulation" %in% names(vista_tbl)) as.character(vista_tbl$regulation) else character(nrow(vista_tbl))
  regulation_identical <- rows_matched && identical(vista_reg[idx], as.character(ref_tbl$regulation))
  up_identical <- setequal(vista_ids[vista_reg %in% "Up"], ref_tbl$gene_id[ref_tbl$regulation %in% "Up"])
  down_identical <- setequal(vista_ids[vista_reg %in% "Down"], ref_tbl$gene_id[ref_tbl$regulation %in% "Down"])
  gene_sets <- gene_sets_identical && rows_matched

  summary <- tibble::tibble(
    method = method,
    comparison = comparison,
    n_genes = nrow(ref_tbl),
    gene_sets_identical = gene_sets,
    up_genes_identical = up_identical,
    down_genes_identical = down_identical,
    deg_sets_identical = up_identical && down_identical,
    regulation_identical = regulation_identical,
    norm_counts_identical = norm_check$identical,
    baseMean_within_tolerance = within("baseMean"),
    log2fc_within_tolerance = within("log2fc"),
    pvalue_within_tolerance = within("pvalue"),
    padj_within_tolerance = within("padj"),
    all_columns_within_tolerance = all_columns,
    columns_compared = paste(value_cols, collapse = ","),
    columns_not_in_reference = paste(unchecked, collapse = ","),
    max_abs_norm_counts_diff = norm_check$max_abs_diff,
    max_abs_baseMean_diff = max_diff("baseMean"),
    max_abs_log2fc_diff = max_diff("log2fc"),
    max_abs_pvalue_diff = max_diff("pvalue"),
    max_abs_padj_diff = max_diff("padj"),
    all_checks_pass = gene_sets &&
      regulation_identical &&
      up_identical &&
      down_identical &&
      norm_check$identical &&
      all_columns
  )

  list(summary = summary)
}

#' @keywords internal
.compare_visual_payload <- function(method,
                                    comparison,
                                    visual,
                                    vista_df,
                                    ref_df,
                                    key,
                                    tolerance = 1e-8,
                                    allow_sign_flip = FALSE) {
  # A column the reference draws that VISTA's payload lacks is a mismatch, not
  # something to skip. Columns only VISTA carries are not compared.
  missing_cols <- setdiff(names(ref_df), names(vista_df))
  if (length(missing_cols)) {
    return(tibble::tibble(
      method = method,
      comparison = if (is.na(comparison)) NA_character_ else comparison,
      visual = visual,
      pass = FALSE,
      detail = paste0("missing from VISTA: ", paste(missing_cols, collapse = ", "))
    ))
  }
  key_of <- function(df) do.call(paste, c(lapply(df[key], as.character), sep = "\r"))
  vista_keys <- key_of(vista_df)
  ref_keys <- key_of(ref_df)
  # Keys must identify rows on both sides; a duplicate could hide a row of the
  # other side from the comparison.
  idx <- match(vista_keys, ref_keys)
  if (anyNA(idx) || nrow(vista_df) != nrow(ref_df) ||
      anyDuplicated(vista_keys) || anyDuplicated(ref_keys)) {
    return(tibble::tibble(
      method = method,
      comparison = if (is.na(comparison)) NA_character_ else comparison,
      visual = visual,
      pass = FALSE,
      detail = "rows differ between VISTA and reference"
    ))
  }
  ref_df <- ref_df[idx, , drop = FALSE]

  if (allow_sign_flip) {
    ref_df <- .align_pca_signs(vista_df, ref_df)
  }

  shared_cols <- intersect(names(vista_df), names(ref_df))
  vista_df <- vista_df[, shared_cols, drop = FALSE]
  ref_df <- ref_df[, shared_cols, drop = FALSE]

  numeric_cols <- shared_cols[vapply(vista_df, is.numeric, logical(1))]
  other_cols <- setdiff(shared_cols, numeric_cols)
  numeric_checks <- lapply(numeric_cols, function(col) {
    .compare_numeric_vectors(vista_df[[col]], ref_df[[col]], tolerance = tolerance)
  })
  names(numeric_checks) <- numeric_cols

  non_numeric_identical <- TRUE
  if (length(other_cols)) {
    non_numeric_identical <- all(vapply(other_cols, function(col) {
      identical(as.character(vista_df[[col]]), as.character(ref_df[[col]]))
    }, logical(1)))
  }

  numeric_identical <- all(vapply(numeric_checks, `[[`, logical(1), "identical"))
  max_abs_diff <- 0
  if (length(numeric_checks)) {
    max_abs_diff <- max(vapply(numeric_checks, `[[`, numeric(1), "max_abs_diff"))
  }

  detail <- if (numeric_identical && non_numeric_identical) {
    "identical"
  } else {
    paste0(
      "non_numeric_identical=", non_numeric_identical,
      ", numeric_identical=", numeric_identical,
      ", max_abs_diff=", signif(max_abs_diff, 4)
    )
  }

  tibble::tibble(
    method = method,
    comparison = if (is.na(comparison)) NA_character_ else comparison,
    visual = visual,
    pass = numeric_identical && non_numeric_identical,
    detail = detail
  )
}

#' @keywords internal
.compare_numeric_vectors <- function(x, y, tolerance = 1e-8) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) != length(y)) {
    return(list(identical = FALSE, max_abs_diff = Inf))
  }
  both_na <- is.na(x) & is.na(y)
  same <- !is.na(x) & !is.na(y) & x == y
  close <- is.finite(x) & is.finite(y) & abs(x - y) <= tolerance
  equal <- both_na | same | close

  # NA against a value, or an infinity against anything else, is reported as
  # an infinite difference rather than hidden as 0.
  diff <- ifelse(both_na | same, 0, abs(x - y))
  diff[is.na(diff)] <- Inf

  list(
    identical = all(equal),
    max_abs_diff = if (length(diff)) max(diff) else 0
  )
}

#' @keywords internal
.compare_numeric_matrix <- function(x, y, tolerance = 1e-8) {
  # Matched by dimnames: equal values in a different gene order must pass, and
  # equal-looking values under the wrong labels must not.
  labelled <- function(m) {
    !is.null(rownames(m)) && !is.null(colnames(m)) &&
      !anyDuplicated(rownames(m)) && !anyDuplicated(colnames(m))
  }
  if (!labelled(x) || !labelled(y)) {
    return(list(identical = FALSE, max_abs_diff = Inf))
  }
  ri <- match(rownames(y), rownames(x))
  ci <- match(colnames(y), colnames(x))
  if (anyNA(ri) || anyNA(ci) || nrow(x) != nrow(y) || ncol(x) != ncol(y)) {
    return(list(identical = FALSE, max_abs_diff = Inf))
  }
  .compare_numeric_vectors(
    as.vector(x[ri, ci, drop = FALSE]),
    as.vector(y),
    tolerance = tolerance
  )
}

#' @keywords internal
.validate_vista_self_consistency <- function(x,
                                             p_value_type = "padj",
                                             log2fc_cutoff = 1,
                                             pval_cutoff = 0.05) {
  comps <- comparisons(x)
  summaries <- deg_summary(x)
  ref_rn <- rownames(norm_counts(x))
  issues <- character()

  for (comp_name in names(comps)) {
    tbl <- .align_de_to_counts(
      df = .tidy_de_results(comps[[comp_name]], rowname_col = "gene_id"),
      ref_rn = ref_rn,
      warn_missing = FALSE
    )
    p_col <- .p_col_for_cutoff(tbl, p_value_type = p_value_type)
    p_vals <- if (!is.na(p_col)) .safe_numeric(tbl[[p_col]], default = 1) else rep(1, nrow(tbl))
    expected_reg <- ifelse(
      !is.na(tbl$log2fc) & tbl$log2fc >= log2fc_cutoff & p_vals <= pval_cutoff,
      "Up",
      ifelse(
        !is.na(tbl$log2fc) & tbl$log2fc <= -log2fc_cutoff & p_vals <= pval_cutoff,
        "Down",
        "Other"
      )
    )

    if (!identical(as.character(tbl$regulation), expected_reg)) {
      issues <- c(issues, sprintf("Comparison '%s' has a regulation column inconsistent with the configured thresholds.", comp_name))
    }

    expected_summary <- tbl |>
      dplyr::filter(.data$regulation %in% c("Up", "Down")) |>
      dplyr::count(.data$regulation, name = "n") |>
      dplyr::mutate(regulation = factor(.data$regulation, levels = c("Up", "Down"))) |>
      dplyr::arrange(.data$regulation)

    observed_summary <- summaries[[comp_name]] |>
      as.data.frame(stringsAsFactors = FALSE) |>
      dplyr::filter(.data$regulation %in% c("Up", "Down")) |>
      dplyr::mutate(regulation = factor(.data$regulation, levels = c("Up", "Down"))) |>
      dplyr::arrange(.data$regulation)

    expected_summary <- tidyr::complete(
      expected_summary,
      regulation = factor(c("Up", "Down"), levels = c("Up", "Down")),
      fill = list(n = 0)
    ) |>
      dplyr::arrange(.data$regulation)

    observed_summary <- tidyr::complete(
      observed_summary,
      regulation = factor(c("Up", "Down"), levels = c("Up", "Down")),
      fill = list(n = 0)
    ) |>
      dplyr::arrange(.data$regulation)

    if (!identical(as.integer(expected_summary$n), as.integer(observed_summary$n))) {
      issues <- c(issues, sprintf("Comparison '%s' has a DEG summary that does not match its regulation calls.", comp_name))
    }
  }

  list(valid = length(issues) == 0, issues = unique(issues))
}

#' @keywords internal
.ma_payload_from_vista <- function(x, comparison) {
  # Read from what get_ma_plot() draws, not re-derived from the stored table,
  # so a defect in the plot function's own data handling is caught too.
  pd <- as.data.frame(get_ma_plot(x, sample_comparison = comparison)$data, stringsAsFactors = FALSE)
  data.frame(
    gene_id = as.character(pd$gene_id),
    baseMean = pd$baseMean,
    log10_baseMean = pd$log10_baseMean,
    log2fc = pd$log2fc,
    regulation = as.character(pd$regulation),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.ma_payload_from_reference <- function(reference, comparison) {
  tbl <- reference$comparisons[[comparison]]
  base_mean <- unname(rowMeans(reference$norm_counts)[tbl$gene_id])
  data.frame(
    gene_id = tbl$gene_id,
    baseMean = base_mean,
    log10_baseMean = log10(base_mean + 1),
    log2fc = tbl$log2fc,
    regulation = tbl$regulation,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.volcano_payload_from_vista <- function(x,
                                        comparison,
                                        p_value_type = "padj",
                                        log2fc_cutoff = 1,
                                        pval_cutoff = 0.05) {
  tbl <- .align_de_to_counts(
    df = .tidy_de_results(comparisons(x)[[comparison]], rowname_col = "gene_id"),
    ref_rn = rownames(norm_counts(x)),
    warn_missing = FALSE
  )
  .volcano_payload_from_table(tbl, p_value_type, log2fc_cutoff, pval_cutoff)
}

#' @keywords internal
.volcano_payload_from_reference <- function(reference,
                                            comparison,
                                            p_value_type = "padj",
                                            log2fc_cutoff = 1,
                                            pval_cutoff = 0.05) {
  tbl <- reference$comparisons[[comparison]]
  p_vals <- tbl[[p_value_type]]
  # The volcano draws a gene the engine did not test at p = 1.
  p_vals[!is.finite(p_vals)] <- 1
  data.frame(
    gene_id = tbl$gene_id,
    log2fc = tbl$log2fc,
    pvalue = p_vals,
    neg_log10_p = -log10(p_vals + 1e-300),
    significant = abs(tbl$log2fc) >= log2fc_cutoff & p_vals <= pval_cutoff,
    regulation = tbl$regulation,
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
.volcano_payload_from_table <- function(tbl,
                                        p_value_type = "padj",
                                        log2fc_cutoff = 1,
                                        pval_cutoff = 0.05) {
  p_col <- .p_col_for_cutoff(tbl, p_value_type = p_value_type)
  p_vals <- if (!is.na(p_col)) .safe_numeric(tbl[[p_col]], default = 1) else rep(1, nrow(tbl))
  neg_log10 <- -log10(p_vals + 1e-300)
  tbl |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      log2fc = .data$log2fc,
      pvalue = p_vals,
      neg_log10_p = neg_log10,
      significant = abs(.data$log2fc) >= log2fc_cutoff & p_vals <= pval_cutoff,
      regulation = as.character(.data$regulation)
    )
}

#' @keywords internal
.deg_count_payload_from_vista <- function(x) {
  df <- suppressWarnings(.collect_deg_count_data(x))
  .standardize_deg_count_payload(df)
}

#' @keywords internal
.deg_count_payload_from_reference <- function(reference) {
  # Built here rather than through .standardize_deg_count_payload(), which
  # formats VISTA's side: shared formatting could drop the same rows on both.
  rows <- lapply(names(reference$deg_summary), function(comp_name) {
    tbl <- reference$deg_summary[[comp_name]]
    data.frame(
      sample_comparisons = comp_name,
      regulation = c("Up", "Down"),
      n = c(sum(tbl$n[tbl$regulation == "Up"]), sum(tbl$n[tbl$regulation == "Down"])),
      stringsAsFactors = FALSE
    )
  })
  if (!length(rows)) {
    return(data.frame(sample_comparisons = character(), regulation = character(), n = numeric()))
  }
  do.call(rbind, rows)
}

#' @keywords internal
.standardize_deg_count_payload <- function(df) {
  template <- tibble::tibble(
    sample_comparisons = character(),
    regulation = character(),
    n = numeric()
  )

  if (is.null(df) || !nrow(df)) {
    return(template)
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE)
  required <- c("sample_comparisons", "regulation", "n")
  if (!all(required %in% names(df))) {
    return(template)
  }

  df <- df |>
    dplyr::transmute(
      sample_comparisons = as.character(.data$sample_comparisons),
      regulation = as.character(.data$regulation),
      n = as.numeric(.data$n)
    )

  if (!nrow(df)) {
    return(template)
  }

  df$sample_comparisons <- factor(df$sample_comparisons, levels = unique(df$sample_comparisons))
  df$regulation <- factor(df$regulation, levels = c("Up", "Down"))
  tidyr::complete(
    df,
    sample_comparisons,
    regulation,
    fill = list(n = 0)
  ) |>
    dplyr::arrange(.data$sample_comparisons, .data$regulation)
}

#' @keywords internal
.pca_payload_from_vista <- function(x) {
  # Read from what get_pca_plot() draws with its defaults.
  pd <- as.data.frame(get_pca_plot(x)$data, stringsAsFactors = FALSE)
  pd[order(pd$sample), , drop = FALSE]
}

#' @keywords internal
.pca_payload_from_reference <- function(reference, vista_obj) {
  # As get_pca_plot() documents: normalized counts, genes with any variance,
  # centred and scaled.
  mat <- as.matrix(reference$norm_counts)
  mat <- mat[apply(mat, 1, stats::var) > 0, , drop = FALSE]
  pca <- stats::prcomp(t(mat), center = TRUE, scale. = TRUE)
  samples <- rownames(pca$x)
  pca_df <- data.frame(
    sample = samples,
    PC1 = unname(pca$x[, 1]),
    PC2 = unname(pca$x[, 2]),
    stringsAsFactors = FALSE
  )
  # Every sample-sheet column, not only the group: the plots annotate and
  # colour by them, so they are compared against VISTA's colData too.
  meta <- as.data.frame(reference$sample_info, stringsAsFactors = FALSE)[samples, , drop = FALSE]
  for (col in setdiff(names(meta), names(pca_df))) {
    pca_df[[col]] <- as.character(meta[[col]])
  }
  pca_df[order(pca_df$sample), , drop = FALSE]
}

#' @keywords internal
.align_pca_signs <- function(vista_df, ref_df) {
  out <- ref_df
  for (pc in c("PC1", "PC2")) {
    if (!pc %in% names(vista_df) || !pc %in% names(ref_df)) {
      next
    }
    corr <- suppressWarnings(stats::cor(vista_df[[pc]], ref_df[[pc]], use = "pairwise.complete.obs"))
    if (is.finite(corr) && corr < 0) {
      out[[pc]] <- -out[[pc]]
    }
  }
  out
}

#' @keywords internal
.plot_reference_ma <- function(df,
                               colors = c(Up = "#a40000", Down = "#16317d", Other = "gray70")) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$log10_baseMean, y = .data$log2fc, color = .data$regulation)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.2) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "log10(baseMean + 1)", y = "log2 fold change", color = "Regulation")
}

#' @keywords internal
.plot_reference_volcano <- function(df,
                                    colors = c(Up = "#a40000", Down = "#007e2f", Other = "grey70")) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$log2fc, y = .data$neg_log10_p, color = .data$regulation)) +
    ggplot2::geom_point(size = 1, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "log2 fold change", y = "-log10(p-value)", color = "Regulation")
}

#' @keywords internal
.plot_reference_deg_count <- function(df,
                                      colors = c(Up = "red4", Down = "blue4")) {
  ggplot2::ggplot(df, ggplot2::aes(x = .data$sample_comparisons, y = .data$n, fill = .data$regulation)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(x = NULL, y = "Gene Count", fill = "Regulation")
}

#' @keywords internal
.plot_reference_pca <- function(df, vista_obj) {
  group_col <- tryCatch(S4Vectors::metadata(vista_obj)$group$column, error = function(e) NULL)
  if (is.null(group_col) || !group_col %in% names(df)) {
    group_col <- names(df)[!names(df) %in% c("sample", "PC1", "PC2")][1]
  }
  color_vals <- group_colors(vista_obj)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data[[group_col]])) +
    ggplot2::geom_point(size = 10) +
    ggplot2::scale_color_manual(values = color_vals) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "PCA Plot", x = "PC1", y = "PC2", color = group_col)
}
