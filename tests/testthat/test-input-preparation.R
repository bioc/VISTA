test_that("read_vista_counts standardizes matrix-like inputs", {
  data("count_data", package = "VISTA", envir = environment())

  counts_in <- count_data[seq_len(12), ]
  cnt <- read_vista_counts(
    counts_in,
    format = "matrix",
    gene_id_column = "gene_id",
    verbose = FALSE
  )

  expect_named(cnt, c("counts", "row_data", "column_geneid", "sample_names", "sample_name_map", "input_format", "report"))
  expect_identical(cnt$column_geneid, "gene_id")
  expect_identical(cnt$input_format, "matrix")
  expect_identical(cnt$counts$gene_id, counts_in$gene_id)
  expect_identical(cnt$sample_names, colnames(counts_in)[-1])
  expect_true(all(vapply(cnt$counts[-1], is.numeric, logical(1))))
})

test_that("read_vista_counts repairs file-path sample columns conservatively", {
  data("count_data", package = "VISTA", envir = environment())

  sample_cols <- colnames(count_data)[colnames(count_data) != "gene_id"][seq_len(4)]
  counts_in <- count_data[seq_len(6), c("gene_id", sample_cols), drop = FALSE]
  names(counts_in)[2:5] <- c(
    "/proj/run/03_alignment/HET_1_U/HET_1_U_star_alignAligned.sortedByCoord.out.bam",
    "/proj/run/03_alignment/HET_1_ovary/HET_1_ovary_star_alignAligned.sortedByCoord.out.bam",
    "/proj/run/03_alignment/WT_1_U/WT_1_U_star_alignAligned.sortedByCoord.out.bam",
    "/proj/run/03_alignment/WT_1_ovary/WT_1_ovary_star_alignAligned.sortedByCoord.out.bam"
  )

  cnt <- read_vista_counts(
    counts_in,
    format = "matrix",
    gene_id_column = "gene_id",
    verbose = FALSE
  )

  expect_identical(
    cnt$sample_names,
    c("HET_1_U", "HET_1_ovary", "WT_1_U", "WT_1_ovary")
  )
  expect_true(all(cnt$sample_name_map$original != cnt$sample_name_map$repaired))
})

test_that("read_vista_counts repairs generic quantification file names via parent directory", {
  data("count_data", package = "VISTA", envir = environment())

  sample_cols <- colnames(count_data)[colnames(count_data) != "gene_id"][seq_len(2)]
  counts_in <- count_data[seq_len(4), c("gene_id", sample_cols), drop = FALSE]
  names(counts_in)[2:3] <- c(
    "/proj/salmon/sample_A/quant.sf",
    "/proj/salmon/sample_B/quant.sf"
  )

  cnt <- read_vista_counts(
    counts_in,
    format = "matrix",
    gene_id_column = "gene_id",
    verbose = FALSE
  )

  expect_identical(cnt$sample_names, c("sample_A", "sample_B"))
})

test_that("read_vista_counts strips lane and read suffixes when unique", {
  data("count_data", package = "VISTA", envir = environment())

  sample_cols <- colnames(count_data)[colnames(count_data) != "gene_id"][seq_len(2)]
  counts_in <- count_data[seq_len(4), c("gene_id", sample_cols), drop = FALSE]
  names(counts_in)[2:3] <- c(
    "sampleA_S1_L001_R1_001.fastq.gz",
    "sampleB_S2_L002_R1_001.fastq.gz"
  )

  cnt <- read_vista_counts(
    counts_in,
    format = "matrix",
    gene_id_column = "gene_id",
    verbose = FALSE
  )

  expect_identical(cnt$sample_names, c("sampleA", "sampleB"))
})

test_that("read_vista_counts parses featureCounts-style tables", {
  data("count_data", package = "VISTA", envir = environment())

  fc <- data.frame(
    Geneid = count_data$gene_id[seq_len(8)],
    Chr = rep("chr1", 8),
    Start = seq(100, by = 100, length.out = 8),
    End = seq(150, by = 100, length.out = 8),
    Strand = "+",
    Length = 50,
    sampleA = unname(count_data[[2]][seq_len(8)]),
    sampleB = unname(count_data[[3]][seq_len(8)]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  tmp <- tempfile(fileext = ".txt")
  writeLines(c(
    "# Program:featureCounts",
    paste(colnames(fc), collapse = "\t"),
    apply(fc, 1, paste, collapse = "\t")
  ), con = tmp)

  cnt <- read_vista_counts(tmp, format = "auto", verbose = FALSE)

  expect_identical(cnt$input_format, "featurecounts")
  expect_identical(colnames(cnt$counts), c("gene_id", "sampleA", "sampleB"))
  expect_true(all(c("Chr", "Start", "End", "Strand", "Length") %in% colnames(cnt$row_data)))
})

test_that("read_vista_counts parses STAR per-sample files", {
  star_lines <- function(scale = 1) {
    c(
      "N_unmapped\t1\t1\t1",
      "N_multimapping\t2\t2\t2",
      "N_noFeature\t3\t3\t3",
      "N_ambiguous\t4\t4\t4",
      paste("geneA", 10 * scale, 11 * scale, 12 * scale, sep = "\t"),
      paste("geneB", 20 * scale, 21 * scale, 22 * scale, sep = "\t")
    )
  }

  f1 <- tempfile(pattern = "sampleA_", fileext = ".ReadsPerGene.out.tab")
  f2 <- tempfile(pattern = "sampleB_", fileext = ".ReadsPerGene.out.tab")
  writeLines(star_lines(1), f1)
  writeLines(star_lines(2), f2)

  cnt <- read_vista_counts(
    c(f1, f2),
    format = "star",
    count_column = "unstranded",
    verbose = FALSE
  )

  expect_identical(cnt$input_format, "star")
  expect_identical(cnt$counts$gene_id, c("geneA", "geneB"))
  expect_equal(as.numeric(cnt$counts[1, -1]), c(10, 20))
})

test_that("read_vista_counts parses HTSeq per-sample files", {
  htseq_lines <- function(scale = 1) {
    c(
      paste("geneA", 5 * scale, sep = "\t"),
      paste("geneB", 7 * scale, sep = "\t"),
      paste("__no_feature", 9 * scale, sep = "\t")
    )
  }

  f1 <- tempfile(pattern = "sampleA_", fileext = ".txt")
  f2 <- tempfile(pattern = "sampleB_", fileext = ".txt")
  writeLines(htseq_lines(1), f1)
  writeLines(htseq_lines(2), f2)

  cnt <- read_vista_counts(c(f1, f2), format = "htseq", verbose = FALSE)

  expect_identical(cnt$input_format, "htseq")
  expect_identical(cnt$counts$gene_id, c("geneA", "geneB"))
  expect_equal(as.numeric(cnt$counts[2, -1]), c(7, 14))
})

test_that("read_vista_counts parses tximport-like inputs with tx2gene", {
  tx_counts <- matrix(
    c(10, 20, 30, 40),
    nrow = 2,
    dimnames = list(c("tx1", "tx2"), c("sampleA", "sampleB"))
  )
  tx_list <- list(
    counts = tx_counts,
    abundance = tx_counts / 10,
    length = tx_counts * 2
  )
  tx2gene <- data.frame(
    transcript_id = c("tx1", "tx2"),
    gene_id = c("geneA", "geneA"),
    stringsAsFactors = FALSE
  )

  cnt <- read_vista_counts(
    tx_list,
    format = "tximport",
    tx2gene = tx2gene,
    verbose = FALSE
  )

  expect_identical(cnt$input_format, "tximport")
  expect_identical(cnt$counts$gene_id, "geneA")
  expect_equal(as.numeric(cnt$counts[1, -1]), c(30, 70))
})

test_that("read_vista_counts parses RSEM gene result files", {
  rsem_df <- data.frame(
    gene_id = c("geneA", "geneB"),
    transcript_id = c("tx1", "tx2"),
    expected_count = c(12, 24),
    TPM = c(1.2, 2.4),
    FPKM = c(0.6, 1.2),
    stringsAsFactors = FALSE
  )
  f1 <- tempfile(pattern = "sampleA_", fileext = ".genes.results")
  write.table(rsem_df, file = f1, sep = "\t", quote = FALSE, row.names = FALSE)

  cnt <- read_vista_counts(f1, format = "rsem", verbose = FALSE)

  expect_identical(cnt$input_format, "rsem")
  expect_equal(ncol(cnt$counts), 2)
  expect_identical(colnames(cnt$counts)[1], "gene_id")
  expect_true("transcript_id" %in% colnames(cnt$row_data))
})

test_that("read_vista_metadata infers sample_names and match_vista_inputs reorders samples", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  counts_in <- count_data[seq_len(10), c("gene_id", sample_metadata$sample_names[seq_len(4)]), drop = FALSE]
  sample_in <- sample_metadata[seq_len(4), c("cond_long", "sample_names"), drop = FALSE]
  sample_in <- sample_in[c(4, 2, 1, 3), , drop = FALSE]
  rownames(sample_in) <- NULL
  names(sample_in)[names(sample_in) == "sample_names"] <- "sample"

  cnt <- read_vista_counts(counts_in, format = "matrix", gene_id_column = "gene_id", verbose = FALSE)
  si <- read_vista_metadata(sample_in, verbose = FALSE)
  matched <- match_vista_inputs(cnt, si, verbose = FALSE)

  expect_identical(matched$sample_info$sample_names, colnames(matched$counts)[-1])
  expect_identical(rownames(matched$sample_info), matched$sample_info$sample_names)
})

test_that("match_vista_inputs drops unmatched samples when requested", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())

  counts_in <- count_data[seq_len(10), c("gene_id", sample_metadata$sample_names[seq_len(4)]), drop = FALSE]
  sample_in <- sample_metadata[seq_len(5), c("sample_names", "cond_long"), drop = FALSE]

  matched <- match_vista_inputs(
    counts = counts_in,
    sample_info = sample_in,
    column_geneid = "gene_id",
    drop_unmatched = TRUE,
    verbose = FALSE
  )

  expect_equal(length(matched$sample_names), 4)
  expect_true(all(matched$sample_names %in% sample_in$sample_names))
})

test_that("derive_vista_metadata parses split-based sample names", {
  meta <- derive_vista_metadata(
    counts = NULL,
    sample_names = c("HET_1_U", "HET_2_ovary", "WT_1_U"),
    parser = "split",
    split = "_",
    fields = c("genotype", "replicate", "tissue"),
    verbose = FALSE
  )

  expect_identical(colnames(meta), c("sample_names", "genotype", "replicate", "tissue"))
  expect_identical(meta$genotype, c("HET", "HET", "WT"))
  expect_identical(meta$tissue, c("U", "ovary", "U"))
})

test_that("derive_vista_metadata parses regex-based sample names", {
  meta <- derive_vista_metadata(
    counts = NULL,
    sample_names = c("SRR1039508", "SRR1039509"),
    parser = "regex",
    pattern = "SRR(\\d+)",
    fields = "run_id",
    verbose = FALSE
  )

  expect_identical(meta$run_id, c("1039508", "1039509"))
})

test_that("derive_vista_metadata auto-detects consistent split sample names", {
  meta <- derive_vista_metadata(
    counts = NULL,
    sample_names = c("WT-1-ovary", "WT-2-uterus"),
    parser = "auto",
    fields = c("genotype", "replicate", "tissue"),
    verbose = FALSE
  )

  expect_identical(meta$genotype, c("WT", "WT"))
  expect_identical(meta$replicate, c("1", "2"))
})

test_that("derive_vista_metadata can derive from repaired count sample names", {
  data("count_data", package = "VISTA", envir = environment())

  sample_cols <- colnames(count_data)[colnames(count_data) != "gene_id"][seq_len(2)]
  counts_in <- count_data[seq_len(4), c("gene_id", sample_cols), drop = FALSE]
  names(counts_in)[2:3] <- c(
    "/proj/run/03_alignment/HET_1_U/HET_1_U_star_alignAligned.sortedByCoord.out.bam",
    "/proj/run/03_alignment/WT_2_ovary/WT_2_ovary_star_alignAligned.sortedByCoord.out.bam"
  )

  meta <- derive_vista_metadata(
    counts = counts_in,
    column_geneid = "gene_id",
    parser = "split",
    fields = c("genotype", "replicate", "tissue"),
    verbose = FALSE
  )

  expect_identical(meta$sample_names, c("HET_1_U", "WT_2_ovary"))
  expect_identical(meta$replicate, c("1", "2"))
  expect_true(is.data.frame(attr(meta, "sample_name_map")))
})

test_that("derive_vista_metadata template adds placeholder columns", {
  meta <- derive_vista_metadata(
    counts = NULL,
    sample_names = c("sampleA", "sampleB"),
    parser = "none",
    return_type = "template",
    verbose = FALSE
  )

  expect_true(all(c("group", "batch") %in% colnames(meta)))
  expect_true(all(is.na(meta$group)))
  expect_true(all(is.na(meta$batch)))
})

# --- B6: multi-file importers must bind by gene identifier, not by position ---

test_that("STAR files with shuffled gene order are matched by identifier", {
  star_file <- function(genes, values) {
    lines <- c(
      "N_unmapped\t1\t1\t1",
      "N_multimapping\t2\t2\t2",
      "N_noFeature\t3\t3\t3",
      "N_ambiguous\t4\t4\t4",
      paste(genes, values, values, values, sep = "\t")
    )
    f <- tempfile(fileext = ".ReadsPerGene.out.tab")
    writeLines(lines, f)
    f
  }

  genes <- c("geneA", "geneB", "geneC")
  f1 <- star_file(genes, c(10, 20, 30))
  # Same genes, different row order, distinct values.
  f2 <- star_file(rev(genes), c(300, 200, 100))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  cnt <- read_vista_counts(
    c(f1, f2), format = "star", sample_names = c("s1", "s2"), verbose = FALSE
  )

  expect_identical(cnt$counts$gene_id, genes)
  expect_equal(cnt$counts$s1, c(10, 20, 30))
  # Positional binding would have produced 300, 200, 100 here.
  expect_equal(cnt$counts$s2, c(100, 200, 300))
})

test_that("HTSeq files with shuffled gene order are matched by identifier", {
  htseq_file <- function(genes, values) {
    f <- tempfile(fileext = ".txt")
    writeLines(c(paste(genes, values, sep = "\t"), "__no_feature\t99"), f)
    f
  }

  genes <- c("geneA", "geneB", "geneC")
  f1 <- htseq_file(genes, c(1, 2, 3))
  f2 <- htseq_file(rev(genes), c(30, 20, 10))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  cnt <- read_vista_counts(
    c(f1, f2), format = "htseq", sample_names = c("s1", "s2"), verbose = FALSE
  )

  expect_identical(cnt$counts$gene_id, genes)
  expect_equal(cnt$counts$s1, c(1, 2, 3))
  expect_equal(cnt$counts$s2, c(10, 20, 30))
})

test_that("RSEM files with shuffled gene order are matched by identifier", {
  rsem_file <- function(genes, values) {
    f <- tempfile(fileext = ".genes.results")
    utils::write.table(
      data.frame(
        gene_id = genes, transcript_id = paste0("tx_", genes),
        expected_count = values, stringsAsFactors = FALSE
      ),
      file = f, sep = "\t", quote = FALSE, row.names = FALSE
    )
    f
  }

  genes <- c("geneA", "geneB", "geneC")
  f1 <- rsem_file(genes, c(5, 6, 7))
  f2 <- rsem_file(rev(genes), c(70, 60, 50))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  cnt <- read_vista_counts(
    c(f1, f2), format = "rsem", sample_names = c("s1", "s2"), verbose = FALSE
  )

  expect_identical(cnt$counts$gene_id, genes)
  expect_equal(cnt$counts$s1, c(5, 6, 7))
  expect_equal(cnt$counts$s2, c(50, 60, 70))
  # row_data is taken from the first file and must stay aligned.
  expect_identical(as.character(cnt$row_data$gene_id), genes)
})

test_that("a file that does not cover the reference gene set is an error", {
  htseq_file <- function(genes, values) {
    f <- tempfile(fileext = ".txt")
    writeLines(paste(genes, values, sep = "\t"), f)
    f
  }

  f1 <- htseq_file(c("geneA", "geneB", "geneC"), c(1, 2, 3))
  f2 <- htseq_file(c("geneA", "geneB", "geneZ"), c(1, 2, 3))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  # Equal row counts, so the old positional bind accepted this silently.
  expect_error(
    read_vista_counts(
      c(f1, f2), format = "htseq", sample_names = c("s1", "s2"), verbose = FALSE
    ),
    "do not cover the reference gene set"
  )

  msg <- tryCatch(
    read_vista_counts(
      c(f1, f2), format = "htseq", sample_names = c("s1", "s2"), verbose = FALSE
    ),
    error = function(e) conditionMessage(e)
  )
  expect_match(msg, basename(f2), fixed = TRUE)
  expect_match(msg, "geneC", fixed = TRUE)
})

test_that("extra identifiers in a later file are dropped with a message", {
  htseq_file <- function(genes, values) {
    f <- tempfile(fileext = ".txt")
    writeLines(paste(genes, values, sep = "\t"), f)
    f
  }

  f1 <- htseq_file(c("geneA", "geneB"), c(1, 2))
  f2 <- htseq_file(c("geneB", "geneA", "geneEXTRA"), c(20, 10, 999))
  on.exit(unlink(c(f1, f2)), add = TRUE)

  expect_message(
    cnt <- read_vista_counts(
      c(f1, f2), format = "htseq", sample_names = c("s1", "s2"), verbose = FALSE
    ),
    "Dropped 1 identifier"
  )
  expect_identical(cnt$counts$gene_id, c("geneA", "geneB"))
  expect_equal(cnt$counts$s2, c(10, 20))
})

test_that("identically ordered files still bind positionally", {
  # Fast path: identical id vectors, including repeated identifiers, must work
  # without requiring uniqueness.
  tabs <- list(
    data.frame(id = c("g1", "g1", "g2"), n = c(1, 2, 3), stringsAsFactors = FALSE),
    data.frame(id = c("g1", "g1", "g2"), n = c(4, 5, 6), stringsAsFactors = FALSE)
  )
  out <- VISTA:::.bind_count_files(
    tabs, id_col = "id", value_col = "n", sample_names = c("a", "b")
  )
  expect_identical(out$gene_id, c("g1", "g1", "g2"))
  expect_equal(out$a, c(1, 2, 3))
  expect_equal(out$b, c(4, 5, 6))
})

test_that("duplicated identifiers with differing order are rejected", {
  tabs <- list(
    data.frame(id = c("g1", "g1", "g2"), n = c(1, 2, 3), stringsAsFactors = FALSE),
    data.frame(id = c("g2", "g1", "g1"), n = c(3, 1, 2), stringsAsFactors = FALSE)
  )
  expect_error(
    VISTA:::.bind_count_files(
      tabs, id_col = "id", value_col = "n", sample_names = c("a", "b")
    ),
    "duplicated gene identifiers"
  )
})

test_that("tximport abundance and length are refused as counts", {
  tx <- list(
    counts = matrix(c(10, 20, 30, 40), 2, dimnames = list(c("g1", "g2"), c("a", "b"))),
    abundance = matrix(1, 2, 2, dimnames = list(c("g1", "g2"), c("a", "b"))),
    length = matrix(1000, 2, 2, dimnames = list(c("g1", "g2"), c("a", "b")))
  )
  expect_error(read_vista_counts(tx, format = "tximport", counts_from = "abundance", verbose = FALSE), "not a count matrix")
  expect_error(read_vista_counts(tx, format = "tximport", counts_from = "length", verbose = FALSE), "not a count matrix")
})

test_that("unscaled tximport counts become length-scaled TPM counts", {
  ids <- list(c("g1", "g2", "g3"), c("a", "b"))
  counts <- matrix(c(100, 50, 10, 300, 40, 20), 3, dimnames = ids)
  abundance <- matrix(c(10, 2, 1, 20, 1, 2), 3, dimnames = ids)
  len <- matrix(c(1000, 4000, 2000, 1200, 3000, 2000), 3, dimnames = ids)
  tx <- list(counts = counts, abundance = abundance, length = len, countsFromAbundance = "no")

  expect_message(
    cnt <- read_vista_counts(tx, format = "tximport", verbose = FALSE),
    "lengthScaledTPM"
  )
  got <- as.matrix(cnt$counts[, c("a", "b")])
  rownames(got) <- cnt$counts$gene_id

  # By hand: abundance times the gene's mean length, rescaled to each sample's
  # original estimated library size.
  scaled <- abundance * rowMeans(len)
  expected <- sweep(scaled, 2, colSums(counts) / colSums(scaled), "*")
  expect_equal(got, expected[rownames(got), ])
  expect_equal(colSums(got), colSums(counts))
})

test_that("length-scaled tximport counts are used as they are", {
  ids <- list(c("g1", "g2"), c("a", "b"))
  counts <- matrix(c(10.5, 20.25, 30, 40.75), 2, dimnames = ids)
  tx <- list(counts = counts, abundance = counts, length = counts * 0 + 1000, countsFromAbundance = "lengthScaledTPM")
  cnt <- read_vista_counts(tx, format = "tximport", verbose = FALSE)
  expect_equal(unname(as.matrix(cnt$counts[, c("a", "b")])), unname(counts))
})

test_that("unscaled tximport counts without abundance or length warn", {
  ids <- list(c("g1", "g2"), c("a", "b"))
  tx <- list(counts = matrix(c(10, 20, 30, 40), 2, dimnames = ids), countsFromAbundance = "no")
  expect_warning(read_vista_counts(tx, format = "tximport", verbose = FALSE), "without transcript-length correction")
})

test_that("transcript-level tximport input is summarised to genes as tximport does", {
  ids <- list(c("t1", "t2", "t3"), c("a", "b"))
  counts <- matrix(c(10, 30, 5, 0, 60, 8), 3, dimnames = ids)
  abundance <- matrix(c(2, 3, 1, 0, 5, 2), 3, dimnames = ids)
  len <- matrix(c(1000, 2000, 1500, 1100, 2100, 1500), 3, dimnames = ids)
  tx2gene <- data.frame(tx = c("t1", "t2", "t3"), gene = c("G1", "G1", "G2"))
  out <- VISTA:::.tximport_summarize_to_gene(counts, abundance, len, gene_ids = c("G1", "G1", "G2"))

  # Gene length is the abundance-weighted mean of its transcripts' lengths.
  expect_equal(out$length["G1", "a"], (2 * 1000 + 3 * 2000) / 5)
  expect_equal(out$length["G1", "b"], (0 * 1100 + 5 * 2100) / 5)
  expect_equal(out$counts["G1", ], c(a = 40, b = 60))

  # A gene with no abundance in a sample takes the geometric mean of its
  # lengths in the other samples.
  abundance["t3", "b"] <- 0
  out <- VISTA:::.tximport_summarize_to_gene(counts, abundance, len, gene_ids = c("G1", "G1", "G2"))
  expect_equal(out$length["G2", "b"], 1500)
})

test_that("RSEM TPM and FPKM columns are refused as counts", {
  rs <- data.frame(gene_id = c("a", "b"), expected_count = c(10.5, 3), TPM = c(1, 2), FPKM = c(1, 2))
  expect_error(read_vista_counts(rs, format = "rsem", count_column = "TPM", verbose = FALSE), "not one of RSEM's count columns")
  expect_error(read_vista_counts(rs, format = "rsem", count_column = "FPKM", verbose = FALSE), "not one of RSEM's count columns")
})

test_that("DESeq2 runs on estimated (non-integer) counts by rounding them", {
  data("count_data", package = "VISTA", envir = environment())
  data("sample_metadata", package = "VISTA", envir = environment())
  cnt <- count_data[seq_len(300), c("gene_id", sample_metadata$sample_names)]
  cnt[, -1] <- cnt[, -1] * 1.13
  expect_message(
    v <- create_vista(cnt, sample_metadata, "gene_id", "dex", "trt", "untrt", method = "deseq2", min_counts = 5),
    "rounding the non-integer"
  )
  expect_true(all(counts(v) == round(counts(v))))
})

test_that("tximport abundance and length are matched to counts by name", {
  ids <- list(c("g1", "g2", "g3"), c("a", "b"))
  counts <- matrix(c(100, 50, 10, 300, 40, 20), 3, dimnames = ids)
  abundance <- matrix(c(10, 2, 1, 20, 1, 2), 3, dimnames = ids)
  len <- matrix(c(1000, 4000, 2000, 1200, 3000, 2000), 3, dimnames = ids)
  ordered <- list(counts = counts, abundance = abundance, length = len, countsFromAbundance = "no")
  shuffled <- ordered
  shuffled$abundance <- abundance[c(3, 1, 2), c(2, 1)]
  shuffled$length <- len[c(2, 3, 1), ]
  a <- suppressMessages(read_vista_counts(ordered, format = "tximport", verbose = FALSE))$counts
  b <- suppressMessages(read_vista_counts(shuffled, format = "tximport", verbose = FALSE))$counts
  expect_equal(b[match(a$gene_id, b$gene_id), c("a", "b")], a[, c("a", "b")], ignore_attr = TRUE)

  renamed <- ordered
  rownames(renamed$abundance)[1] <- "other"
  expect_error(read_vista_counts(renamed, format = "tximport", verbose = FALSE), "do not match")
})

test_that("tximport inputs that cannot be scaled are refused", {
  ids <- list(c("g1", "g2"), c("a", "b"))
  base <- list(
    counts = matrix(c(10, 20, 30, 40), 2, dimnames = ids),
    abundance = matrix(c(1, 2, 3, 4), 2, dimnames = ids),
    length = matrix(1000, 2, 2, dimnames = ids),
    countsFromAbundance = "no"
  )
  with_na <- base; with_na$abundance[1, 2] <- NA
  expect_error(read_vista_counts(with_na, format = "tximport", verbose = FALSE), "missing, non-finite or negative")
  zero <- base; zero$abundance[, 2] <- 0
  expect_error(read_vista_counts(zero, format = "tximport", verbose = FALSE), "total abundance is zero")
})

test_that("a tx2gene that matches no transcript is a clear error", {
  ids <- list(c("tx1", "tx2"), c("a", "b"))
  tx <- list(
    counts = matrix(c(10, 20, 30, 40), 2, dimnames = ids),
    abundance = matrix(c(1, 2, 3, 4), 2, dimnames = ids),
    length = matrix(1000, 2, 2, dimnames = ids),
    countsFromAbundance = "no"
  )
  t2g <- data.frame(tx = c("tx1.1", "tx2.1"), gene = c("G", "G"))
  expect_error(read_vista_counts(tx, format = "tximport", tx2gene = t2g, verbose = FALSE), "None of the transcript IDs")
})

test_that("scaled transcript-level counts are rebuilt at gene level, as summarizeToGene() does", {
  ids <- list(c("t1", "t2", "t3"), c("a", "b"))
  counts <- matrix(c(10, 30, 5, 2, 60, 8), 3, dimnames = ids)
  abundance <- matrix(c(2, 3, 1, 1, 5, 2), 3, dimnames = ids)
  len <- matrix(c(1000, 2000, 1500, 1100, 2100, 1500), 3, dimnames = ids)
  tx <- list(counts = counts, abundance = abundance, length = len, countsFromAbundance = "lengthScaledTPM")
  t2g <- data.frame(tx = c("t1", "t2", "t3"), gene = c("G1", "G1", "G2"))
  got <- suppressMessages(read_vista_counts(tx, format = "tximport", tx2gene = t2g, verbose = FALSE))$counts
  got <- as.matrix(got[, c("a", "b")]); rownames(got) <- c("G1", "G2")

  g <- VISTA:::.tximport_summarize_to_gene(counts, abundance, len, c("G1", "G1", "G2"))
  expected <- VISTA:::.tximport_length_scaled_counts(g$counts, g$abundance, g$length)
  expect_equal(got, expected[rownames(got), ], ignore_attr = TRUE)
  # Not the plain sum of the scaled transcript counts.
  expect_false(isTRUE(all.equal(got["G1", ], colSums(counts[1:2, ]), check.attributes = FALSE)))
})

test_that("only RSEM count columns are accepted", {
  rs <- data.frame(gene_id = c("a", "b"), expected_count = c(10.5, 3), posterior_mean_count = c(10, 3),
                   pme_TPM = c(1, 2), TPM = c(1, 2))
  expect_error(read_vista_counts(rs, format = "rsem", count_column = "pme_TPM", verbose = FALSE), "not one of RSEM's count columns")
  expect_silent(read_vista_counts(rs, format = "rsem", count_column = "posterior_mean_count", verbose = FALSE))
})

test_that("an unrecognised countsFromAbundance value is refused", {
  ids <- list(c("g1", "g2"), c("a", "b"))
  tx <- list(counts = matrix(c(10, 20, 30, 40), 2, dimnames = ids), countsFromAbundance = "lengthscaledtpm")
  expect_error(read_vista_counts(tx, format = "tximport", verbose = FALSE), "must be one of")
})

test_that("abundance and length named only by row are matched by row name", {
  ids <- list(c("g1", "g2", "g3"), c("a", "b"))
  counts <- matrix(c(100, 50, 10, 300, 40, 20), 3, dimnames = ids)
  abundance <- matrix(c(10, 2, 1, 20, 1, 2), 3, dimnames = ids)
  len <- matrix(c(1000, 4000, 2000, 1200, 3000, 2000), 3, dimnames = ids)
  ref <- suppressMessages(read_vista_counts(list(counts = counts, abundance = abundance, length = len, countsFromAbundance = "no"),
                                            format = "tximport", verbose = FALSE))$counts
  shuffled_ab <- abundance[c(3, 1, 2), ]; colnames(shuffled_ab) <- NULL
  shuffled_len <- len[c(2, 3, 1), ]; colnames(shuffled_len) <- NULL
  got <- suppressMessages(read_vista_counts(list(counts = counts, abundance = shuffled_ab, length = shuffled_len, countsFromAbundance = "no"),
                                            format = "tximport", verbose = FALSE))$counts
  expect_equal(got, ref)
})

test_that("missing values only matter for transcripts that tx2gene keeps", {
  ids <- list(c("t1", "t2", "t3"), c("a", "b"))
  counts <- matrix(c(10, 30, 5, 2, 60, 8), 3, dimnames = ids)
  abundance <- matrix(c(2, 3, NA, 1, 5, 2), 3, dimnames = ids)
  len <- matrix(c(1000, 2000, 1500, 1100, 2100, 1500), 3, dimnames = ids)
  tx <- list(counts = counts, abundance = abundance, length = len, countsFromAbundance = "no")
  t2g <- data.frame(tx = c("t1", "t2"), gene = c("G1", "G1"))
  cnt <- suppressMessages(read_vista_counts(tx, format = "tximport", tx2gene = t2g, verbose = FALSE))$counts
  expect_identical(cnt$gene_id, "G1")
  t2g_all <- data.frame(tx = c("t1", "t2", "t3"), gene = c("G1", "G1", "G2"))
  expect_error(suppressMessages(read_vista_counts(tx, format = "tximport", tx2gene = t2g_all, verbose = FALSE)), "missing, non-finite")
})

test_that("a merged RSEM table is read from its named sample columns", {
  merged <- data.frame(gene_id = c("a", "b"), s1 = c(10.5, 3), s2 = c(7, 1.25))
  cnt <- read_vista_counts(merged, format = "rsem", sample_names = c("s1", "s2"), count_column = "s1", verbose = FALSE)
  expect_equal(unname(as.matrix(cnt$counts[, c("s1", "s2")])), unname(as.matrix(merged[, c("s1", "s2")])))
})
