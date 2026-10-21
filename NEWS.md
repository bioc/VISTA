# VISTA 1.1.6

## New features

- **edgeR's quasi-likelihood F-test.** `create_vista()`,
  `run_edger_analysis()`, `benchmark_vista_equivalence()` and
  `validate_vista_deep()` gain `edger_test = c("LRT", "QL")`. `"QL"` runs
  `glmQLFit()`/`glmQLFTest()` (edgeR v4 pipeline), which the edgeR authors
  recommend for bulk RNA-seq; the likelihood-ratio test treats the dispersion
  as known and is liberal with few replicates. The test used is recorded in
  `cutoffs(v)$edger_test`, and the report config accepts `edger_test:`; for an
  edgeR object made before 1.1.6, which does not record the test, the report
  states the likelihood-ratio test and marks it as inferred. The equivalence
  check runs its own QL reference. The vignettes now set `edger_test`.

## Default changing in VISTA 1.4.0

- **`edger_test` will default to `"QL"`.** It stays `"LRT"` in this release
  so that no results change without warning. A call that runs edgeR
  (`method = "edger"` or `"both"`) without setting `edger_test`, or with
  `edger_test = NULL`, prints a message, classed `vista_default_change` and
  throttled like the deprecation warnings, saying so. Set
  `edger_test = "LRT"` to keep current results after the change, or
  `edger_test = "QL"` to use the new default now; either silences the message.
  DESeq2 and limma-voom analyses are not affected. See `?VISTA-deprecated`.

## Bug fixes

- **tximport abundance or lengths could be analysed as counts.**
  `read_vista_counts(format = "tximport", counts_from = "abundance")` imported
  TPM, and `counts_from = "length"` transcript lengths, as if they were counts,
  and edgeR and limma tested them without complaint. Both now raise an error.
- **tximport counts were analysed without transcript-length correction.**
  tximport's default `countsFromAbundance = "no"` leaves estimated counts
  unscaled and expects the DE method to use length offsets, which VISTA does
  not. Such input that carries `abundance` and `length` (also when it does not
  say how it was scaled) is now converted to length-scaled TPM counts
  (tximport's `countsFromAbundance = "lengthScaledTPM"`, recommended when the
  method does not use offsets), with a message. With `tx2gene`, transcript-level
  input is summarised to genes as tximport's `summarizeToGene()` does, which
  rebuilds gene counts from abundance and length even when the transcript
  counts were already scaled; summing them, as VISTA did, undid the gene-level
  length correction. The results agree with tximport's own to within
  numerical precision (relative differences below 1e-14). `abundance` and
  `length` are matched to `counts` by name, and missing or negative values, a
  sample with zero total abundance, or a `tx2gene` that matches no transcript
  are clear errors, as is an unrecognised `countsFromAbundance` value. Without
  `abundance`/`length`, VISTA warns and uses the counts as given. *Were you
  affected?* If you imported tximport output made with the default
  `countsFromAbundance = "no"`, a hand-built list without that field, or
  transcript-level `"lengthScaledTPM"`/`"dtuScaledTPM"` output with `tx2gene`
  (scaledTPM sums to the same gene counts); re-run those analyses.
- **DESeq2 stopped on estimated counts.** Non-integer counts from tximport or
  RSEM made `create_vista(method = "deseq2")` fail with "some values in assay
  are not integers". DESeq2 now receives them rounded, with a message, as
  `DESeqDataSetFromTximport()` does. Filtering is applied to the counts as
  given, as for edgeR and limma, so `method = "both"` keeps one set of genes,
  and its `counts(v)` holds DESeq2's rounded counts whichever source is
  active.
- **Values that are not counts could be analysed without comment.** All three
  engines now refuse negative values, and warn (class
  `vista_not_counts_warning`) about non-integer values with an unmistakable
  sign of not being counts: every sample summing to exactly one million (TPM,
  CPM) or small, almost all fractional values (log scale). These are warnings
  because no rule on the values alone separates counts from normalised data:
  estimated counts from libraries of about a million reads must not be
  refused, and TMM-CPM, FPKM or filtered TPM cannot be recognised. Supply raw
  or estimated counts. Integer counts are unaffected.
- **RSEM TPM or FPKM could be imported as counts.** Where `count_column`
  selects the values read from RSEM output, only RSEM's count columns,
  `expected_count` (the default) and `posterior_mean_count`, are now accepted.
  A merged table whose sample columns are named in `sample_names` is read as
  given.

- **Untested genes were stored with `pvalue = 1`.** DESeq2 reports `NA` for
  genes it does not test (Cook's-distance count outliers); VISTA replaced that
  with 1 when it stored the comparison, beside a `padj` that stayed `NA`, so an
  untested gene read as tested and the row disagreed with itself. On airway
  this affected 51 genes. Tables from `comparisons()` and those written by
  `export_vista_assets()` now show `NA`, as the engine reported, and so do the
  per-backend `pvalue_deseq2`/`pvalue_edger` columns of a consensus table.
  Plots and rankings already substituted their own value for a missing
  p-value, so figures are unchanged.
- **Filtering a consensus table on `padj` returned the union, not the
  intersection.** With `consensus_mode = "intersection"`, a gene called by only
  one backend is not called, but its `padj` was that backend's significant
  value, contrary to the documented "less-significant of the two". Selecting
  `padj <= 0.05` on the table, as `get_foldchange_scatter()` and the volcano
  thresholds do, therefore picked up every single-backend gene. Such genes now
  carry the less-significant `pvalue` and `padj` of the two backends, and their
  `log2fc` follows `consensus_log2fc` (the mean of the two by default) instead
  of the calling backend's estimate, as for every other uncalled gene. On
  airway this changes `log2fc` for all 119 single-backend genes (by up to 0.36)
  and the p-values of most of them. In union mode, where they are called, they
  keep the calling backend's values, as before. `regulation` was always correct
  and is unchanged.
- **DEG-count plots dropped comparisons with no DEGs.**
  `get_deg_count_barplot()`, `get_deg_count_pieplot()` and
  `get_deg_count_donutplot()` left out a comparison whose summary was all
  "Other", asking for it by name aborted with "not found in DEG summary", and
  when no requested comparison had any DEG they returned `NULL` with a warning.
  Such comparisons are now drawn at zero. Comparison order is unchanged. A
  summary that uses none of the labels `Up`, `Down` or `Other` (possible with
  `as_vista()`) is left out with a warning naming the labels found, rather than
  drawn as zero.
- **limma MA plots put the wrong quantity on the x axis.** For limma the
  table's `baseMean` is `AveExpr`, already log2-CPM, and `get_ma_plot()`
  preferred it to the object's mean expression, so the axis showed
  `log10(log2-CPM + 1)`. It now uses the object's mean expression for every
  engine, as its own comment intended.
- **An explicit intercept was not removed for edgeR and limma.** For
  `design_formula = "~ 1 + cell + group"`, VISTA prefixed `0 +`, which `+ 1`
  then undid, so the model kept its intercept. VISTA now reads the intercept
  from the formula's terms. The fitted model spans the same space, so results
  differ only within numerical precision.

## Verification

- **`benchmark_vista_equivalence()` no longer shares code with what it
  checks.** Its reference pipelines used VISTA's own helpers for sample
  matching, design construction, contrasts, DEG calls and row alignment, so a
  defect in any of them was reproduced on both sides and passed -- which is
  how the first and third defects above went unnoticed. The reference now
  calls only base R, `stats` and the engines, re-implements the documented
  contract from scratch (see `?benchmark_vista_equivalence`), derives contrasts
  from the design itself, reads VISTA's tables through its exported accessors
  and matches rows by identifier. It compares every numeric column the engine
  reports, the retained gene set and the set of comparisons; the summary gains
  `gene_sets_identical`, `comparison_sets_identical`,
  `all_columns_within_tolerance` and `columns_compared`.
  `baseMean_within_tolerance` is `NA` for edgeR, which reports no `baseMean`
  (it was a vacuous `TRUE`); a numeric column VISTA stores but the reference
  does not report (`columns_not_in_reference`) counts as a mismatch; and a
  mismatch between a value and `NA` is reported as an infinite difference
  instead of 0. The MA, PCA and DEG-count checks now read the data
  `get_ma_plot()`, `get_pca_plot()` and the DEG-count plots draw, and the PCA
  check compares every sample-sheet column; the volcano check still rebuilds
  its inputs from the stored table.
- Each method's `reference` element has a new shape. Its `comparisons` are
  data frames in VISTA's column names (`gene_id`, `log2fc`, `stat`, `pvalue`,
  `padj`, plus `baseMean` and `lfcSE` for DESeq2, `logCPM` for edgeR, and
  `baseMean` (limma's `AveExpr`), `B`, `P.Value` and `adj.P.Val` for limma),
  not the engines' own tables. `row_data` is removed (use
  `rowMeans(reference$norm_counts)`), `group_column` is added, and each
  `deg_summary` table has exactly two rows, `Up` and `Down`.
- New tests check that no reference function uses VISTA code -- as a call,
  as a function value, through `VISTA::`/`VISTA:::`, or by name in a string. They inject defects into VISTA's own helpers that the
  benchmark must report: a flipped contrast sign, swapped Up/Down calls,
  values under the wrong gene labels, samples matched to the wrong metadata
  rows, a dropped covariate, untested genes stored as `p = 1`, and a
  comparison dropped from DEG-count plot data. The 1.1.5 check caught two of
  the first five; those five were chosen to target the helpers it shared, so
  this shows the new check can fail where the old one could not, not a
  detection rate.

## Other changes

- Objects record the VISTA and R versions and the installed DESeq2, edgeR and
  limma versions in `metadata(x)$provenance`, as `?VISTA-class` already
  described. The stored object behind `example_vista()` was regenerated to
  carry them.
- `?create_vista` documents how a sample sheet without `sample_names` is
  matched to count columns, and what the consensus `padj` is: a conservative
  ranking score, not an FDR-adjusted p-value of the consensus. In simulations
  the intersection lowered, but did not control, the single backends' false
  discovery rate, and the union raised it.
- The DESCRIPTION now names limma-voom among the engines.
- The code-economy vignette now computes every tally from its own chunks at
  build time, states its counting rule, attaches what each workflow needs, and
  checks that both workflows give the same DESeq2 results.
- `min_counts` and `min_replicates` in `benchmark_vista_equivalence()` and
  `validate_vista_deep()` are now described as they behave: `min_replicates`
  counts samples across the whole experiment, not per group.

# VISTA 1.1.5

## New features

- `get_foldchange_matrix()` gains `display_id`, `display_from` and
  `display_orgdb`, matching [get_foldchange_heatmap()]. `genes` may now be
  given as display labels, and the returned matrix is labelled with them.
  Because symbols are not unique, duplicated labels are made unique with a
  warning rather than left to collide -- duplicated rownames would make
  `m["SYM", ]` return whichever row came first and silently hide the other.
  Rows keep the object's identifiers unless `display_id` is supplied, so
  existing calls are unaffected.

## Bug fixes

- **`get_foldchange_heatmap(annotate_columns = TRUE)` annotated nothing.** The
  columns of a fold-change heatmap are comparisons, but the annotation was
  built by filtering `colData` for samples named in those columns -- an
  intersection that is always empty -- so the track carried no levels and no
  colours. It now annotates by comparison, coloured from the object's
  comparison palette, and accepts `column_anno_colors` like the expression
  heatmap.
- **A `display_id` supplied without a `display_orgdb` failed with
  `missing value where TRUE/FALSE needed`.** The internal comparison
  `to_type == from_type` is `NA` when `from_type` is `NA` and `logical(0)`
  when it is `NULL`, and `||` can evaluate neither. It now names the missing
  argument.
- **`get_expression_lineplot()` accepted `display_from`/`display_orgdb` and
  ignored them** -- its private copy of the identifier mapping handled only the
  `rowData()` path. All gene-taking plots now share one resolver.

# VISTA 1.1.4

## Internal

- Plot-stability snapshots for PCA and MDS no longer depend on eigenvector
  sign. `prcomp()` and `cmdscale()` may return an axis or its negation for
  identical input, decided by the LAPACK/BLAS build, so the snapshots could
  fail on one machine and pass on another with nothing in the package
  changed. Each axis is now canonicalised before comparison. No user-visible
  behaviour changes; only the test suite is affected.

# VISTA 1.1.3

## Bug fixes

- **Statistical annotations failed under ggpubr 1.0.0.** Passing
  `stats_group = TRUE` produced a plot that constructed normally and then
  died when drawn, with `could not find function "create_p_label"`. ggpubr
  1.0.0 builds `p.label = "p.format"` into an expression calling that
  package-internal helper, which is not visible from the environment of the
  mapping VISTA supplies, so the annotation layer could not be evaluated.
  VISTA now names the statistic's own computed columns via
  `ggplot2::after_stat()`. The rendered label text is unchanged on every
  ggpubr version, and an unrecognised `p.label` still falls through to
  ggpubr rather than silently dropping the annotation.

  Affects `get_expression_boxplot()`, `get_expression_violinplot()`,
  `get_expression_barplot()`, `get_expression_raincloud()`,
  `get_foldchange_boxplot()` and `get_foldchange_raincloud()` — anywhere
  `stats_group = TRUE` was combined with the default `p.label = "p.format"`.
  Users on ggpubr < 1.0.0 were unaffected.

## Documentation

- Figure cropping is disabled (`crop = NULL`). BiocStyle's crop hook shells
  out to `magick`; without it every figure emitted a warning.

# VISTA 1.1.2

Targeted at Bioconductor 3.24 (VISTA 1.2.0). This release fixes several
defects that produced silently wrong results, makes VISTA behave like a
proper Bioconductor object, and begins a deprecation cycle for a number of
inconsistent argument names. **No existing script stops working**: every
renamed argument still functions and warns.

## Data-integrity fixes

Check whether you were affected:

- **`get_expression_matrix(genes = ..., summarise = TRUE)` mislabelled rows.**
  Rows were subset in assay order but labelled in the order `genes` was
  supplied, so whenever those differed every row carried another gene's
  identifier and values.

  *Scope:* no other VISTA function called `get_expression_matrix()`, so no
  plot, heatmap, DEG table, export or generated report was affected. Internal
  group averaging goes through `norm_counts(summarise = TRUE)`, which labels
  rows from the matrix itself and was never affected by this. Only code that
  called `get_expression_matrix()` directly could be wrong.

  *Were you affected?* Only if you passed `genes` **and** `summarise = TRUE`
  **and** your gene vector was in some order other than the object's own row
  order. Passing genes already in `rownames(x)` order (or omitting `genes`)
  always produced correct output. Note that ranked inputs — for example
  `get_genes_by_regulation(top_n = ...)`, which sorts by absolute fold change —
  are normally *not* in row order, so those calls were affected. Re-run them.

  The bug dates back to the first release, so it is present in 1.0.0 and every
  0.99.x.
- **GSEA ranked vectors slid scores onto the wrong genes.** `get_gsea()`
  reassigned `names(rank_vec)` from a mapping helper that drops unmapped and
  duplicate identifiers, so the shorter vector was padded with `NA` and every
  score after the first gap moved to a different gene. The resulting NES
  values and leading-edge gene lists were meaningless.

  *Scope:* only `set_type = "kegg"` was affected in practice. That path always
  converts to ENTREZID, and any real dataset contains identifiers with no
  ENTREZ mapping, so the vector always shortened. `set_type = "msigdb"`
  converted an identifier type to *itself*, which short-circuits before
  querying the annotation database and cannot shorten. `set_type = "go"`
  never reassigned the names at all. Verified on a 22-gene input containing
  two identifiers absent from the OrgDb: `"kegg"` returned 20 names for 22
  scores, `"msigdb"` returned 22.

  No other VISTA function called `get_gsea()`, so nothing else in the package
  consumed a corrupted ranking. Re-run any `get_gsea(set_type = "kegg")`
  results; MSigDB and GO results are unaffected.
- **`get_corr_heatmap()` drew a staircase, not a triangle.** The `triangle`
  mask was computed against input order while the axes were re-levelled to
  clustered order, under the default `cluster_by = "correlation"`.
- **Multi-file count importers bound columns positionally.** The STAR, HTSeq
  and RSEM importers took gene identifiers from the first file only and never
  checked that row *i* described the same gene in the others, so files with
  equal row counts but different ordering were silently misaligned and every
  count for those samples was attributed to the wrong gene.

  *Scope:* only `read_vista_counts()` given a **vector of two or more file
  paths** with `format` of `"star"`, `"htseq"` or `"rsem"`, and only when
  those files did not already list the same genes in the same order. A single
  file, a data frame, a matrix, `"featurecounts"` and `"tximport"` all take a
  name-keyed path and were never affected. Files written by one pipeline run
  normally do share an order, so a uniform set of inputs was fine; the risk is
  mixing files produced at different times or by different pipeline versions,
  where HTSeq's trailing `__no_feature`/`__ambiguous` rows or STAR's four
  leading `N_*` rows may differ.

  Imports now match by gene identifier and **error**, naming the file and up
  to five missing identifiers, when a file does not cover the reference gene
  set. This can surface an error where a script previously appeared to work;
  it was producing wrong numbers. Anything built from a multi-file STAR,
  HTSeq or RSEM import should be re-imported and re-analysed.
- **`run_vista_report()` reported its own defaults as your parameters.** When a
  prebuilt object was supplied via `vista_rds`, the report printed
  deseq2 / 1 / 0.05 / padj regardless of how the object was built.
- **`run_cell_deconvolution()`'s SYMBOL retry told xCell2 the wrong identifier
  type.** When the first scoring attempt failed, VISTA retried with a
  SYMBOL-keyed matrix built from `rowData`, but carried the original
  `gene_id_type` (for example `"ensembl"`) into the retry -- handing xCell2
  symbols while declaring them Ensembl IDs. The retry therefore either failed
  with the real error masked behind a generic message, or scored against a
  near-empty gene overlap.
- **Deconvolution was broken for `ENSG...:SYMBOL` rownames** — the common case
  reached automatically by `gene_id_type = "auto"`.
- Fixed `get_ma_plot()` reading a non-existent metadata key for its threshold
  fallback, `get_expression_lollipop()` failing to resolve symbols through an
  OrgDb, comparison colours collapsing to `character(0)` and crashing
  `get_foldchange_lineplot()`, `NaN` columns when summarising over a factor
  with unused levels, and YAML/filename escaping in report and asset output.
- **Hardened the DE-table alignment helper**, which reindexes every comparison
  onto the counts matrix and so underpins the whole package. It now rejects a
  duplicated reference gene list instead of emitting rows whose names R had
  silently de-duplicated (`g2` -> `g2.1`) away from the `gene_id` column they
  came from, coerces tibble input to a plain data frame — a tibble ignores
  `rownames<-`, so every gene read as absent and was replaced by an `NA` row
  with a perfectly correct label — and asserts as a post-condition that the
  returned rownames and `gene_id` column both match the reference exactly.
  No shipped code path reached the tibble case; both are closed against future
  callers, because a wrong answer here is invisible in the output's shape.

## Reproducibility

- **Pathway gene capping is now deterministic.** `get_pathway_heatmap(max_genes)`
  and `get_enrichment_chord(max_genes)` truncated the pathway gene set
  positionally. Pathway membership is a *set*, and clusterProfiler's ordering of
  it is arbitrary and, verified here, differs between R sessions -- so the same
  call could plot a different subset of genes on different runs. The heatmap now
  ranks candidates by expression variance (most informative rows first) and the
  chord by pathway participation, both breaking ties on the gene identifier.
- **Deconvolution scores are no longer aligned to samples by position when
  their own labels disagree.** `run_cell_deconvolution()` accepted scores whose
  row count matched the sample count and stamped the object's sample names onto
  them, even when the scores carried different labels of their own -- which
  would attach every sample's cell fractions to the wrong sample. Alignment is
  now by name whenever names are available in either orientation; contradictory
  labels are an error; and position is used only when the scores carry no
  labels at all, with a message saying so.
- The airway vignette sets a seed. GSEA estimates p-values by permutation and
  ComplexHeatmap's row k-means is unseeded, so the rendered document previously
  changed on every build for reasons unrelated to the code.

## Analysis changes

- The edgeR and limma backends now filter with the same predicate as DESeq2
  (a gene is kept when at least `min_replicates` samples reach `min_counts`);
  they previously ignored `min_counts` and used a hardcoded `cpm > 1`. They
  also filter *before* `calcNormFactors()`, per the edgeR user's guide. edgeR
  and limma results shift slightly relative to 1.1.1.
- `min_replicates` is documented accurately: it counts samples across the
  whole experiment, never within each group.

## New features

- **Raw counts are retained** as a `counts` assay, reachable with `counts(v)`.
  `as_deseq_dataset()` converts a VISTA object back into a `DESeqDataSet`.
  Roughly doubles the assay footprint; opt out with
  `create_vista(keep_raw_counts = FALSE)`.
- **`[` method.** Subsetting now keeps metadata consistent: row subsetting
  reindexes every DE table and recounts DEG summaries, column subsetting
  prunes orphaned colour maps and warns when a comparison loses its samples.
  Previously `v[1:10, ]` left the DE tables describing all original genes and
  `validObject()` still returned `TRUE`.
- **`show()` method.** Objects now display as `class: VISTA` with their
  grouping, comparisons, active DE source, cutoffs and schema version.
- **`updateObject()`** migrates objects from older metadata schemas.
- **`?VISTA`** resolves, and there is a package overview page.
- Unknown arguments passed through `...` are now rejected with a did-you-mean
  suggestion. A typo such as `get_expression_heatmap(v, gene = my_genes)`
  previously plotted the default gene set instead of erroring.
- `get_foldchange_lineplot()` gains `return_type`. Its default stays `"both"`
  (the list it has always returned); pass `"plot"` for a bare ggplot.
- `get_enrichment_chord()` gains `return_type`. Its default stays `"data"`
  because it draws to the active device and has always returned its table
  invisibly; `"plot"` returns a recorded plot that `save_vista_plot()` can now
  write to a file. VISTA enables the device display list before drawing, since
  `png()` and `pdf()` leave it off and would otherwise yield an empty
  recording that silently saved nothing.

## Deprecations

All still work and warn once per session. See `?VISTA-deprecated` for the
full table and timelines; each becomes defunct in 1.4.0.

| Function | Deprecated | Use instead |
|---|---|---|
| `get_deg_count_pieplot()`, `get_deg_count_donutplot()` | `label` | `label_type` |
| `get_corr_heatmap()` | `cluster_by` | `order_by` |
| `get_corr_heatmap()` | `show_corr_values`, `col_corr_values` | `label`, `label_color` |
| `get_celltype_group_dotplot()` | `error` | `errorbar` |
| four `get_expression_*()` plots | `comparisons` | `stat_comparisons` |
| `get_pca_plot()`, `get_mds_plot()`, `get_umap_plot()` | `top_n_genes` | `top_n` |
| `get_expression_lollipop()`, `get_foldchange_lollipop()` | `line_size` | `linewidth` |
| `get_volcano_plot()` | `col_up`, `col_down`, `col_other(s)`, `lab_size` | `colors`, `label_size` |
| `get_pca_plot()` | `sample.seed` | removed; it never had any effect |
| `get_expression_barplot()`, `get_expression_lollipop()` | `facet_scale` | `facet_scales` |
| `get_expression_violinplot()`, `get_expression_lineplot()` | `value_transform` | `log_transform` |

`return_type` values are unified on `c("plot", "data", "both")`. The legacy
spellings still work and warn: `"heatmap"`/`"clusters"` (both heatmaps),
`"matrix"` (`get_celltype_heatmap()`), `"genes"` (`get_pathway_heatmap()`).
Arguments that select the *shape* of a purely tabular result --
`get_pathway_genes()`, `get_genes_by_regulation()`, `read_vista_counts()`,
`derive_vista_metadata()` -- keep their own vocabularies, because they have no
plot component.

Undocumented gene caps (20 for embeddings, 15 for lollipop, 25 for barplot)
are now the `max_genes` argument, with the existing values as defaults.

## Testing

- Added label/value integrity tests covering every exported function that
  accepts a gene set. They assert two properties the 1.0.0 defect violated:
  results do not depend on the order genes were supplied in, and the value
  reported for a gene is that gene's own value, checked against an independent
  recomputation from the object. The previous suite asserted shape --
  `is.matrix()`, `nrow()`, class -- all of which the broken output satisfied.
- Verified by mutation testing: reintroducing the exact 1.0.0 defect now fails
  five tests across three files, and replacing a name-keyed lookup with a
  positional one elsewhere fails nine.
- The DE-table alignment helper is now covered end to end: identifier resolution
  and its fallbacks, DESeq2/edgeR/limma column canonicalisation, type-preserving
  `NA` padding across numeric, integer, character, logical, factor, `Date` and
  `POSIXct` columns, and the reordering guarantee itself. Both faults fixed above
  were found by these tests.

## Internal

- Accessor generics dispatch only on the object rather than on every formal.
- `Depends: R (>= 4.6.0)`; `viridis` dropped from Imports (unused);
  `BiocGenerics` added. `grDevices`, `graphics`, `stats`, `tools` and `utils`
  are now declared — all five were already used via `::` without being listed.
  The maintainer's ORCID iD is recorded in `Authors@R`.
- `print.vista` removed — no object could carry that class.
- Plot layer data is pinned by snapshot tests so API renames cannot silently
  change output.

# VISTA 1.1.1

## Bug fixes

- `get_volcano_plot()` now inherits `log2fc_cutoff`, `pval_cutoff`, and the new
  `p_value_type` argument from `cutoffs(x)` instead of hard-coding the raw
  `pvalue` column at 0.05. Previously the volcano could colour a different set
  of genes than `deg_summary()`, `get_deg_count_barplot()`, and
  `get_genes_by_regulation()` reported for the same object. Pass the arguments
  explicitly to restore the old behaviour
  (`p_value_type = "pvalue"`, `pval_cutoff = 0.05`). The y-axis label now
  reflects which p-value column is plotted.

- Consensus DE tables (`method = "both"`) no longer overwrite `pvalue` and
  `padj` with `1` for every gene that is not called by both backends. Each gene
  now carries the less-significant of the two backend values, so the columns
  stay continuous and usable for volcano/MA plots and ranking. Per-backend
  values remain available in `pvalue_deseq2`/`pvalue_edger`/`padj_deseq2`/
  `padj_edger`, and DEG calls are unchanged.

- `get_ma_plot()` read its threshold fallback from the non-existent
  `metadata(x)$cutoffs` key, silently using package defaults instead of the
  object's stored cutoffs.

## Reproducibility

- **Pathway gene capping is now deterministic.** `get_pathway_heatmap(max_genes)`
  and `get_enrichment_chord(max_genes)` truncated the pathway gene set
  positionally. Pathway membership is a *set*, and clusterProfiler's ordering of
  it is arbitrary and, verified here, differs between R sessions -- so the same
  call could plot a different subset of genes on different runs. The heatmap now
  ranks candidates by expression variance (most informative rows first) and the
  chord by pathway participation, both breaking ties on the gene identifier.
- **Deconvolution scores are no longer aligned to samples by position when
  their own labels disagree.** `run_cell_deconvolution()` accepted scores whose
  row count matched the sample count and stamped the object's sample names onto
  them, even when the scores carried different labels of their own -- which
  would attach every sample's cell fractions to the wrong sample. Alignment is
  now by name whenever names are available in either orientation; contradictory
  labels are an error; and position is used only when the scores carry no
  labels at all, with a message saying so.
- The airway vignette sets a seed. GSEA estimates p-values by permutation and
  ComplexHeatmap's row k-means is unseeded, so the rendered document previously
  changed on every build for reasons unrelated to the code.

## Analysis changes

- The edgeR and limma-voom backends now filter genes with the same predicate as
  the DESeq2 backend — a gene is kept when at least `min_replicates` samples
  reach `min_counts`. They previously ignored `min_counts` at this step and used
  a hard-coded `cpm > 1` threshold, so the same arguments produced different
  feature sets across backends.

- The edgeR and limma-voom backends now filter *before* calling
  `edgeR::calcNormFactors()`, per the edgeR user's guide; normalization factors
  were previously estimated from genes that were about to be discarded. Both
  changes shift edgeR/limma results slightly relative to 1.1.0.

## Documentation

- `min_replicates` is documented accurately: it counts samples across the whole
  experiment, not within each group, which is what all backends have always
  implemented.
- `create_vista()` now documents the consensus table's column semantics,
  including the `support` column.
- Removed a reference to `get_enrichment_network()` from the README; no such
  function exists. Corrected the citation block and listed the missing
  expression-plot entries.

## Testing

- Added label/value integrity tests covering every exported function that
  accepts a gene set. They assert two properties the 1.0.0 defect violated:
  results do not depend on the order genes were supplied in, and the value
  reported for a gene is that gene's own value, checked against an independent
  recomputation from the object. The previous suite asserted shape --
  `is.matrix()`, `nrow()`, class -- all of which the broken output satisfied.
- Verified by mutation testing: reintroducing the exact 1.0.0 defect now fails
  five tests across three files, and replacing a name-keyed lookup with a
  positional one elsewhere fails nine.

## Internal

- Removed duplicate definitions of `.prepare_sample_metadata()` and
  `.filter_genes()`, which existed in both `utils-internal.R` and
  `viz_related.R`, with collation order silently deciding which one ran.
- Removed unreachable helpers `.plot_pca()`, `.plot_mds()`,
  `.prepare_mds_dataframe()`, `.prepare_corr_matrix()`, `.plot_corr_heatmap()`,
  and `.cluster_log2fc_matrix()` along with their generated man pages.
- pkgdown-only articles under `vignettes/guides/` and `vignettes/workflows/`
  are no longer shipped in the source tarball.

# VISTA 1.0.0

- First Bioconductor release (Bioconductor 3.23).

# VISTA 0.99.8

# VISTA 0.99.7

# VISTA 0.99.6

# VISTA 0.99.5

- Added `read_vista_counts()`, `read_vista_metadata()`, and
  `match_vista_inputs()` to standardize common RNA-seq input formats without
  changing the existing `create_vista()` API.
- Added `derive_vista_metadata()` to bootstrap starter sample metadata from
  count-derived sample names using split- or regex-based parsing.
- Added lightweight import support for plain count tables, featureCounts,
  STAR gene counts, HTSeq-count, tximport-like inputs, and RSEM gene result
  files.

# VISTA 0.99.4

# VISTA 0.99.3

# VISTA 0.99.2

# VISTA 0.99.1

- `example_vista()` now uses a precomputed default object to reduce example,
  test, and package-check runtime while preserving the existing API.

# VISTA 0.99.0

*Submitted to Bioconductor 2026-02-11*

## Overview

VISTA (Visualization Toolkit for Transcriptomic Analysis) provides a unified S4-based framework for differential expression analysis of RNA-seq data, wrapping DESeq2 and edgeR workflows with consistent metadata management and rich visualization capabilities.

## Key Features

### Core Infrastructure
- S4 `VISTA` class extending `SummarizedExperiment` for standardized data management
- Unified differential expression workflow supporting DESeq2 and edgeR backends
- Consistent metadata structure for comparisons, cutoffs, and group information
- Flexible color palette system for visualizations

### Visualization Suite (28+ functions)
- **Dimension reduction**: PCA, MDS plots with customizable aesthetics
- **DE results**: Volcano plots, MA plots, DEG count barplots
- **Expression**: Barplots, boxplots, violin plots, density plots, joyplots, heatmaps
- **Comparisons**: Venn diagrams, alluvial plots, correlation heatmaps, pairwise plots
- **Fold-change**: Scatter plots, barplots, matrix visualizations, chromosome plots

### Functional Analysis
- MSigDB enrichment with flexible ID mapping (SYMBOL, ENSEMBL, ENTREZID)
- GO enrichment analysis (BP, MF, CC ontologies)
- KEGG pathway enrichment
- GSEA support with customizable gene sets
- Integrated visualization functions for enrichment results

### Optional Features
- Cell-type deconvolution via xCell2 integration
- Automated report generation with Quarto
- Accessor functions for all metadata components

## Implementation Details

- Comprehensive input validation and edge case handling
- Extensive test suite (>70% coverage)
- Complete roxygen2 documentation with runnable examples
- BiocStyle vignettes demonstrating complete workflows
- Proper namespace management and import declarations

## Bug Fixes

- Fixed contradictory roxygen documentation markers in internal utilities
- Added missing `@importFrom` declarations across all modules
- Improved error messages for invalid inputs
- Enhanced edge case handling in visualization functions
- Heatmap utilities now validate non-character `genes` input explicitly, support
  minimal-call defaults, and allow custom colours for multi-column annotations
- DEG count pie/donut plots now optionally include non-DE genes as an `"Other"`
  slice and support configurable label text colour
- `get_genes_by_regulation()` now supports top-gene ranking by `abs(log2fc)`
  and optional annotated table output
- PCA/MDS/UMAP plots now accept the standardized `use_group_colors` argument
  while keeping `use_vista_colors` as a deprecated compatibility alias

---
