#!/usr/bin/env Rscript
# =============================================================================
# ATAC-seq differential accessibility: code to reproduce the manuscript results
#
#   1. DiffBind differential accessibility (4 comparisons)
#   2. ChIPseeker peak annotation
#   3. Figures for a user-supplied gene list:
#        - log2 fold-change waterfall barplot
#        - heatmap of row z-scored normalized accessibility
#
# Usage:
#   Rscript ATAC_DiffBind_manuscript.R                    # uses settings below
#   Rscript ATAC_DiffBind_manuscript.R genes.txt          # override the gene list
#
# Gene list: plain text, one gene symbol per line. The heatmap rows follow
# the order of the file.
#
# Sample sheet: standard DiffBind CSV (SampleID, Factor, Replicate, bamReads,
# Peaks, PeakCaller). The values in `Factor` must match the group names
# used in `contrasts` below.
# =============================================================================

suppressPackageStartupMessages({
  library(DiffBind)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(forcats)
  library(pheatmap)
  library(viridisLite)
})

options(stringsAsFactors = FALSE)

# ---------------------------- Settings ---------------------------------------
sample_sheet <- "./sample_sheet.csv"
output_dir   <- "./results"
gene_list    <- "./genes.txt"

# Name = output prefix; c(group1, group2) -> log2FC = group1 / group2
contrasts <- list(
  A_v_WT  = c("A",  "WT"),
  D_v_WT  = c("D",  "WT"),
  A_v_D   = c("A",  "D"),
  WT_v_EV = c("WT", "EV")
)

n_top_peaks <- 50000     # consensus peaks kept for counting (ranked by mean signal)
fdr_filter  <- NULL      # e.g. 0.05 to plot only significant genes; NULL = plot all listed genes

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) gene_list <- args[1]

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
for (d in c("DE", "annotation", "figures")) {
  dir.create(file.path(output_dir, d), showWarnings = FALSE)
}

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


# =============================================================================
# 1. DiffBind differential accessibility
# =============================================================================
run_diffbind <- function(samples, contrasts, output_dir) {

  db <- dba(sampleSheet = samples)

  # PCA (peak occupancy)
  pdf(file.path(output_dir, "PCA_plot.pdf"), height = 8, width = 8)
  dba.plotPCA(db, DBA_REPLICATE, label = DBA_ID)
  dev.off()

  # Consensus peak set: rank merged peaks by mean peak score across samples
  # and keep the top n_top_peaks
  cons_df    <- dba.peakset(db, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
  count_cols <- which(!(colnames(cons_df) %in% c("CHR", "START", "END")))
  cons_df$avg_signal <- rowMeans(cons_df[, count_cols], na.rm = TRUE)
  cons_df    <- cons_df[order(cons_df$avg_signal, decreasing = TRUE), ]
  cons_df    <- cons_df[seq_len(min(n_top_peaks, nrow(cons_df))), ]

  consensus_filtered <- GRanges(
    seqnames = cons_df$CHR,
    ranges   = IRanges(cons_df$START, cons_df$END)
  )

  # Count reads in consensus peaks, normalize (DiffBind defaults)
  db <- dba.count(db, peaks = consensus_filtered, score = DBA_SCORE_READS)
  db <- dba.normalize(db)

  # Contrasts
  for (con in contrasts) {
    db <- dba.contrast(
      db,
      minMembers = 2,
      group1 = db$masks[[con[1]]],
      group2 = db$masks[[con[2]]],
      name1  = con[1],
      name2  = con[2]
    )
  }

  db <- dba.analyze(db, method = DBA_ALL_METHODS)
  db
}

# All peaks for one contrast (th = 1), edgeR statistics, normalized counts per sample
extract_diffbind_results <- function(db, contrast_index, outname) {

  res <- dba.report(
    db,
    contrast    = contrast_index,
    th          = 1,
    method      = DBA_EDGER,
    bCounts     = TRUE,
    bNormalized = TRUE,
    DataType    = DBA_DATA_FRAME
  )

  if (is.null(res) || nrow(res) == 0) {
    warning("No results for contrast ", contrast_index)
    return(NULL)
  }

  df <- res[complete.cases(res[, c("Start", "End")]), ]
  df <- df[df$End > df$Start, ]
  df <- df[!is.na(df$Fold), ]

  write.csv(df, outname, row.names = FALSE)
  df
}


# =============================================================================
# 2. Peak annotation (ChIPseeker, hg38 knownGene, TSS +/- 2 kb)
# =============================================================================
annotate_peak_df <- function(df, name, output_dir) {

  gr <- makeGRangesFromDataFrame(df)

  anno <- annotatePeak(
    gr,
    TxDb      = txdb,
    tssRegion = c(-2000, 2000),
    annoDb    = "org.Hs.eg.db"
  )

  anno_df <- as.data.frame(anno)

  anno_df$annotation <- ifelse(
    startsWith(anno_df$annotation, "Intron"), "Intron",
    ifelse(startsWith(anno_df$annotation, "Exon"), "Exon", anno_df$annotation)
  )

  ## Genome distribution pie chart
  p1 <- as.data.frame(table(anno_df$annotation)) %>%
    mutate(
      csum = rev(cumsum(rev(Freq))),
      pos  = Freq / 2 + dplyr::lead(csum, 1),
      pos  = if_else(is.na(pos), Freq / 2, pos),
      per  = round((Freq / sum(Freq)) * 100, 2)
    )

  pie <- ggplot(p1, aes(x = "", y = Freq, fill = fct_inorder(Var1))) +
    geom_col(width = 1, color = 1) +
    coord_polar(theta = "y") +
    geom_label_repel(aes(y = pos, label = paste0(per, "%")),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    theme_void() +
    theme(legend.title = element_blank())

  pdf(file.path(output_dir, paste0(name, "_genome_distribution.pdf")), height = 8, width = 8)
  print(pie)
  dev.off()

  final <- merge(df, anno_df, by = c(1, 2, 3))
  if ("FDR" %in% colnames(final)) final <- final[order(final$FDR), ]

  write.csv(final, file.path(output_dir, paste0(name, "_annotated.csv")), row.names = FALSE)
  final
}


# =============================================================================
# 3. Gene-list figures
#    One value per gene = the gene's most significant peak (lowest FDR)
# =============================================================================
select_gene_peaks <- function(anno_df, genes, fdr_filter = NULL) {

  sel <- anno_df %>%
    filter(SYMBOL %in% genes, !is.na(FDR)) %>%
    group_by(SYMBOL) %>%
    slice_min(FDR, n = 1, with_ties = FALSE) %>%
    ungroup()

  if (!is.null(fdr_filter)) sel <- sel %>% filter(FDR < fdr_filter)

  missing <- setdiff(genes, sel$SYMBOL)
  if (length(missing) > 0) {
    message("  Not plotted (no annotated peak", if (!is.null(fdr_filter)) " passing FDR", "): ",
            paste(missing, collapse = ", "))
  }

  # keep gene-list order
  sel[order(match(sel$SYMBOL, genes)), ]
}

plot_fc_waterfall <- function(sel, group1, group2, title = "") {

  df <- sel %>%
    mutate(
      direction = ifelse(Fold > 0, paste("Higher in", group1), paste("Higher in", group2)),
      direction = factor(direction, levels = c(paste("Higher in", group2),
                                               paste("Higher in", group1))),
      SYMBOL    = reorder(SYMBOL, Fold)
    )

  ggplot(df, aes(x = Fold, y = SYMBOL, fill = direction)) +
    geom_col(width = 0.8) +
    geom_vline(xintercept = 0, color = "black") +
    scale_fill_manual(values = setNames(c("#1F1F8F", "#E41A1C"), levels(df$direction)),
                      drop = FALSE) +
    labs(x = expression(log[2] ~ "fold change"), y = NULL, fill = NULL, title = title) +
    theme_classic() +
    theme(axis.text.y     = element_text(face = "italic"),
          legend.position = "top",
          plot.title      = element_text(hjust = 0.5))
}

plot_gene_heatmap <- function(sel, sample_cols, out_pdf, title = "") {

  mat <- as.matrix(sel[, sample_cols])
  rownames(mat) <- sel$SYMBOL

  # log2 normalized counts, then z-score each gene across samples
  mat <- log2(mat + 1)
  z   <- t(scale(t(mat)))
  z   <- z[stats::complete.cases(z), , drop = FALSE]   # drops genes with zero variance

  pdf(out_pdf, width = 5, height = max(4, 0.25 * nrow(z) + 2))
  pheatmap(
    z,
    color          = viridis(100),
    cluster_rows   = FALSE,
    cluster_cols   = FALSE,
    border_color   = NA,
    fontsize_row   = 9,
    main           = title
  )
  dev.off()

  z
}


# =============================================================================
# Run
# =============================================================================
samples <- read.csv(sample_sheet)

# --- 1. Differential accessibility -------------------------------------------
db <- run_diffbind(samples, contrasts, output_dir)
saveRDS(db, file.path(output_dir, "dba_object.rds"))
capture.output(dba.show(db), dba.show(db, bContrasts = TRUE),
               file = file.path(output_dir, "dba_summary.txt"))

results_list <- list()
for (i in seq_along(contrasts)) {
  nm <- names(contrasts)[i]
  results_list[[nm]] <- extract_diffbind_results(
    db, contrast_index = i,
    outname = file.path(output_dir, "DE", paste0(nm, "_DE_peaks.csv"))
  )
}

# Counts in all consensus peaks for every sample
counts_df <- dba.peakset(db, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)
write.csv(counts_df, file.path(output_dir, "all_peak_counts.csv"), row.names = FALSE)

# --- 2. Annotation -------------------------------------------------------------
anno_list <- list()
for (nm in names(results_list)) {
  df <- results_list[[nm]]
  if (is.null(df)) next
  colnames(df)[1:3] <- c("chr", "start", "end")
  anno_list[[nm]] <- annotate_peak_df(df, paste0(nm, "_DE_peaks"),
                                      file.path(output_dir, "annotation"))
}

# Summary of significant peaks per contrast
summary_tab <- do.call(rbind, lapply(names(results_list), function(nm) {
  r <- results_list[[nm]]
  data.frame(contrast = nm, total_peaks = nrow(r),
             FDR05_up   = sum(r$FDR < 0.05 & r$Fold > 0),
             FDR05_down = sum(r$FDR < 0.05 & r$Fold < 0))
}))
write.csv(summary_tab, file.path(output_dir, "DE_summary.csv"), row.names = FALSE)
print(summary_tab)

# --- 3. Gene-list figures --------------------------------------------------------
if (file.exists(gene_list)) {
  genes <- unique(trimws(readLines(gene_list)))
  genes <- genes[genes != ""]
  message("Gene list: ", length(genes), " genes")

  for (nm in names(anno_list)) {
    message(nm)
    g1 <- contrasts[[nm]][1]
    g2 <- contrasts[[nm]][2]
    sel <- select_gene_peaks(anno_list[[nm]], genes, fdr_filter)
    if (nrow(sel) == 0) next

    # Samples in this comparison, group1 first
    sample_cols <- c(samples$SampleID[samples$Factor == g1],
                     samples$SampleID[samples$Factor == g2])
    sample_cols <- intersect(sample_cols, colnames(sel))

    # Table behind both figures
    write.csv(sel %>% dplyr::select(SYMBOL, chr, start, end, annotation, distanceToTSS,
                                    Conc, Fold, dplyr::any_of(c("p-value", "p.value")),
                                    FDR, dplyr::all_of(sample_cols)),
              file.path(output_dir, "figures", paste0(nm, "_genelist_peaks.csv")),
              row.names = FALSE)

    p <- plot_fc_waterfall(sel, g1, g2, title = nm)
    ggsave(file.path(output_dir, "figures", paste0(nm, "_FC_waterfall.pdf")), p,
           width = 4.5, height = max(4, 0.22 * nrow(sel) + 1.5))

    plot_gene_heatmap(sel, sample_cols,
                      file.path(output_dir, "figures", paste0(nm, "_genelist_heatmap.pdf")),
                      title = nm)
  }
} else {
  message("Gene list not found (", gene_list, "): skipping figures.")
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
