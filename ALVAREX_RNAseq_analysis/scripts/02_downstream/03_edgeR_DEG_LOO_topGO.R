############################################################
##
## RNAseq Protection Pipeline — v7 topGO + improved viz
## ---------------------------------------------------------
## Main objective:
## Identify transcriptomic signatures associated with the
## protected phenotype in treated carrot genotypes.
##
## CHANGELOG v6 -> v7:
## - GO enrichment (Parts 6 & 8) reverted to topGO (weight01 Fisher)
## - Visualizations completely redesigned (clusterProfiler-style):
##     dotplot (GeneRatio × -log10p, size = DEG_Count)
##     barplot (-log10p, fill = FDR)
##     combined bubble (Ontology × Direction / GeneSet)
## - GeneRatio, FDR, DEG_Genes columns added to topGO output
## - Removed clusterProfiler / PlanT2T / OrgDb dependencies
##
############################################################


# ============================================================
# PART 1 — GLOBAL SETTINGS
# ============================================================

## ------------------------------------------------------------
## 1.1 Packages
## ------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_pkgs <- c(
  "ggplot2", "ggrepel", "dplyr", "tibble", "tidyr", "forcats",
  "stringr", "svglite", "openxlsx", "patchwork", "scales",
  "grid", "gridExtra"
)

bioc_pkgs <- c(
  "edgeR", "limma", "pheatmap",
  "topGO", "WGCNA"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

message("✅ Packages loaded")


## ------------------------------------------------------------
## 1.2 Global parameters
## ------------------------------------------------------------
SEED          <- 123
alpha         <- 0.05
lfc_threshold <- 1

set.seed(SEED)
options(stringsAsFactors = FALSE)
options(scipen = 999)
options(error = NULL)   ## reset any recover/debug mode from the session

message("✅ Global parameters set")


## ------------------------------------------------------------
## 1.3 Analysis switches
## ------------------------------------------------------------
RUN_QC          <- TRUE
RUN_PRIMARY_DEG <- TRUE
RUN_GO_PRIMARY  <- TRUE
RUN_LOO         <- TRUE
RUN_WGCNA       <- TRUE

RUN_PAIRWISE    <- FALSE
RUN_CROSSN      <- FALSE
RUN_VENN        <- FALSE
RUN_UPSET       <- FALSE
RUN_CONCORDANCE <- FALSE

message("✅ Analysis switches defined")


## ------------------------------------------------------------
## 1.4 USER CONFIGURATION
## ------------------------------------------------------------
## For security and portability, cluster-specific paths are not
## stored in this public script. Replace the placeholders below.
counts_file   <- "/path/to/input/Count.csv"
meta_file     <- "/path/to/input/MetaData.csv"

## Custom gene-to-GO mapping used by topGO
go_annot_file <- "/path/to/annotation/gene_to_GO.csv"

## topGO parameters
GO_P_CUTOFF   <- 0.05
GO_MIN_GENES  <- 4
TOP_N_TERMS   <- 10    # terms per ontology in plots
ontologies    <- c("BP", "MF", "CC")

message("\u2705 Input paths defined")


## ------------------------------------------------------------
## 1.5 Output structure
## ------------------------------------------------------------
main_dir <- "/path/to/results/edgeR_LOO_topGO"
dir_qc    <- file.path(main_dir, "01_QC")
dir_deg   <- file.path(main_dir, "02_DEG")
dir_go    <- file.path(main_dir, "03_GO")
dir_wgcna <- file.path(main_dir, "04_WGCNA")
dir_supp  <- file.path(main_dir, "05_Supplementary")

for (d in c(main_dir, dir_qc, dir_deg, dir_go, dir_wgcna, dir_supp)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

dir_pairwise    <- file.path(dir_supp, "Pairwise")
dir_crossn      <- file.path(dir_supp, "CrossN")
dir_venn        <- file.path(dir_supp, "Venn")
dir_upset       <- file.path(dir_supp, "UpSet")
dir_concordance <- file.path(dir_supp, "Concordance")

if (RUN_PAIRWISE)    dir.create(dir_pairwise,    recursive = TRUE, showWarnings = FALSE)
if (RUN_CROSSN)      dir.create(dir_crossn,      recursive = TRUE, showWarnings = FALSE)
if (RUN_VENN)        dir.create(dir_venn,        recursive = TRUE, showWarnings = FALSE)
if (RUN_UPSET)       dir.create(dir_upset,       recursive = TRUE, showWarnings = FALSE)
if (RUN_CONCORDANCE) dir.create(dir_concordance, recursive = TRUE, showWarnings = FALSE)

message("✅ Output folders created")


## ------------------------------------------------------------
## 1.6 Color system
## ------------------------------------------------------------
GENOTYPE_PALETTE <- c(
  "Amsterdam_N"      = "#0A9E6E",
  "NantaiseInbred_N" = "#1B4FBF",
  "Orleans_N"        = "#007B8A",
  "Dijon_N"          = "#3949AB",
  "Genevieve_N"      = "#006B3C",
  "Deep_Purple_P"    = "#D62728",
  "Presto_P"         = "#E8820C",
  "Neva_P"           = "#9B2ECA",
  "Robila_P"         = "#C9A800",
  "Genevieve_P"      = "#F57F17",
  "Orleans_P"        = "#5E9E1A",
  "Dijon_P"          = "#C2185B",
  "Oxhella_P"        = "#8B0000"
)

PROTECTION_COLORS <- c(
  "Protected"    = "#C0392B",
  "NonProtected" = "#2471A3"
)

DEG_COLORS <- c(
  "Up"   = "#C0392B",
  "Down" = "#2471A3",
  "NS"   = "#D5D8DC"
)

ONT_COLORS <- c(
  "BP" = "#E76F51",
  "MF" = "#2A9D8F",
  "CC" = "#5C4B8A"
)

message("✅ Color system defined")


## ------------------------------------------------------------
## 1.7 Publication theme
## ------------------------------------------------------------
theme_pub <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.border      = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.8),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.background  = ggplot2::element_rect(fill = "white"),
      axis.line         = ggplot2::element_blank(),
      axis.ticks        = ggplot2::element_line(colour = "black", linewidth = 0.6),
      axis.ticks.length = ggplot2::unit(3, "mm"),
      axis.text         = ggplot2::element_text(colour = "black", size = base_size + 1),
      axis.title        = ggplot2::element_text(colour = "black", size = base_size + 2, face = "bold"),
      plot.title        = ggplot2::element_text(size = base_size + 4, face = "bold", hjust = 0),
      plot.subtitle     = ggplot2::element_text(size = base_size + 1, colour = "grey40"),
      plot.caption      = ggplot2::element_text(size = base_size - 1, colour = "grey55"),
      legend.title      = ggplot2::element_text(size = base_size + 1, face = "bold"),
      legend.text       = ggplot2::element_text(size = base_size),
      legend.key.size   = ggplot2::unit(5, "mm"),
      strip.text        = ggplot2::element_text(face = "bold", size = base_size + 1),
      strip.background  = ggplot2::element_rect(fill = "grey94", colour = "black", linewidth = 0.4)
    )
}

message("✅ Publication theme ready")


## ------------------------------------------------------------
## 1.8 Export helpers
## ------------------------------------------------------------
save_plot_pdf_svg <- function(plot_obj, file_base, width = 7, height = 5, dpi = 600) {
  ## PDF
  ggplot2::ggsave(
    filename = paste0(file_base, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
  
  ## SVG
  svglite::svglite(file = paste0(file_base, ".svg"), width = width, height = height)
  print(plot_obj)
  grDevices::dev.off()
  
  ## PNG 600 dpi
  ggplot2::ggsave(
    filename = paste0(file_base, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = "white"
  )
}
save_base_plot_pdf_svg <- function(plot_fn, file_base, width = 7, height = 5, dpi = 600) {
  ## PDF
  grDevices::pdf(paste0(file_base, ".pdf"), width = width, height = height, useDingbats = FALSE)
  plot_fn()
  grDevices::dev.off()
  
  ## SVG
  svglite::svglite(paste0(file_base, ".svg"), width = width, height = height)
  plot_fn()
  grDevices::dev.off()
  
  ## PNG 600 dpi
  grDevices::png(
    filename = paste0(file_base, ".png"),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white"
  )
  plot_fn()
  grDevices::dev.off()
}
save_pheatmap_pdf_svg <- function(ph_args, file_base, w = 9, h = 8, dpi = 600) {
  ## PDF
  do.call(
    pheatmap::pheatmap,
    c(ph_args, list(filename = paste0(file_base, ".pdf"), width = w, height = h))
  )
  
  ## SVG
  ph_obj <- do.call(pheatmap::pheatmap, c(ph_args, list(silent = TRUE)))
  svglite::svglite(paste0(file_base, ".svg"), width = w, height = h)
  grid::grid.newpage()
  grid::grid.draw(ph_obj$gtable)
  grDevices::dev.off()
  
  ## PNG 600 dpi
  grDevices::png(
    filename = paste0(file_base, ".png"),
    width = w,
    height = h,
    units = "in",
    res = dpi,
    bg = "white"
  )
  grid::grid.newpage()
  grid::grid.draw(ph_obj$gtable)
  grDevices::dev.off()
}

## ------------------------------------------------------------
## 1.9  topGO helpers + publication-quality visualization
## ------------------------------------------------------------

## ── topGO enrichment core ─────────────────────────────────────────────────

## Parse GO annotation file -> named list geneID2GO
load_go_annotation <- function(go_annot_file) {
  go_raw <- read.csv(go_annot_file, header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)
  colnames(go_raw) <- trimws(colnames(go_raw))
  
  gene_col <- intersect(c("Gene_ID", "gene", "Gene"), colnames(go_raw))[1]
  go_col   <- intersect(c("Annotation", "GO", "go", "GO_ID"), colnames(go_raw))[1]
  
  if (is.na(gene_col) || is.na(go_col))
    stop("Cannot find Gene/GO columns in GO annotation file.")
  
  go_df <- go_raw[, c(gene_col, go_col)] |>
    dplyr::rename(gene_id = 1, go_term = 2) |>
    dplyr::mutate(
      gene_id = trimws(as.character(gene_id)),
      go_term = trimws(as.character(go_term))
    ) |>
    tidyr::separate_rows(go_term, sep = "[,;\\s]+") |>
    dplyr::filter(grepl("^GO:\\d+$", go_term)) |>
    dplyr::distinct()
  
  cat("Valid gene-GO pairs:", nrow(go_df), "\n")
  split(go_df$go_term, go_df$gene_id)
}

## Run topGO for one gene set × one ontology
## Returns data.frame with standardised columns
run_topgo_one <- function(genes_oi, gene_universe, geneID2GO, ont,
                          go_p_cutoff = GO_P_CUTOFF,
                          min_genes   = GO_MIN_GENES) {
  
  genes_oi <- intersect(genes_oi, gene_universe)
  if (length(genes_oi) < min_genes) return(NULL)
  
  geneList <- factor(as.integer(gene_universe %in% genes_oi), levels = c(0, 1))
  names(geneList) <- gene_universe
  
  GOdata <- tryCatch(
    new("topGOdata",
        ontology = ont,
        allGenes = geneList,
        annot    = annFUN.gene2GO,
        gene2GO  = geneID2GO),
    error = function(e) NULL
  )
  if (is.null(GOdata)) return(NULL)
  
  res_w01  <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  res_clas <- runTest(GOdata, algorithm = "classic",  statistic = "fisher")
  
  allRes <- GenTable(
    GOdata,
    weight01 = res_w01,
    classic  = res_clas,
    topNodes = length(score(res_w01)),
    numChar  = 120
  )
  
  clean_p <- function(x) suppressWarnings(as.numeric(gsub("< ", "", x)))
  allRes$weight01 <- clean_p(allRes$weight01)
  allRes$classic  <- clean_p(allRes$classic)
  allRes$weight01[is.na(allRes$weight01) | allRes$weight01 == 0] <- 1e-300
  allRes$classic [is.na(allRes$classic)  | allRes$classic  == 0] <- 1e-300
  allRes$FDR      <- p.adjust(allRes$weight01, method = "BH")
  
  allRes$DEG_Count <- vapply(allRes$GO.ID, function(gid) {
    tryCatch(sum(genesInTerm(GOdata, gid)[[1]] %in% genes_oi),
             error = function(e) 0L)
  }, integer(1))
  
  allRes$DEG_Genes <- vapply(allRes$GO.ID, function(gid) {
    tryCatch(paste(sort(intersect(genesInTerm(GOdata, gid)[[1]], genes_oi)),
                   collapse = ";"),
             error = function(e) "")
  }, character(1))
  
  ## GeneRatio mirrors clusterProfiler: DEG_Count / |genes_oi in universe|
  n_query <- length(genes_oi)
  allRes$GeneRatio <- allRes$DEG_Count / n_query
  allRes$negLogP   <- -log10(allRes$weight01)
  
  out <- allRes[allRes$weight01 < go_p_cutoff & allRes$DEG_Count >= min_genes, ]
  if (nrow(out) == 0) return(NULL)
  out
}

## Run topGO for all 3 ontologies for one gene set
run_topgo_all_ont <- function(genes_oi, gene_universe, geneID2GO,
                              ontologies  = c("BP", "MF", "CC"),
                              go_p_cutoff = GO_P_CUTOFF,
                              min_genes   = GO_MIN_GENES) {
  results <- list()
  for (ont in ontologies) {
    res <- run_topgo_one(genes_oi, gene_universe, geneID2GO,
                         ont, go_p_cutoff, min_genes)
    if (!is.null(res) && nrow(res) > 0) {
      res$Ontology <- ont
      results[[ont]] <- res
      cat("  topGO", ont, ":", nrow(res), "significant terms\n")
    } else {
      cat("  topGO", ont, ": no significant terms\n")
    }
  }
  if (length(results) == 0) return(NULL)
  dplyr::bind_rows(results)
}

## ── GO visualization helpers (simple clusterProfiler-like style) ─────────────

plot_go_dot_simple <- function(go_df,
                               title = NULL,
                               top_n = 10,
                               label_max = 55,
                               base_size = 10) {
  
  if (is.null(go_df) || nrow(go_df) == 0) return(NULL)
  
  df <- go_df |>
    dplyr::mutate(
      FDR_plot = dplyr::if_else(!is.na(FDR) & FDR > 0, FDR, weight01),
      score    = -log10(FDR_plot),
      Term_lab = ifelse(
        nchar(Term) > label_max,
        paste0(substr(Term, 1, label_max - 3), "..."),
        Term
      ),
      Ontology = factor(Ontology, levels = c("BP", "MF", "CC"))
    ) |>
    dplyr::group_by(Ontology) |>
    dplyr::slice_min(order_by = weight01, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()
  
  if (nrow(df) == 0) return(NULL)
  
  df <- df |>
    dplyr::group_by(Ontology) |>
    dplyr::arrange(GeneRatio, .by_group = TRUE) |>
    dplyr::mutate(Term_lab = factor(Term_lab, levels = unique(Term_lab))) |>
    dplyr::ungroup()
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = GeneRatio, y = Term_lab, size = DEG_Count, colour = score)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_wrap(~ Ontology, scales = "free_y", ncol = 1) +
    ggplot2::scale_colour_gradient(
      low = "#6BAED6",
      high = "#CB181D",
      name = expression(-log[10](italic(FDR)))
    ) +
    ggplot2::scale_size_continuous(
      range = c(3, 9),
      name = "Gene count"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(accuracy = 0.1),
      expand = ggplot2::expansion(mult = c(0.02, 0.10))
    ) +
    ggplot2::labs(
      title = title,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_pub(base_size = base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = base_size - 1),
      strip.text  = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}


plot_go_dot_grouped <- function(go_df,
                                title = NULL,
                                top_n = 6,
                                label_max = 45,
                                base_size = 9) {
  
  if (is.null(go_df) || nrow(go_df) == 0) return(NULL)
  if (!"GeneSet" %in% colnames(go_df)) stop("Column 'GeneSet' is required.")
  
  df <- go_df |>
    dplyr::mutate(
      FDR_plot = dplyr::if_else(!is.na(FDR) & FDR > 0, FDR, weight01),
      score    = -log10(FDR_plot),
      Term_lab = ifelse(
        nchar(Term) > label_max,
        paste0(substr(Term, 1, label_max - 3), "..."),
        Term
      ),
      Ontology = factor(Ontology, levels = c("BP", "MF", "CC"))
    ) |>
    dplyr::group_by(GeneSet, Ontology) |>
    dplyr::slice_min(order_by = weight01, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()
  
  if (nrow(df) == 0) return(NULL)
  
  df <- df |>
    dplyr::arrange(GeneSet, Ontology, GeneRatio) |>
    dplyr::mutate(
      Term_key = paste(GeneSet, Ontology, Term_lab, sep = "___"),
      Term_key = factor(Term_key, levels = unique(Term_key))
    )
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = GeneRatio, y = Term_key, size = DEG_Count, colour = score)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(Ontology ~ GeneSet, scales = "free_y", space = "free_y") +
    ggplot2::scale_y_discrete(labels = function(x) sub("^.*___.*___", "", x)) +
    ggplot2::scale_colour_gradient(
      low = "#6BAED6",
      high = "#CB181D",
      name = expression(-log[10](italic(FDR)))
    ) +
    ggplot2::scale_size_continuous(
      range = c(2.5, 8),
      name = "Gene count"
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_percent(accuracy = 0.1),
      expand = ggplot2::expansion(mult = c(0.02, 0.10))
    ) +
    ggplot2::labs(
      title = title,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_pub(base_size = base_size) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = base_size - 1),
      strip.text  = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

message("\u2705 topGO helpers + publication GO visualization functions ready")


cat("\n----------------------------------------\n")
cat("PART 1 completed\n")
cat("----------------------------------------\n")
cat("Main directory :", main_dir, "\n")
cat("QC directory   :", dir_qc, "\n")
cat("DEG directory  :", dir_deg, "\n")
cat("GO directory   :", dir_go, "\n")
cat("WGCNA directory:", dir_wgcna, "\n")
cat("----------------------------------------\n\n")



# ============================================================
# PART 2 — DATA IMPORT AND CLEANING
# ============================================================

cat("---- PART 2: data import and cleaning ----\n")


## ------------------------------------------------------------
## 2.1 Helper functions
## ------------------------------------------------------------
normalize_sample_name <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^X", "", x)
  x <- gsub("\"", "", x)
  x <- gsub("\\s+", "", x)
  x <- gsub("\r|\n", "", x)
  x <- gsub("_F-([0-9]{4}clean_mapping)$", "_F_\\1", x)
  x
}

fix_encoding <- function(x) {
  x2 <- iconv(as.character(x), from = "", to = "UTF-8", sub = "")
  x2[is.na(x2)] <- x[is.na(x2)]
  x2
}

resolve_genotype_colors <- function(genotype_levels) {
  known   <- genotype_levels[genotype_levels %in% names(GENOTYPE_PALETTE)]
  unknown <- genotype_levels[!genotype_levels %in% names(GENOTYPE_PALETTE)]
  
  if (length(unknown) > 0) {
    hues <- seq(15, 375, length.out = length(unknown) + 1)[seq_len(length(unknown))]
    fallback <- setNames(grDevices::hsv(hues / 360, s = 0.82, v = 0.75), unknown)
    warning(length(unknown),
            " genotype(s) absent from GENOTYPE_PALETTE -> auto-colors assigned: ",
            paste(unknown, collapse = ", "),
            call. = FALSE)
    return(c(GENOTYPE_PALETTE[known], fallback)[genotype_levels])
  }
  GENOTYPE_PALETTE[genotype_levels]
}

message("✅ Helper functions loaded")


## ------------------------------------------------------------
## 2.2 Check input files
## ------------------------------------------------------------
stopifnot(file.exists(counts_file))
stopifnot(file.exists(meta_file))

if (RUN_GO_PRIMARY && !file.exists(go_annot_file)) {
  warning("GO annotation file not found: ", go_annot_file,
          "\nGO analyses will be skipped.", call. = FALSE)
  RUN_GO_PRIMARY <- FALSE
}

message("✅ Input files found")


## ------------------------------------------------------------
## 2.3 Import count matrix
## ------------------------------------------------------------
counts_data <- read.csv(counts_file, header = TRUE, check.names = FALSE)

if (!"Geneid" %in% colnames(counts_data)) {
  stop("The count matrix must contain a column named 'Geneid'.")
}

rownames(counts_data) <- counts_data$Geneid
counts_data$Geneid <- NULL
counts_data <- as.matrix(counts_data)
mode(counts_data) <- "numeric"
colnames(counts_data) <- normalize_sample_name(colnames(counts_data))

cat("Count matrix dimensions:", nrow(counts_data), "genes x", ncol(counts_data), "samples\n")


## ------------------------------------------------------------
## 2.4 Import metadata
## ------------------------------------------------------------
sample_info <- read.csv(meta_file, header = TRUE, check.names = FALSE)

required_meta_cols <- c("Sample", "Genotype", "Replicate")
missing_meta_cols  <- setdiff(required_meta_cols, colnames(sample_info))
if (length(missing_meta_cols) > 0) {
  stop("Missing required metadata column(s): ", paste(missing_meta_cols, collapse = ", "))
}

sample_info <- sample_info |>
  dplyr::mutate(
    Sample    = normalize_sample_name(Sample),
    Genotype  = fix_encoding(trimws(Genotype)),
    Replicate = as.integer(Replicate)
  )

cat("Metadata rows:", nrow(sample_info), "\n")


## ------------------------------------------------------------
## 2.5 Match counts and metadata
## ------------------------------------------------------------
common_samples <- intersect(colnames(counts_data), sample_info$Sample)

cat("Matched samples:", length(common_samples), "\n")

if (length(common_samples) == 0) {
  stop("No common sample names found between count matrix and metadata.")
}

sample_info <- sample_info |> dplyr::filter(Sample %in% common_samples)
counts_data <- counts_data[, common_samples, drop = FALSE]
sample_info <- sample_info |>
  dplyr::arrange(match(Sample, colnames(counts_data)))
rownames(sample_info) <- sample_info$Sample

stopifnot(all(colnames(counts_data) == sample_info$Sample))
message("✅ Count matrix and metadata aligned")


## ------------------------------------------------------------
## 2.6 Derive biological variables
## ------------------------------------------------------------
sample_info <- sample_info |>
  dplyr::mutate(
    Protection = dplyr::case_when(
      grepl("_P$", Genotype) ~ "Protected",
      grepl("_N$", Genotype) ~ "NonProtected",
      TRUE ~ NA_character_
    )
  )

if (any(is.na(sample_info$Protection))) {
  warning("Samples without _P or _N suffix detected and removed.", call. = FALSE)
  sample_info <- sample_info |> dplyr::filter(!is.na(Protection))
  counts_data <- counts_data[, sample_info$Sample, drop = FALSE]
  rownames(sample_info) <- sample_info$Sample
}

sample_info <- sample_info |>
  dplyr::mutate(
    Protection = factor(Protection, levels = c("NonProtected", "Protected")),
    Genotype   = factor(Genotype,   levels = sort(unique(Genotype))),
    Replicate  = factor(Replicate,  levels = sort(unique(Replicate)))
  )

conditions_P <- levels(sample_info$Genotype)[grepl("_P$", levels(sample_info$Genotype))]
conditions_N <- levels(sample_info$Genotype)[grepl("_N$", levels(sample_info$Genotype))]
n_P <- length(conditions_P)
n_N <- length(conditions_N)

if (n_P == 0 || n_N == 0) {
  stop("Both protected (_P) and non-protected (_N) genotypes are required.")
}

cat("\nProtected genotypes (_P):\n");    print(conditions_P)
cat("\nNon-protected genotypes (_N):\n"); print(conditions_N)
cat("\nProtection distribution:\n");     print(table(sample_info$Protection))
cat("\nGenotype distribution:\n");       print(table(sample_info$Genotype))


## ------------------------------------------------------------
## 2.7 Resolve active genotype colors
## ------------------------------------------------------------
genotype_colors <- resolve_genotype_colors(levels(sample_info$Genotype))

write.csv(
  data.frame(Genotype = names(genotype_colors), Color = unname(genotype_colors),
             row.names = NULL),
  file.path(main_dir, "Active_Genotype_Colors.csv"),
  row.names = FALSE
)

message("✅ Active genotype colors resolved")


## ------------------------------------------------------------
## 2.8 Export cleaned metadata
## ------------------------------------------------------------
write.csv(sample_info, file.path(main_dir, "Cleaned_Metadata.csv"), row.names = TRUE)
message("✅ Cleaned metadata exported")


## ------------------------------------------------------------
## 2.9 Final sanity checks
## ------------------------------------------------------------
stopifnot(ncol(counts_data) == nrow(sample_info))
stopifnot(all(colnames(counts_data) == rownames(sample_info)))
stopifnot(all(!is.na(sample_info$Protection)))
stopifnot(all(!is.na(sample_info$Genotype)))
stopifnot(all(!is.na(sample_info$Replicate)))

cat("\n----------------------------------------\n")
cat("PART 2 completed\n")
cat("----------------------------------------\n")
cat("Genes   :", nrow(counts_data), "\n")
cat("Samples :", ncol(counts_data), "\n")
cat("Protected genotypes    :", n_P, "\n")
cat("Non-protected genotypes:", n_N, "\n")
cat("----------------------------------------\n\n")



# ============================================================
# PART 3 — edgeR MAIN MODEL
# ============================================================

cat("---- PART 3: edgeR main model ----\n")


## ------------------------------------------------------------
## 3.1 Build DGEList
## ------------------------------------------------------------
dge <- edgeR::DGEList(counts = round(counts_data))

cat("Initial DGEList:\n")
cat("  Genes   :", nrow(dge), "\n")
cat("  Samples :", ncol(dge), "\n\n")


## ------------------------------------------------------------
## 3.2 Filter lowly expressed genes
## ------------------------------------------------------------
design_for_filter <- model.matrix(~ Protection + Replicate, data = sample_info)
keep <- edgeR::filterByExpr(dge, design = design_for_filter)

filter_summary <- data.frame(
  Step         = c("Before filtering", "After filtering"),
  N_genes      = c(nrow(dge), sum(keep)),
  Percent_genes = c(100, round(100 * sum(keep) / nrow(dge), 2))
)

dge <- dge[keep, , keep.lib.sizes = FALSE]

write.csv(filter_summary,
          file.path(dir_deg, "Gene_Filtering_Summary.csv"),
          row.names = FALSE)

cat("Filtering summary:\n")
print(filter_summary)
cat("\n")


## ------------------------------------------------------------
## 3.3 Library size normalization
## ------------------------------------------------------------
dge <- edgeR::calcNormFactors(dge, method = "TMMwsp")

norm_factors_df <- data.frame(
  Sample        = colnames(dge),
  LibrarySize   = dge$samples$lib.size,
  NormFactor    = dge$samples$norm.factors,
  EffectiveSize = dge$samples$lib.size * dge$samples$norm.factors,
  Genotype      = sample_info[colnames(dge), "Genotype"],
  Protection    = sample_info[colnames(dge), "Protection"],
  Replicate     = sample_info[colnames(dge), "Replicate"],
  row.names     = NULL
)

write.csv(norm_factors_df,
          file.path(dir_deg, "Normalization_Factors.csv"),
          row.names = FALSE)

cat("Normalization factors:\n")
print(norm_factors_df[, c("Sample", "NormFactor")])
cat("\n")


## ------------------------------------------------------------
## 3.4 Main design matrix
## ------------------------------------------------------------
design <- model.matrix(~ 0 + Genotype + Replicate, data = sample_info)
colnames(design) <- gsub("^Genotype", "", colnames(design))
colnames(design) <- make.names(colnames(design))

geno_map <- setNames(
  make.names(levels(sample_info$Genotype)),
  levels(sample_info$Genotype)
)

design_df <- as.data.frame(design)
design_df$Sample <- rownames(sample_info)
write.csv(design_df, file.path(dir_deg, "Design_Matrix.csv"), row.names = FALSE)

cat("Design columns:\n")
print(colnames(design))
cat("\n")


## ------------------------------------------------------------
## 3.5 Estimate dispersions
## ------------------------------------------------------------
dge <- edgeR::estimateDisp(dge, design, robust = TRUE)

dispersion_summary <- data.frame(
  CommonDispersion  = dge$common.dispersion,
  CommonBCV         = sqrt(dge$common.dispersion),
  AveTagwiseDisp    = mean(dge$tagwise.dispersion, na.rm = TRUE),
  MedianTagwiseDisp = median(dge$tagwise.dispersion, na.rm = TRUE)
)

write.csv(dispersion_summary,
          file.path(dir_deg, "Dispersion_Summary.csv"),
          row.names = FALSE)

cat("Dispersion summary:\n")
print(dispersion_summary)
cat("\n")


## ------------------------------------------------------------
## 3.6 Fit quasi-likelihood model
## ------------------------------------------------------------
fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
message("✅ Quasi-likelihood model fitted")


## ------------------------------------------------------------
## 3.7 Precompute useful matrices
## ------------------------------------------------------------
logcpm_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
avg_log_cpm <- edgeR::aveLogCPM(dge)

write.csv(
  data.frame(Gene = rownames(dge), aveLogCPM = avg_log_cpm),
  file.path(dir_deg, "Average_logCPM.csv"),
  row.names = FALSE
)

message("✅ logCPM matrix and aveLogCPM computed")


## ------------------------------------------------------------
## 3.8 Save compact model summary
## ------------------------------------------------------------
model_summary <- list(
  n_genes_after_filtering  = nrow(dge),
  n_samples                = ncol(dge),
  n_protected_genotypes    = n_P,
  n_nonprotected_genotypes = n_N,
  design_columns           = colnames(design),
  common_bcv               = sqrt(dge$common.dispersion)
)

capture.output(model_summary,
               file = file.path(dir_deg, "Model_Summary.txt"))

cat("----------------------------------------\n")
cat("PART 3 completed\n")
cat("----------------------------------------\n")
cat("Genes after filtering  :", nrow(dge), "\n")
cat("Samples                :", ncol(dge), "\n")
cat("Protected genotypes    :", n_P, "\n")
cat("Non-protected genotypes:", n_N, "\n")
cat("Common BCV             :", round(sqrt(dge$common.dispersion), 4), "\n")
cat("----------------------------------------\n\n")



# ============================================================
# PART 4 — QUALITY CONTROL
# ============================================================

if (!RUN_QC) {
  message("QC skipped")
} else {
  
  cat("---- PART 4: QC ----\n")
  
  
  ## 4.1 PCA — colored by genotype
  pca_res <- prcomp(t(logcpm_mat), scale. = FALSE)
  
  pca_df <- as.data.frame(pca_res$x[, 1:2])
  pca_df$Genotype   <- sample_info$Genotype
  pca_df$Protection <- sample_info$Protection
  pca_df$Replicate  <- sample_info$Replicate
  
  percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)
  
  rep_shapes <- setNames(
    c(16, 17, 15, 18, 8, 3)[seq_along(levels(sample_info$Replicate))],
    levels(sample_info$Replicate)
  )
  
  p_pca_genotype <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(PC1, PC2, colour = Genotype, shape = Replicate)
  ) +
    ggplot2::geom_point(size = 4, stroke = 0.5, alpha = 0.92) +
    ggplot2::scale_colour_manual(values = genotype_colors) +
    ggplot2::scale_shape_manual(values = rep_shapes) +
    ggplot2::xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
    ggplot2::ylab(paste0("PC2 (", percentVar[2], "% variance)")) +
    ggplot2::labs(
      title    = "PCA — logCPM",
      subtitle = "Colour = genotype | Shape = replicate",
      caption  = paste0("n = ", nrow(dge), " genes | TMM normalization")
    ) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "right")
  
  save_plot_pdf_svg(p_pca_genotype,
                    file.path(dir_qc, "PCA_Genotype"),
                    width = 8, height = 5.5)
  message("✅ PCA by genotype saved")
  
  
  ## 4.2 PCA — colored by protection
  p_pca_protection <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(PC1, PC2, colour = Protection, shape = Replicate)
  ) +
    ggplot2::geom_point(size = 4, stroke = 0.5, alpha = 0.92) +
    ggplot2::scale_colour_manual(values = PROTECTION_COLORS) +
    ggplot2::scale_shape_manual(values = rep_shapes) +
    ggplot2::xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
    ggplot2::ylab(paste0("PC2 (", percentVar[2], "% variance)")) +
    ggplot2::labs(
      title    = "PCA — logCPM",
      subtitle = "Colour = protection status | Shape = replicate",
      caption  = paste0("n = ", nrow(dge), " genes | TMM normalization")
    ) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "right")
  
  save_plot_pdf_svg(p_pca_protection,
                    file.path(dir_qc, "PCA_Protection"),
                    width = 8, height = 5.5)
  message("✅ PCA by protection saved")
  
  
  ## 4.3 MDS — colored by genotype
  save_base_plot_pdf_svg(
    function() {
      col_mds <- genotype_colors[as.character(sample_info$Genotype)]
      pch_mds <- c(16, 17, 15, 18, 8, 3)[as.integer(sample_info$Replicate)]
      limma::plotMDS(dge, col = col_mds, pch = pch_mds, cex = 1.6,
                     main = "MDS plot — leading logFC distances")
      legend("topright",
             legend = levels(sample_info$Genotype),
             col    = genotype_colors[levels(sample_info$Genotype)],
             pch    = 16, cex = 0.75, bty = "n", title = "Genotype")
    },
    file.path(dir_qc, "MDS_Genotype"),
    width = 8, height = 6
  )
  message("✅ MDS by genotype saved")
  
  
  ## 4.4 MDS — colored by protection
  save_base_plot_pdf_svg(
    function() {
      col_mds <- PROTECTION_COLORS[as.character(sample_info$Protection)]
      pch_mds <- c(16, 17, 15, 18, 8, 3)[as.integer(sample_info$Replicate)]
      limma::plotMDS(dge, col = col_mds, pch = pch_mds, cex = 1.6,
                     main = "MDS plot — leading logFC distances")
      legend("topright",
             legend = levels(sample_info$Protection),
             col    = PROTECTION_COLORS[levels(sample_info$Protection)],
             pch    = 16, cex = 0.85, bty = "n", title = "Protection")
    },
    file.path(dir_qc, "MDS_Protection"),
    width = 8, height = 6
  )
  message("✅ MDS by protection saved")
  
  
  ## 4.5 Sample distance heatmap
  sampleDists <- stats::dist(t(logcpm_mat))
  sdm <- as.matrix(sampleDists)
  hm_labels <- sample_info$Sample
  rownames(sdm) <- hm_labels
  colnames(sdm) <- hm_labels
  
  ann_qc <- data.frame(
    Protection = sample_info$Protection,
    Genotype   = sample_info$Genotype,
    row.names  = hm_labels
  )
  
  save_pheatmap_pdf_svg(
    list(
      mat                      = sdm,
      clustering_distance_rows = sampleDists,
      clustering_distance_cols = sampleDists,
      annotation_row           = ann_qc,
      annotation_col           = ann_qc,
      annotation_colors        = list(Protection = PROTECTION_COLORS,
                                      Genotype   = genotype_colors),
      border_color             = NA,
      fontsize                 = 8,
      treeheight_row           = 20,
      treeheight_col           = 20,
      main                     = "Sample-to-sample distances | logCPM TMM"
    ),
    file.path(dir_qc, "SampleDistanceHeatmap"),
    w = 9, h = 8
  )
  message("✅ Sample distance heatmap saved")
  
  
  ## 4.6 BCV plot
  save_base_plot_pdf_svg(
    function() {
      edgeR::plotBCV(
        dge,
        main = paste0("Biological coefficient of variation\nCommon BCV = ",
                      round(sqrt(dge$common.dispersion), 3))
      )
    },
    file.path(dir_qc, "BCV_plot"),
    width = 7, height = 5.5
  )
  message("✅ BCV plot saved")
  
  
  ## 4.7 Mean-variance relationship
  save_base_plot_pdf_svg(
    function() {
      edgeR::plotMeanVar(dge, show.raw.vars = TRUE, show.tagwise.vars = TRUE,
                         NBline = TRUE, main = "Mean-variance relationship")
    },
    file.path(dir_qc, "MeanVariance_plot"),
    width = 7, height = 5.5
  )
  message("✅ Mean-variance plot saved")
  
  
  cat("----------------------------------------\n")
  cat("PART 4 completed — QC plots saved\n")
  cat("----------------------------------------\n\n")
}



# ============================================================
# PART 5 — PRIMARY DIFFERENTIAL EXPRESSION: Protection_mean
# ============================================================

if (!RUN_PRIMARY_DEG) {
  message("Primary DEG analysis skipped")
} else {
  
  cat("---- PART 5: Protection_mean differential expression ----\n")
  
  
  ## 5.1 Build Protection_mean contrast
  p_design_names <- geno_map[conditions_P]
  n_design_names <- geno_map[conditions_N]
  
  contrast_mean_str <- paste0(
    "(", paste(p_design_names, collapse = " + "), ") / ", length(p_design_names),
    " - (", paste(n_design_names, collapse = " + "), ") / ", length(n_design_names)
  )
  
  cat("Protection_mean contrast:\n  ", contrast_mean_str, "\n\n")
  
  
  ## 5.2 Run glmQLFTest
  contrast_mean <- limma::makeContrasts(contrasts = contrast_mean_str, levels = design)
  
  qlf_mean <- edgeR::glmQLFTest(fit, contrast = contrast_mean)
  
  res_mean <- edgeR::topTags(qlf_mean, n = Inf, sort.by = "PValue")$table |>
    tibble::rownames_to_column("Gene") |>
    dplyr::rename(log2FoldChange = logFC,
                  pvalue         = PValue,
                  padj           = FDR,
                  aveLogCPM      = logCPM) |>
    dplyr::mutate(
      padj      = ifelse(is.na(padj), 1, padj),
      Contrast  = "Protection_mean",
      Status    = dplyr::case_when(
        padj < alpha & log2FoldChange >  lfc_threshold ~ "Up",
        padj < alpha & log2FoldChange < -lfc_threshold ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  message("✅ Protection_mean contrast computed")
  
  
  ## 5.3 Export DEG tables
  write.csv(res_mean,
            file.path(dir_deg, "DEG_Protection_mean.csv"),
            row.names = FALSE)
  
  res_mean_sig <- res_mean |>
    dplyr::filter(Status != "NS") |>
    dplyr::arrange(padj)
  
  write.csv(res_mean_sig,
            file.path(dir_deg, "DEG_Protection_mean_significant.csv"),
            row.names = FALSE)
  message("✅ DEG tables exported")
  
  
  ## 5.4 DEG summary
  n_up   <- sum(res_mean$Status == "Up")
  n_down <- sum(res_mean$Status == "Down")
  n_deg  <- n_up + n_down
  
  deg_summary <- data.frame(
    Contrast      = "Protection_mean",
    N_tested      = nrow(res_mean),
    N_up          = n_up,
    N_down        = n_down,
    N_total_DEG   = n_deg,
    log2FC_cutoff = lfc_threshold,
    FDR_cutoff    = alpha
  )
  
  write.csv(deg_summary,
            file.path(dir_deg, "DEG_Protection_mean_summary.csv"),
            row.names = FALSE)
  
  cat("DEG summary:\n")
  print(deg_summary)
  cat("\n")
  
  
  ## 5.5 Volcano plot
  lfc_lim  <- max(abs(res_mean$log2FoldChange), na.rm = TRUE) * 1.05
  pval_safe <- res_mean$padj[res_mean$padj > 0 & !is.na(res_mean$padj)]
  pval_lim  <- min(max(-log10(pval_safe), na.rm = TRUE) * 1.05, 50)
  
  top_label <- res_mean |>
    dplyr::filter(Status != "NS") |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = 15)
  
  p_volcano <- ggplot2::ggplot(
    res_mean,
    ggplot2::aes(x = log2FoldChange,
                 y = -log10(pmax(padj, 10^(-pval_lim))),
                 colour = Status)
  ) +
    ggplot2::geom_point(data = dplyr::filter(res_mean, Status == "NS"),
                        colour = DEG_COLORS["NS"], alpha = 0.35, size = 0.8) +
    ggplot2::geom_point(data = dplyr::filter(res_mean, Status != "NS"),
                        alpha = 0.80, size = 1.4) +
    ggplot2::scale_colour_manual(
      values = DEG_COLORS,
      labels = c(Up   = paste0("Up (n = ", n_up, ")"),
                 Down = paste0("Down (n = ", n_down, ")"),
                 NS   = "Not significant")
    ) +
    ggrepel::geom_text_repel(
      data = top_label,
      ggplot2::aes(label = Gene),
      size = 2.6, max.overlaps = 18, segment.linewidth = 0.3,
      segment.colour = "grey55", colour = "black", fontface = "italic",
      box.padding = 0.35, min.segment.length = 0.2
    ) +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::scale_x_continuous(limits = c(-lfc_lim, lfc_lim)) +
    ggplot2::scale_y_continuous(limits = c(0, pval_lim),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title    = "Volcano plot — Protection_mean",
      subtitle = "Protected vs NonProtected",
      x        = expression(log[2]~"fold change"),
      y        = expression(-log[10]~"(FDR)"),
      colour   = NULL
    ) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  
  save_plot_pdf_svg(p_volcano,
                    file.path(dir_deg, "Volcano_Protection_mean"),
                    width = 6.5, height = 5.5)
  message("✅ Volcano plot saved")
  
  
  ## 5.6 MA plot
  p_ma <- ggplot2::ggplot(
    res_mean,
    ggplot2::aes(x = aveLogCPM, y = log2FoldChange)
  ) +
    ggplot2::geom_point(data = dplyr::filter(res_mean, Status == "NS"),
                        colour = DEG_COLORS["NS"], alpha = 0.30, size = 0.6) +
    ggplot2::geom_point(data = dplyr::filter(res_mean, Status == "Up"),
                        colour = DEG_COLORS["Up"],   alpha = 0.70, size = 0.9) +
    ggplot2::geom_point(data = dplyr::filter(res_mean, Status == "Down"),
                        colour = DEG_COLORS["Down"], alpha = 0.70, size = 0.9) +
    ggplot2::geom_hline(yintercept = 0,
                        colour = "black", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_smooth(colour = "#C0392B", se = TRUE, method = "loess",
                         method.args = list(span = 0.4, family = "symmetric"),
                         linewidth = 0.8, fill = "#F1948A", alpha = 0.20,
                         na.rm = TRUE) +
    ggplot2::labs(
      title    = "MA plot — Protection_mean",
      subtitle = "Protected vs NonProtected",
      x        = "Average logCPM",
      y        = expression(log[2]~"fold change")
    ) +
    theme_pub(base_size = 11)
  
  save_plot_pdf_svg(p_ma,
                    file.path(dir_deg, "MA_Protection_mean"),
                    width = 6.5, height = 5)
  message("✅ MA plot saved")
  
  ## 5.7 Save main results object
  
  ## 5.7 Heatmap publication-style des gènes DE
  dir_deg_hm <- file.path(dir_deg, "Heatmaps")
  dir.create(dir_deg_hm, recursive = TRUE, showWarnings = FALSE)
  
  ## Sélection des gènes significatifs
  heatmap_genes <- res_mean |>
    dplyr::filter(Status != "NS") |>
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) |>
    dplyr::pull(Gene)
  
  ## Option de sécurité : si trop de gènes, on garde les plus informatifs
  MAX_HEATMAP_GENES <- 300
  
  if (length(heatmap_genes) > MAX_HEATMAP_GENES) {
    heatmap_genes <- res_mean |>
      dplyr::filter(Status != "NS") |>
      dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) |>
      dplyr::slice_head(n = MAX_HEATMAP_GENES) |>
      dplyr::pull(Gene)
  }
  
  if (length(heatmap_genes) >= 2) {
    
    hm_mat <- logcpm_mat[heatmap_genes, , drop = FALSE]
    
    ## Z-score par gène
    hm_mat_z <- t(scale(t(hm_mat)))
    hm_mat_z[is.na(hm_mat_z)] <- 0
    
    ann_col <- data.frame(
      Protection = sample_info$Protection,
      Genotype   = sample_info$Genotype,
      Replicate  = sample_info$Replicate,
      row.names  = rownames(sample_info)
    )
    
    rep_levels <- levels(sample_info$Replicate)
    rep_cols <- setNames(
      grDevices::hcl.colors(length(rep_levels), palette = "Dark 3"),
      rep_levels
    )
    
    ann_colors <- list(
      Protection = PROTECTION_COLORS,
      Genotype   = genotype_colors,
      Replicate  = rep_cols
    )
    
    save_pheatmap_pdf_svg(
      list(
        mat               = hm_mat_z,
        annotation_col    = ann_col,
        annotation_colors = ann_colors,
        cluster_rows      = TRUE,
        cluster_cols      = TRUE,
        show_rownames     = FALSE,
        show_colnames     = FALSE,
        border_color      = NA,
        fontsize          = 8,
        treeheight_row    = 30,
        treeheight_col    = 35,
        color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        breaks            = seq(-3, 3, length.out = 101),
        main = paste0(
          "Heatmap — Protection_mean DE genes",
          "\nProtected vs NonProtected | n = ", nrow(hm_mat_z), " genes"
        )
      ),
      file.path(dir_deg_hm, "Heatmap_Protection_mean_DEG"),
      w = 9,
      h = max(7, min(18, 4 + nrow(hm_mat_z) * 0.02))
    )
    
    message("✅ Publication-style DEG heatmap saved")
  } else {
    message("⚠️ Not enough DE genes for heatmap")
  }
  
  ## 5.7 Save main results object
  primary_deg_results <- res_mean
  
  cat("----------------------------------------\n")
  cat("PART 5 completed\n")
  cat("----------------------------------------\n")
  cat("Contrast  :", "Protection_mean", "\n")
  cat("Up DEGs   :", n_up, "\n")
  cat("Down DEGs :", n_down, "\n")
  cat("Total DEGs:", n_deg, "\n")
  cat("----------------------------------------\n\n")
}





# ============================================================
# ============================================================
# PART 6 — GO ENRICHMENT: Protection_mean  (topGO)
# ============================================================
#
#  6.1  Load GO annotation + build universe
#  6.2  Define DEG gene sets (Up / Down)
#  6.3  Run topGO  (weight01 Fisher, 3 ontologies × 2 directions)
#  6.4  Export raw tables
#  6.5  Visualizations per direction:
#         Dotplot  (GeneRatio × -log10p, size = DEG_Count)
#         Barplot  (-log10p bars, fill = FDR)
#  6.6  Combined bubble plot (Ontology × Direction)
#  6.7  Summary table
#
# ============================================================

if (!RUN_GO_PRIMARY) {
  message("GO analysis skipped")
} else {
  
  cat("---- PART 6: GO enrichment (topGO) — Protection_mean ----\n")
  
  
  ## 6.1 Load GO annotation + universe ----------------------
  if (!file.exists(go_annot_file))
    stop("GO annotation file not found: ", go_annot_file)
  
  geneID2GO     <- load_go_annotation(go_annot_file)
  genes_tested  <- rownames(dge)
  gene_universe <- intersect(genes_tested, names(geneID2GO))
  
  cat("Genes tested  :", length(genes_tested), "\n")
  cat("GO universe   :", length(gene_universe), "\n")
  
  if (length(gene_universe) == 0)
    stop("No overlap between tested genes and GO annotation.")
  
  message("\u2705 GO annotation loaded")
  
  
  ## 6.2 Define DEG gene sets --------------------------------
  genes_up   <- res_mean$Gene[res_mean$Status == "Up"]
  genes_down <- res_mean$Gene[res_mean$Status == "Down"]
  
  cat("Up DEGs  :", length(genes_up),   "\n")
  cat("Down DEGs:", length(genes_down), "\n\n")
  
  
  ## 6.3 Run topGO -------------------------------------------
  go_results <- list()
  
  for (direction in c("Up", "Down")) {
    genes_oi <- if (direction == "Up") genes_up else genes_down
    cat("Running topGO —", direction, "\n")
    
    res_all_ont <- run_topgo_all_ont(
      genes_oi      = genes_oi,
      gene_universe = gene_universe,
      geneID2GO     = geneID2GO
    )
    
    if (!is.null(res_all_ont) && nrow(res_all_ont) > 0) {
      res_all_ont$Direction <- direction
      res_all_ont$GeneSet   <- "Protection_mean"
      go_results[[direction]] <- res_all_ont
    }
  }
  
  
  ## 6.4 Export raw tables -----------------------------------
  for (direction in names(go_results)) {
    df <- go_results[[direction]]
    for (ont in unique(df$Ontology)) {
      write.csv(
        df[df$Ontology == ont, ],
        file.path(dir_go,
                  paste0("GO_Protection_mean_", direction, "_", ont, ".csv")),
        row.names = FALSE
      )
    }
  }
  
  if (length(go_results) > 0) {
    go_combined <- dplyr::bind_rows(go_results) |>
      dplyr::arrange(Direction, Ontology, weight01)
    
    write.csv(go_combined,
              file.path(dir_go, "GO_Protection_mean_combined.csv"),
              row.names = FALSE)
    message("\u2705 Combined GO table exported")
  } else {
    warning("No significant GO terms found for Protection_mean.", call. = FALSE)
    cat("----------------------------------------\n")
    cat("PART 6 completed (no significant terms)\n")
    cat("----------------------------------------\n\n")
  }
  
  
  ## 6.5 Publication figure: grouped view --------------------------------
  ## 6.5 Publication figure: same style as Tier1_both --------------------
  if (length(go_results) > 0) {
    
    ## Fusionner Up + Down si présents
    go_main_simple <- go_combined |>
      dplyr::mutate(
        GeneSet = "Protection_mean"
      )
    
    p_pub <- plot_go_dot_simple(
      go_main_simple,
      title = "GO enrichment — Protection_mean",
      top_n = TOP_N_TERMS
    )
    
    save_plot_pdf_svg(
      p_pub,
      file.path(dir_go, "GO_Protection_mean_Publication"),
      width  = 8,
      height = 7.5
    )
    message("✅ Publication GO figure saved: GO_Protection_mean_Publication")
    
    
    ## 6.6 Per-direction figures
    for (direction in names(go_results)) {
      df_dir <- go_results[[direction]]
      
      p_dir <- plot_go_dot_simple(
        df_dir,
        title = paste0("GO enrichment — Protection_mean_", direction),
        top_n = TOP_N_TERMS
      )
      
      save_plot_pdf_svg(
        p_dir,
        file.path(dir_go, paste0("GO_Protection_mean_", direction, "_Supp")),
        width  = 8,
        height = 7
      )
    }
    message("✅ Supplementary per-direction GO figures saved")
    
    
    ## 6.7 Summary table ----------------------------------------
    go_summary <- go_combined |>
      dplyr::group_by(Direction, Ontology) |>
      dplyr::summarise(N_significant = dplyr::n(),
                       Top_term      = Term[which.min(weight01)],
                       .groups = "drop") |>
      dplyr::arrange(Direction, Ontology)
    
    write.csv(go_summary,
              file.path(dir_go, "GO_Protection_mean_summary.csv"),
              row.names = FALSE)
    
    cat("\nGO enrichment summary:\n")
    print(go_summary)
  }
  
  
  cat("----------------------------------------\n")
  cat("PART 6 completed\n")
  cat("----------------------------------------\n")
  cat("GO results in:", dir_go, "\n")
  cat("----------------------------------------\n\n")
}
# PART 7 — LOO VALIDATION: SYMMETRIC + COMPARATIVE ANALYSIS
# ============================================================
#
#  BLOCK A — LOO-P
#  BLOCK B — LOO-N
#  BLOCK C — COMPARATIVE ANALYSIS
#
# ============================================================

if (!RUN_LOO) {
  message("LOO analysis skipped")
} else {
  
  cat("========================================================\n")
  cat("PART 7 — LOO symmetric validation\n")
  cat("========================================================\n\n")
  
  
  ## 7.0 Helper functions ------------------------------------
  
  run_loo_contrast <- function(contrast_str, contrast_name, removed_genotype, side) {
    cvec <- tryCatch(
      limma::makeContrasts(contrasts = contrast_str, levels = design),
      error = function(e) {
        warning("makeContrasts failed: ", contrast_name, " — ", e$message,
                call. = FALSE)
        NULL
      }
    )
    if (is.null(cvec)) return(NULL)
    
    qlf <- edgeR::glmQLFTest(fit, contrast = cvec)
    
    res <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table |>
      tibble::rownames_to_column("Gene") |>
      dplyr::rename(log2FoldChange = logFC,
                    pvalue         = PValue,
                    padj           = FDR,
                    aveLogCPM      = logCPM) |>
      dplyr::mutate(
        padj     = ifelse(is.na(padj), 1, padj),
        Contrast = contrast_name,
        Removed  = removed_genotype,
        Side     = side,
        Status   = dplyr::case_when(
          padj < alpha & log2FoldChange >  lfc_threshold ~ "Up",
          padj < alpha & log2FoldChange < -lfc_threshold ~ "Down",
          TRUE ~ "NS"
        )
      )
    res
  }
  
  primary_sig <- res_mean |>
    dplyr::filter(Status != "NS") |>
    dplyr::select(Gene, log2FoldChange, padj, aveLogCPM, Status)
  
  score_genes <- function(gene_vec, direction, loo_list) {
    if (length(gene_vec) == 0 || length(loo_list) == 0) {
      return(data.frame(Gene = gene_vec, LOO_Score = integer(length(gene_vec))))
    }
    scores <- sapply(gene_vec, function(g) {
      sum(sapply(loo_list, function(df) {
        idx <- match(g, df$Gene)
        if (is.na(idx)) return(0L)
        as.integer(df$Status[idx] == direction)
      }))
    })
    data.frame(Gene = gene_vec, LOO_Score = as.integer(scores))
  }
  
  message("✅ LOO helper functions ready")
  
  
  # ============================================================
  # BLOCK A — LOO-P
  # ============================================================
  
  cat("\n-------- BLOCK A: LOO-P --------\n")
  
  ## 7A.1 Build LOO-P contrasts
  loo_contrasts_P <- lapply(conditions_P, function(gP) {
    remaining <- p_design_names[p_design_names != geno_map[gP]]
    if (length(remaining) == 0) {
      warning("LOO-P: only one P genotype — cannot remove ", gP, call. = FALSE)
      return(NULL)
    }
    ct_str <- paste0(
      "(", paste(remaining, collapse = " + "), ") / ", length(remaining),
      " - (", paste(n_design_names, collapse = " + "), ") / ", length(n_design_names)
    )
    list(name    = paste0("LOOP_without_", gP),
         removed = gP, side = "P", str = ct_str)
  })
  loo_contrasts_P <- Filter(Negate(is.null), loo_contrasts_P)
  n_loo_P <- length(loo_contrasts_P)
  
  cat("LOO-P contrasts:", n_loo_P, "\n")
  cat("  Removed one by one:",
      paste(sapply(loo_contrasts_P, `[[`, "removed"), collapse = ", "), "\n\n")
  
  
  ## 7A.2 Run LOO-P contrasts
  loo_results_P <- list()
  loo_summary_P <- list()
  
  for (item in loo_contrasts_P) {
    cat("  Running:", item$name, "\n")
    res <- run_loo_contrast(item$str, item$name, item$removed, "P")
    if (!is.null(res)) {
      loo_results_P[[item$name]] <- res
      loo_summary_P[[item$name]] <- data.frame(
        Contrast = item$name, Removed = item$removed, Side = "P",
        N_up = sum(res$Status == "Up"), N_down = sum(res$Status == "Down"),
        N_total = sum(res$Status != "NS")
      )
      write.csv(res, file.path(dir_deg, paste0("DEG_", item$name, ".csv")),
                row.names = FALSE)
    }
  }
  
  loo_summary_P_df <- dplyr::bind_rows(loo_summary_P)
  write.csv(loo_summary_P_df, file.path(dir_deg, "LOO_P_summary.csv"), row.names = FALSE)
  message("✅ LOO-P contrasts done")
  
  
  ## 7A.3 LOO-P robustness scores
  score_P_up   <- score_genes(primary_sig$Gene[primary_sig$Status == "Up"],
                              "Up",   loo_results_P)
  score_P_down <- score_genes(primary_sig$Gene[primary_sig$Status == "Down"],
                              "Down", loo_results_P)
  
  loo_score_P <- primary_sig |>
    dplyr::left_join(
      dplyr::bind_rows(score_P_up, score_P_down) |>
        dplyr::rename(LOO_Score_P = LOO_Score),
      by = "Gene"
    ) |>
    dplyr::mutate(
      LOO_Score_P = dplyr::coalesce(LOO_Score_P, 0L),
      LOO_Max_P   = n_loo_P,
      LOO_Norm_P  = ifelse(n_loo_P > 0, LOO_Score_P / n_loo_P, NA_real_)
    )
  
  write.csv(loo_score_P, file.path(dir_deg, "LOO_P_gene_scores.csv"), row.names = FALSE)
  message("✅ LOO-P scores computed")
  
  
  ## 7A.4 Tier classification LOO-P
  tier_P_df <- loo_score_P |>
    dplyr::mutate(
      Tier_P = dplyr::case_when(
        LOO_Score_P == LOO_Max_P              ~ "Tier1",
        LOO_Score_P >= ceiling(LOO_Max_P / 2) ~ "Tier2",
        TRUE                                  ~ "Tier3"
      ),
      Tier_P = factor(Tier_P, levels = c("Tier1", "Tier2", "Tier3"))
    ) |>
    dplyr::arrange(Tier_P, padj)
  
  write.csv(tier_P_df, file.path(dir_deg, "DEG_Tiers_LOOP.csv"), row.names = FALSE)
  
  tier_P_summary <- tier_P_df |>
    dplyr::count(Tier_P, Status, name = "N_genes") |>
    dplyr::arrange(Tier_P, Status)
  
  write.csv(tier_P_summary, file.path(dir_deg, "DEG_TierSummary_LOOP.csv"), row.names = FALSE)
  
  cat("\nTier classification — LOO-P:\n")
  print(tier_P_summary)
  cat("\n")
  message("✅ LOO-P tier classification done")
  
  
  ## 7A.5 Plots LOO-P
  p_bar_P <- ggplot2::ggplot(
    loo_score_P |>
      dplyr::count(Status, LOO_Score_P) |>
      dplyr::mutate(LOO_Score_P = factor(LOO_Score_P, levels = 0:n_loo_P)),
    ggplot2::aes(x = LOO_Score_P, y = n, fill = Status)
  ) +
    ggplot2::geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = DEG_COLORS[c("Up", "Down")]) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(title = "LOO-P robustness score",
                  subtitle = paste0("Each protected genotype removed once | max score = ", n_loo_P),
                  x = "LOO-P score", y = "Number of genes", fill = NULL) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  
  save_plot_pdf_svg(p_bar_P, file.path(dir_deg, "LOO_P_Score_Barplot"),
                    width = 7, height = 5)
  
  p_tier_P <- ggplot2::ggplot(tier_P_df, ggplot2::aes(x = Tier_P, fill = Tier_P)) +
    ggplot2::geom_bar(colour = "white", linewidth = 0.4) +
    ggplot2::facet_wrap(~ Status, scales = "free_y") +
    ggplot2::geom_text(stat = "count",
                       ggplot2::aes(label = ggplot2::after_stat(count)),
                       vjust = -0.4, fontface = "bold", size = 4) +
    ggplot2::scale_fill_manual(
      values = c("Tier1" = "#1A237E", "Tier2" = "#F57F17", "Tier3" = "#78909C")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(title    = "Tier classification — LOO-P only",
                  subtitle = paste0("Tier1 = robust across all ", n_loo_P, " LOO-P contrasts"),
                  x = NULL, y = "Number of genes", fill = NULL) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "none")
  
  save_plot_pdf_svg(p_tier_P, file.path(dir_deg, "Tier_Barplot_LOOP"),
                    width = 8, height = 5.5)
  
  tv_P <- res_mean |>
    dplyr::left_join(tier_P_df |> dplyr::select(Gene, Tier_P), by = "Gene") |>
    dplyr::mutate(Tier_P = as.character(dplyr::coalesce(as.character(Tier_P), "NS")))
  
  pv_lim  <- min(max(-log10(tv_P$padj[tv_P$padj > 0 & !is.na(tv_P$padj)]),
                     na.rm = TRUE) * 1.05, 50)
  lfc_lim <- max(abs(tv_P$log2FoldChange), na.rm = TRUE) * 1.05
  
  p_volc_P <- ggplot2::ggplot(
    tv_P,
    ggplot2::aes(x = log2FoldChange,
                 y = -log10(pmax(padj, 10^(-pv_lim))),
                 colour = Tier_P)
  ) +
    ggplot2::geom_point(data = dplyr::filter(tv_P, Tier_P == "NS"),
                        colour = "#E0E0E0", alpha = 0.35, size = 0.8) +
    ggplot2::geom_point(data = dplyr::filter(tv_P, Tier_P != "NS"),
                        alpha = 0.85, size = 1.5) +
    ggplot2::scale_colour_manual(
      values = c("Tier1" = "#1A237E", "Tier2" = "#F57F17",
                 "Tier3" = "#78909C", "NS" = "#E0E0E0")) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(tv_P, Tier_P == "Tier1") |>
        dplyr::arrange(padj) |> dplyr::slice_head(n = 15),
      ggplot2::aes(label = Gene),
      size = 2.6, max.overlaps = 20, segment.linewidth = 0.3,
      segment.colour = "grey55", colour = "black", fontface = "italic",
      box.padding = 0.35
    ) +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::scale_x_continuous(limits = c(-lfc_lim, lfc_lim)) +
    ggplot2::scale_y_continuous(limits = c(0, pv_lim),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(title    = "Volcano — LOO-P tiers",
                  subtitle = paste0("Tier1 labels | LOO-P max score = ", n_loo_P),
                  x = expression(log[2]~"fold change"),
                  y = expression(-log[10]~"(FDR)"),
                  colour = NULL) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  
  save_plot_pdf_svg(p_volc_P, file.path(dir_deg, "Volcano_Tiers_LOOP"),
                    width = 7, height = 6)
  message("✅ LOO-P plots saved")
  
  
  # ============================================================
  # BLOCK B — LOO-N
  # ============================================================
  
  cat("\n-------- BLOCK B: LOO-N --------\n")
  
  ## 7B.1 Build LOO-N contrasts
  loo_contrasts_N <- lapply(conditions_N, function(gN) {
    remaining <- n_design_names[n_design_names != geno_map[gN]]
    if (length(remaining) == 0) {
      warning("LOO-N: only one N genotype — cannot remove ", gN, call. = FALSE)
      return(NULL)
    }
    ct_str <- paste0(
      "(", paste(p_design_names, collapse = " + "), ") / ", length(p_design_names),
      " - (", paste(remaining, collapse = " + "), ") / ", length(remaining)
    )
    list(name    = paste0("LOON_without_", gN),
         removed = gN, side = "N", str = ct_str)
  })
  loo_contrasts_N <- Filter(Negate(is.null), loo_contrasts_N)
  n_loo_N <- length(loo_contrasts_N)
  
  cat("LOO-N contrasts:", n_loo_N, "\n")
  cat("  Removed one by one:",
      paste(sapply(loo_contrasts_N, `[[`, "removed"), collapse = ", "), "\n\n")
  
  
  ## 7B.2 Run LOO-N contrasts
  loo_results_N <- list()
  loo_summary_N <- list()
  
  for (item in loo_contrasts_N) {
    cat("  Running:", item$name, "\n")
    res <- run_loo_contrast(item$str, item$name, item$removed, "N")
    if (!is.null(res)) {
      loo_results_N[[item$name]] <- res
      loo_summary_N[[item$name]] <- data.frame(
        Contrast = item$name, Removed = item$removed, Side = "N",
        N_up = sum(res$Status == "Up"), N_down = sum(res$Status == "Down"),
        N_total = sum(res$Status != "NS")
      )
      write.csv(res, file.path(dir_deg, paste0("DEG_", item$name, ".csv")),
                row.names = FALSE)
    }
  }
  
  loo_summary_N_df <- dplyr::bind_rows(loo_summary_N)
  write.csv(loo_summary_N_df, file.path(dir_deg, "LOO_N_summary.csv"), row.names = FALSE)
  message("✅ LOO-N contrasts done")
  
  
  ## 7B.3 LOO-N robustness scores
  score_N_up   <- score_genes(primary_sig$Gene[primary_sig$Status == "Up"],
                              "Up",   loo_results_N)
  score_N_down <- score_genes(primary_sig$Gene[primary_sig$Status == "Down"],
                              "Down", loo_results_N)
  
  loo_score_N <- primary_sig |>
    dplyr::left_join(
      dplyr::bind_rows(score_N_up, score_N_down) |>
        dplyr::rename(LOO_Score_N = LOO_Score),
      by = "Gene"
    ) |>
    dplyr::mutate(
      LOO_Score_N = dplyr::coalesce(LOO_Score_N, 0L),
      LOO_Max_N   = n_loo_N,
      LOO_Norm_N  = ifelse(n_loo_N > 0, LOO_Score_N / n_loo_N, NA_real_)
    )
  
  write.csv(loo_score_N, file.path(dir_deg, "LOO_N_gene_scores.csv"), row.names = FALSE)
  message("✅ LOO-N scores computed")
  
  
  ## 7B.4 Tier classification LOO-N
  tier_N_df <- loo_score_N |>
    dplyr::mutate(
      Tier_N = dplyr::case_when(
        LOO_Score_N == LOO_Max_N              ~ "Tier1",
        LOO_Score_N >= ceiling(LOO_Max_N / 2) ~ "Tier2",
        TRUE                                  ~ "Tier3"
      ),
      Tier_N = factor(Tier_N, levels = c("Tier1", "Tier2", "Tier3"))
    ) |>
    dplyr::arrange(Tier_N, padj)
  
  write.csv(tier_N_df, file.path(dir_deg, "DEG_Tiers_LOON.csv"), row.names = FALSE)
  
  tier_N_summary <- tier_N_df |>
    dplyr::count(Tier_N, Status, name = "N_genes") |>
    dplyr::arrange(Tier_N, Status)
  
  write.csv(tier_N_summary, file.path(dir_deg, "DEG_TierSummary_LOON.csv"), row.names = FALSE)
  
  cat("\nTier classification — LOO-N:\n")
  print(tier_N_summary)
  cat("\n")
  message("✅ LOO-N tier classification done")
  
  
  ## 7B.5 Plots LOO-N
  p_bar_N <- ggplot2::ggplot(
    loo_score_N |>
      dplyr::count(Status, LOO_Score_N) |>
      dplyr::mutate(LOO_Score_N = factor(LOO_Score_N, levels = 0:n_loo_N)),
    ggplot2::aes(x = LOO_Score_N, y = n, fill = Status)
  ) +
    ggplot2::geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = DEG_COLORS[c("Up", "Down")]) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(title = "LOO-N robustness score",
                  subtitle = paste0("Each non-protected genotype removed once | max score = ", n_loo_N),
                  x = "LOO-N score", y = "Number of genes", fill = NULL) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  
  save_plot_pdf_svg(p_bar_N, file.path(dir_deg, "LOO_N_Score_Barplot"),
                    width = 7, height = 5)
  
  p_tier_N <- ggplot2::ggplot(tier_N_df, ggplot2::aes(x = Tier_N, fill = Tier_N)) +
    ggplot2::geom_bar(colour = "white", linewidth = 0.4) +
    ggplot2::facet_wrap(~ Status, scales = "free_y") +
    ggplot2::geom_text(stat = "count",
                       ggplot2::aes(label = ggplot2::after_stat(count)),
                       vjust = -0.4, fontface = "bold", size = 4) +
    ggplot2::scale_fill_manual(
      values = c("Tier1" = "#1A237E", "Tier2" = "#F57F17", "Tier3" = "#78909C")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(title    = "Tier classification — LOO-N only",
                  subtitle = paste0("Tier1 = robust across all ", n_loo_N, " LOO-N contrasts"),
                  x = NULL, y = "Number of genes", fill = NULL) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "none")
  
  save_plot_pdf_svg(p_tier_N, file.path(dir_deg, "Tier_Barplot_LOON"),
                    width = 8, height = 5.5)
  
  tv_N <- res_mean |>
    dplyr::left_join(tier_N_df |> dplyr::select(Gene, Tier_N), by = "Gene") |>
    dplyr::mutate(Tier_N = as.character(dplyr::coalesce(as.character(Tier_N), "NS")))
  
  pv_lim_N  <- min(max(-log10(tv_N$padj[tv_N$padj > 0 & !is.na(tv_N$padj)]),
                       na.rm = TRUE) * 1.05, 50)
  lfc_lim_N <- max(abs(tv_N$log2FoldChange), na.rm = TRUE) * 1.05
  
  p_volc_N <- ggplot2::ggplot(
    tv_N,
    ggplot2::aes(x = log2FoldChange,
                 y = -log10(pmax(padj, 10^(-pv_lim_N))),
                 colour = Tier_N)
  ) +
    ggplot2::geom_point(data = dplyr::filter(tv_N, Tier_N == "NS"),
                        colour = "#E0E0E0", alpha = 0.35, size = 0.8) +
    ggplot2::geom_point(data = dplyr::filter(tv_N, Tier_N != "NS"),
                        alpha = 0.85, size = 1.5) +
    ggplot2::scale_colour_manual(
      values = c("Tier1" = "#1A237E", "Tier2" = "#F57F17",
                 "Tier3" = "#78909C", "NS" = "#E0E0E0")) +
    ggrepel::geom_text_repel(
      data = dplyr::filter(tv_N, Tier_N == "Tier1") |>
        dplyr::arrange(padj) |> dplyr::slice_head(n = 15),
      ggplot2::aes(label = Gene),
      size = 2.6, max.overlaps = 20, segment.linewidth = 0.3,
      segment.colour = "grey55", colour = "black", fontface = "italic",
      box.padding = 0.35
    ) +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::scale_x_continuous(limits = c(-lfc_lim_N, lfc_lim_N)) +
    ggplot2::scale_y_continuous(limits = c(0, pv_lim_N),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(title    = "Volcano — LOO-N tiers",
                  subtitle = paste0("Tier1 labels | LOO-N max score = ", n_loo_N),
                  x = expression(log[2]~"fold change"),
                  y = expression(-log[10]~"(FDR)"),
                  colour = NULL) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  
  save_plot_pdf_svg(p_volc_N, file.path(dir_deg, "Volcano_Tiers_LOON"),
                    width = 7, height = 6)
  message("✅ LOO-N plots saved")
  
  
  # ============================================================
  # BLOCK C — COMPARATIVE ANALYSIS
  # ============================================================
  
  cat("\n-------- BLOCK C: LOO-P vs LOO-N comparative --------\n")
  
  ## 7C.1 Merge tier tables
  tier_compare_df <- primary_sig |>
    dplyr::left_join(
      tier_P_df |> dplyr::select(Gene, Tier_P, LOO_Score_P, LOO_Max_P, LOO_Norm_P),
      by = "Gene"
    ) |>
    dplyr::left_join(
      tier_N_df |> dplyr::select(Gene, Tier_N, LOO_Score_N, LOO_Max_N, LOO_Norm_N),
      by = "Gene"
    ) |>
    dplyr::mutate(
      Tier_P = dplyr::coalesce(as.character(Tier_P), "Tier3"),
      Tier_N = dplyr::coalesce(as.character(Tier_N), "Tier3"),
      Tier_Global = dplyr::case_when(
        Tier_P == "Tier1" & Tier_N == "Tier1" ~ "Tier1_both",
        Tier_P == "Tier1" & Tier_N != "Tier1" ~ "Tier1_P_only",
        Tier_P != "Tier1" & Tier_N == "Tier1" ~ "Tier1_N_only",
        Tier_P == "Tier2" | Tier_N == "Tier2" ~ "Tier2",
        TRUE                                  ~ "Tier3"
      ),
      Tier_Global = factor(Tier_Global,
                           levels = c("Tier1_both", "Tier1_P_only",
                                      "Tier1_N_only", "Tier2", "Tier3")),
      Tier_Changed    = Tier_P != Tier_N,
      Tier_Change_Dir = dplyr::case_when(
        !Tier_Changed   ~ "Stable",
        Tier_P < Tier_N ~ "Better_in_P",
        Tier_P > Tier_N ~ "Better_in_N",
        TRUE            ~ "Stable"
      )
    ) |>
    dplyr::arrange(Tier_Global, padj)
  
  write.csv(tier_compare_df,
            file.path(dir_deg, "DEG_Tiers_Comparative.csv"),
            row.names = FALSE)
  message("✅ Comparative tier table built")
  
  
  ## 7C.2 Concordance table
  concordance_tbl <- table(
    "Tier_LOO_P" = tier_compare_df$Tier_P,
    "Tier_LOO_N" = tier_compare_df$Tier_N
  )
  write.csv(as.data.frame.matrix(concordance_tbl),
            file.path(dir_deg, "LOO_P_vs_N_Concordance_Table.csv"))
  
  cat("\nConcordance table — Tier-P vs Tier-N:\n")
  print(concordance_tbl)
  cat("\n")
  
  conc_df <- as.data.frame(concordance_tbl) |>
    dplyr::rename(Tier_P = Tier_LOO_P, Tier_N = Tier_LOO_N, Count = Freq) |>
    dplyr::mutate(
      Tier_P = factor(Tier_P, levels = c("Tier1", "Tier2", "Tier3")),
      Tier_N = factor(Tier_N, levels = c("Tier1", "Tier2", "Tier3"))
    )
  
  p_concordance <- ggplot2::ggplot(
    conc_df,
    ggplot2::aes(x = Tier_N, y = Tier_P, fill = Count)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = Count),
                       fontface = "bold", size = 5, colour = "white") +
    ggplot2::scale_fill_gradient(low = "#BDC3C7", high = "#1A237E", name = "N genes") +
    ggplot2::labs(title    = "Concordance: LOO-P tier vs LOO-N tier",
                  subtitle = paste0("Diagonal = genes with same tier on both sides\n",
                                    "Off-diagonal = genes whose robustness differs by side"),
                  x = "Tier from LOO-N", y = "Tier from LOO-P") +
    theme_pub(base_size = 12) +
    ggplot2::theme(legend.position = "right")
  
  save_plot_pdf_svg(p_concordance, file.path(dir_deg, "LOO_Concordance_Heatmap"),
                    width = 6, height = 5)
  message("✅ Concordance heatmap saved")
  
  
  ## 7C.3 Gene sub-lists
  tier1_both   <- tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_both")
  tier1_P_only <- tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_P_only")
  tier1_N_only <- tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_N_only")
  tier_changed <- tier_compare_df |> dplyr::filter(Tier_Changed)
  
  cat("Genes Tier1 on BOTH sides     :", nrow(tier1_both),   "\n")
  cat("Genes Tier1 on LOO-P side only:", nrow(tier1_P_only), "\n")
  cat("Genes Tier1 on LOO-N side only:", nrow(tier1_N_only), "\n")
  cat("Genes with any tier change     :", nrow(tier_changed), "\n\n")
  
  write.csv(tier1_both,   file.path(dir_deg, "DEG_Tier1_Both_Sides.csv"),    row.names = FALSE)
  write.csv(tier1_P_only, file.path(dir_deg, "DEG_Tier1_P_only.csv"),        row.names = FALSE)
  write.csv(tier1_N_only, file.path(dir_deg, "DEG_Tier1_N_only.csv"),        row.names = FALSE)
  write.csv(tier_changed, file.path(dir_deg, "DEG_Tier_Changed_P_vs_N.csv"), row.names = FALSE)
  message("✅ Discordance gene lists exported")
  
  
  ## 7C.4 Scatter LOO-P vs LOO-N score
  p_scatter <- ggplot2::ggplot(
    tier_compare_df,
    ggplot2::aes(x = LOO_Score_P, y = LOO_Score_N,
                 colour = Tier_Global, shape = Status)
  ) +
    ggplot2::geom_jitter(width = 0.15, height = 0.12, size = 2.2, alpha = 0.82) +
    ggplot2::scale_colour_manual(
      values = c("Tier1_both"   = "#1A237E", "Tier1_P_only" = "#E53935",
                 "Tier1_N_only" = "#00897B", "Tier2"        = "#F57F17",
                 "Tier3"        = "#B0BEC5"),
      labels = c("Tier1_both"   = "Tier1 both sides",
                 "Tier1_P_only" = "Tier1 LOO-P only",
                 "Tier1_N_only" = "Tier1 LOO-N only",
                 "Tier2"        = "Tier2",
                 "Tier3"        = "Tier3")
    ) +
    ggplot2::scale_shape_manual(values = c("Up" = 16, "Down" = 17)) +
    ggplot2::scale_x_continuous(breaks = 0:n_loo_P,
                                limits = c(-0.5, n_loo_P + 0.5),
                                name   = paste0("LOO-P score (max = ", n_loo_P, ")")) +
    ggplot2::scale_y_continuous(breaks = 0:n_loo_N,
                                limits = c(-0.5, n_loo_N + 0.5),
                                name   = paste0("LOO-N score (max = ", n_loo_N, ")")) +
    ggplot2::annotate("text", x = n_loo_P, y = n_loo_N, label = "Tier1\nboth sides",
                      hjust = 1.1, vjust = 1.3, size = 3, colour = "#1A237E", fontface = "bold") +
    ggplot2::annotate("text", x = 0, y = n_loo_N, label = "Fragile P\nrobust N",
                      hjust = -0.05, vjust = 1.3, size = 3, colour = "#00897B") +
    ggplot2::annotate("text", x = n_loo_P, y = 0, label = "Robust P\nfragile N",
                      hjust = 1.1, vjust = -0.3, size = 3, colour = "#E53935") +
    ggplot2::labs(title    = "LOO-P score vs LOO-N score",
                  subtitle = "Each point = one significant gene | top-right = most robust",
                  colour   = "Tier (global)", shape = "Direction") +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "right")
  
  save_plot_pdf_svg(p_scatter, file.path(dir_deg, "LOO_P_vs_N_Scatter"),
                    width = 7.5, height = 6)
  message("✅ LOO-P vs LOO-N scatter saved")
  
  
  ## 7C.5 Consistency heatmap (Tier1_both genes)
  genes_hm <- tier_compare_df |>
    dplyr::filter(Tier_Global == "Tier1_both") |>
    dplyr::pull(Gene)
  
  all_loo_results <- c(loo_results_P, loo_results_N)
  
  if (length(genes_hm) > 0 && length(all_loo_results) > 0) {
    
    loo_mat <- sapply(all_loo_results, function(df) {
      sapply(genes_hm, function(g) {
        idx   <- match(g, df$Gene)
        if (is.na(idx)) return(0L)
        dir_g <- tier_compare_df$Status[tier_compare_df$Gene == g]
        as.integer(df$Status[idx] == dir_g)
      })
    })
    rownames(loo_mat) <- genes_hm
    colnames(loo_mat) <- names(all_loo_results)
    
    ann_col <- data.frame(
      LOO_Side = sapply(all_loo_results, function(df) unique(df$Side)),
      row.names = names(all_loo_results)
    )
    ann_row <- data.frame(
      Direction = tier_compare_df$Status[match(genes_hm, tier_compare_df$Gene)],
      row.names = genes_hm
    )
    
    save_pheatmap_pdf_svg(
      list(
        mat               = loo_mat,
        color             = c("#F5F5F5", "#1A237E"),
        breaks            = c(-0.5, 0.5, 1.5),
        annotation_col    = ann_col,
        annotation_row    = ann_row,
        annotation_colors = list(LOO_Side  = c("P" = "#C0392B", "N" = "#2471A3"),
                                 Direction = DEG_COLORS[c("Up", "Down")]),
        cluster_rows      = TRUE,
        cluster_cols      = FALSE,
        show_rownames     = length(genes_hm) <= 80,
        show_colnames     = TRUE,
        fontsize          = 7,
        fontsize_row      = 6,
        border_color      = "grey90",
        legend_breaks     = c(0, 1),
        legend_labels     = c("Not DE", "DE"),
        gaps_col          = n_loo_P,
        main              = paste0(
          "LOO consistency — Tier1 (both sides) genes\n",
          "Left = LOO-P (n=", n_loo_P, ") | Right = LOO-N (n=", n_loo_N, ")"
        )
      ),
      file.path(dir_deg, "LOO_Consistency_Heatmap_Tier1both"),
      w = max(8, (n_loo_P + n_loo_N) * 0.9 + 3),
      h = max(7, min(length(genes_hm) * 0.16 + 3, 24))
    )
    message("✅ LOO consistency heatmap saved")
    
  } else {
    message("No Tier1-both genes — heatmap skipped")
  }
  
  
  ## 7C.6 Interpretation summary
  interp_df <- data.frame(
    Category = c("Tier1_both", "Tier1_P_only", "Tier1_N_only", "Tier2", "Tier3"),
    N_genes = c(
      nrow(tier1_both),
      nrow(tier1_P_only),
      nrow(tier1_N_only),
      sum(tier_compare_df$Tier_Global == "Tier2"),
      sum(tier_compare_df$Tier_Global == "Tier3")
    ),
    Interpretation = c(
      "Robust regardless of which genotype (P or N) is removed — strongest candidates",
      "Robust to P-side perturbations but sensitive to N composition",
      "Robust to N-side perturbations but sensitive to P composition",
      "Partially robust — consistent in at least half of all LOO tests",
      "Weak LOO support — signal is genotype-dependent, interpret with caution"
    )
  )
  
  write.csv(interp_df,
            file.path(dir_deg, "LOO_Interpretation_Summary.csv"),
            row.names = FALSE)
  
  cat("\n=== LOO Comparative Interpretation ===\n")
  print(interp_df, right = FALSE)
  
  ## 7C.7 Split tiers by direction (Up / Down)
  dir_tier_split <- file.path(dir_deg, "Tier_Split_By_Direction")
  dir.create(dir_tier_split, recursive = TRUE, showWarnings = FALSE)
  
  tier_split_lists <- list(
    Tier1_both_Up      = tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_both",    Status == "Up"),
    Tier1_both_Down    = tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_both",    Status == "Down"),
    Tier1_P_only_Up    = tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_P_only",  Status == "Up"),
    Tier1_P_only_Down  = tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_P_only",  Status == "Down"),
    Tier1_N_only_Up    = tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_N_only",  Status == "Up"),
    Tier1_N_only_Down  = tier_compare_df |> dplyr::filter(Tier_Global == "Tier1_N_only",  Status == "Down"),
    Tier1P_Tier2N_Up   = tier_compare_df |> dplyr::filter(Tier_P == "Tier1", Tier_N == "Tier2", Status == "Up"),
    Tier1P_Tier2N_Down = tier_compare_df |> dplyr::filter(Tier_P == "Tier1", Tier_N == "Tier2", Status == "Down")
  )
  
  for (nm in names(tier_split_lists)) {
    write.csv(
      tier_split_lists[[nm]],
      file.path(dir_tier_split, paste0(nm, ".csv")),
      row.names = FALSE
    )
  }
  
  tier_split_summary <- data.frame(
    GeneSet = names(tier_split_lists),
    N_genes = sapply(tier_split_lists, nrow),
    row.names = NULL
  )
  
  write.csv(
    tier_split_summary,
    file.path(dir_tier_split, "Tier_Split_By_Direction_Summary.csv"),
    row.names = FALSE
  )
  
  cat("\nTier split by direction:\n")
  print(tier_split_summary)
  cat("\n")
  message("✅ Tier split by direction exported")
  
  ## 7.FINAL
  loo_results_list <- all_loo_results
  tier_results     <- tier_compare_df
  
  cat("\n========================================================\n")
  cat("PART 7 completed\n")
  cat("========================================================\n")
  cat("  LOO-P contrasts     :", n_loo_P, "\n")
  cat("  LOO-N contrasts     :", n_loo_N, "\n")
  cat("  Total LOO contrasts :", n_loo_P + n_loo_N, "\n")
  cat("  ---\n")
  cat("  Tier1 BOTH sides    :", nrow(tier1_both),   "\n")
  cat("  Tier1 LOO-P only    :", nrow(tier1_P_only), "\n")
  cat("  Tier1 LOO-N only    :", nrow(tier1_N_only), "\n")
  cat("  Tier2               :", sum(tier_compare_df$Tier_Global == "Tier2"), "\n")
  cat("  Tier3               :", sum(tier_compare_df$Tier_Global == "Tier3"), "\n")
  cat("========================================================\n\n")
}



# ============================================================
# ============================================================
# PART 8 — GO ENRICHMENT ON LOO-BASED TIERS  (topGO)
# ============================================================
#
#  Gene sets:
#    8A  Tier1_both    — robust on both sides
#    8B  Tier1_P_only  — robust on LOO-P side only
#    8C  Tier1_N_only  — robust on LOO-N side only
#    8D  Tier1P_Tier2N — Tier1-P but Tier2-N
#
#  For each gene set:
#    topGO weight01 Fisher (BP, MF, CC)
#    Dotplot + Barplot per gene set
#  Combined tier comparative dotplot (GeneSet × Ontology)
#
# ============================================================

if (!RUN_GO_PRIMARY) {
  message("Tier-based GO analysis skipped")
} else {
  
  cat("---- PART 8: GO enrichment on LOO-based tiers (topGO) ----\n")
  
  
  ## 8.1 Check required objects
  if (!exists("tier_compare_df"))
    stop("tier_compare_df not found. Run PART 7 first.")
  
  ## Reload GO annotation if needed (Part 6 may have been skipped)
  if (!exists("geneID2GO") || is.null(geneID2GO)) {
    if (!file.exists(go_annot_file))
      stop("GO annotation file not found: ", go_annot_file)
    geneID2GO     <- load_go_annotation(go_annot_file)
    genes_tested  <- rownames(dge)
    gene_universe <- intersect(genes_tested, names(geneID2GO))
  }
  
  dir_go_tiers <- file.path(dir_go, "Tier_based_GO")
  dir.create(dir_go_tiers, recursive = TRUE, showWarnings = FALSE)
  
  message("\u2705 Tier GO output folder ready")
  
  
  ## 8.2 Define gene sets (split by direction: Up / Down) ----
  tier_compare_df <- tier_compare_df |>
    dplyr::mutate(
      Tier_P      = as.character(Tier_P),
      Tier_N      = as.character(Tier_N),
      Tier_Global = as.character(Tier_Global),
      Status      = as.character(Status)
    )
  
  gene_sets <- list(
    Tier1_both_Up      = tier_compare_df |>
      dplyr::filter(Tier_Global == "Tier1_both",   Status == "Up")   |> dplyr::pull(Gene),
    Tier1_both_Down    = tier_compare_df |>
      dplyr::filter(Tier_Global == "Tier1_both",   Status == "Down") |> dplyr::pull(Gene),
    
    Tier1_P_only_Up    = tier_compare_df |>
      dplyr::filter(Tier_Global == "Tier1_P_only", Status == "Up")   |> dplyr::pull(Gene),
    Tier1_P_only_Down  = tier_compare_df |>
      dplyr::filter(Tier_Global == "Tier1_P_only", Status == "Down") |> dplyr::pull(Gene),
    
    Tier1_N_only_Up    = tier_compare_df |>
      dplyr::filter(Tier_Global == "Tier1_N_only", Status == "Up")   |> dplyr::pull(Gene),
    Tier1_N_only_Down  = tier_compare_df |>
      dplyr::filter(Tier_Global == "Tier1_N_only", Status == "Down") |> dplyr::pull(Gene),
    
    Tier1P_Tier2N_Up   = tier_compare_df |>
      dplyr::filter(Tier_P == "Tier1", Tier_N == "Tier2", Status == "Up")   |> dplyr::pull(Gene),
    Tier1P_Tier2N_Down = tier_compare_df |>
      dplyr::filter(Tier_P == "Tier1", Tier_N == "Tier2", Status == "Down") |> dplyr::pull(Gene)
  )
  
  gene_set_summary <- data.frame(
    GeneSet = names(gene_sets),
    N_genes = sapply(gene_sets, length),
    row.names = NULL
  )
  
  write.csv(gene_set_summary,
            file.path(dir_go_tiers, "Tier_GO_GeneSet_Summary.csv"),
            row.names = FALSE)
  
  cat("Gene sets for GO enrichment:\n")
  print(gene_set_summary)
  cat("\n")
  
  
  ## 8.3 Run topGO for each gene set -------------------------
  all_tier_go <- list()
  
  for (gs_name in names(gene_sets)) {
    
    genes_oi <- unique(gene_sets[[gs_name]])
    cat("\n--- Gene set:", gs_name, "(n =", length(genes_oi), ") ---\n")
    
    if (length(genes_oi) < GO_MIN_GENES) {
      warning("Gene set too small for topGO: ", gs_name,
              " (n=", length(genes_oi), ")", call. = FALSE)
      next
    }
    
    res_gs <- run_topgo_all_ont(
      genes_oi      = genes_oi,
      gene_universe = gene_universe,
      geneID2GO     = geneID2GO
    )
    
    if (is.null(res_gs) || nrow(res_gs) == 0) {
      message("  No significant GO terms for: ", gs_name)
      next
    }
    
    res_gs$GeneSet <- gs_name
    all_tier_go[[gs_name]] <- res_gs
    
    ## Export per gene set
    gs_dir <- file.path(dir_go_tiers, gs_name)
    dir.create(gs_dir, recursive = TRUE, showWarnings = FALSE)
    
    write.csv(res_gs,
              file.path(gs_dir, paste0("GO_", gs_name, "_all_ont.csv")),
              row.names = FALSE)
    
    for (ont in unique(res_gs$Ontology)) {
      write.csv(res_gs[res_gs$Ontology == ont, ],
                file.path(gs_dir, paste0("GO_", gs_name, "_", ont, ".csv")),
                row.names = FALSE)
    }
    
    ## Publication dotplot per gene set
    p_gs_pub <- plot_go_dot_simple(
      res_gs,
      title = paste0("GO enrichment — ", gs_name),
      top_n = TOP_N_TERMS
    )
    
    save_plot_pdf_svg(
      p_gs_pub,
      file.path(gs_dir, paste0("GO_", gs_name, "_Publication")),
      width  = 8,
      height = max(5, min(14, nrow(res_gs) * 0.22))
    )
    message("✅ Publication GO figure saved: ", gs_name)
  }
  
  
  ## 8.4 Combine + comparative plot --------------------------
  if (length(all_tier_go) == 0) {
    
    warning("No significant GO terms found for any tier gene set.", call. = FALSE)
    
  } else {
    
    go_tier_combined <- dplyr::bind_rows(all_tier_go) |>
      dplyr::arrange(GeneSet, Ontology, weight01)
    
    write.csv(go_tier_combined,
              file.path(dir_go_tiers, "GO_Tier_Based_Combined.csv"),
              row.names = FALSE)
    
    message("\u2705 Combined tier-based GO table exported")
    
    ## Comparative publication figure across all tiers
    n_gs <- length(unique(go_tier_combined$GeneSet))
    
    p_tier_pub <- plot_go_dot_grouped(
      go_tier_combined,
      title = "GO enrichment — LOO-based tier categories",
      top_n = 6
    )
    
    if (!is.null(p_tier_pub)) {
      save_plot_pdf_svg(
        p_tier_pub,
        file.path(dir_go_tiers, "GO_Tier_Comparative_Publication"),
        width  = max(10, n_gs * 2.8 + 2),
        height = max(6, min(20, nrow(go_tier_combined) * 0.16))
      )
      message("✅ Tier-comparative publication GO figure saved")
    }
  }
  
  
  cat("----------------------------------------\n")
  cat("PART 8 completed\n")
  cat("----------------------------------------\n")
  cat("topGO enrichment done for:\n")
  for (gs in names(gene_sets)) {
    n_terms <- if (gs %in% names(all_tier_go)) nrow(all_tier_go[[gs]]) else 0
    cat("  -", gs, "(n_genes =", length(gene_sets[[gs]]),
        "| n_terms =", n_terms, ")\n")
  }
  cat("Results in:", dir_go_tiers, "\n")
  cat("----------------------------------------\n\n")
}

## ============================================================
## SESSION INFO
## ============================================================
sink(file.path(main_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("✅ Pipeline v6 completed — session info saved")