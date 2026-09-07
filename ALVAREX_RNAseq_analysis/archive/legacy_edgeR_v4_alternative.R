############################################################
##
##  edgeR — Differential Expression Analysis
##  ─────────────────────────────────────────────────────────
##  VERSION 4.0 — ALTERNATIVE STRATEGY
##
##  RATIONALE :
##    Le design ne contient QUE des génotypes traités.
##    Une comparaison _P vs _P directe capture l'effet génotype,
##    pas l'effet protection. Cette version implémente une
##    stratégie en 3 niveaux pour isoler un signal Protection
##    robuste malgré l'absence de contrôle non-traité.
##
##  STRATÉGIE PRINCIPALE :
##    [S1] Contraste Protection_mean  (~ 0 + Genotype + Replicate)
##         → Test formel de l'effet Protection moyen (≡ design ~ Protection)
##         → Contrôle de l'effet génotype via le modèle (pas par intersection)
##         → C'est l'ANALYSE PRIMAIRE
##
##    [S2] Leave-One-Out sur les _P (LOO)
##         → Retirer un _P à la fois, recalculer Protection_mean
##         → Identifier les gènes STABLES = robustes au fond génétique
##         → Score de robustesse LOO par gène (0 → n_P)
##
##    [S3] Cross-N intersection (formalisation de ta stratégie initiale)
##         → Pour chaque _N : intersection des DEGs de tous les _P vs ce _N
##         → Intersection ENTRE les deux _N de ce core
##         → Ces gènes sont DEGs indépendamment du N de référence
##         → Interprétation : signal Protection non artefactuel de la basale
##
##    [S4] Robustness Score global
##         → Pour chaque gène : sur combien de contrastes pairwise est-il DE ?
##         → Combiné avec LOO → Tier 1 (très robuste) / Tier 2 / Tier 3
##
##  CORRECTIONS TECHNIQUES (héritées v4) :
##   [C2]  filterByExpr()
##   [C3]  BCV + QL dispersion plots
##   [C4]  MDS plot
##   [C5]  MA plots par contraste
##   [C6]  Mean-variance plot
##   [C7]  pval_lim cappé
##   [C8]  FDR par contraste BH
##   [C9]  GO par contraste + intersections
##  [C10]  separate_longer_delim
##  [C13]  Concordance log2FC Spearman
##   [B1]  write.csv
##   [B2]  avg_log_cpm précalculé
##   [B5]  coalesce(Status_a, "NS")
##   [B6]  bind_rows une seule fois
##
##  NOUVEAU v4-ALT :
##   [A1]  LOO (Leave-One-Out) sur les _P
##   [A2]  Score de robustesse LOO par gène
##   [A3]  Cross-N intersection formalisée
##   [A4]  Tier system (Tier1/2/3) pour priorisation des gènes
##   [A5]  PCA coloré par Protection + Génotype (vérification séparation)
##   [A6]  Volcano annoté par Tier
##   [A7]  GO sur Tier1 uniquement
##
##  Références :
##    Robinson et al. 2010 Bioinformatics
##    Chen et al. 2016 F1000Research
##    McCarthy et al. 2012 Nucleic Acids Research
##    Law et al. 2016 F1000Research
##    Conway et al. 2017 Bioinformatics
##
############################################################


# ============================================================
# SECTION 0 — PACKAGES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

cran_pkgs <- c(
  "ggplot2", "ggrepel", "RColorBrewer", "dplyr", "tibble",
  "tidyr", "forcats", "stringr", "svglite", "openxlsx",
  "UpSetR", "VennDiagram", "ggVennDiagram",
  "scales", "patchwork", "grid", "gridExtra", "corrplot"
)
bioc_pkgs <- c("edgeR", "limma", "pheatmap", "topGO")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg, ask = FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
message("✅ All packages loaded")


# ============================================================
# SECTION 1 — MASTER COLOR SYSTEM
# ============================================================

GENOTYPE_PALETTE <- c(
  ## _N genotypes — cool tones
  "Amsterdam_N"      = "#0A9E6E",
  "NantaiseInbred_N" = "#1B4FBF",
  "Orleans_N"        = "#007B8A",
  "Dijon_N"          = "#3949AB",
  "Genevieve_N"      = "#006B3C",
  ## _P genotypes — warm tones
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

## [A4] Tier colors
TIER_COLORS <- c(
  "Tier1" = "#1A237E",   ## Robuste LOO + Cross-N → très haute confiance
  "Tier2" = "#F57F17",   ## Robuste dans l'une des deux validations
  "Tier3" = "#78909C",   ## Significatif dans Protection_mean seulement
  "NS"    = "#ECEFF1"
)

ONT_COLORS <- c(
  "BP" = "#E76F51",
  "MF" = "#2A9D8F",
  "CC" = "#5C4B8A"
)

HEATMAP_DIVERGING <- grDevices::colorRampPalette(
  c("#2166AC","#4393C3","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B")
)(100)

HEATMAP_SEQUENTIAL <- grDevices::colorRampPalette(
  c("#FFFFFF","#FED8B1","#FF7F00","#8B1A0A")
)(100)

HEATMAP_DISTANCE <- grDevices::colorRampPalette(
  c("#1A1A2E","#16213E","#0F3460","#A8C7FA","#E8F4F8","#FFFFFF")
)(255)

resolve_genotype_colors <- function(genotype_levels) {
  known   <- genotype_levels[genotype_levels %in% names(GENOTYPE_PALETTE)]
  unknown <- genotype_levels[!genotype_levels %in% names(GENOTYPE_PALETTE)]
  if (length(unknown) > 0) {
    hues     <- seq(15, 375, length.out = length(unknown) + 1)[seq_len(length(unknown))]
    fallback <- setNames(grDevices::hsv(hues / 360, s = 0.82, v = 0.75), unknown)
    warning(length(unknown), " genotype(s) absent de GENOTYPE_PALETTE → auto-color: ",
            paste(unknown, collapse = ", "), call. = FALSE)
    return(c(GENOTYPE_PALETTE[known], fallback)[genotype_levels])
  }
  GENOTYPE_PALETTE[genotype_levels]
}

message("✅ Master color system defined")


# ============================================================
# SECTION 2 — PUBLICATION THEME
# ============================================================

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

message("✅ Publication theme defined")


# ============================================================
# SECTION 3 — PARAMÈTRES GLOBAUX & CHEMINS
# ============================================================

SEED          <- 123
alpha         <- 0.05
lfc_threshold <- 1

set.seed(SEED)
options(stringsAsFactors = FALSE)

## Fichiers d'entrée
counts_file   <- "/path/to/input/Count.csv"
meta_file     <- "/path/to/input/MetaData.csv"
go_annot_file <- "/path/to/annotation/gene_to_GO.csv"

## Répertoires de sortie
main_dir     <- "/path/to/results/legacy_edgeR_v4_alternative"
qc_dir       <- file.path(main_dir, "QC_Plots")
diag_dir     <- file.path(main_dir, "QC_ModelDiagnostics")
deg_dir      <- file.path(main_dir, "DEG_Results")
volcano_dir  <- file.path(main_dir, "Volcano_Plots")
ma_dir       <- file.path(main_dir, "MA_Plots")
heatmap_dir  <- file.path(main_dir, "Heatmaps")
summary_dir  <- file.path(main_dir, "Summaries")
venn_dir     <- file.path(main_dir, "Venn_Plots")
upset_dir    <- file.path(main_dir, "UpSet_Plots")
go_dir       <- file.path(main_dir, "GO_Analysis")
concord_dir  <- file.path(main_dir, "Concordance")
loo_dir      <- file.path(main_dir, "LOO_Validation")       ## [A1]
crossN_dir   <- file.path(main_dir, "CrossN_Intersection")  ## [A3]
robust_dir   <- file.path(main_dir, "Robustness_Score")     ## [A4]

for (d in c(main_dir, qc_dir, diag_dir, deg_dir, volcano_dir, ma_dir,
            heatmap_dir, summary_dir, venn_dir, upset_dir, go_dir,
            concord_dir, loo_dir, crossN_dir, robust_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("✅ Directories created")


# ============================================================
# SECTION 4 — FONCTIONS UTILITAIRES
# ============================================================

normalize_sample_name <- function(x) {
  x <- trimws(as.character(x))
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

save_figure <- function(p, base_path, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(base_path, ".pdf"), p,
                  width = width, height = height, bg = "white")
  ggplot2::ggsave(paste0(base_path, ".png"), p,
                  width = width, height = height, dpi = dpi, bg = "white")
  svglite::svglite(paste0(base_path, ".svg"), width = width, height = height)
  print(p)
  grDevices::dev.off()
  message("   Saved: ", basename(base_path), "  [PDF / PNG / SVG]")
}

save_base_plot <- function(plot_fn, base_path, width = 8, height = 6) {
  grDevices::pdf(paste0(base_path, ".pdf"),
                 width = width, height = height, useDingbats = FALSE)
  plot_fn(); grDevices::dev.off()
  grDevices::png(paste0(base_path, ".png"),
                 width = width * 220, height = height * 220, res = 220)
  plot_fn(); grDevices::dev.off()
  svglite::svglite(paste0(base_path, ".svg"), width = width, height = height)
  plot_fn(); grDevices::dev.off()
  message("   Saved: ", basename(base_path), "  [PDF / PNG / SVG]")
}

save_pheatmap <- function(ph_args, base_path, w = 12, h = 10) {
  do.call(pheatmap::pheatmap,
          c(ph_args, list(filename = paste0(base_path, ".pdf"), width = w, height = h)))
  do.call(pheatmap::pheatmap,
          c(ph_args, list(filename = paste0(base_path, ".png"), width = w, height = h)))
  ph_obj <- do.call(pheatmap::pheatmap, c(ph_args, list(silent = TRUE)))
  svglite::svglite(paste0(base_path, ".svg"), width = w, height = h)
  grid::grid.newpage(); grid::grid.draw(ph_obj$gtable); grDevices::dev.off()
  message("✅ Heatmap: ", basename(base_path), "  [PDF / PNG / SVG]")
}

get_set_colors <- function(set_names) {
  vapply(set_names, function(nm) {
    if (nm %in% names(genotype_colors)) return(unname(genotype_colors[nm]))
    hit <- names(genotype_colors)[startsWith(names(genotype_colors), nm)]
    if (length(hit) > 0) return(unname(genotype_colors[hit[1]]))
    "#999999"
  }, character(1))
}

export_membership_csv <- function(gene_lists, file_path) {
  all_g   <- unique(unlist(gene_lists))
  if (length(all_g) == 0) return(invisible(NULL))
  bin_mat <- as.data.frame(
    lapply(gene_lists, function(v) as.integer(all_g %in% v))
  )
  rownames(bin_mat) <- all_g
  bin_mat$N_sets  <- rowSums(bin_mat)
  bin_mat$Robust  <- ifelse(bin_mat$N_sets == length(gene_lists),
                            "Robust", "Specific")
  write.csv(tibble::rownames_to_column(bin_mat, "Gene"), file_path, row.names = FALSE)
  cat("    Core:", sum(bin_mat$Robust == "Robust"),
      "| Specific:", sum(bin_mat$Robust == "Specific"), "\n")
}

message("✅ Helper functions defined")


# ============================================================
# SECTION 5 — IMPORT DES DONNÉES
# ============================================================
cat("---- IMPORT ----\n")

counts_data           <- read.csv(counts_file, header = TRUE, check.names = FALSE)
rownames(counts_data) <- counts_data$Geneid
counts_data$Geneid    <- NULL
counts_data           <- as.matrix(counts_data)
mode(counts_data)     <- "numeric"
colnames(counts_data) <- normalize_sample_name(colnames(counts_data))

sample_info_raw <- read.csv(meta_file, header = TRUE, check.names = FALSE)
sample_info <- sample_info_raw |>
  dplyr::mutate(
    Sample    = normalize_sample_name(Sample),
    Genotype  = fix_encoding(trimws(Genotype)),
    Replicate = as.integer(Replicate)
  )

common_samples <- intersect(colnames(counts_data), sample_info$Sample)
cat("Matched samples:", length(common_samples), "\n")
stopifnot(length(common_samples) > 0)

sample_info <- sample_info |>
  dplyr::filter(Sample %in% common_samples) |>
  dplyr::arrange(match(Sample, colnames(counts_data)))

counts_data           <- counts_data[, sample_info$Sample, drop = FALSE]
rownames(sample_info) <- sample_info$Sample
stopifnot(all(colnames(counts_data) == sample_info$Sample))

sample_info <- sample_info |>
  dplyr::mutate(
    Protection = factor(
      dplyr::case_when(
        grepl("_P$", Genotype) ~ "Protected",
        grepl("_N$", Genotype) ~ "NonProtected",
        TRUE ~ NA_character_
      ),
      levels = c("NonProtected", "Protected")
    ),
    Genotype  = factor(Genotype,  levels = sort(unique(Genotype))),
    Replicate = factor(Replicate, levels = sort(unique(as.integer(Replicate))))
  )

if (any(is.na(sample_info$Protection))) {
  warning("Samples sans suffixe _P ou _N détectés — exclus.", call. = FALSE)
  sample_info <- dplyr::filter(sample_info, !is.na(Protection))
  counts_data <- counts_data[, sample_info$Sample, drop = FALSE]
}

conditions_P <- levels(sample_info$Genotype)[grepl("_P$", levels(sample_info$Genotype))]
conditions_N <- levels(sample_info$Genotype)[grepl("_N$", levels(sample_info$Genotype))]
n_P          <- length(conditions_P)
n_N          <- length(conditions_N)

cat("\n_P:", paste(conditions_P, collapse = ", "), "\n")
cat("_N:", paste(conditions_N, collapse = ", "), "\n")
cat("Distribution Protection:\n"); print(table(sample_info$Protection))

genotype_colors <- resolve_genotype_colors(levels(sample_info$Genotype))
message("✅ Data imported")


# ============================================================
# SECTION 6 — MODÈLE edgeR (FULL MODEL)
# ============================================================
cat("---- EDGER MODEL ----\n")

dge <- edgeR::DGEList(
  counts = round(counts_data),
  group  = sample_info$Genotype
)

design_for_filter <- model.matrix(~ Protection + Replicate, data = sample_info)
keep <- edgeR::filterByExpr(dge, design = design_for_filter)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
cat("Genes retenus (filterByExpr):", sum(keep), "/", length(keep),
    "(", round(100 * sum(keep) / length(keep), 1), "%)\n")

dge <- edgeR::calcNormFactors(dge, method = "TMMwsp")
cat("TMM normalization factors:\n"); print(round(dge$samples$norm.factors, 4))

## Design ~ 0 + Genotype + Replicate
## Le contraste Protection_mean = (∑_P/n_P) − (∑_N/n_N)
## est équivalent au coefficient Protection dans un modèle
## ~ Protection + Replicate, mais permet aussi les pairwise.
design <- model.matrix(~ 0 + Genotype + Replicate, data = sample_info)
colnames(design) <- gsub("^Genotype", "", colnames(design))
colnames(design) <- make.names(colnames(design))
cat("\nDesign columns:\n"); print(colnames(design))

geno_map <- setNames(make.names(levels(sample_info$Genotype)),
                     levels(sample_info$Genotype))

dge <- edgeR::estimateDisp(dge, design, robust = TRUE)
cat("\nCommon BCV:", round(sqrt(dge$common.dispersion), 4), "\n")

fit <- edgeR::glmQLFit(dge, design, robust = TRUE)
avg_log_cpm <- edgeR::aveLogCPM(dge)

message("✅ edgeR model fitted")


# ============================================================
# SECTION 7 — QC + DIAGNOSTICS (ÉTENDU) [A5]
#
##  Ajout : PCA coloré par Protection (vérifier séparation)
##  Si PC1 capture le statut Protection → signal interprétable
##  Si PC1 capture le génotype → bruit génétique dominant
# ============================================================
cat("---- QC ----\n")

logcpm_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

## ── PCA par Génotype ────────────────────────────────────────────────
pca_res  <- prcomp(t(logcpm_mat), scale. = FALSE)
pca_df   <- as.data.frame(pca_res$x[, 1:2])
pca_df$Genotype   <- sample_info$Genotype
pca_df$Replicate  <- sample_info$Replicate
pca_df$Protection <- sample_info$Protection
percentVar <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2))

rep_shapes <- setNames(
  c(16, 17, 15, 18, 8, 3)[seq_along(levels(sample_info$Replicate))],
  levels(sample_info$Replicate)
)

p_pca_geno <- ggplot2::ggplot(pca_df,
                              ggplot2::aes(PC1, PC2, colour = Genotype, shape = Replicate)) +
  ggplot2::geom_point(size = 4, stroke = 0.5, alpha = 0.92) +
  ggplot2::scale_colour_manual(values = genotype_colors) +
  ggplot2::scale_shape_manual(values = rep_shapes) +
  ggplot2::xlab(paste0("PC1  (", percentVar[1], "% variance)")) +
  ggplot2::ylab(paste0("PC2  (", percentVar[2], "% variance)")) +
  ggplot2::labs(title    = "PCA — logCPM TMM | Couleur = Génotype",
                subtitle = "Vérifier si les _P et _N se séparent globalement",
                caption  = paste0("n = ", nrow(dge), " genes  |  TMM")) +
  theme_pub(base_size = 11)
save_figure(p_pca_geno, file.path(qc_dir, "QC_PCA_Genotype"), width = 8, height = 5.5)

## ── [A5] PCA coloré par Protection ──────────────────────────────────
## LECTURE CLEF : si PC1 sépare Protected vs NonProtected → signal fort
## Si seul PC2 ou PC3 le fait → signal faible, interpréter avec prudence
p_pca_prot <- ggplot2::ggplot(pca_df,
                              ggplot2::aes(PC1, PC2, colour = Protection, shape = Replicate,
                                           label = Genotype)) +
  ggplot2::geom_point(size = 5, stroke = 0.7, alpha = 0.90) +
  ggrepel::geom_text_repel(size = 3, fontface = "italic",
                           max.overlaps = 20, segment.linewidth = 0.3,
                           box.padding = 0.4, colour = "grey30") +
  ggplot2::scale_colour_manual(values = PROTECTION_COLORS,
                               guide = ggplot2::guide_legend(
                                 title = "Protection status",
                                 override.aes = list(size = 5, shape = 16))) +
  ggplot2::scale_shape_manual(values = rep_shapes) +
  ggplot2::xlab(paste0("PC1  (", percentVar[1], "% variance)")) +
  ggplot2::ylab(paste0("PC2  (", percentVar[2], "% variance)")) +
  ggplot2::labs(
    title    = "PCA — Protection status  [A5]",
    subtitle = paste0(
      "Si PC1 sépare Protected vs NonProtected → signal Protection interprétable\n",
      "Si PC1 ≈ génotype → bruit génétique dominant, interpréter avec prudence"),
    caption = paste0("n = ", nrow(dge), " genes  |  filterByExpr  |  TMM")
  ) +
  theme_pub(base_size = 11) +
  ggplot2::theme(legend.position = "right")
save_figure(p_pca_prot, file.path(qc_dir, "QC_PCA_Protection_STATUS"), width = 9, height = 6)

## ── [C4] MDS plot ────────────────────────────────────────────────────
save_base_plot(
  function() {
    col_mds <- genotype_colors[as.character(sample_info$Genotype)]
    pch_mds <- c(16, 17, 15, 18, 8, 3)[as.integer(sample_info$Replicate)]
    limma::plotMDS(dge, col = col_mds, pch = pch_mds, cex = 1.6,
                   main = "MDS plot — leading logFC distances  |  edgeR QC")
    legend("topright", legend = levels(sample_info$Genotype),
           col = genotype_colors[levels(sample_info$Genotype)],
           pch = 16, cex = 0.75, bty = "n", title = "Genotype")
  },
  file.path(qc_dir, "QC_MDS"), width = 8, height = 6
)

## ── Distance heatmap ─────────────────────────────────────────────────
sampleDists <- stats::dist(t(logcpm_mat))
sdm         <- as.matrix(sampleDists)
hm_labels   <- paste(sample_info$Genotype, sample_info$Replicate, sep = "_")
rownames(sdm) <- colnames(sdm) <- hm_labels

ann_qc <- data.frame(Protection = sample_info$Protection,
                     Genotype   = sample_info$Genotype,
                     row.names  = hm_labels)
save_pheatmap(
  list(mat = sdm, color = HEATMAP_DISTANCE, border_color = NA,
       clustering_distance_rows = sampleDists,
       clustering_distance_cols = sampleDists,
       annotation_row    = ann_qc, annotation_col = ann_qc,
       annotation_colors = list(Protection = PROTECTION_COLORS,
                                Genotype   = genotype_colors),
       fontsize = 8, treeheight_row = 20, treeheight_col = 20,
       main = "Sample-to-sample distances  |  logCPM TMM"),
  file.path(qc_dir, "QC_SampleDistances"), w = 9, h = 8
)

## ── [C3] BCV + QL dispersion ─────────────────────────────────────────
save_base_plot(
  function() {
    edgeR::plotBCV(dge,
                   main = paste0("BCV  |  Common BCV = ",
                                 round(sqrt(dge$common.dispersion), 3)))
    abline(h = sqrt(dge$common.dispersion), col = "#C0392B", lty = 2, lwd = 1.5)
  },
  file.path(diag_dir, "Diag_BCV"), width = 7, height = 5.5
)

save_base_plot(
  function() edgeR::plotQLDisp(fit,
                               main = "QL dispersion  |  quasi-likelihood F-test  |  robust = TRUE"),
  file.path(diag_dir, "Diag_QLDisp"), width = 7, height = 5.5
)

save_base_plot(
  function() edgeR::plotMeanVar(dge,
                                show.raw.vars = TRUE, show.tagwise.vars = TRUE, NBline = TRUE,
                                main = "Mean-variance relationship  |  TMM"),
  file.path(diag_dir, "Diag_MeanVar"), width = 7, height = 5.5
)

message("✅ QC + diagnostics saved")


# ============================================================
# SECTION 8 — CONSTRUCTION DES CONTRASTES
#
##  HIÉRARCHIE DES ANALYSES :
##
##  [PRIMARY] Protection_mean
##    → Effet moyen Protection corrigé par le modèle (~ 0 + Genotype)
##    → Équivalent formel à ~ Protection + Genotype
##    → Produit les DEGs avec le plus grand pouvoir statistique
##
##  [VALIDATION 1] LOO (Section 9b) — vérification robustesse sur _P
##  [VALIDATION 2] Cross-N intersection (Section 9c) — vérification
##    que le signal n'est pas artefact de la basale choisie
##
##  [SUPPORT] Pairwise (Section 9a) — alimentation des Venns/UpSets
# ============================================================
cat("---- CONTRASTES ----\n")

p_design_names <- geno_map[conditions_P]
n_design_names <- geno_map[conditions_N]

contrast_mean_str <- paste0(
  "(", paste(p_design_names, collapse = " + "), ") / ", n_P,
  " - (",
  paste(n_design_names, collapse = " + "), ") / ", n_N
)
cat("Contraste PRIMARY Protection_mean:\n  ", contrast_mean_str, "\n\n")

contrast_pairs <- expand.grid(
  Condition_P = conditions_P,
  Condition_N = conditions_N,
  stringsAsFactors = FALSE
) |> dplyr::arrange(Condition_P, Condition_N)

cat("Contrastes pairwise (support):", nrow(contrast_pairs), "\n")

## Contrastes LOO : pour chaque _P retiré, recalculer la moyenne
## sans lui → (∑_{P\P_i} / (n_P-1)) - (∑_N / n_N)
loo_contrasts <- lapply(conditions_P, function(gP_removed) {
  remaining_P <- p_design_names[p_design_names != geno_map[gP_removed]]
  if (length(remaining_P) == 0) return(NULL)
  ct_str <- paste0(
    "(", paste(remaining_P, collapse = " + "), ") / ", length(remaining_P),
    " - (",
    paste(n_design_names, collapse = " + "), ") / ", n_N
  )
  list(name = paste0("LOO_without_", gP_removed),
       str  = ct_str,
       removed = gP_removed)
})
loo_contrasts <- Filter(Negate(is.null), loo_contrasts)
cat("Contrastes LOO:", length(loo_contrasts), "\n\n")


# ============================================================
# SECTION 9a — ANALYSE DIFFÉRENTIELLE (CONTRASTES PRIMAIRES + PAIRWISE)
# ============================================================
cat("---- ANALYSE DIFFÉRENTIELLE ----\n")

run_contrast <- function(contrast_str, contrast_name, gP = NA, gN = NA) {
  cvec <- tryCatch(
    limma::makeContrasts(contrasts = contrast_str, levels = design),
    error = function(e) {
      message("⚠️  makeContrasts échoué: ", contrast_name, " — ", e$message)
      NULL
    }
  )
  if (is.null(cvec)) return(NULL)
  
  qlf <- edgeR::glmQLFTest(fit, contrast = cvec)
  
  res <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table |>
    tibble::rownames_to_column("Gene") |>
    dplyr::rename(log2FoldChange = logFC, pvalue = PValue, padj = FDR) |>
    dplyr::mutate(
      padj          = ifelse(is.na(padj), 1, padj),
      Contrast      = contrast_name,
      Condition_P   = gP,
      Condition_N   = gN,
      aveLogCPM     = avg_log_cpm[Gene],
      Status = dplyr::case_when(
        padj < alpha & log2FoldChange >  lfc_threshold ~ "Up",
        padj < alpha & log2FoldChange < -lfc_threshold ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  n_up   <- sum(res$Status == "Up")
  n_down <- sum(res$Status == "Down")
  cat(sprintf("  %-55s  UP: %4d | DOWN: %4d\n", contrast_name, n_up, n_down))
  
  write.csv(res,
            file.path(deg_dir, paste0("DEG_", contrast_name, ".csv")),
            row.names = FALSE)
  res
}

all_results  <- list()
summary_list <- list()

## [PRIMARY] Contraste Protection_mean
cat("---- [PRIMARY] Protection_mean ----\n")
res_mean <- run_contrast(contrast_mean_str, "Protection_mean", "All_P", "All_N")
if (!is.null(res_mean)) {
  all_results[["Protection_mean"]] <- res_mean
  summary_list[["Protection_mean"]] <- data.frame(
    Contrast = "Protection_mean", Type = "Primary",
    Condition_P = "All_P", Condition_N = "All_N",
    N_up   = sum(res_mean$Status == "Up"),
    N_down = sum(res_mean$Status == "Down"),
    N_total = sum(res_mean$Status != "NS"))
}

## Contrastes pairwise (support pour Venn/UpSet)
cat("---- Contrastes pairwise (support) ----\n")
for (i in seq_len(nrow(contrast_pairs))) {
  gP <- contrast_pairs$Condition_P[i]
  gN <- contrast_pairs$Condition_N[i]
  ct_name <- paste0(gP, "_vs_", gN)
  ct_str  <- paste0(geno_map[gP], " - ", geno_map[gN])
  
  res <- run_contrast(ct_str, ct_name, gP, gN)
  if (!is.null(res)) {
    all_results[[ct_name]] <- res
    summary_list[[ct_name]] <- data.frame(
      Contrast    = ct_name, Type = "Pairwise",
      Condition_P = gP, Condition_N = gN,
      N_up        = sum(res$Status == "Up"),
      N_down      = sum(res$Status == "Down"),
      N_total     = sum(res$Status != "NS"))
  }
}

all_res_df <- dplyr::bind_rows(all_results)
summary_df <- dplyr::bind_rows(summary_list)

write.csv(all_res_df, file.path(deg_dir, "All_DEG_results_combined.csv"), row.names = FALSE)
write.csv(summary_df, file.path(summary_dir, "DEG_summary_counts.csv"),   row.names = FALSE)

pairwise_names <- names(all_results)[names(all_results) != "Protection_mean"]
cat("\nRésumé:\n"); print(summary_df); cat("\n")
message("✅ Primary + pairwise contrasts completed")


# ============================================================
# SECTION 9b — [A1][A2] LEAVE-ONE-OUT VALIDATION
#
##  OBJECTIF : Vérifier que les DEGs de Protection_mean ne sont pas
##  portés par un seul génotype _P atypique.
##
##  Pour chaque _P retiré :
##    → Recalculer le contraste mean sans lui
##    → Compter les DEGs restants + leur concordance avec Protection_mean
##
##  SCORE LOO par gène = nombre de sous-ensembles où le gène reste DE
##  dans le même sens (Up ou Down) → entre 0 et n_P
##
##  Gène avec score LOO = n_P (max) → robuste à la sortie de n'importe quel _P
##  Gène avec score LOO = 0       → n'est DE que grâce à un seul _P
# ============================================================
cat("---- [A1] LOO LEAVE-ONE-OUT VALIDATION ----\n")

loo_results  <- list()
loo_summary  <- list()

for (loo_item in loo_contrasts) {
  ct_name <- loo_item$name
  ct_str  <- loo_item$str
  removed <- loo_item$removed
  
  cat("  LOO: retirer", removed, "\n")
  res_loo <- run_contrast(ct_str, ct_name, paste0("All_P_except_", removed), "All_N")
  
  if (!is.null(res_loo)) {
    loo_results[[ct_name]] <- res_loo
    loo_summary[[ct_name]] <- data.frame(
      LOO_name    = ct_name,
      P_removed   = removed,
      N_up        = sum(res_loo$Status == "Up"),
      N_down      = sum(res_loo$Status == "Down"),
      N_total     = sum(res_loo$Status != "NS"))
    
    write.csv(res_loo,
              file.path(loo_dir, paste0("LOO_", ct_name, ".csv")),
              row.names = FALSE)
  }
}

loo_summary_df <- dplyr::bind_rows(loo_summary)
write.csv(loo_summary_df,
          file.path(loo_dir, "LOO_summary_counts.csv"),
          row.names = FALSE)

## ── [A2] Score de robustesse LOO par gène ───────────────────────────
## Pour chaque gène dans Protection_mean : dans combien de LOO
## reste-t-il DE dans le même sens ?

if (!is.null(res_mean) && length(loo_results) > 0) {
  genes_primary_up   <- res_mean$Gene[res_mean$Status == "Up"]
  genes_primary_down <- res_mean$Gene[res_mean$Status == "Down"]
  
  compute_loo_score <- function(gene_set, direction) {
    if (length(gene_set) == 0) return(data.frame(Gene = character(0), LOO_Score = integer(0)))
    scores <- sapply(gene_set, function(g) {
      sum(sapply(loo_results, function(r) {
        idx <- match(g, r$Gene)
        if (is.na(idx)) return(0L)
        as.integer(r$Status[idx] == direction)
      }))
    })
    data.frame(Gene = gene_set, LOO_Score = as.integer(scores))
  }
  
  loo_score_up   <- compute_loo_score(genes_primary_up,   "Up")
  loo_score_down <- compute_loo_score(genes_primary_down, "Down")
  
  loo_score_all <- dplyr::bind_rows(
    dplyr::mutate(loo_score_up,   Direction = "Up"),
    dplyr::mutate(loo_score_down, Direction = "Down")
  ) |>
    dplyr::mutate(
      LOO_Max        = length(loo_results),
      LOO_Robust     = (LOO_Score == LOO_Max),
      LOO_Score_Norm = LOO_Score / LOO_Max
    ) |>
    dplyr::arrange(dplyr::desc(LOO_Score))
  
  write.csv(loo_score_all,
            file.path(loo_dir, "LOO_robustness_scores.csv"),
            row.names = FALSE)
  
  ## Barplot résumé LOO
  loo_plot_df <- loo_score_all |>
    dplyr::count(Direction, LOO_Score) |>
    dplyr::mutate(
      LOO_Score = factor(LOO_Score, levels = seq(0, length(loo_results))),
      Label = paste0("Score = ", LOO_Score, "/", length(loo_results))
    )
  
  p_loo <- ggplot2::ggplot(loo_plot_df,
                           ggplot2::aes(x = LOO_Score, y = n, fill = Direction)) +
    ggplot2::geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = DEG_COLORS[c("Up", "Down")]) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(
      title    = "Score de robustesse LOO  [A2]",
      subtitle = paste0(
        "Score = nb de sous-ensembles LOO où le gène reste DE dans le même sens\n",
        "Score max = ", length(loo_results), " (robuste à la sortie de tout _P)"),
      x = paste0("Score LOO  (0 → ", length(loo_results), ")"),
      y = "Nombre de gènes",
      fill = NULL,
      caption = paste0("Basé sur ", length(loo_results), " contrastes LOO  |  FDR < ",
                       alpha, "  |  |log2FC| ≥ ", lfc_threshold)
    ) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  
  save_figure(p_loo, file.path(loo_dir, "LOO_Score_Barplot"), width = 7, height = 5)
  
  ## Scatter : LOO Score vs log2FC (Protection_mean)
  loo_scatter_df <- dplyr::left_join(
    res_mean |>
      dplyr::filter(Status != "NS") |>
      dplyr::select(Gene, log2FoldChange, padj, Status),
    loo_score_all |> dplyr::select(Gene, LOO_Score, LOO_Robust),
    by = "Gene"
  ) |> dplyr::filter(!is.na(LOO_Score))
  
  p_loo_scatter <- ggplot2::ggplot(loo_scatter_df,
                                   ggplot2::aes(x = log2FoldChange, y = LOO_Score, colour = Status,
                                                shape = LOO_Robust)) +
    ggplot2::geom_jitter(width = 0, height = 0.15, size = 2.5, alpha = 0.75) +
    ggplot2::scale_colour_manual(values = DEG_COLORS[c("Up", "Down")]) +
    ggplot2::scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = paste0("Robuste (score = ", length(loo_results), ")"),
                 "FALSE" = "Instable (score < max)")) +
    ggplot2::scale_y_continuous(breaks = seq(0, length(loo_results))) +
    ggplot2::labs(
      title    = "LOO Score vs log2FC  [A2]",
      subtitle = "Gènes avec score max = robustes à la sortie de n'importe quel _P",
      x        = expression(log[2]~"fold change  (Protection_mean)"),
      y        = paste0("Score LOO  (/ ", length(loo_results), ")"),
      colour   = NULL, shape = "Robustesse LOO"
    ) +
    theme_pub(base_size = 11)
  
  save_figure(p_loo_scatter, file.path(loo_dir, "LOO_Score_vs_LFC"),
              width = 7, height = 5)
  
  cat("✅ LOO: ", sum(loo_score_all$LOO_Robust), "gènes robustes /",
      nrow(loo_score_all), "DEGs primaires\n")
}

message("✅ LOO validation completed")


# ============================================================
# SECTION 9c — [A3] CROSS-N INTERSECTION
#
##  OBJECTIF : Formaliser ta stratégie initiale d'intersection
##  Cette analyse répond à : "Le signal observé est-il artefact
##  du génotype non-protégé choisi comme référence ?"
##
##  Pour chaque _N :
##    Core_N = ∩ { DEGs de tous les _P vs ce même _N }
##    (= "gènes DE par rapport à ce N, quel que soit le _P")
##
##  Cross-N = Core_N1 ∩ Core_N2
##    (= gènes DE indépendamment du N de référence)
##    → Signal Protection le plus propre possible
##
##  LECTURE :
##    Si Cross-N grand → signal Protection robuste à la basale choisie
##    Si Cross-N ≈ 0  → le signal change selon la basale → problème
# ============================================================
cat("---- [A3] CROSS-N INTERSECTION ----\n")

core_by_N <- list()   ## Pour chaque _N : intersection de tous les _P

for (gN in conditions_N) {
  cts_gN <- pairwise_names[grepl(paste0("_vs_", gN, "$"), pairwise_names)]
  if (length(cts_gN) == 0) next
  res_list <- all_results[cts_gN]
  
  core_up   <- Reduce(intersect, lapply(res_list, function(d) d$Gene[d$Status == "Up"]))
  core_down <- Reduce(intersect, lapply(res_list, function(d) d$Gene[d$Status == "Down"]))
  union_all <- unique(unlist(lapply(res_list, function(d) d$Gene[d$Status != "NS"])))
  
  core_by_N[[gN]] <- list(Up = core_up, Down = core_down, All = union_all)
  
  cat(sprintf("  Core_%s : Up=%d | Down=%d | Union=%d\n",
              gN, length(core_up), length(core_down), length(union_all)))
  
  write.csv(
    data.frame(
      Gene = c(core_up, core_down),
      Direction = c(rep("Up", length(core_up)), rep("Down", length(core_down))),
      Core_N = gN,
      N_contrasts_DE = sapply(c(core_up, core_down), function(g) {
        sum(sapply(res_list, function(d) {
          idx <- match(g, d$Gene); if (is.na(idx)) 0L else as.integer(d$Status[idx] != "NS")
        }))
      })
    ),
    file.path(crossN_dir, paste0("Core_DEGs_", gN, ".csv")),
    row.names = FALSE
  )
}

## Cross-N : intersection des cores de TOUS les _N
crossN_up   <- Reduce(intersect, lapply(core_by_N, `[[`, "Up"))
crossN_down <- Reduce(intersect, lapply(core_by_N, `[[`, "Down"))
crossN_all  <- c(crossN_up, crossN_down)

cat(sprintf("\n  ╔══════════════════════════════════╗\n"))
cat(sprintf("  ║ CROSS-N Up   : %5d gènes        ║\n", length(crossN_up)))
cat(sprintf("  ║ CROSS-N Down : %5d gènes        ║\n", length(crossN_down)))
cat(sprintf("  ║ CROSS-N Total: %5d gènes        ║\n", length(crossN_all)))
cat(sprintf("  ╚══════════════════════════════════╝\n\n"))

crossN_df <- data.frame(
  Gene          = c(crossN_up, crossN_down),
  Direction     = c(rep("Up", length(crossN_up)), rep("Down", length(crossN_down))),
  CrossN_Signal = TRUE
)
write.csv(crossN_df,
          file.path(crossN_dir, "CrossN_Intersection_genes.csv"),
          row.names = FALSE)

## Venn par _N si exactement 2 N (cas typique)
if (length(conditions_N) == 2) {
  gN1 <- conditions_N[1]; gN2 <- conditions_N[2]
  for (dir_tag in c("Up", "Down")) {
    gl <- list(
      core_by_N[[gN1]][[dir_tag]],
      core_by_N[[gN2]][[dir_tag]]
    )
    names(gl) <- c(gN1, gN2)
    gl <- gl[sapply(gl, length) > 0]
    
    if (length(gl) < 2) next
    
    col_v <- get_set_colors(names(gl))
    vp <- VennDiagram::venn.diagram(
      x = gl, filename = NULL,
      fill = col_v, alpha = 0.5, lwd = 2, col = "grey20",
      cex = 1.4, fontface = "bold",
      cat.cex = 1.1, cat.fontface = "bold", cat.col = col_v,
      main = paste0("Cross-N Intersection — Core ", dir_tag, "regulated  [A3]"),
      sub = paste0("Intersection = signal Protection indépendant du N\n",
                   "Zone exclusive = signal lié à la basale choisie"),
      main.cex = 1.2, sub.cex = 0.9, margin = 0.12
    )
    for (ext in c("pdf", "png", "svg")) {
      fpath <- file.path(crossN_dir, paste0("Venn_CrossN_", dir_tag, ".", ext))
      if      (ext == "pdf") grDevices::pdf(fpath, width = 9, height = 8, useDingbats = FALSE)
      else if (ext == "png") grDevices::png(fpath, width = 9 * 220, height = 8 * 220, res = 220)
      else                   svglite::svglite(fpath, width = 9, height = 8)
      grid::grid.draw(vp); grDevices::dev.off()
    }
    cat("  Cross-N Venn", dir_tag, ": saved\n")
  }
}

message("✅ Cross-N intersection completed")


# ============================================================
# SECTION 9d — [A4] TIER SYSTEM
#
##  Combiner les trois niveaux de preuve en un système de priorité :
##
##  TIER 1 (haute confiance) :
##    → DE dans Protection_mean (analyse primaire)
##    → ET LOO_Score = max (robuste à tout _P retiré)
##    → ET présent dans Cross-N (robuste à la basale choisie)
##
##  TIER 2 (confiance moyenne) :
##    → DE dans Protection_mean
##    → ET (LOO_Score = max OU présent dans Cross-N)
##    → Un des deux critères de validation seulement
##
##  TIER 3 (confiance limitée) :
##    → DE dans Protection_mean seulement
##    → Ni LOO ni Cross-N ne confirment
##    → Peut refléter un effet génotype résiduel
# ============================================================
cat("---- [A4] TIER SYSTEM ----\n")

if (!is.null(res_mean)) {
  primary_de_genes <- res_mean |>
    dplyr::filter(Status != "NS") |>
    dplyr::select(Gene, log2FoldChange, padj, Status, aveLogCPM)
  
  ## Jointure LOO score
  if (exists("loo_score_all")) {
    primary_de_genes <- dplyr::left_join(
      primary_de_genes,
      loo_score_all |> dplyr::select(Gene, LOO_Score, LOO_Robust),
      by = "Gene"
    )
  } else {
    primary_de_genes$LOO_Score  <- NA_integer_
    primary_de_genes$LOO_Robust <- NA
  }
  
  ## Jointure Cross-N signal
  primary_de_genes <- dplyr::left_join(
    primary_de_genes,
    crossN_df |> dplyr::select(Gene, CrossN_Signal),
    by = "Gene"
  ) |>
    dplyr::mutate(
      LOO_Robust     = dplyr::coalesce(LOO_Robust, FALSE),
      CrossN_Signal  = dplyr::coalesce(CrossN_Signal, FALSE),
      Tier = dplyr::case_when(
        LOO_Robust & CrossN_Signal ~ "Tier1",
        LOO_Robust | CrossN_Signal ~ "Tier2",
        TRUE                       ~ "Tier3"
      )
    ) |>
    dplyr::arrange(Tier, padj)
  
  write.csv(primary_de_genes,
            file.path(robust_dir, "DEG_Tier_classification.csv"),
            row.names = FALSE)
  
  ## Résumé Tier
  tier_summary <- primary_de_genes |>
    dplyr::count(Tier, Status, name = "n_genes") |>
    dplyr::arrange(Tier, Status)
  
  cat("\n  Tier summary:\n"); print(tier_summary); cat("\n")
  write.csv(tier_summary,
            file.path(robust_dir, "Tier_summary.csv"),
            row.names = FALSE)
  
  ## Barplot Tier
  p_tier <- ggplot2::ggplot(primary_de_genes,
                            ggplot2::aes(x = Tier, fill = Tier)) +
    ggplot2::geom_bar(colour = "white", linewidth = 0.4) +
    ggplot2::facet_wrap(~ Status, scales = "free_y") +
    ggplot2::geom_text(stat = "count",
                       ggplot2::aes(label = ggplot2::after_stat(count)),
                       vjust = -0.4, fontface = "bold", size = 4) +
    ggplot2::scale_fill_manual(values = TIER_COLORS[c("Tier1","Tier2","Tier3")]) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(
      title    = "Tier classification des DEGs  [A4]",
      subtitle = paste0(
        "Tier1 = LOO robuste + Cross-N | Tier2 = l'un des deux | Tier3 = Protection_mean seul"),
      x = NULL, y = "Nombre de gènes", fill = NULL
    ) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "none")
  
  save_figure(p_tier, file.path(robust_dir, "Tier_Barplot"), width = 8, height = 5)
  
  message("✅ Tier classification: ",
          sum(primary_de_genes$Tier == "Tier1"), " Tier1 | ",
          sum(primary_de_genes$Tier == "Tier2"), " Tier2 | ",
          sum(primary_de_genes$Tier == "Tier3"), " Tier3")
}


# ============================================================
# SECTION 10 — DEG SUMMARY BARPLOTS
# ============================================================
cat("---- SUMMARY PLOTS ----\n")

summary_long <- summary_df |>
  tidyr::pivot_longer(cols = c("N_up", "N_down"),
                      names_to = "Direction", values_to = "n") |>
  dplyr::mutate(
    Direction = dplyr::recode(Direction, "N_up" = "Up", "N_down" = "Down"),
    n_signed  = ifelse(Direction == "Down", -n, n),
    Contrast  = forcats::fct_reorder(as.character(Contrast), N_total, .fun = max)
  )

max_n_abs <- max(abs(summary_long$n_signed), na.rm = TRUE)
sub_str   <- paste0("|log2FC| ≥ ", lfc_threshold, "  |  FDR < ", alpha,
                    "  |  edgeR QL F-test")
cap_str   <- paste0("n = ", nrow(dge), " genes  |  filterByExpr  |  TMM  |  FDR BH")

summary_totals <- summary_df |>
  dplyr::mutate(Contrast = factor(as.character(Contrast),
                                  levels = levels(summary_long$Contrast)))

p_div <- ggplot2::ggplot(summary_long,
                         ggplot2::aes(x = Contrast, y = n_signed, fill = Direction)) +
  ggplot2::geom_col(width = 0.72, colour = "white", linewidth = 0.25) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.55, colour = "black") +
  ggplot2::scale_fill_manual(values = DEG_COLORS[c("Up", "Down")]) +
  ggplot2::scale_y_continuous(
    labels = function(x) formatC(abs(x), format = "d", big.mark = ","),
    expand = ggplot2::expansion(mult = c(0.14, 0.14))) +
  ggplot2::labs(title = "DEGs par contraste", subtitle = sub_str,
                x = NULL, y = "Nombre de DEGs", fill = NULL, caption = cap_str) +
  ggplot2::coord_flip() +
  theme_pub(base_size = 11) +
  ggplot2::theme(legend.position = "bottom",
                 axis.text.y = ggplot2::element_text(size = 8.5))
save_figure(p_div, file.path(summary_dir, "DEG_Barplot_Divergent"),
            width = 10, height = max(5, nrow(summary_df) * 0.55 + 2))

message("✅ Summary barplots saved")


# ============================================================
# SECTION 11 — [A6] VOLCANO PLOTS — ANNOTÉS PAR TIER
#
##  Volcano principal sur Protection_mean, coloré par Tier
##  Volcano par contraste pairwise (standard)
# ============================================================
cat("---- VOLCANO PLOTS ----\n")

pval_safe <- all_res_df$padj[all_res_df$padj > 0 & !is.na(all_res_df$padj)]
lfc_lim   <- max(abs(all_res_df$log2FoldChange), na.rm = TRUE) * 1.05
pval_lim  <- min(max(-log10(pval_safe), na.rm = TRUE) * 1.05, 50)

## ── [A6] Volcano Tier annoté (Protection_mean seulement) ─────────────
if (!is.null(res_mean) && exists("primary_de_genes")) {
  res_tier <- dplyr::left_join(
    res_mean |>
      dplyr::mutate(neg_log_padj = -log10(pmax(padj, 10^(-pval_lim)))),
    primary_de_genes |> dplyr::select(Gene, Tier),
    by = "Gene"
  ) |>
    dplyr::mutate(
      Tier_plot = dplyr::case_when(
        !is.na(Tier) & Tier == "Tier1" ~ "Tier1",
        !is.na(Tier) & Tier == "Tier2" ~ "Tier2",
        !is.na(Tier) & Tier == "Tier3" ~ "Tier3",
        TRUE ~ "NS"
      ),
      Tier_plot = factor(Tier_plot, levels = c("Tier1", "Tier2", "Tier3", "NS"))
    )
  
  top_t1 <- res_tier |>
    dplyr::filter(Tier_plot == "Tier1") |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = 15)
  
  p_vol_tier <- ggplot2::ggplot(
    res_tier,
    ggplot2::aes(x = log2FoldChange,
                 y = neg_log_padj,
                 colour = Tier_plot,
                 size   = Tier_plot)) +
    ggplot2::geom_point(data = dplyr::filter(res_tier, Tier_plot == "NS"),
                        colour = TIER_COLORS["NS"], alpha = 0.3, size = 0.7) +
    ggplot2::geom_point(data = dplyr::filter(res_tier, Tier_plot %in% c("Tier3","Tier2","Tier1")),
                        alpha = 0.85) +
    ggplot2::scale_colour_manual(
      values = TIER_COLORS,
      labels = c(
        "Tier1" = paste0("Tier1 — LOO+CrossN (n=",  sum(primary_de_genes$Tier=="Tier1"),")"),
        "Tier2" = paste0("Tier2 — LOO ou CrossN (n=",sum(primary_de_genes$Tier=="Tier2"),")"),
        "Tier3" = paste0("Tier3 — Mean seul (n=",    sum(primary_de_genes$Tier=="Tier3"),")"),
        "NS"    = "Not significant")) +
    ggplot2::scale_size_manual(values = c("Tier1" = 2.5, "Tier2" = 1.8, "Tier3" = 1.2, "NS" = 0.5),
                               guide = "none") +
    ggrepel::geom_text_repel(
      data = top_t1,
      ggplot2::aes(label = Gene),
      size = 2.7, max.overlaps = 20, segment.linewidth = 0.3,
      segment.colour = "grey55", colour = "black", fontface = "italic",
      box.padding = 0.35) +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::scale_x_continuous(limits = c(-lfc_lim, lfc_lim)) +
    ggplot2::scale_y_continuous(limits = c(0, pval_lim),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title    = "Volcano — Protection_mean  |  Tier Classification  [A6]",
      subtitle = "Tier1 = double validation LOO + Cross-N  |  Labels = Tier1 les plus significatifs",
      x        = expression(log[2]~"fold change  (Protected vs NonProtected)"),
      y        = expression(-log[10]~"(FDR)"),
      colour   = NULL,
      caption  = "Analyse primaire : contraste moyen Protection  |  FDR BH"
    ) +
    theme_pub(base_size = 11) +
    ggplot2::theme(legend.position = "top",
                   legend.text = ggplot2::element_text(size = 9))
  
  save_figure(p_vol_tier,
              file.path(volcano_dir, "Volcano_Protection_mean_TIER"),
              width = 8, height = 6.5)
}

## Volcanos standards (tous contrastes)
for (ct_name in names(all_results)) {
  res_df  <- all_results[[ct_name]]
  n_up    <- sum(res_df$Status == "Up")
  n_down  <- sum(res_df$Status == "Down")
  top_lab <- res_df |> dplyr::filter(Status != "NS") |>
    dplyr::arrange(padj) |> dplyr::slice_head(n = 15)
  
  p_v <- ggplot2::ggplot(
    res_df,
    ggplot2::aes(x = log2FoldChange,
                 y = -log10(pmax(padj, 10^(-pval_lim))),
                 colour = Status)) +
    ggplot2::geom_point(data = dplyr::filter(res_df, Status == "NS"),
                        colour = DEG_COLORS["NS"], alpha = 0.35, size = 0.8) +
    ggplot2::geom_point(data = dplyr::filter(res_df, Status != "NS"),
                        alpha = 0.80, size = 1.4) +
    ggplot2::scale_colour_manual(values = DEG_COLORS,
                                 labels = c("Up"   = paste0("Up (n=",   n_up,   ")"),
                                            "Down" = paste0("Down (n=", n_down, ")"),
                                            "NS"   = "NS")) +
    ggrepel::geom_text_repel(data = top_lab,
                             ggplot2::aes(label = Gene),
                             size = 2.6, max.overlaps = 18, segment.linewidth = 0.3,
                             segment.colour = "grey55", colour = "black", fontface = "italic") +
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    ggplot2::scale_x_continuous(limits = c(-lfc_lim, lfc_lim)) +
    ggplot2::scale_y_continuous(limits = c(0, pval_lim),
                                expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(title = ct_name,
                  x = expression(log[2]~"fold change"),
                  y = expression(-log[10]~"(FDR)"),
                  colour = NULL) +
    theme_pub(base_size = 11)
  
  save_figure(p_v,
              file.path(volcano_dir, paste0("Volcano_", ct_name)),
              width = 6.5, height = 5.5)
}
message("✅ Volcano plots saved")


# ============================================================
# SECTION 12 — MA PLOTS
# ============================================================
cat("---- MA PLOTS ----\n")

for (ct_name in names(all_results)) {
  res_df <- all_results[[ct_name]]
  p_ma <- ggplot2::ggplot(res_df, ggplot2::aes(x = aveLogCPM, y = log2FoldChange)) +
    ggplot2::geom_point(data = dplyr::filter(res_df, Status == "NS"),
                        colour = DEG_COLORS["NS"], alpha = 0.3, size = 0.6) +
    ggplot2::geom_point(data = dplyr::filter(res_df, Status == "Up"),
                        colour = DEG_COLORS["Up"], alpha = 0.7, size = 0.9) +
    ggplot2::geom_point(data = dplyr::filter(res_df, Status == "Down"),
                        colour = DEG_COLORS["Down"], alpha = 0.7, size = 0.9) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(-lfc_threshold, lfc_threshold),
                        linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    ggplot2::geom_smooth(colour = "#C0392B", se = TRUE, method = "loess",
                         linewidth = 0.8, fill = "#F1948A", alpha = 0.25) +
    ggplot2::labs(title    = paste0("MA plot — ", ct_name),
                  x        = "Average logCPM",
                  y        = expression(log[2]~"fold change")) +
    theme_pub(base_size = 11)
  save_figure(p_ma, file.path(ma_dir, paste0("MA_", ct_name)), width = 6.5, height = 5)
}
message("✅ MA plots saved")


# ============================================================
# SECTION 13 — HEATMAP UNIQUE
#
##  Heatmap en 2 versions :
##    A. Union de tous les DEGs (héritage v4)
##    B. Tier1 uniquement (nouveau) — signal haute confiance
# ============================================================
cat("---- HEATMAP ----\n")

col_labels    <- paste(sample_info$Genotype, sample_info$Replicate, sep = "_")
ann_col_hm    <- data.frame(Protection = sample_info$Protection,
                            Genotype   = sample_info$Genotype,
                            row.names  = col_labels)
ann_colors_hm <- list(Protection = PROTECTION_COLORS, Genotype = genotype_colors)

genes_de_union <- unique(all_res_df$Gene[all_res_df$Status != "NS"])
cat("DEGs union (≥ 1 contraste) :", length(genes_de_union), "\n")

make_heatmap <- function(gene_set, file_prefix, title_str) {
  if (length(gene_set) < 2) {
    message("⚠️  Heatmap ignorée — pas assez de gènes (n=", length(gene_set), ")")
    return(invisible(NULL))
  }
  mat_de   <- logcpm_mat[rownames(logcpm_mat) %in% gene_set, , drop = FALSE]
  mat_z    <- t(scale(t(mat_de)))
  mat_z    <- pmax(pmin(mat_z, 3), -3)
  colnames(mat_z) <- col_labels
  
  show_rn  <- nrow(mat_z) <= 80
  fsr      <- ifelse(nrow(mat_z) > 300, 4,
                     ifelse(nrow(mat_z) > 150, 5,
                            ifelse(nrow(mat_z) > 80,  6, 8)))
  N_CLUST  <- min(6, nrow(mat_z))
  
  hclust_rows <- hclust(dist(mat_z), method = "ward.D2")
  cluster_ids <- cutree(hclust_rows, k = N_CLUST)
  cluster_pal <- setNames(
    grDevices::colorRampPalette(
      c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628"))(N_CLUST),
    paste0("C", seq_len(N_CLUST))
  )
  ann_row_hm    <- data.frame(Cluster = factor(paste0("C", cluster_ids)),
                              row.names = names(cluster_ids))
  ann_colors_full <- c(ann_colors_hm, list(Cluster = cluster_pal))
  
  write.csv(data.frame(Gene = names(cluster_ids),
                       Cluster = paste0("C", cluster_ids)) |>
              dplyr::arrange(Cluster),
            file.path(heatmap_dir, paste0(file_prefix, "_clusters.csv")),
            row.names = FALSE)
  
  save_pheatmap(
    list(mat               = mat_z,
         color             = HEATMAP_DIVERGING,
         breaks            = seq(-3, 3, length.out = 101),
         annotation_col    = ann_col_hm,
         annotation_row    = ann_row_hm,
         annotation_colors = ann_colors_full,
         cluster_rows      = TRUE, cluster_cols = TRUE,
         clustering_method = "ward.D2",
         cutree_rows       = N_CLUST,
         show_rownames     = show_rn,
         show_colnames     = TRUE,
         fontsize_row      = fsr, fontsize_col = 9,
         fontsize          = 9, cellheight = NA,
         border_color      = NA,
         treeheight_row    = 20, treeheight_col = 20,
         main = title_str),
    file.path(heatmap_dir, file_prefix), w = 12, h = 14
  )
}

## A. Heatmap union
make_heatmap(
  genes_de_union,
  "Heatmap_All_DE_union",
  paste0("Tous DEGs — union >= 1 contraste  |  n = ", length(genes_de_union),
         "  |  Z-score logCPM  |  Ward.D2")
)

## B. Heatmap Tier1 seulement — signal haute confiance
if (exists("primary_de_genes")) {
  tier1_genes <- primary_de_genes$Gene[primary_de_genes$Tier == "Tier1"]
  if (length(tier1_genes) >= 2) {
    make_heatmap(
      tier1_genes,
      "Heatmap_Tier1_HighConfidence",
      paste0("Tier1 DEGs — LOO + Cross-N robustes  |  n = ", length(tier1_genes),
             "  |  Z-score logCPM  |  Ward.D2")
    )
    cat("✅ Heatmap Tier1:", length(tier1_genes), "gènes\n")
  } else {
    message("⚠️  Tier1 < 2 gènes — heatmap Tier1 ignorée")
  }
}

message("✅ Heatmaps saved")


# ============================================================
# SECTIONS 14–16 — VENN / UPSET / CONCORDANCE
# (héritées de v4, inchangées — voir code v4 original)
# ============================================================

## Pour garder ce script autonome et lisible, les sections Venn,
## UpSet et Concordance sont identiques à la v4 et peuvent être
## copiées directement depuis RNAseq_edgeR_v4.R (sections 14-16).
## Elles opèrent sur all_results et pairwise_names qui sont
## définis ci-dessus de façon identique.

message("ℹ️  Sections Venn/UpSet/Concordance : reporter depuis v4 (identiques)")


# ============================================================
# SECTION 17 — [A7] GO ENRICHMENT — PRIORITÉ TIER1
#
##  GO sur 4 niveaux (héritage v4) + niveau Tier1 ajouté :
##    0. Tier1 seul  ← NOUVEAU : gènes haute confiance uniquement
##    1. Par contraste individuel
##    2. Par _N : intersection (héritage)
##    3. Par _P : intersection (héritage)
##    4. Cross-N intersection (≡ ta stratégie initiale formalisée)
# ============================================================
cat("---- GO ENRICHMENT ----\n")

if (!file.exists(go_annot_file)) {
  warning("⚠️  GO annotation file not found: ", go_annot_file)
} else {
  
  go_raw        <- read.csv(go_annot_file, header = TRUE, stringsAsFactors = FALSE)
  colnames(go_raw) <- trimws(colnames(go_raw))
  gene_col <- intersect(c("Gene_ID", "gene", "Gene"), colnames(go_raw))[1]
  go_col   <- intersect(c("Annotation", "GO", "go", "GO_ID"), colnames(go_raw))[1]
  if (is.na(gene_col) || is.na(go_col))
    stop("Colonnes Gene/GO introuvables dans: ", go_annot_file)
  
  go_df <- go_raw[, c(gene_col, go_col)] |>
    dplyr::rename(gene_id = 1, go_term = 2) |>
    dplyr::mutate(gene_id = trimws(as.character(gene_id)),
                  go_term = trimws(as.character(go_term))) |>
    tidyr::separate_longer_delim(go_term, delim = stringr::regex("[,;\\s]+")) |>
    dplyr::filter(grepl("^GO:\\d+$", go_term)) |>
    dplyr::distinct()
  
  cat("Valid gene-GO pairs:", nrow(go_df), "\n")
  genes_tested  <- rownames(dge)
  geneID2GO     <- split(go_df$go_term, go_df$gene_id)
  gene_universe <- intersect(genes_tested, names(geneID2GO))
  cat("Universe:", length(gene_universe), "\n\n")
  
  if (length(gene_universe) == 0)
    stop("Aucun gène commun entre edgeR et annotation GO.")
  
  GO_P_CUTOFF  <- 0.05
  GO_MIN_GENES <- 4
  TOP_N_TERMS  <- 8
  ontologies   <- c("BP", "MF", "CC")
  
  run_topgo <- function(genes_oi, ont) {
    genes_oi <- intersect(genes_oi, gene_universe)
    if (length(genes_oi) < GO_MIN_GENES) return(NULL)
    geneList <- factor(as.integer(gene_universe %in% genes_oi), levels = c(0, 1))
    names(geneList) <- gene_universe
    GOdata <- tryCatch(
      new("topGOdata", ontology = ont, allGenes = geneList,
          annot = annFUN.gene2GO, gene2GO = geneID2GO),
      error = function(e) NULL)
    if (is.null(GOdata)) return(NULL)
    res_w01  <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    res_clas <- runTest(GOdata, algorithm = "classic",  statistic = "fisher")
    allRes   <- GenTable(GOdata, weight01 = res_w01, classic = res_clas,
                         topNodes = length(score(res_w01)), numChar = 120)
    clean_p  <- function(x) suppressWarnings(as.numeric(gsub("< ", "", x)))
    allRes$weight01 <- clean_p(allRes$weight01)
    allRes$classic  <- clean_p(allRes$classic)
    allRes$weight01[is.na(allRes$weight01) | allRes$weight01 == 0] <- 1e-300
    allRes$classic [is.na(allRes$classic)  | allRes$classic  == 0] <- 1e-300
    allRes$FDR <- p.adjust(allRes$weight01, method = "BH")
    allRes$DEG_Count <- vapply(allRes$GO.ID, function(gid)
      tryCatch(sum(genesInTerm(GOdata, gid)[[1]] %in% genes_oi),
               error = function(e) 0L), integer(1))
    allRes$DEG_Genes <- vapply(allRes$GO.ID, function(gid)
      tryCatch(paste(sort(intersect(genesInTerm(GOdata, gid)[[1]], genes_oi)),
                     collapse = ";"), error = function(e) ""), character(1))
    out <- allRes[allRes$weight01 < GO_P_CUTOFF & allRes$DEG_Count >= GO_MIN_GENES, ]
    if (nrow(out) == 0) return(NULL)
    out$negLogP <- -log10(out$weight01); out
  }
  
  ## [Bubble plot + Excel export — identiques à v4, reporter depuis v4]
  
  run_go_set <- function(gene_lists, set_id, set_label, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    go_collector <- list()
    for (dir_tag in names(gene_lists)) {
      genes_oi <- gene_lists[[dir_tag]]
      if (length(intersect(genes_oi, gene_universe)) < GO_MIN_GENES) next
      for (ont in ontologies) {
        res_go <- run_topgo(genes_oi, ont)
        if (!is.null(res_go) && nrow(res_go) > 0) {
          res_go$Ontology  <- ont
          res_go$Direction <- dir_tag
          res_go$GeneSet   <- set_id
          write.csv(res_go,
                    file.path(out_dir, paste0("GO_", set_id, "_", dir_tag, "_", ont, ".csv")),
                    row.names = FALSE)
          go_collector[[paste(dir_tag, ont, sep = "_")]] <- res_go
          cat("    ", dir_tag, ont, ":", nrow(res_go), "terms\n")
        }
      }
    }
    if (length(go_collector) > 0) {
      all_go_df <- dplyr::bind_rows(go_collector) |> dplyr::arrange(weight01)
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "All_Terms")
      openxlsx::writeData(wb, "All_Terms", all_go_df)
      openxlsx::saveWorkbook(wb, file.path(out_dir, paste0("GO_", set_id, ".xlsx")),
                             overwrite = TRUE)
      message("✅ GO saved: ", set_id)
    } else {
      cat("  ⚠️  Aucun terme significatif:", set_id, "\n")
    }
    invisible(go_collector)
  }
  
  ## ── [A7] GO Tier1 ────────────────────────────────────────────────
  cat("\nGO Tier1 (signal haute confiance)...\n")
  if (exists("primary_de_genes")) {
    tier1_up   <- primary_de_genes$Gene[primary_de_genes$Tier == "Tier1" &
                                          primary_de_genes$Status == "Up"]
    tier1_down <- primary_de_genes$Gene[primary_de_genes$Tier == "Tier1" &
                                          primary_de_genes$Status == "Down"]
    run_go_set(
      list(Up = tier1_up, Down = tier1_down),
      "Tier1_HighConfidence",
      "Tier1 — LOO + Cross-N robuste",
      file.path(go_dir, "00_Tier1_HighConfidence")
    )
  }
  
  ## ── GO Cross-N ────────────────────────────────────────────────────
  cat("\nGO Cross-N intersection (ta stratégie formalisée [A3])...\n")
  run_go_set(
    list(Up   = crossN_up,
         Down = crossN_down),
    "CrossN_intersection",
    "Cross-N — signal Protection indépendant de la basale",
    file.path(go_dir, "01_CrossN_Intersection")
  )
  
  ## ── GO Protection_mean (analyse primaire) ─────────────────────────
  cat("\nGO Protection_mean (analyse primaire)...\n")
  if (!is.null(res_mean)) {
    run_go_set(
      list(Up   = res_mean$Gene[res_mean$Status == "Up"],
           Down = res_mean$Gene[res_mean$Status == "Down"]),
      "Protection_mean",
      "Protection_mean — contraste primaire",
      file.path(go_dir, "02_Protection_mean")
    )
  }
  
  message("✅ GO analysis done — résultats dans: ", go_dir)
}


# ============================================================
# SECTION 18 — SESSION SUMMARY + MATÉRIELS & MÉTHODES
# ============================================================

cat("\n====================================================\n")
cat("  edgeR v4-ALT — Analysis complete\n")
cat("====================================================\n")
cat("  Stratégie     : 3 niveaux de validation\n")
cat("  [PRIMARY]     : Protection_mean (contraste moyen)\n")
cat("  [VALIDATION1] : LOO Leave-One-Out sur _P\n")
cat("  [VALIDATION2] : Cross-N intersection\n")
cat("  [TIER]        : Tier1/2/3 par combinaison des validations\n")
cat("  Contrastes    :", length(all_results), "(primary + pairwise)\n")
cat("  LOO           :", length(loo_results), "contrastes\n")
cat("  Cross-N Up    :", length(crossN_up), "gènes\n")
cat("  Cross-N Down  :", length(crossN_down), "gènes\n")
if (exists("primary_de_genes")) {
  cat("  Tier1 gènes   :", sum(primary_de_genes$Tier == "Tier1"), "\n")
  cat("  Tier2 gènes   :", sum(primary_de_genes$Tier == "Tier2"), "\n")
  cat("  Tier3 gènes   :", sum(primary_de_genes$Tier == "Tier3"), "\n")
}
cat("  Seuil         : |log2FC| >=", lfc_threshold, "& FDR <", alpha, "\n")
cat("  Output        :", main_dir, "\n")
cat("====================================================\n\n")

mm_text <- c(
  "================================================",
  "MATERIALS AND METHODS — edgeR v4-ALT",
  "================================================",
  "",
  "Experimental design and limitation:",
  "This RNAseq experiment includes only treated plant samples. Direct pairwise",
  "comparison between protected (_P) and non-protected (_N) genotypes without",
  "untreated controls risks confounding the treatment response with constitutive",
  "genetic differences. The following three-level validation strategy was adopted",
  "to isolate a protection-related transcriptomic signal.",
  "",
  "Data filtering and normalization (identical to v4):",
  "Raw counts were filtered using edgeR::filterByExpr() (Chen et al. 2016) and",
  "normalized by TMM (Robinson & Oshlack 2010). A ~ 0 + Genotype + Replicate design",
  paste0("was used, retaining ", nrow(dge), " genes."),
  "",
  "Statistical model and PRIMARY contrast:",
  "The primary analysis used the Protection_mean contrast:",
  paste0("  ", contrast_mean_str),
  "This contrast is formally equivalent to the Protection coefficient in a",
  "~ Protection + Genotype + Replicate model: the genotype term absorbs genetic",
  "background variation, and the contrast directly tests the average effect of",
  "the protection treatment across all genotypes. Differential expression thresholds:",
  paste0("|log2FC| >= ", lfc_threshold, " and FDR < ", alpha, " (Benjamini-Hochberg)."),
  "",
  paste0("VALIDATION 1 — Leave-One-Out (LOO) on _P genotypes (", length(loo_results), " contrasts):"),
  "To verify that the primary signal is not driven by a single atypical protected",
  "genotype, the Protection_mean contrast was recomputed ", length(loo_results), " times,",
  "each time excluding one _P genotype. A LOO robustness score (0 to n_P) was assigned",
  "to each DEG reflecting how many LOO subsets it remained differentially expressed",
  "(same direction). Genes with maximum LOO score are robust to the removal of any",
  "single _P genotype.",
  "",
  paste0("VALIDATION 2 — Cross-N intersection (", length(conditions_N), " non-protected baselines):"),
  "For each _N genotype, the intersection of DEGs from all _P vs that _N was computed",
  "(Core_N). The cross-N intersection (∩ Core_N across all _N genotypes) identifies",
  "genes that are differentially expressed regardless of which non-protected reference",
  paste0("genotype is used (Cross-N Up: ", length(crossN_up),
         " genes; Cross-N Down: ", length(crossN_down), " genes)."),
  "This formalizes and validates the rationale of comparing each _P genotype against",
  "multiple baselines.",
  "",
  "GENE TIER CLASSIFICATION:",
  "DEGs from the primary contrast were classified into three confidence tiers:",
  "  Tier 1 (high confidence): significant in Protection_mean AND LOO-robust (max score)",
  "          AND present in Cross-N intersection.",
  "  Tier 2 (medium confidence): significant in Protection_mean AND either LOO-robust",
  "          OR present in Cross-N.",
  "  Tier 3 (limited confidence): significant in Protection_mean only; neither",
  "          validation confirms.",
  if (exists("primary_de_genes")) {
    paste0("  Results: Tier1=", sum(primary_de_genes$Tier=="Tier1"),
           " | Tier2=", sum(primary_de_genes$Tier=="Tier2"),
           " | Tier3=", sum(primary_de_genes$Tier=="Tier3"), " genes.")
  } else "",
  "",
  "GO enrichment was performed with topGO weight01 Fisher (Alexa & Rahnenfuhrer 2023)",
  "on Tier1 genes, Cross-N intersection genes, and Protection_mean DEGs.",
  "",
  "References:",
  "  Robinson MD & Oshlack A (2010) Genome Biology 11:R25.",
  "  Robinson MD et al. (2010) Bioinformatics 26:139-140.",
  "  McCarthy DJ et al. (2012) Nucleic Acids Res 40:4288-4297.",
  "  Chen Y et al. (2016) F1000Research 5:1408.",
  "  Conway JR et al. (2017) Bioinformatics 33:2750-2752.",
  "  Alexa A & Rahnenfuhrer J (2023) topGO R package v2.x.",
  "  R Core Team (2024) R v4.x. R Foundation.",
  "================================================"
)
writeLines(mm_text, file.path(main_dir, "Materials_and_Methods_v4ALT.txt"))
writeLines(c("R session information", "=====================",
             capture.output(sessionInfo())),
           file.path(main_dir, "session_info.txt"))

message("✅ edgeR v4-ALT COMPLETE")

