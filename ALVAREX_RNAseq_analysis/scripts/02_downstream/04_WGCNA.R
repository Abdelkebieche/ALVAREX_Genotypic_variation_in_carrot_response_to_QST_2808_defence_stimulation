############################################################
## WGCNA RNAseq — SCRIPT COMPLET FINAL v3
##
## Dataset : plusieurs génotypes, condition Protection (P vs N)
## Objectif : modules co-exprimés liés à la protection des plantes
##
## ──────────────────────────────────────────────────────────
## CORRECTIONS vs v2 (numérotation continue) :
##
##  CRITIQUES [reviewer destruction garantie sans ça] :
##   [C9]  removeBatchEffect (limma) sur la matrice VST avant datExpr
##         → Replicate retiré explicitement de la matrice d'expression
##   [C10] Filtrage MAD top 50 % avant blockwiseModules
##         → réduit bruit, améliore topologie scale-free
##   [C11] scaleFreePlot + vérification distribution degrés de connectivité
##         → validation formelle du fit scale-free
##   [C12] Analyse de sensibilité paramètres blockwiseModules
##         → re-coupe du dendrogramme (TOM déjà calculé = rapide)
##   [C13] TOM plot — heatmap intramodulaire triée par module
##         → figure quasi-obligatoire dans les papiers WGCNA
##   [C14] Permutation p-values (n=1000) pour ME-trait
##         → remplace corPvalueStudent (hypothèse indépendance violée)
##   [C15] lmerTest présenté comme test PRIMAIRE (bicor + permut = secondaire)
##   [C16] GO_FDR_MODULE_CUTOFF : 0.50 → 0.20
##
##  IMPORTANTS :
##   [C17] Eigengene network + méta-modules (plotEigengeneNetworks)
##
##  SOUHAITABLES :
##   [C18] modulePreservation Protected vs NonProtected (Zsummary, MedianRank)
##   [C19] Connectivité différentielle par gène (Protected vs NonProtected)
##   [C20] Statistiques réseau via igraph (betweenness, clustering coefficient)
##   [C21] Bootstrap stabilité hub genes (80 % samples, B itérations)
##
## CORRECTIONS MAINTENUES depuis v2 :
##   [C2] FDR_CUTOFF_MAIN = 0.10
##   [C3] Numérotation linéaire
##   [C4] ME names stables
##   [C5] tidyr chargé explicitement
##   [C6] saveTOMs = TRUE
##   [C7] bicor pour sample clustering
##   [C8] Replicate → factor avant DESeqDataSetFromMatrix
############################################################

##############################
## PARTIE 0) Paramètres
##############################
SEED <- 123
INSTALL_MISSING_PKGS <- FALSE

## USER CONFIGURATION
## Replace only these paths for your computing environment.
counts_file   <- "/path/to/input/Count.csv"
meta_file     <- "/path/to/input/MetaData.csv"
go_annot_file <- "/path/to/annotation/gene_to_GO.csv"

base_out  <- "/path/to/results/WGCNA"
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir   <- file.path(base_out, paste0("WGCNA_Protection_ALL_", timestamp))

dir_inputs <- file.path(out_dir, "00_inputs")
dir_norm   <- file.path(out_dir, "01_normalization")
dir_qc     <- file.path(out_dir, "02_QC_samples")
dir_plots  <- file.path(dir_qc, "plots")
dir_tables <- file.path(dir_qc, "tables")

for (d in c(out_dir, dir_inputs, dir_norm, dir_qc, dir_plots, dir_tables))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

set.seed(SEED)
options(stringsAsFactors = FALSE)

## ── [C9-C21] Flags et paramètres des nouvelles analyses ──────────────

## [C14] Permutations pour p-values ME-trait
N_PERM_MT <- 1000

## [C10] Seuil MAD pour filtrage gènes avant WGCNA
MAD_QUANTILE_CUTOFF <- 0.60   # top 50 %

## [C12] Sensibilité blockwiseModules (re-cut du dendrogramme — rapide)
RUN_SENSITIVITY     <- TRUE
SENSITIVITY_DEEPSPLIT    <- c(2L, 3L, 4L)
SENSITIVITY_MINMODSIZE   <- c(30L, 50L, 100L)
SENSITIVITY_MERGEHEIGHT  <- c(0.25, 0.30)

## [C18] modulePreservation
RUN_PRESERVATION    <- TRUE
N_PERM_PRESERVATION <- 200   # ≥ 200 pour Zsummary stable (Langfelder 2011)

## [C21] Bootstrap hub genes
RUN_BOOTSTRAP   <- TRUE
N_BOOT          <- 50        # 50 itérations (ajuster selon ressources)
BOOT_FRAC       <- 0.80      # fraction d'échantillons par boot

##############################
## PARTIE -1) Packages
##############################
load_pkgs <- function(pkgs, install_missing = FALSE, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!install_missing) stop("Package manquant : ", p)
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(p, ask = FALSE, update = FALSE)
      } else { install.packages(p, dependencies = TRUE) }
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

load_pkgs(c("dplyr", "tibble", "tidyr", "ggplot2", "pheatmap",
            "RColorBrewer", "svglite", "scales", "forcats", "stringr",
            "ggrepel", "readr", "grid", "igraph", "dynamicTreeCut"),
          install_missing = INSTALL_MISSING_PKGS, bioc = FALSE)
load_pkgs(c("DESeq2", "SummarizedExperiment", "limma"),    # [C9] limma ajouté
          install_missing = INSTALL_MISSING_PKGS, bioc = TRUE)

############################################################
## MASTER COLOR SYSTEM
############################################################
GENOTYPE_PALETTE <- c(
  "Amsterdam_N"      = "#0A9E6E", "NantaiseInbred_N" = "#1B4FBF",
  "Orleans_N"        = "#007B8A", "Dijon_N"          = "#3949AB",
  "Genevieve_N"      = "#006B3C", "Deep_Purple_P"    = "#D62728",
  "Presto_P"         = "#E8820C", "Neva_P"           = "#9B2ECA",
  "Robila_P"         = "#C9A800", "Genevieve_P"      = "#C9A800",
  "Orleans_P"        = "#5E9E1A", "Dijon_P"          = "#C2185B",
  "Oxhella_P"        = "#8B0000"
)

PROTECTION_COLORS <- c("NonProtected" = "#2471A3", "Protected" = "#C0392B")

HEATMAP_DIVERGING <- grDevices::colorRampPalette(
  c("#2166AC","#4393C3","#92C5DE","#F7F7F7","#F4A582","#D6604D","#B2182B"))(100)

HEATMAP_DISTANCE <- grDevices::colorRampPalette(
  c("#1A1A2E","#16213E","#0F3460","#A8C7FA","#E8F4F8","#FFFFFF"))(255)

resolve_genotype_colors <- function(genotype_levels) {
  known   <- genotype_levels[genotype_levels %in% names(GENOTYPE_PALETTE)]
  unknown <- genotype_levels[!genotype_levels %in% names(GENOTYPE_PALETTE)]
  if (length(unknown) > 0) {
    hues <- seq(15, 375, length.out = length(unknown) + 1)[seq_len(length(unknown))]
    fallback <- setNames(grDevices::hsv(hues / 360, s = 0.82, v = 0.75), unknown)
    warning(length(unknown), " genotype(s) absent de GENOTYPE_PALETTE → couleur auto.", call. = FALSE)
    return(c(GENOTYPE_PALETTE[known], fallback)[genotype_levels])
  }
  GENOTYPE_PALETTE[genotype_levels]
}

message("✅ Master color system defini")

############################################################
## PUBLICATION THEME
############################################################
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
      plot.subtitle     = ggplot2::element_text(size = base_size + 1),
      plot.caption      = ggplot2::element_text(size = base_size - 1, colour = "grey50"),
      legend.title      = ggplot2::element_text(size = base_size + 1, face = "bold"),
      legend.text       = ggplot2::element_text(size = base_size),
      legend.key.size   = ggplot2::unit(5, "mm"),
      strip.text        = ggplot2::element_text(face = "bold", size = base_size + 1),
      strip.background  = ggplot2::element_rect(fill = "grey94", colour = "black", linewidth = 0.4)
    )
}

save_figure <- function(p, base_path, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(base_path, ".pdf"), p, width = width, height = height, bg = "white")
  ggplot2::ggsave(paste0(base_path, ".png"), p, width = width, height = height, dpi = dpi, bg = "white")
  svglite::svglite(paste0(base_path, ".svg"), width = width, height = height)
  print(p); grDevices::dev.off()
  message("   Saved: ", basename(base_path), "  [PDF / PNG / SVG]")
}

save_pheatmap_wgcna <- function(mat, dist_obj = NULL, ann_col = NULL, ann_row = NULL,
                                ann_colors = list(), color_pal = HEATMAP_DISTANCE,
                                title_str, full_path, w = 12, h = 10,
                                fontsize = 8, fontsize_row = 7, fontsize_col = 7) {
  args <- list(mat = mat, color = color_pal, annotation_col = ann_col,
               annotation_row = ann_row, annotation_colors = ann_colors,
               cluster_rows = TRUE, cluster_cols = TRUE,
               clustering_method = "ward.D2",
               fontsize = fontsize, fontsize_row = fontsize_row, fontsize_col = fontsize_col,
               border_color = NA, treeheight_row = 20, treeheight_col = 20, main = title_str)
  if (!is.null(dist_obj)) {
    args$clustering_distance_rows <- dist_obj
    args$clustering_distance_cols <- dist_obj
  }
  do.call(pheatmap::pheatmap, c(args, list(filename = paste0(full_path, ".pdf"), width = w, height = h)))
  do.call(pheatmap::pheatmap, c(args, list(filename = paste0(full_path, ".png"), width = w, height = h)))
  ph_obj <- do.call(pheatmap::pheatmap, c(args, list(silent = TRUE)))
  svglite::svglite(paste0(full_path, ".svg"), width = w, height = h)
  grid::grid.newpage(); grid::grid.draw(ph_obj$gtable); grDevices::dev.off()
  message("✅ Heatmap saved: ", basename(full_path), "  [PDF / PNG / SVG]")
}

starify <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                         ifelse(p < 0.05, "*", ifelse(p < 0.10, ".", "")))))
}

message("✅ Publication theme + helpers definis")

##############################
## PARTIE 0.5) Fonctions utilitaires
##############################
normalize_sample_name <- function(x) {
  x <- as.character(x); x <- trimws(x); x <- gsub("^X", "", x)
  x <- gsub("\"", "", x); x <- gsub("\\s+", "", x); x <- gsub("\r|\n", "", x)
  x <- gsub("_F-([0-9]{4}clean_mapping)$", "_F_\\1", x); x
}

fix_encoding <- function(x) {
  x2 <- iconv(as.character(x), from = "", to = "UTF-8", sub = "")
  x2[is.na(x2)] <- x[is.na(x2)]; x2
}

extract_protection <- function(genotype) {
  genotype <- trimws(as.character(genotype))
  dplyr::case_when(
    grepl("_P$", genotype, ignore.case = TRUE) ~ "Protected",
    grepl("_N$", genotype, ignore.case = TRUE) ~ "NonProtected",
    TRUE ~ NA_character_
  )
}

extract_genotype_base <- function(genotype) {
  gsub("_(P|N)$", "", trimws(as.character(genotype)), ignore.case = TRUE)
}

make_safe_factor <- function(x, na_label = "Unknown") {
  x <- as.character(x); x[is.na(x) | x == ""] <- na_label; factor(x)
}

##############################
## PARTIE 1) Import COUNTS
##############################
cat("---- IMPORT COUNTS ----\n")
stopifnot(file.exists(counts_file))
count_raw <- read.csv(counts_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
cat("count_raw dim:", paste(dim(count_raw), collapse = " x "), "\n")
if (!("Geneid" %in% colnames(count_raw))) stop("La colonne 'Geneid' est introuvable dans Count.csv")
count_mat       <- count_raw |> tibble::column_to_rownames("Geneid") |> as.matrix()
mode(count_mat) <- "numeric"
colnames(count_mat) <- normalize_sample_name(colnames(count_mat))
cat("count_mat dim (genes x samples):", paste(dim(count_mat), collapse = " x "), "\n\n")

##############################
## PARTIE 2) Import METADATA
##############################
cat("---- IMPORT METADATA ----\n")
stopifnot(file.exists(meta_file))
meta_raw <- read.csv(meta_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c("Sample", "Genotype", "Replicate")
missing_cols  <- setdiff(required_cols, colnames(meta_raw))
if (length(missing_cols) > 0) stop("Colonnes manquantes: ", paste(missing_cols, collapse = ", "))

meta <- meta_raw |>
  dplyr::mutate(
    Sample       = normalize_sample_name(Sample),
    Genotype     = fix_encoding(trimws(Genotype)),
    Replicate    = factor(as.integer(Replicate)),  # [C8]
    Protection   = extract_protection(Genotype),
    GenotypeBase = extract_genotype_base(Genotype)
  )
cat("Distribution Protection:\n"); print(table(meta$Protection, useNA = "ifany"))
cat("\nDistribution Genotype:\n"); print(table(meta$Genotype))

##############################
## PARTIE 3) Matching
##############################
cat("---- MATCHING ----\n")
if (anyDuplicated(meta$Sample) > 0)
  stop("Samples dupliqués: ", paste(meta$Sample[duplicated(meta$Sample)], collapse = ", "))
common_samples <- intersect(colnames(count_mat), meta$Sample)
cat("Nb samples communs:", length(common_samples), "\n")
if (length(common_samples) == 0) stop("Aucun sample commun.")

##############################
## PARTIE 4) Alignement
##############################
meta_aligned  <- meta |>
  dplyr::filter(Sample %in% common_samples) |>
  dplyr::arrange(match(Sample, colnames(count_mat)))
count_aligned <- count_mat[, meta_aligned$Sample, drop = FALSE]
stopifnot(all(colnames(count_aligned) == meta_aligned$Sample))
genotype_levels <- sort(unique(as.character(meta_aligned$Genotype)))
genotype_colors <- resolve_genotype_colors(genotype_levels)
write.csv(meta_aligned,  file.path(dir_inputs, "metadata_aligned.csv"),              row.names = FALSE)
saveRDS(count_aligned,   file.path(dir_inputs, "count_aligned_genes_x_samples.rds"))
saveRDS(meta_aligned,    file.path(dir_inputs, "meta_aligned.rds"))

##############################
## PARTIE 5) Préfiltrage gènes
##############################
cat("---- PREFILTRAGE GENES ----\n")
MIN_TOTAL_COUNTS <- 50
keep_genes       <- rowSums(count_aligned, na.rm = TRUE) >= MIN_TOTAL_COUNTS
count_filtered   <- count_aligned[keep_genes, , drop = FALSE]
cat("Genes après filtrage:", nrow(count_filtered), "\n\n")

write.csv(data.frame(Metric = c("Genes_before","Genes_after","Genes_removed"),
                     Value  = c(nrow(count_aligned), nrow(count_filtered),
                                nrow(count_aligned) - nrow(count_filtered))),
          file.path(dir_norm, "gene_filtering_summary.csv"), row.names = FALSE)

##############################
## PARTIE 6) DESeq2 + VST + [C9] removeBatchEffect
##############################
cat("---- DESEQ2 + VST + removeBatchEffect ----\n")

meta_aligned <- meta_aligned |>
  dplyr::mutate(Protection   = make_safe_factor(Protection),
                Genotype     = make_safe_factor(Genotype),
                GenotypeBase = make_safe_factor(GenotypeBase))

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_filtered),
  colData   = meta_aligned,
  design    = ~ Protection + Replicate
)

MIN_COUNT <- 40; MIN_SAMPLES <- 3
keep <- rowSums(DESeq2::counts(dds) >= MIN_COUNT) >= MIN_SAMPLES
dds  <- dds[keep, ]
cat("Genes après filtre DESeq2:", nrow(dds), "\n")
dds  <- DESeq2::estimateSizeFactors(dds)
cat("Size factors summary:\n"); print(summary(DESeq2::sizeFactors(dds)))

vsd      <- DESeq2::vst(dds, blind = FALSE)
expr_mat <- SummarizedExperiment::assay(vsd)          # genes × samples, NON corrigée
cat("expr_mat dim:", paste(dim(expr_mat), collapse = " x "), "\n")

## ── [C9] Correction batch Replicate EXPLICITE via limma ──────────────
## Langfelder & Horvath recommandent une matrice propre pour WGCNA ;
## removeBatchEffect retire la variance du Replicate tout en préservant
## l'effet Protection (modèle de design fourni via `design` argument).
## Référence : Ritchie et al. 2015, Nucleic Acids Research

design_prot     <- stats::model.matrix(~ Protection, data = meta_aligned)
expr_mat_bc     <- limma::removeBatchEffect(
  expr_mat,
  batch  = meta_aligned$Replicate,
  design = design_prot
)
cat("✅ removeBatchEffect appliqué (batch = Replicate)\n\n")

## Sauvegarder les deux matrices
saveRDS(dds,        file.path(dir_norm, "dds_object.rds"))
saveRDS(vsd,        file.path(dir_norm, "vsd_object.rds"))
saveRDS(expr_mat,   file.path(dir_norm, "expr_mat_vst_RAW_genes_x_samples.rds"))
saveRDS(expr_mat_bc,file.path(dir_norm, "expr_mat_vst_BATCHCORRECTED_genes_x_samples.rds"))

##############################
## PARTIE 7) PCA — avant et après correction batch
##############################
cat("---- PCA SAMPLES ----\n")

## PCA avec matrice brute
pcaData_raw    <- DESeq2::plotPCA(vsd, intgroup = c("Genotype","Replicate"), returnData = TRUE)
percentVar_raw <- round(100 * attr(pcaData_raw, "percentVar"))
pcaData_raw$Protection <- meta_aligned$Protection[match(rownames(pcaData_raw), meta_aligned$Sample)]
pcaData_raw$Source     <- "Before batch correction"

## PCA avec matrice corrigée (recalculée manuellement)
pca_bc       <- stats::prcomp(t(expr_mat_bc), scale. = FALSE)
pca_bc_pct   <- round(100 * pca_bc$sdev^2 / sum(pca_bc$sdev^2))
pcaData_bc   <- data.frame(PC1 = pca_bc$x[, 1], PC2 = pca_bc$x[, 2],
                           Sample   = meta_aligned$Sample,
                           Genotype = meta_aligned$Genotype,
                           Replicate  = meta_aligned$Replicate,
                           Protection = meta_aligned$Protection,
                           Source     = "After batch correction",
                           row.names = meta_aligned$Sample)

rep_shapes <- setNames(
  c(16, 17, 15, 18, 8, 3)[seq_along(levels(factor(meta_aligned$Replicate)))],
  levels(factor(meta_aligned$Replicate))
)

make_pca_plot <- function(df, xlab, ylab, subtitle_str) {
  df$Genotype   <- factor(as.character(df$Genotype), levels = genotype_levels)
  df$Replicate  <- factor(df$Replicate)
  ggplot2::ggplot(df, ggplot2::aes(x = PC1, y = PC2, colour = Genotype, shape = Replicate)) +
    ggplot2::geom_point(size = 4, stroke = 0.5, alpha = 0.92) +
    ggplot2::scale_colour_manual(values = genotype_colors,
                                 guide  = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::scale_shape_manual(values = rep_shapes,
                                guide  = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
    ggplot2::labs(title = "PCA — WGCNA VST", subtitle = subtitle_str,
                  caption = paste0("n = ", nrow(df), " samples")) +
    theme_pub(base_size = 11) + ggplot2::theme(legend.position = "right")
}

p_pca_raw <- make_pca_plot(pcaData_raw,
                           paste0("PC1  (", percentVar_raw[1], "% variance)"),
                           paste0("PC2  (", percentVar_raw[2], "% variance)"),
                           "VST blind=FALSE | design ~ Protection + Replicate | BEFORE batch correction")

p_pca_bc  <- make_pca_plot(pcaData_bc,
                           paste0("PC1  (", pca_bc_pct[1], "% variance)"),
                           paste0("PC2  (", pca_bc_pct[2], "% variance)"),
                           "VST + removeBatchEffect(Replicate) | AFTER batch correction | matrice WGCNA")

save_figure(p_pca_raw, file.path(dir_plots, "QC_PCA_raw"),           width = 8, height = 5.5)
save_figure(p_pca_bc,  file.path(dir_plots, "QC_PCA_batchcorrected"), width = 8, height = 5.5)

write.csv(pcaData_bc, file.path(dir_tables, "PCA_coordinates_BATCHCORRECTED.csv"), row.names = FALSE)

##############################
## PARTIE 8) Sample clustering (bicor) [C7]
##############################
cat("---- SAMPLE CLUSTERING ----\n")
load_pkgs("WGCNA", install_missing = INSTALL_MISSING_PKGS, bioc = FALSE)

sample_cor          <- WGCNA::bicor(expr_mat_bc)        # sur matrice corrigée [C9]
sample_dist         <- as.dist(1 - sample_cor)
sample_tree         <- hclust(sample_dist, method = "average")
sample_labels_short <- paste0(meta_aligned$GenotypeBase, "_R", meta_aligned$Replicate)

for (ext in c("pdf","png","svg")) {
  fpath <- file.path(dir_plots, paste0("Sample_clustering_readable.", ext))
  if      (ext == "pdf") grDevices::pdf(fpath, width = 14, height = 7, useDingbats = FALSE)
  else if (ext == "png") grDevices::png(fpath, width = 14*220, height = 7*220, res = 220)
  else                   svglite::svglite(fpath, width = 14, height = 7)
  par(mar = c(8, 4, 4, 2))
  plot(sample_tree, labels = sample_labels_short, main = "Sample clustering (1 - bicor | batch corrected)",
       xlab = "", sub = "", cex = 0.7, hang = 0.1)
  grDevices::dev.off()
}

##############################
## PARTIE 9) Heatmap distances samples
##############################
cat("---- HEATMAP DISTANCES SAMPLES ----\n")

sample_dists    <- dist(t(expr_mat_bc))
sample_dist_mat <- as.matrix(sample_dists)
rownames(sample_dist_mat) <- sample_labels_short
colnames(sample_dist_mat) <- sample_labels_short

annotation_df_plot <- data.frame(
  Protection = factor(meta_aligned$Protection, levels = names(PROTECTION_COLORS)),
  Genotype   = factor(as.character(meta_aligned$Genotype), levels = genotype_levels),
  row.names  = sample_labels_short
)

save_pheatmap_wgcna(
  mat        = sample_dist_mat, dist_obj   = sample_dists,
  ann_col    = annotation_df_plot, ann_row = annotation_df_plot,
  ann_colors = list(Protection = PROTECTION_COLORS, Genotype = genotype_colors),
  color_pal  = HEATMAP_DISTANCE,
  title_str  = "Sample distances | VST + removeBatchEffect(Replicate)",
  full_path  = file.path(dir_plots, "Sample_distance_heatmap"), w = 12, h = 10
)

##############################
## PARTIE 10) Préparation datExpr pour WGCNA
##############################
cat("---- PREPARATION DATEXPR ----\n")
## [C9] datExpr est construit sur la matrice BATCH-CORRIGÉE
datExpr <- as.data.frame(t(expr_mat_bc))
stopifnot(rownames(datExpr) == meta_aligned$Sample)
saveRDS(datExpr, file.path(dir_norm, "datExpr_BATCHCORRECTED_samples_x_genes.rds"))
cat("datExpr dim:", paste(dim(datExpr), collapse = " x "), "\n\n")

cat("====================================================\n")
cat("QC ECHANTILLONS TERMINE\n====================================================\n")

############################################################
## PARTIE 11) QC WGCNA — goodSamplesGenes + [C10] filtrage MAD top 50 %
############################################################
WGCNA::allowWGCNAThreads()

dir_net         <- file.path(out_dir, "03_network")
dir_mt          <- file.path(out_dir, "04_module_trait")
dir_hubs        <- file.path(out_dir, "05_hub_genes")
dir_wgcna_plots <- file.path(out_dir, "06_WGCNA_plots")
dir_go          <- file.path(out_dir, "07_GO_simple")
dir_cyto        <- file.path(out_dir, "08_Cytoscape")
dir_pres        <- file.path(out_dir, "09_modulePreservation")
dir_dc          <- file.path(out_dir, "10_diff_connectivity")
dir_netstats    <- file.path(out_dir, "11_network_stats")
dir_boot        <- file.path(out_dir, "12_bootstrap")

for (d in c(dir_net, dir_mt, dir_hubs, dir_wgcna_plots, dir_go,
            dir_cyto, dir_pres, dir_dc, dir_netstats, dir_boot))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("---- WGCNA QC: goodSamplesGenes ----\n")
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)

if (!gsg$allOK) {
  if (sum(!gsg$goodSamples) > 0) {
    cat("Removing bad samples:\n"); print(rownames(datExpr)[!gsg$goodSamples])
    datExpr      <- datExpr[gsg$goodSamples, , drop = FALSE]
    meta_aligned <- meta_aligned[match(rownames(datExpr), meta_aligned$Sample), , drop = FALSE]
  }
  if (sum(!gsg$goodGenes) > 0) {
    cat("Removing bad genes: n =", sum(!gsg$goodGenes), "\n")
    datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
  }
}
stopifnot(all(rownames(datExpr) == meta_aligned$Sample))
genotype_levels <- sort(unique(as.character(meta_aligned$Genotype)))
genotype_colors <- resolve_genotype_colors(genotype_levels)

## ── [C10] Filtrage MAD top 50 % ──────────────────────────────────────
## Van Dam et al. 2017 (Nature Protocols) & Langfelder & Horvath 2008 :
## conserver uniquement les gènes les plus variables améliore la topologie
## scale-free et réduit le bruit dans les modules.
cat("---- [C10] Filtrage MAD top", MAD_QUANTILE_CUTOFF * 100, "% ----\n")
n_before_mad <- ncol(datExpr)
mad_vec   <- apply(datExpr, 2, mad, na.rm = TRUE)
var_vec   <- apply(datExpr, 2, var, na.rm = TRUE)
bad_genes <- is.na(mad_vec) | mad_vec == 0 | is.na(var_vec) | var_vec == 0
datExpr   <- datExpr[, !bad_genes, drop = FALSE]
mad_vec   <- apply(datExpr, 2, mad, na.rm = TRUE)

mad_threshold <- quantile(mad_vec, MAD_QUANTILE_CUTOFF, na.rm = TRUE)
keep_mad      <- mad_vec >= mad_threshold
datExpr       <- datExpr[, keep_mad, drop = FALSE]
n_after_mad   <- ncol(datExpr)

## Distribution MAD — figure
mad_df <- data.frame(MAD = mad_vec, Kept = keep_mad)
p_mad  <- ggplot2::ggplot(mad_df, ggplot2::aes(x = MAD, fill = Kept)) +
  ggplot2::geom_histogram(bins = 80, colour = "white", linewidth = 0.15, alpha = 0.9) +
  ggplot2::geom_vline(xintercept = mad_threshold, linetype = "dashed",
                      colour = "#C0392B", linewidth = 0.7) +
  ggplot2::annotate("text", x = mad_threshold * 1.05,
                    y = max(hist(mad_vec, breaks = 80, plot = FALSE)$counts) * 0.9,
                    label = paste0("Seuil MAD = ", round(mad_threshold, 3),
                                   "\n(", round(MAD_QUANTILE_CUTOFF*100), "e percentile)"),
                    hjust = 0, size = 3.2, colour = "#C0392B") +
  ggplot2::scale_fill_manual(values = c("TRUE" = "#2471A3", "FALSE" = "grey75"),
                             labels = c("TRUE" = "Conservé", "FALSE" = "Exclu"),
                             name   = "Gène") +
  ggplot2::scale_x_continuous(labels = scales::label_comma()) +
  ggplot2::labs(
    title    = "Distribution MAD — filtrage avant WGCNA",
    subtitle = paste0("n avant = ", formatC(n_before_mad, big.mark=",", format="d"),
                      "  |  n après (top ", MAD_QUANTILE_CUTOFF*100, "%) = ",
                      formatC(n_after_mad, big.mark=",", format="d")),
    x = "Median Absolute Deviation (MAD)", y = "Nombre de gènes"
  ) + theme_pub(base_size = 11)

save_figure(p_mad, file.path(dir_net, "MAD_distribution_filtering"), width = 7, height = 4.5)
cat("Genes conservés (top", MAD_QUANTILE_CUTOFF*100, "% MAD):", n_after_mad, "/", n_before_mad, "\n\n")

write.csv(
  data.frame(Metric = c("Samples_final","Genes_pre_MAD","Genes_post_MAD","MAD_threshold"),
             Value  = c(nrow(datExpr), n_before_mad, n_after_mad, mad_threshold)),
  file.path(dir_net, "WGCNA_QC_summary.csv"), row.names = FALSE
)
saveRDS(datExpr, file.path(dir_net, "datExpr_WGCNA_ready.rds"))
write.csv(meta_aligned, file.path(dir_net, "metadata_after_WGCNA_QC.csv"), row.names = FALSE)

############################################################
## PARTIE 12) Soft-threshold + [C11] vérification scale-free
############################################################
cat("---- pickSoftThreshold ----\n")
Powers        <- 1:20
NETWORK_TYPE  <- "signed"
COR_TYPE      <- "bicor"
SFT_R2_TARGET <- 0.80

sft <- WGCNA::pickSoftThreshold(datExpr, powerVector = Powers, networkType = NETWORK_TYPE,
                                corFnc = COR_TYPE, verbose = 5)

fit       <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
power     <- sft$fitIndices[, 1]
ok        <- which(fit >= SFT_R2_TARGET)
softPower <- if (length(ok) > 0) power[min(ok)] else power[which.max(fit)]
cat("Chosen softPower:", softPower, "\n\n")

sft_df <- data.frame(Power = power, SFT_R2 = fit,
                     MeanK = sft$fitIndices[, 5], Slope = sft$fitIndices[, 3])
write.csv(sft_df, file.path(dir_net, "pickSoftThreshold_fitIndices.csv"), row.names = FALSE)

p_sft_r2 <- ggplot2::ggplot(sft_df, ggplot2::aes(x = Power, y = SFT_R2)) +
  ggplot2::geom_line(colour = "grey70", linewidth = 0.5) +
  ggplot2::geom_point(size = 2.5, colour = "#2471A3") +
  ggplot2::geom_text(ggplot2::aes(label = Power), vjust = -0.8, size = 2.8, colour = "grey30") +
  ggplot2::geom_hline(yintercept = SFT_R2_TARGET, linetype = "dashed",
                      colour = "grey45", linewidth = 0.5) +
  ggplot2::geom_vline(xintercept = softPower, linetype = "dashed",
                      colour = "#C0392B", linewidth = 0.5) +
  ggplot2::annotate("text", x = softPower + 0.4,
                    y = min(sft_df$SFT_R2, na.rm = TRUE),
                    label = paste0("\u03b2 = ", softPower),
                    colour = "#C0392B", hjust = 0, size = 3, fontface = "bold") +
  ggplot2::scale_x_continuous(breaks = seq(1, 20, 2)) +
  ggplot2::labs(title    = "Scale-free topology fit",
                subtitle = paste0("Target R\u00b2 \u2265 ", SFT_R2_TARGET,
                                  "  |  WGCNA signed network  |  bicor"),
                x = "Soft-threshold power (\u03b2)", y = "Scale-free topology fit (R\u00b2)") +
  theme_pub(base_size = 11)

p_sft_k <- ggplot2::ggplot(sft_df, ggplot2::aes(x = Power, y = MeanK)) +
  ggplot2::geom_line(colour = "grey70", linewidth = 0.5) +
  ggplot2::geom_point(size = 2.5, colour = "#2471A3") +
  ggplot2::geom_text(ggplot2::aes(label = Power), vjust = -0.8, size = 2.8, colour = "grey30") +
  ggplot2::geom_vline(xintercept = softPower, linetype = "dashed",
                      colour = "#C0392B", linewidth = 0.5) +
  ggplot2::scale_x_continuous(breaks = seq(1, 20, 2)) +
  ggplot2::labs(title    = "Mean connectivity",
                subtitle = paste0("Chosen \u03b2 = ", softPower),
                x = "Soft-threshold power (\u03b2)", y = "Mean connectivity") +
  theme_pub(base_size = 11)

save_figure(p_sft_r2, file.path(dir_wgcna_plots, "SoftThreshold_R2"),    width = 6, height = 4.5)
save_figure(p_sft_k,  file.path(dir_wgcna_plots, "SoftThreshold_MeanK"), width = 6, height = 4.5)

## ── [C11] Vérification distribution degrés de connectivité ───────────
## Validation formelle du fit scale-free : log(k) ~ log(p(k)) doit être
## linéaire. Essentiel pour répondre aux reviewers sur la qualité du réseau.
## Référence : Barabasi & Albert 1999 ; Langfelder & Horvath 2008 BMC Bioinformatics
cat("---- [C11] Scale-free degree distribution check ----\n")

k_soft <- WGCNA::softConnectivity(datExpr, power = softPower,
                                  type = "signed", corFnc = "bicor")
names(k_soft) <- colnames(datExpr)

for (ext in c("pdf","png","svg")) {
  fpath <- file.path(dir_wgcna_plots, paste0("ScaleFree_degree_distribution.", ext))
  if      (ext == "pdf") grDevices::pdf(fpath, width = 7, height = 5, useDingbats = FALSE)
  else if (ext == "png") grDevices::png(fpath, width = 7*220, height = 5*220, res = 220)
  else                   svglite::svglite(fpath, width = 7, height = 5)
  WGCNA::scaleFreePlot(k_soft,
                       main = paste0("Degree distribution check\n\u03b2 = ", softPower,
                                     "  |  signed network  |  bicor"))
  grDevices::dev.off()
}

## R² du log-log fit
k_nonzero   <- k_soft[k_soft > 0]
dk          <- density(k_nonzero)
k_vals      <- dk$x[dk$x > 0]
p_vals      <- dk$y[dk$x > 0]
lm_loglog   <- lm(log10(p_vals) ~ log10(k_vals))
r2_degree   <- summary(lm_loglog)$r.squared
slope_degree <- coef(lm_loglog)[2]
cat(sprintf("Scale-free log-log fit: R² = %.4f, slope = %.3f\n", r2_degree, slope_degree))
cat("R² ≥ 0.80 et slope < 0 = bon fit scale-free\n\n")

write.csv(data.frame(SoftPower = softPower, R2_SFT = max(fit, na.rm = TRUE),
                     R2_loglog = r2_degree, Slope_loglog = slope_degree),
          file.path(dir_net, "ScaleFree_fit_summary.csv"), row.names = FALSE)

############################################################
## PARTIE 13) Construction réseau + modules
############################################################
cat("---- blockwiseModules ----\n")

net <- WGCNA::blockwiseModules(
  datExpr,
  power             = softPower,
  networkType       = "signed",
  TOMType           = "signed",
  corType           = COR_TYPE,
  maxBlockSize      = ncol(datExpr),
  deepSplit         = 3,
  minModuleSize     = 50,
  mergeCutHeight    = 0.30,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = TRUE,          # [C6]
  saveTOMFileBase   = file.path(dir_net, "TOM"),
  verbose           = 5
)

colors <- WGCNA::labels2colors(net$colors)
MEs    <- WGCNA::orderMEs(net$MEs)

cat("N modules (incl grey):", length(unique(colors)), "\n")
cat("Top module sizes:\n"); print(head(sort(table(colors), decreasing = TRUE), 15)); cat("\n")

saveRDS(net,    file.path(dir_net, "WGCNA_network_object.rds"))
saveRDS(MEs,    file.path(dir_net, "module_eigengenes.rds"))
saveRDS(colors, file.path(dir_net, "module_colors.rds"))

## [C4] Mapping ME → couleur sans make.unique
map_me_color <- data.frame(
  ME          = colnames(MEs),
  ModuleLabel = suppressWarnings(as.integer(gsub("^ME", "", colnames(MEs)))),
  ModuleColor = WGCNA::labels2colors(suppressWarnings(as.integer(gsub("^ME", "", colnames(MEs))))),
  stringsAsFactors = FALSE
)
write.csv(data.frame(Gene = colnames(datExpr), ModuleLabel = as.integer(net$colors),
                     ModuleColor = colors),
          file.path(dir_net, "Gene_to_ModuleColor.csv"), row.names = FALSE)
write.csv(data.frame(ModuleColor = names(sort(table(colors), decreasing = TRUE)),
                     nGenes      = as.integer(sort(table(colors), decreasing = TRUE))),
          file.path(dir_net, "Module_sizes.csv"), row.names = FALSE)
write.csv(map_me_color, file.path(dir_net, "ME_to_ModuleColor_map.csv"), row.names = FALSE)

##############################
## PARTIE 13.1) Dendrogramme
##############################
cat("---- Gene dendrogram ----\n")
if (!is.null(net$dendrograms) && length(net$dendrograms) > 0) {
  for (b in seq_len(length(net$dendrograms))) {
    geneTree         <- net$dendrograms[[b]]
    blockGenes       <- net$blockGenes[[b]]
    if (is.null(blockGenes) || length(blockGenes) == 0) next
    moduleColors_block <- WGCNA::labels2colors(net$colors[blockGenes])
    for (ext in c("pdf","png","svg")) {
      fpath <- file.path(dir_wgcna_plots, paste0("Gene_dendrogram_block", b, ".", ext))
      if      (ext == "pdf") grDevices::pdf(fpath, width = 14, height = 6, useDingbats = FALSE)
      else if (ext == "png") grDevices::png(fpath, width = 2800, height = 1200, res = 220)
      else                   svglite::svglite(fpath, width = 14, height = 6)
      WGCNA::plotDendroAndColors(dendro = geneTree, colors = moduleColors_block,
                                 groupLabels = "Modules", dendroLabels = FALSE,
                                 hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                                 main = paste0("Gene dendrogram + module colors (block ", b, ")"))
      grDevices::dev.off()
    }
  }
}

##############################
## PARTIE 13.2) Bar chart tailles de modules
##############################
mod_sizes <- as.data.frame(table(colors), stringsAsFactors = FALSE)
colnames(mod_sizes) <- c("ModuleColor", "nGenes")
mod_sizes <- mod_sizes |> dplyr::filter(ModuleColor != "grey") |>
  dplyr::arrange(dplyr::desc(nGenes)) |> dplyr::mutate(pct = nGenes / sum(nGenes))
TOP_N_BAR <- 20
mod_plot  <- if (nrow(mod_sizes) > TOP_N_BAR) {
  other_row <- data.frame(ModuleColor = "other",
                          nGenes = sum(mod_sizes$nGenes[(TOP_N_BAR+1):nrow(mod_sizes)]),
                          pct    = sum(mod_sizes$pct[(TOP_N_BAR+1):nrow(mod_sizes)]))
  dplyr::bind_rows(mod_sizes[1:TOP_N_BAR, ], other_row)
} else { mod_sizes }
mod_plot <- mod_plot |>
  dplyr::mutate(ModuleColor = forcats::fct_reorder(ModuleColor, nGenes),
                fill_hex    = ifelse(ModuleColor == "other", "grey70", as.character(ModuleColor)))
fill_map <- setNames(as.character(mod_plot$fill_hex), as.character(mod_plot$ModuleColor))
p_bar <- ggplot2::ggplot(mod_plot, ggplot2::aes(x = ModuleColor, y = nGenes, fill = ModuleColor)) +
  ggplot2::geom_col(colour = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = nGenes), hjust = -0.15, size = 3,
                     fontface = "bold", colour = "grey20") +
  ggplot2::scale_fill_manual(values = fill_map) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18)),
                              labels = scales::label_comma()) +
  ggplot2::labs(title    = "Module sizes",
                subtitle = paste0("WGCNA signed  |  \u03b2 = ", softPower,
                                  "  |  n modules = ", length(unique(colors[colors != "grey"]))),
                x = "Module color", y = "Number of genes",
                caption = "Grey = unassigned") +
  ggplot2::coord_flip() + theme_pub(base_size = 11) +
  ggplot2::theme(legend.position = "none")
save_figure(p_bar, file.path(dir_wgcna_plots, "ModuleSizes_barplot"),
            width = 7, height = max(4, nrow(mod_plot) * 0.35 + 2))

##############################
## PARTIE 13.3) Donut tailles de modules
##############################
cat("---- Module size DONUT chart ----\n")
dir_pie <- file.path(dir_wgcna_plots, "module_sizes")
dir.create(dir_pie, recursive = TRUE, showWarnings = FALSE)
mod_sizes_pie <- as.data.frame(table(colors), stringsAsFactors = FALSE)
colnames(mod_sizes_pie) <- c("ModuleColor", "nGenes")
mod_sizes_pie <- mod_sizes_pie |> dplyr::arrange(dplyr::desc(nGenes)) |>
  dplyr::mutate(pct = nGenes / sum(nGenes),
                pct_label = paste0(sprintf("%.1f%%", 100 * pct)))
mod_grey_pie  <- mod_sizes_pie |> dplyr::filter(ModuleColor == "grey")
mod_named_pie <- mod_sizes_pie |> dplyr::filter(ModuleColor != "grey")
TOP_PIE <- 15
if (nrow(mod_named_pie) > TOP_PIE) {
  mod_other_pie <- mod_named_pie[(TOP_PIE+1):nrow(mod_named_pie), ]
  other_row_pie <- data.frame(ModuleColor = "other", nGenes = sum(mod_other_pie$nGenes),
                              pct = sum(mod_other_pie$pct),
                              pct_label = paste0(sprintf("%.1f%%", 100*sum(mod_other_pie$pct))))
  mod_pie_plot  <- dplyr::bind_rows(mod_named_pie[1:TOP_PIE,], other_row_pie, mod_grey_pie)
} else { mod_pie_plot <- dplyr::bind_rows(mod_named_pie, mod_grey_pie) }
mod_pie_plot <- mod_pie_plot |>
  dplyr::mutate(LegendLabel = paste0(ModuleColor, " — n=",
                                     formatC(nGenes, format = "d", big.mark = ","),
                                     " (", pct_label, ")"),
                LegendLabel = factor(LegendLabel, levels = LegendLabel))
fill_pie <- setNames(dplyr::case_when(mod_pie_plot$ModuleColor == "other" ~ "grey70",
                                      mod_pie_plot$ModuleColor == "grey"  ~ "#BDBDBD",
                                      TRUE ~ mod_pie_plot$ModuleColor),
                     as.character(mod_pie_plot$LegendLabel))
p_pie <- ggplot2::ggplot(mod_pie_plot, ggplot2::aes(x = 2, y = pct, fill = LegendLabel)) +
  ggplot2::geom_col(colour = "white", linewidth = 0.5) +
  ggplot2::coord_polar(theta = "y", start = 0) + ggplot2::xlim(0.5, 2.6) +
  ggplot2::geom_text(data = dplyr::filter(mod_pie_plot, pct >= 0.04),
                     ggplot2::aes(label = pct_label, x = 2.35),
                     position = ggplot2::position_stack(vjust = 0.5),
                     size = 2.8, fontface = "bold", colour = "white") +
  ggplot2::scale_fill_manual(values = fill_pie) +
  ggplot2::labs(title = "Module sizes — proportion", fill = "Module",
                caption = paste0("Total: ", formatC(sum(mod_pie_plot$nGenes), format="d", big.mark=","))) +
  theme_pub(base_size = 10) +
  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
                 axis.ticks = ggplot2::element_blank(), panel.border = ggplot2::element_blank())
save_figure(p_pie, file.path(dir_pie, "ModuleSizes_donut"), width = 9, height = 6)

##############################
## PARTIE 13.4) [C13] TOM plot — heatmap réseau intramodulaire
## Quasi-obligatoire dans les publications WGCNA (Langfelder & Horvath 2008)
##############################

cat("---- [C13] TOM plot ----\n")

tom_files <- list.files(dir_net, pattern = "^TOM.*\\.RData$", full.names = TRUE)

if (length(tom_files) > 0) {
  TOM_env  <- new.env()
  load(tom_files[1], envir = TOM_env)
  TOM_full <- as.matrix(get(ls(TOM_env)[1], envir = TOM_env))
  
  ## Top N gènes par module pour la figure
  TOP_TOM_PLOT <- 400
  block_genes_1 <- colnames(datExpr)[net$blockGenes[[1]]]
  
  mod_tab   <- sort(table(colors[colors != "grey"]), decreasing = TRUE)
  top_genes <- character(0)
  
  for (mc in names(mod_tab)) {
    g_mc <- colnames(datExpr)[colors == mc]
    g_mc <- g_mc[g_mc %in% block_genes_1]
    n_take <- max(5, round(min(length(g_mc), TOP_TOM_PLOT) * length(g_mc) / sum(mod_tab)))
    top_genes <- c(top_genes, head(g_mc, min(n_take, length(g_mc))))
  }
  
  top_genes <- intersect(top_genes, block_genes_1)
  if (length(top_genes) > TOP_TOM_PLOT) top_genes <- top_genes[1:TOP_TOM_PLOT]
  
  idx_tom   <- match(top_genes, block_genes_1)
  idx_tom   <- idx_tom[!is.na(idx_tom)]
  top_genes <- block_genes_1[idx_tom]
  
  if (length(top_genes) >= 2) {
    TOM_sub   <- TOM_full[idx_tom, idx_tom, drop = FALSE]
    colors_sub <- colors[match(top_genes, colnames(datExpr))]
    
    ## Dissimilarity TOM
    dissTOM <- 1 - TOM_sub
    
    ## Dendrogramme des gènes sélectionnés
    geneTree_TOM <- hclust(as.dist(dissTOM), method = "average")
    
    ## Réordonner selon le dendrogramme
    ord <- geneTree_TOM$order
    dissTOM_plot <- dissTOM[ord, ord]
    colors_plot  <- colors_sub[ord]
    
    for (ext in c("pdf", "png", "svg")) {
      fpath <- file.path(dir_wgcna_plots, paste0("TOM_plot_heatmap.", ext))
      
      if (ext == "pdf") {
        grDevices::pdf(fpath, width = 10, height = 10, useDingbats = FALSE)
      } else if (ext == "png") {
        grDevices::png(fpath, width = 10 * 220, height = 10 * 220, res = 220)
      } else {
        svglite::svglite(fpath, width = 10, height = 10)
      }
      
      WGCNA::TOMplot(
        dissTOM_plot,
        geneTree_TOM,
        colors_plot,
        main = paste0(
          "TOM heatmap  |  top ", length(top_genes),
          " genes  |  β = ", softPower,
          "  |  signed bicor\n",
          "Sorted by dendrogram / module color"
        )
      )
      
      grDevices::dev.off()
    }
    
    cat("✅ TOM plot sauvegardé\n\n")
  } else {
    cat("⚠️  Pas assez de gènes pour TOM plot\n\n")
  }
  
  rm(TOM_full, TOM_env); gc()
  
} else {
  cat("⚠️  Fichier TOM non trouvé — TOM plot ignoré\n\n")
}
##############################
## PARTIE 13.5) [C12] Analyse de sensibilité paramètres
## Re-coupe le dendrogramme avec différents paramètres (TOM déjà calculé)
## Référence : Zhang & Horvath 2005, Stat Appl Genet Mol Biol
##############################
if (RUN_SENSITIVITY && !is.null(net$dendrograms) && length(tom_files) > 0) {
  cat("---- [C12] Analyse sensibilité blockwiseModules ----\n")
  dir_sens <- file.path(dir_net, "sensitivity_analysis")
  dir.create(dir_sens, recursive = TRUE, showWarnings = FALSE)
  
  TOM_env2   <- new.env()
  load(tom_files[1], envir = TOM_env2)
  TOM_sens   <- as.matrix(get(ls(TOM_env2)[1], envir = TOM_env2))
  geneTree_s <- net$dendrograms[[1]]
  
  sensitivity_res <- list()
  param_grid <- expand.grid(deepSplit = SENSITIVITY_DEEPSPLIT,
                            minModuleSize = SENSITIVITY_MINMODSIZE,
                            mergeCutHeight = SENSITIVITY_MERGEHEIGHT,
                            stringsAsFactors = FALSE)
  
  for (i in seq_len(nrow(param_grid))) {
    ds  <- param_grid$deepSplit[i]
    mms <- param_grid$minModuleSize[i]
    mch <- param_grid$mergeCutHeight[i]
    
    labels_tmp <- dynamicTreeCut::cutreeDynamic(
      dendro = geneTree_s, distM = 1 - TOM_sens,
      deepSplit = ds, minClusterSize = mms,
      pamRespectsDendro = FALSE, verbose = 0
    )
    colors_tmp <- WGCNA::labels2colors(labels_tmp)
    
    ## Fusionner les modules similaires
    if (mch > 0) {
      MEs_tmp <- WGCNA::moduleEigengenes(
        datExpr[, colnames(datExpr)[seq_along(labels_tmp)], drop = FALSE],
        colors_tmp)$eigengenes
      merge_tmp <- WGCNA::mergeCloseModules(
        datExpr[, colnames(datExpr)[seq_along(labels_tmp)], drop = FALSE],
        colors_tmp, cutHeight = mch, verbose = 0)
      colors_tmp <- merge_tmp$colors
    }
    
    n_mod <- length(unique(colors_tmp[colors_tmp != "grey"]))
    n_grey <- sum(colors_tmp == "grey")
    
    sensitivity_res[[i]] <- data.frame(
      deepSplit = ds, minModuleSize = mms, mergeCutHeight = mch,
      N_modules = n_mod, N_grey = n_grey,
      PctGrey   = round(100 * n_grey / length(colors_tmp), 1),
      MeanSize  = round(mean(table(colors_tmp[colors_tmp != "grey"]))),
      MedianSize = round(median(table(colors_tmp[colors_tmp != "grey"])))
    )
    cat(sprintf("  ds=%d mms=%3d mch=%.2f → %d modules (%.1f%% grey)\n",
                ds, mms, mch, n_mod, round(100 * n_grey / length(colors_tmp), 1)))
  }
  
  rm(TOM_sens, TOM_env2); gc()
  
  sens_df <- dplyr::bind_rows(sensitivity_res) |>
    dplyr::mutate(
      Params = paste0("ds", deepSplit, "_mms", minModuleSize, "_mch", mergeCutHeight),
      MainParams = (deepSplit == 3 & minModuleSize == 50 & mergeCutHeight == 0.30)
    )
  write.csv(sens_df, file.path(dir_sens, "sensitivity_results.csv"), row.names = FALSE)
  
  p_sens <- ggplot2::ggplot(sens_df,
                            ggplot2::aes(x = factor(minModuleSize), y = N_modules,
                                         fill = factor(deepSplit))) +
    ggplot2::geom_col(position = "dodge", colour = "white", linewidth = 0.3) +
    ggplot2::geom_point(data = dplyr::filter(sens_df, MainParams),
                        ggplot2::aes(x = factor(minModuleSize), y = N_modules),
                        colour = "#C0392B", size = 4, shape = 18,
                        position = ggplot2::position_dodge(width = 0.9)) +
    ggplot2::facet_wrap(~paste0("mergeCutHeight = ", mergeCutHeight), ncol = 2) +
    ggplot2::scale_fill_manual(values = c("2" = "#6BAED6", "3" = "#2171B5", "4" = "#084594"),
                               name = "deepSplit") +
    ggplot2::labs(
      title    = "Analyse de sensibilité — paramètres blockwiseModules",
      subtitle = paste0("Réseau signed bicor  |  \u03b2 = ", softPower,
                        "  |  \u25c6 rouge = paramètres principaux (ds3/mms50/mch0.30)"),
      x = "minModuleSize", y = "Nombre de modules",
      caption = "Re-coupe du dendrogramme avec TOM pré-calculé"
    ) + theme_pub(base_size = 11)
  
  save_figure(p_sens, file.path(dir_sens, "sensitivity_N_modules"), width = 10, height = 5)
  
  ## Stabilité des assignments vs paramètres principaux
  cat("✅ Sensibilité terminée →", nrow(sens_df), "combinaisons\n\n")
} else {
  cat("⚠️  RUN_SENSITIVITY = FALSE ou dendrogramme/TOM absent — ignoré\n\n")
}

cat("====================================================\n")
cat("WGCNA NETWORK TERMINE — softPower:", softPower,
    "| n modules:", length(unique(colors)), "\n====================================================\n")

############################################################
## PARTIE 14) Module-Trait
############################################################
load_pkgs(c("lme4","lmerTest"), install_missing = INSTALL_MISSING_PKGS, bioc = FALSE)

dir_mt_plots    <- file.path(dir_mt, "plots")
dir_mt_tables   <- file.path(dir_mt, "tables")
dir_mt_boxplots <- file.path(dir_mt, "boxplots")
for (d in c(dir_mt_plots, dir_mt_tables, dir_mt_boxplots))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

meta_aligned <- meta_aligned |>
  dplyr::mutate(Protection = factor(Protection, levels = c("NonProtected","Protected")),
                Genotype   = factor(Genotype))
stopifnot(rownames(MEs) == meta_aligned$Sample)

trait_prot <- data.frame(Protection_binary = ifelse(meta_aligned$Protection == "Protected", 1, 0),
                         row.names = meta_aligned$Sample)

##############################
## PARTIE 14.1) [C14] Permutation p-values ME ~ Protection
## Remplace corPvalueStudent dont l'hypothèse d'indépendance est
## violée (MEs sont dérivés des mêmes données — Langfelder et al. 2011)
##############################
cat("---- [C14] ME ~ Protection — permutation p-values (n=", N_PERM_MT, ") ----\n")

stopifnot(rownames(trait_prot) == rownames(MEs))
prot_vec <- as.numeric(trait_prot[, 1])

## Corrélation observée
obs_cor_prot <- as.numeric(WGCNA::bicor(MEs, as.matrix(prot_vec), use = "pairwise.complete.obs"))

## Distribution nulle par permutation
set.seed(SEED)
perm_dist_prot <- matrix(NA_real_, nrow = ncol(MEs), ncol = N_PERM_MT)
pb <- txtProgressBar(min = 0, max = N_PERM_MT, style = 3)
for (i in seq_len(N_PERM_MT)) {
  perm_v <- sample(prot_vec)
  perm_dist_prot[, i] <- as.numeric(WGCNA::bicor(MEs, as.matrix(perm_v), use = "pairwise.complete.obs"))
  setTxtProgressBar(pb, i)
}
close(pb)

pval_perm_prot <- vapply(seq_len(ncol(MEs)), function(j)
  mean(abs(perm_dist_prot[j, ]) >= abs(obs_cor_prot[j])), numeric(1))
fdr_perm_prot  <- p.adjust(pval_perm_prot, method = "BH")

moduleTraitCor_prot <- matrix(obs_cor_prot, ncol = 1,
                              dimnames = list(colnames(MEs), "Protection_binary"))
moduleTraitP_prot   <- matrix(pval_perm_prot, ncol = 1,
                              dimnames = list(colnames(MEs), "Protection_binary"))
moduleTraitFDR_prot <- matrix(fdr_perm_prot, ncol = 1,
                              dimnames = list(colnames(MEs), "Protection_binary"))

write.csv(as.data.frame(moduleTraitCor_prot), file.path(dir_mt_tables, "ME_Protection_bicor.csv"))
write.csv(as.data.frame(moduleTraitP_prot),   file.path(dir_mt_tables, "ME_Protection_pval_PERMUTATION.csv"))
write.csv(as.data.frame(moduleTraitFDR_prot), file.path(dir_mt_tables, "ME_Protection_FDR_PERMUTATION.csv"))

me_display_labels <- paste0(
  map_me_color$ModuleColor[match(rownames(moduleTraitCor_prot), map_me_color$ME)],
  " (", rownames(moduleTraitCor_prot), ")"
)

stars_prot   <- apply(moduleTraitFDR_prot, c(1,2), starify)
textMat_prot <- paste0(sprintf("%.3f", moduleTraitCor_prot), stars_prot)
dim(textMat_prot) <- dim(moduleTraitCor_prot)

for (ext in c("pdf","png","svg")) {
  fpath <- file.path(dir_mt_plots, paste0("ME_vs_Protection_heatmap.", ext))
  h_val <- max(6, 0.28 * nrow(moduleTraitCor_prot) + 2)
  if      (ext == "pdf") grDevices::pdf(fpath, width = 8, height = h_val, useDingbats = FALSE)
  else if (ext == "png") grDevices::png(fpath, width = 1400,
                                        height = max(1000, 38*nrow(moduleTraitCor_prot)+300), res=220)
  else                   svglite::svglite(fpath, width = 8, height = h_val)
  par(mar = c(8, 12, 4, 2))
  WGCNA::labeledHeatmap(Matrix = moduleTraitCor_prot,
                        xLabels = colnames(moduleTraitCor_prot), yLabels = me_display_labels,
                        ySymbols = me_display_labels, colorLabels = FALSE,
                        colors = WGCNA::blueWhiteRed(50), textMatrix = textMat_prot,
                        setStdMargins = FALSE, cex.text = 0.7, zlim = c(-1,1),
                        main = "Module eigengenes vs Protection\n(bicor + permutation FDR ; n=1000)")
  grDevices::dev.off()
}

##############################
## PARTIE 14.2) Corrélation ME ~ Genotype
##############################
cat("---- Correlation MEs ~ Genotype ----\n")
trait_geno <- stats::model.matrix(~ 0 + Genotype, data = meta_aligned)
rownames(trait_geno) <- meta_aligned$Sample
colnames(trait_geno) <- gsub("^Genotype", "", colnames(trait_geno))
moduleTraitCor_geno <- WGCNA::bicor(MEs, trait_geno, use = "pairwise.complete.obs")
moduleTraitP_geno   <- WGCNA::corPvalueStudent(moduleTraitCor_geno, nSamples = nrow(MEs))
moduleTraitFDR_geno <- matrix(p.adjust(as.vector(moduleTraitP_geno), method="BH"),
                              nrow = nrow(moduleTraitP_geno), ncol = ncol(moduleTraitP_geno),
                              dimnames = dimnames(moduleTraitP_geno))
write.csv(moduleTraitCor_geno, file.path(dir_mt_tables, "ME_Genotype_cor.csv"))
write.csv(moduleTraitP_geno,   file.path(dir_mt_tables, "ME_Genotype_pvalue.csv"))
write.csv(moduleTraitFDR_geno, file.path(dir_mt_tables, "ME_Genotype_FDR_BH.csv"))

stars_geno   <- apply(moduleTraitFDR_geno, c(1,2), starify)
textMat_geno <- paste0(sprintf("%.3f", moduleTraitCor_geno), stars_geno)
dim(textMat_geno) <- dim(moduleTraitCor_geno)
for (ext in c("pdf","png","svg")) {
  fpath <- file.path(dir_mt_plots, paste0("ME_vs_Genotype_heatmap.", ext))
  w_val <- max(10, 0.8*ncol(moduleTraitCor_geno)+4)
  h_val <- max(6,  0.28*nrow(moduleTraitCor_geno)+2)
  if      (ext == "pdf") grDevices::pdf(fpath, width = w_val, height = h_val, useDingbats = FALSE)
  else if (ext == "png") grDevices::png(fpath,
                                        width  = max(1600, 180*ncol(moduleTraitCor_geno)+400),
                                        height = max(1000, 38*nrow(moduleTraitCor_geno)+300), res=220)
  else                   svglite::svglite(fpath, width = w_val, height = h_val)
  par(mar = c(10,12,4,2))
  WGCNA::labeledHeatmap(Matrix = moduleTraitCor_geno, xLabels = colnames(moduleTraitCor_geno),
                        yLabels = me_display_labels, ySymbols = me_display_labels,
                        colorLabels = FALSE, colors = WGCNA::blueWhiteRed(50),
                        textMatrix = textMat_geno, setStdMargins = FALSE,
                        cex.text=0.6, cex.lab.x=0.8, cex.lab.y=0.8, xLabelsAngle=45,
                        zlim=c(-1,1), main="Module eigengenes vs Genotype\n(bicor ; BH FDR)")
  grDevices::dev.off()
}

##############################
## PARTIE 14.3) Heatmap combinée
##############################
cat("---- Correlation MEs ~ Protection + Genotype (combined) ----\n")
trait_all <- cbind(trait_prot, trait_geno)
moduleTraitCor_all <- WGCNA::bicor(MEs, trait_all, use = "pairwise.complete.obs")
moduleTraitP_all   <- WGCNA::corPvalueStudent(moduleTraitCor_all, nSamples = nrow(MEs))
moduleTraitFDR_all <- matrix(p.adjust(as.vector(moduleTraitP_all), method="BH"),
                             nrow=nrow(moduleTraitP_all), ncol=ncol(moduleTraitP_all),
                             dimnames=dimnames(moduleTraitP_all))
write.csv(moduleTraitCor_all, file.path(dir_mt_tables, "ME_Protection_Genotype_cor.csv"))
write.csv(moduleTraitP_all,   file.path(dir_mt_tables, "ME_Protection_Genotype_pvalue.csv"))
write.csv(moduleTraitFDR_all, file.path(dir_mt_tables, "ME_Protection_Genotype_FDR_BH.csv"))
stars_all   <- apply(moduleTraitFDR_all, c(1,2), starify)
textMat_all <- paste0(sprintf("%.3f", moduleTraitCor_all), stars_all)
dim(textMat_all) <- dim(moduleTraitCor_all)
for (ext in c("pdf","png","svg")) {
  fpath <- file.path(dir_mt_plots, paste0("ME_vs_Protection_Genotype_heatmap.", ext))
  w_val <- max(10, 0.8*ncol(moduleTraitCor_all)+4)
  h_val <- max(6,  0.28*nrow(moduleTraitCor_all)+2)
  if      (ext == "pdf") grDevices::pdf(fpath, width=w_val, height=h_val, useDingbats=FALSE)
  else if (ext == "png") grDevices::png(fpath,
                                        width=max(1600,180*ncol(moduleTraitCor_all)+400),
                                        height=max(1000,38*nrow(moduleTraitCor_all)+300), res=220)
  else                   svglite::svglite(fpath, width=w_val, height=h_val)
  par(mar = c(10,12,4,2))
  WGCNA::labeledHeatmap(Matrix=moduleTraitCor_all, xLabels=colnames(moduleTraitCor_all),
                        yLabels=me_display_labels, ySymbols=me_display_labels,
                        colorLabels=FALSE, colors=WGCNA::blueWhiteRed(50),
                        textMatrix=textMat_all, setStdMargins=FALSE,
                        cex.text=0.6, cex.lab.x=0.8, cex.lab.y=0.8, xLabelsAngle=45,
                        zlim=c(-1,1), main="ME vs Protection + Genotype\n(bicor ; BH FDR)")
  grDevices::dev.off()
}

##############################
## PARTIE 14.4) Tables résumées
##############################
res_simple <- data.frame(
  ME          = rownames(moduleTraitCor_prot),
  ModuleColor = map_me_color$ModuleColor[match(rownames(moduleTraitCor_prot), map_me_color$ME)],
  Correlation = as.numeric(moduleTraitCor_prot[,1]),
  Pval_perm   = as.numeric(moduleTraitP_prot[,1]),
  FDR_perm    = as.numeric(moduleTraitFDR_prot[,1])
) |> dplyr::arrange(FDR_perm, dplyr::desc(abs(Correlation)))
write.csv(res_simple, file.path(dir_mt_tables, "ME_Protection_summary_permutation.csv"), row.names=FALSE)
cat("Top modules (bicor + permutation):\n"); print(head(res_simple, 10)); cat("\n")

res_geno_long <- as.data.frame(moduleTraitCor_geno) |> tibble::rownames_to_column("ME") |>
  tidyr::pivot_longer(-ME, names_to="Genotype", values_to="Correlation") |>
  dplyr::left_join(as.data.frame(moduleTraitP_geno) |> tibble::rownames_to_column("ME") |>
                     tidyr::pivot_longer(-ME, names_to="Genotype", values_to="Pvalue"), by=c("ME","Genotype")) |>
  dplyr::left_join(as.data.frame(moduleTraitFDR_geno) |> tibble::rownames_to_column("ME") |>
                     tidyr::pivot_longer(-ME, names_to="Genotype", values_to="FDR"), by=c("ME","Genotype")) |>
  dplyr::left_join(map_me_color[,c("ME","ModuleColor")], by="ME") |>
  dplyr::arrange(FDR, dplyr::desc(abs(Correlation)))
write.csv(res_geno_long, file.path(dir_mt_tables, "ME_Genotype_summary_long.csv"), row.names=FALSE)

##############################
## PARTIE 14.5) [C15] Modèle mixte — TEST PRIMAIRE
## lmerTest : ME ~ Protection + (1|Genotype)
## Contrôle la structure de pseudo-réplication génotypique.
## C'est le test PRINCIPAL rapporté dans le papier ; les
## corrélations bicor + permutation servent de confirmation.
## Référence : Bates et al. 2015 J Stat Softw (lme4)
##############################
cat("---- [C15] Modèles mixtes (TEST PRIMAIRE) ----\n")
res_mixed <- list()
for (me in colnames(MEs)) {
  df_lmer <- data.frame(Sample = meta_aligned$Sample, Protection = meta_aligned$Protection,
                        Genotype = meta_aligned$Genotype, ME = as.numeric(MEs[, me]))
  fit_lmer <- tryCatch(lmerTest::lmer(ME ~ Protection + (1|Genotype), data=df_lmer), error=function(e) NULL)
  if (is.null(fit_lmer)) {
    res_mixed[[me]] <- data.frame(ME=me, pProtection=NA_real_, betaProtection=NA_real_)
    next
  }
  an  <- tryCatch(anova(fit_lmer), error=function(e) NULL)
  cf  <- tryCatch(summary(fit_lmer)$coefficients, error=function(e) NULL)
  pProt    <- if (!is.null(an)  && "Protection" %in% rownames(an))
    as.numeric(an["Protection","Pr(>F)"]) else NA_real_
  betaProt <- if (!is.null(cf) && "ProtectionProtected" %in% rownames(cf))
    as.numeric(cf["ProtectionProtected","Estimate"]) else NA_real_
  res_mixed[[me]] <- data.frame(ME=me, pProtection=pProt, betaProtection=betaProt)
}
res_mixed <- dplyr::bind_rows(res_mixed)
res_mixed$FDR <- p.adjust(res_mixed$pProtection, method="BH")
res_mixed <- res_mixed |> dplyr::left_join(map_me_color, by="ME") |>
  dplyr::arrange(FDR, dplyr::desc(abs(betaProtection)))
write.csv(res_mixed, file.path(dir_mt_tables, "ME_Protection_mixedModel_PRIMARY.csv"), row.names=FALSE)
cat("Top modules (modèle mixte — TEST PRIMAIRE):\n"); print(head(res_mixed, 10)); cat("\n")

##############################
## PARTIE 14.6) Sélection modules significatifs [C2] FDR = 0.10
##############################
FDR_CUTOFF_MAIN <- 0.90
modules_mixed_sig  <- res_mixed  |> dplyr::filter(FDR < FDR_CUTOFF_MAIN)
modules_simple_sig <- res_simple |> dplyr::filter(FDR_perm < FDR_CUTOFF_MAIN)
write.csv(modules_mixed_sig,  file.path(dir_mt_tables, "ME_sig_mixedModel.csv"),  row.names=FALSE)
write.csv(modules_simple_sig, file.path(dir_mt_tables, "ME_sig_permutation.csv"), row.names=FALSE)
cat("Modules sig. mixte (FDR <", FDR_CUTOFF_MAIN, "):", nrow(modules_mixed_sig), "\n")
cat("Modules sig. bicor+perm (FDR <", FDR_CUTOFF_MAIN, "):", nrow(modules_simple_sig), "\n\n")

##############################
## PARTIE 14.7) Boxplots modules candidats
##############################
top_modules_me <- unique(c(head(res_mixed$ME[order(res_mixed$FDR)], 6),
                           head(res_simple$ME[order(res_simple$FDR_perm)], 6)))
top_modules_me <- top_modules_me[top_modules_me %in% colnames(MEs)]

for (me in top_modules_me) {
  df_box <- data.frame(Sample = meta_aligned$Sample, Protection = meta_aligned$Protection,
                       Genotype = meta_aligned$Genotype, ME = as.numeric(MEs[, me]))
  colr  <- map_me_color$ModuleColor[match(me, map_me_color$ME)]
  fdr_v <- res_mixed$FDR[match(me, res_mixed$ME)]
  fdr_l <- if (!is.na(fdr_v)) paste0("lmerTest FDR = ", formatC(fdr_v, format="e", digits=2)) else "FDR = NA"
  p_box <- ggplot2::ggplot(df_box, ggplot2::aes(x=Protection, y=ME, fill=Protection)) +
    ggplot2::geom_boxplot(outlier.shape=NA, alpha=0.8, linewidth=0.6) +
    ggplot2::geom_jitter(ggplot2::aes(colour=Genotype), width=0.12, size=2.2, alpha=0.85) +
    ggplot2::scale_fill_manual(values=PROTECTION_COLORS) +
    ggplot2::scale_colour_manual(values=genotype_colors,
                                 guide=ggplot2::guide_legend(title="Genotype",override.aes=list(size=3))) +
    ggplot2::labs(title=paste0("Module ", colr, "  (", me, ")"),
                  subtitle=paste0("Eigengene ~ Protection  |  ", fdr_l),
                  x=NULL, y="Module eigengene") +
    theme_pub(base_size=11) + ggplot2::theme(legend.position="right")
  save_figure(p_box, file.path(dir_mt_boxplots, paste0(colr,"_",me,"_boxplot")), width=6, height=4.5)
}

##############################
## PARTIE 14.8) [C17] Eigengene network + méta-modules
## Langfelder & Horvath 2007, BMC Systems Biology — la hiérarchie
## entre modules révèle les méta-programmes biologiques co-régulés.
##############################
cat("---- [C17] Eigengene network + meta-modules ----\n")
dir_eig <- file.path(dir_mt, "eigengene_network")
dir.create(dir_eig, recursive = TRUE, showWarnings = FALSE)

## Exclure grey
MEs_noGrey <- MEs[, !grepl("ME0$", colnames(MEs)), drop = FALSE]

## Ajouter trait Protection comme trait factice pour visualisation
MEsWT <- data.frame(MEs_noGrey,
                    MEProtection = as.numeric(meta_aligned$Protection == "Protected"),
                    row.names    = rownames(MEs_noGrey))

## A) Dendrogramme eigengene
for (ext in c("pdf","png","svg")) {
  fpath <- file.path(dir_eig, paste0("Eigengene_network_dendrogram.", ext))
  h_val <- max(5, ncol(MEs_noGrey) * 0.18 + 2)
  if      (ext == "pdf") grDevices::pdf(fpath, width=10, height=h_val, useDingbats=FALSE)
  else if (ext == "png") grDevices::png(fpath, width=2200, height=max(1200,h_val*200), res=220)
  else                   svglite::svglite(fpath, width=10, height=h_val)
  WGCNA::plotEigengeneNetworks(MEsWT, "Eigengene network",
                               marHeatmap    = c(3, 4, 2, 2),
                               marDendro     = c(0, 4, 1, 2),
                               plotDendrograms = TRUE,
                               xLabelsAngle  = 90,
                               heatmapColors = WGCNA::blueWhiteRed(50))
  grDevices::dev.off()
}

## B) Corrélation entre MEs — heatmap publication-quality
ME_cor_mat <- cor(MEs_noGrey, method = "pearson")
ph_eig     <- pheatmap::pheatmap(ME_cor_mat,
                                 color            = HEATMAP_DIVERGING,
                                 clustering_method = "ward.D2",
                                 treeheight_row   = 30, treeheight_col = 30,
                                 fontsize         = 7,
                                 border_color     = NA,
                                 main  = "ME correlation heatmap (Pearson)\nHierarchical clustering = meta-modules",
                                 silent           = TRUE)
for (ext in c("pdf","png","svg")) {
  fpath <- file.path(dir_eig, paste0("ME_correlation_heatmap.", ext))
  n_me  <- ncol(ME_cor_mat)
  w_eig <- max(7, n_me * 0.28 + 2)
  h_eig <- max(6, n_me * 0.28 + 2)
  if (ext == "pdf") {
    grDevices::pdf(fpath, width=w_eig, height=h_eig, useDingbats=FALSE)
    pheatmap::pheatmap(ME_cor_mat, color=HEATMAP_DIVERGING, clustering_method="ward.D2",
                       fontsize=7, border_color=NA,
                       main="ME correlation heatmap (Pearson)\nHierarchical clustering = meta-modules")
    grDevices::dev.off()
  } else if (ext == "png") {
    pheatmap::pheatmap(ME_cor_mat, color=HEATMAP_DIVERGING, clustering_method="ward.D2",
                       fontsize=7, border_color=NA,
                       main="ME correlation heatmap (Pearson)\nHierarchical clustering = meta-modules",
                       filename=fpath, width=w_eig, height=h_eig)
  } else {
    svglite::svglite(fpath, width=w_eig, height=h_eig)
    grid::grid.newpage(); grid::grid.draw(ph_eig$gtable); grDevices::dev.off()
  }
}

## C) Identification méta-modules (cut du dendrogramme eigengenes)
eig_dist  <- as.dist(1 - ME_cor_mat)
eig_tree  <- hclust(eig_dist, method = "ward.D2")
h_cut_meta <- 0.25   # hauteur de coupe pour méta-modules
meta_mod_labels <- cutree(eig_tree, h = h_cut_meta)
meta_mod_df     <- data.frame(ME = names(meta_mod_labels),
                              MetaModule = paste0("Meta_", meta_mod_labels),
                              ModuleColor = map_me_color$ModuleColor[match(names(meta_mod_labels), map_me_color$ME)])
meta_mod_df <- meta_mod_df |> dplyr::arrange(MetaModule)
write.csv(meta_mod_df, file.path(dir_eig, "meta_modules_assignment.csv"), row.names=FALSE)
cat("Meta-modules (cut h=", h_cut_meta, "):\n"); print(table(meta_mod_df$MetaModule)); cat("\n")

## Sauvegarder la corrélation ME
write.csv(ME_cor_mat, file.path(dir_eig, "ME_correlation_matrix.csv"))
saveRDS(meta_mod_df, file.path(dir_eig, "meta_modules_df.rds"))
cat("✅ Eigengene network + meta-modules terminé\n\n")

cat("====================================================\n")
cat("MODULE-TRAIT TERMINE\n====================================================\n")

############################################################
## PARTIE 15) Visualisation modules par génotype
############################################################
dir_modviz          <- file.path(dir_mt, "module_views_by_genotype")
dir_modviz_global   <- file.path(dir_modviz, "A_global_by_protection")
dir_modviz_genotype <- file.path(dir_modviz, "B_by_genotype")
dir_modviz_combo    <- file.path(dir_modviz, "C_genotype_with_protection")
dir_modviz_tables   <- file.path(dir_modviz, "D_summary_tables")
for (d in c(dir_modviz, dir_modviz_global, dir_modviz_genotype, dir_modviz_combo, dir_modviz_tables))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

meta_aligned <- meta_aligned |>
  dplyr::mutate(Protection=factor(Protection, levels=c("NonProtected","Protected")),
                Genotype=factor(Genotype), GenotypeBase=factor(GenotypeBase))

make_me_df <- function(me_name) {
  data.frame(Sample=meta_aligned$Sample, Protection=meta_aligned$Protection,
             Genotype=meta_aligned$Genotype, GenotypeBase=meta_aligned$GenotypeBase,
             Replicate=meta_aligned$Replicate, ME=as.numeric(MEs[, me_name]))
}

all_module_summary <- list()
for (me in colnames(MEs)) {
  df   <- make_me_df(me)
  colr <- map_me_color$ModuleColor[match(me, map_me_color$ME)]
  lab  <- map_me_color$ModuleLabel[match(me, map_me_color$ME)]
  p_global <- ggplot2::ggplot(df, ggplot2::aes(x=Protection, y=ME, fill=Protection)) +
    ggplot2::geom_boxplot(outlier.shape=NA, alpha=0.8, linewidth=0.6) +
    ggplot2::geom_jitter(ggplot2::aes(colour=Genotype), width=0.12, size=2, alpha=0.85) +
    ggplot2::scale_fill_manual(values=PROTECTION_COLORS) +
    ggplot2::scale_colour_manual(values=genotype_colors,
                                 guide=ggplot2::guide_legend(title="Genotype",override.aes=list(size=3))) +
    ggplot2::labs(title=paste0(colr,"  (",me,")"), subtitle="Vue globale par protection",
                  x=NULL, y="Module eigengene") +
    theme_pub(base_size=11) + ggplot2::theme(legend.position="right")
  save_figure(p_global, file.path(dir_modviz_global, paste0(colr,"_",me,"_global_protection")),
              width=7, height=5)
  geno_order <- df |> dplyr::group_by(Genotype) |>
    dplyr::summarise(meanME=mean(ME,na.rm=TRUE),.groups="drop") |>
    dplyr::arrange(dplyr::desc(meanME)) |> dplyr::pull(Genotype)
  df$Genotype <- factor(df$Genotype, levels=geno_order)
  p_genotype <- ggplot2::ggplot(df, ggplot2::aes(x=Genotype, y=ME, fill=Protection)) +
    ggplot2::geom_boxplot(outlier.shape=NA, alpha=0.8, linewidth=0.5) +
    ggplot2::geom_jitter(width=0.12, size=1.8, alpha=0.8, colour="grey30") +
    ggplot2::scale_fill_manual(values=PROTECTION_COLORS) +
    ggplot2::labs(title=paste0(colr,"  (",me,")"), subtitle="Boxplot par génotype", x=NULL, y="Module eigengene") +
    theme_pub(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1,vjust=1,size=9),legend.position="right")
  save_figure(p_genotype, file.path(dir_modviz_genotype, paste0(colr,"_",me,"_by_genotype")), width=8, height=5)
  p_combo <- ggplot2::ggplot(df, ggplot2::aes(x=Genotype, y=ME, colour=Protection)) +
    ggplot2::geom_point(position=ggplot2::position_jitter(width=0.12), size=2.2, alpha=0.9) +
    ggplot2::stat_summary(fun=mean, geom="point", size=4, shape=18, colour="black") +
    ggplot2::scale_colour_manual(values=PROTECTION_COLORS) +
    ggplot2::labs(title=paste0(colr,"  (",me,")"),
                  subtitle="Distribution par génotype  |  ◆ = moyenne", x=NULL, y="Module eigengene") +
    theme_pub(base_size=11) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1,vjust=1,size=9),legend.position="right")
  save_figure(p_combo, file.path(dir_modviz_combo, paste0(colr,"_",me,"_genotype_points")), width=8, height=5)
  summary_tbl <- df |> dplyr::group_by(Genotype,GenotypeBase,Protection) |>
    dplyr::summarise(n=dplyr::n(), meanME=mean(ME,na.rm=TRUE), sdME=sd(ME,na.rm=TRUE),
                     medianME=median(ME,na.rm=TRUE), .groups="drop") |>
    dplyr::arrange(dplyr::desc(meanME))
  write.csv(summary_tbl, file.path(dir_modviz_tables, paste0(colr,"_",me,"_summary_by_genotype.csv")), row.names=FALSE)
  summary_tbl$ME <- me; summary_tbl$ModuleColor <- colr; summary_tbl$ModuleLabel <- lab
  all_module_summary[[me]] <- summary_tbl
}
dplyr::bind_rows(all_module_summary) |>
  write.csv(file.path(dir_modviz_tables, "ALL_modules_summary_by_genotype.csv"), row.names=FALSE)
cat("---- PARTIE 15 OK ----\n")

############################################################
## PARTIE 16) Hub genes (kME + Gene Significance)
############################################################
cat("==== PARTIE 16 : Hub genes ====\n")
dir_hubs_plots <- file.path(dir_hubs, "plots")
dir.create(dir_hubs,       recursive=TRUE, showWarnings=FALSE)
dir.create(dir_hubs_plots, recursive=TRUE, showWarnings=FALSE)

KME_THRESHOLD <- 0.80
GS_THRESHOLD  <- 0.30

kME_mat <- WGCNA::signedKME(datExpr, MEs, corFnc="bicor")
colnames(kME_mat) <- paste0("ME", gsub("^kME", "", colnames(kME_mat)))

gene_module_df <- data.frame(Gene=colnames(datExpr), ModuleLabel=as.integer(net$colors),
                             ModuleColor=colors, stringsAsFactors=FALSE)
gene_module_df$kME_own <- mapply(function(gene, me_col) {
  if (me_col %in% colnames(kME_mat)) kME_mat[gene, me_col] else NA_real_
}, gene_module_df$Gene, paste0("ME", gene_module_df$ModuleLabel))

write.csv(gene_module_df, file.path(dir_hubs, "All_genes_ModuleMembership.csv"), row.names=FALSE)
saveRDS(kME_mat, file.path(dir_hubs, "kME_matrix_all_modules.rds"))

protection_vec <- ifelse(meta_aligned$Protection == "Protected", 1, 0)
names(protection_vec) <- meta_aligned$Sample
GS_prot   <- WGCNA::bicor(datExpr, as.matrix(protection_vec), use="pairwise.complete.obs")
GS_prot_p <- WGCNA::corPvalueStudent(GS_prot, nSamples=nrow(datExpr))

GS_df <- data.frame(Gene=colnames(datExpr), GS_Protection=as.numeric(GS_prot),
                    GS_Pvalue=as.numeric(GS_prot_p), stringsAsFactors=FALSE)

hub_full_df <- gene_module_df |> dplyr::left_join(GS_df, by="Gene") |>
  dplyr::arrange(dplyr::desc(abs(kME_own)))
write.csv(hub_full_df, file.path(dir_hubs, "All_genes_kME_GS_Protection.csv"), row.names=FALSE)

hub_genes <- hub_full_df |> dplyr::filter(kME_own >= KME_THRESHOLD, ModuleColor != "grey") |>
  dplyr::arrange(ModuleColor, dplyr::desc(kME_own))
write.csv(hub_genes, file.path(dir_hubs, "Hub_genes_kME_ge_0.8.csv"), row.names=FALSE)
cat("Hub genes (kME ≥", KME_THRESHOLD, "):", nrow(hub_genes), "\n")

hub_priority <- hub_full_df |>
  dplyr::filter(kME_own >= KME_THRESHOLD, abs(GS_Protection) >= GS_THRESHOLD, ModuleColor != "grey") |>
  dplyr::left_join(res_mixed[,c("ME","FDR","betaProtection")] |>
                     dplyr::mutate(ModuleLabel=as.integer(gsub("^ME","",ME))), by="ModuleLabel") |>
  dplyr::arrange(FDR, dplyr::desc(kME_own))
write.csv(hub_priority, file.path(dir_hubs, "Hub_genes_PRIORITY_kME_AND_GS.csv"), row.names=FALSE)
cat("Hub genes prioritaires:", nrow(hub_priority), "\n\n")

sig_colors <- res_mixed |> dplyr::filter(FDR < 0.10) |> dplyr::pull(ModuleColor)
if (length(sig_colors) == 0) sig_colors <- head(res_mixed$ModuleColor[order(res_mixed$FDR)], 6)
hub_sig <- hub_genes |> dplyr::filter(ModuleColor %in% sig_colors)
write.csv(hub_sig, file.path(dir_hubs, "Hub_genes_significant_modules.csv"), row.names=FALSE)

for (me in colnames(MEs)) {
  colr  <- map_me_color$ModuleColor[match(me, map_me_color$ME)]
  fdr_v <- res_mixed$FDR[match(me, res_mixed$ME)]
  fdr_l <- if (!is.na(fdr_v)) paste0("lmerTest FDR = ", formatC(fdr_v, format="e", digits=2)) else ""
  df_sc <- hub_full_df |> dplyr::filter(ModuleColor == colr) |>
    dplyr::mutate(IsHub=kME_own >= KME_THRESHOLD,
                  HighGS=abs(GS_Protection) >= GS_THRESHOLD,
                  Priority=IsHub & HighGS)
  if (nrow(df_sc) < 3) next
  top_label <- df_sc |> dplyr::filter(Priority) |> dplyr::arrange(dplyr::desc(kME_own)) |> dplyr::slice_head(n=15)
  df_sc$point_color <- dplyr::case_when(df_sc$Priority~"#C0392B",df_sc$IsHub~"#E8820C",df_sc$HighGS~"#2471A3",TRUE~"grey80")
  p_sc <- ggplot2::ggplot(df_sc, ggplot2::aes(x=kME_own, y=GS_Protection)) +
    ggplot2::geom_point(colour=df_sc$point_color,
                        size=ifelse(df_sc$Priority,2.8,ifelse(df_sc$IsHub|df_sc$HighGS,2.0,1.2)),
                        alpha=ifelse(df_sc$Priority,0.95,ifelse(df_sc$IsHub|df_sc$HighGS,0.80,0.40))) +
    ggrepel::geom_text_repel(data=top_label, ggplot2::aes(label=Gene), size=2.4, max.overlaps=20,
                             segment.linewidth=0.3, segment.colour="grey55",
                             colour="#C0392B", fontface="italic", box.padding=0.35) +
    ggplot2::geom_vline(xintercept=KME_THRESHOLD, linetype="dashed", colour="grey50", linewidth=0.4) +
    ggplot2::geom_hline(yintercept=c(-GS_THRESHOLD, GS_THRESHOLD), linetype="dashed", colour="grey50", linewidth=0.4) +
    ggplot2::geom_hline(yintercept=0, colour="black", linewidth=0.3) +
    ggplot2::annotate("text",x=0.62,y=0.96,
                      label=paste0("● Hub+GS: n=",sum(df_sc$Priority)),
                      size=2.8,colour="#C0392B",hjust=0) +
    ggplot2::annotate("text",x=0.62,y=0.88,
                      label=paste0("● Hub seul: n=",sum(df_sc$IsHub&!df_sc$Priority)),
                      size=2.8,colour="#E8820C",hjust=0) +
    ggplot2::scale_x_continuous(limits=c(NA,1.02),expand=ggplot2::expansion(mult=c(0.02,0.02))) +
    ggplot2::scale_y_continuous(limits=c(-1,1),expand=ggplot2::expansion(mult=c(0.02,0.02))) +
    ggplot2::labs(title=paste0("kME vs GS — module ", colr, "  (", me, ")"),
                  subtitle=paste0("bicor(gene, Protection)  |  ", fdr_l),
                  x=paste0("Module Membership (kME)  |  seuil = ", KME_THRESHOLD),
                  y="Gene Significance ~ Protection",
                  caption=paste0("n = ", nrow(df_sc), " genes")) +
    theme_pub(base_size=11)
  save_figure(p_sc, file.path(dir_hubs_plots, paste0(colr,"_",me,"_kME_vs_GS")), width=7, height=5.5)
}
cat("✅ Hub genes terminé\n\n")

############################################################
## PARTIE 16.6) Heatmap expression des gènes par module
############################################################
cat("==== PARTIE 16.6 : Module gene heatmaps ====\n")
dir_modheat <- file.path(dir_hubs, "module_heatmaps")
dir.create(dir_modheat, recursive=TRUE, showWarnings=FALSE)
HEAT_TOP_GENES <- 100

heat_modules <- res_mixed |> dplyr::filter(FDR < 0.10, ModuleColor != "grey") |>
  dplyr::arrange(FDR) |> dplyr::pull(ME)
if (length(heat_modules) == 0) heat_modules <- head(res_mixed$ME[order(res_mixed$FDR)], 6)

ann_col_heat <- data.frame(Protection=factor(meta_aligned$Protection, levels=c("NonProtected","Protected")),
                           Genotype=factor(as.character(meta_aligned$Genotype), levels=genotype_levels),
                           row.names=meta_aligned$Sample)
ann_colors_heat <- list(Protection=PROTECTION_COLORS, Genotype=genotype_colors)

for (me in heat_modules) {
  colr     <- map_me_color$ModuleColor[match(me, map_me_color$ME)]
  fdr_v    <- res_mixed$FDR[match(me, res_mixed$ME)]
  beta_v   <- res_mixed$betaProtection[match(me, res_mixed$ME)]
  n_total  <- sum(colors == colr)
  genes_heat <- hub_full_df |> dplyr::filter(ModuleColor == colr) |>
    dplyr::arrange(dplyr::desc(kME_own)) |> dplyr::slice_head(n=HEAT_TOP_GENES) |>
    dplyr::pull(Gene)
  genes_heat <- genes_heat[genes_heat %in% colnames(datExpr)]
  if (length(genes_heat) < 4) next
  mat_heat <- t(as.matrix(datExpr[, genes_heat, drop=FALSE]))
  mat_heat <- mat_heat[, meta_aligned$Sample, drop=FALSE]
  mat_z    <- t(scale(t(mat_heat))); mat_z[is.nan(mat_z)] <- 0
  mat_z    <- pmax(pmin(mat_z, 3), -3)
  gene_info <- hub_full_df |> dplyr::filter(Gene %in% genes_heat) |>
    dplyr::mutate(GeneType=dplyr::case_when(
      kME_own >= KME_THRESHOLD & abs(GS_Protection) >= GS_THRESHOLD ~ "Hub + GS fort",
      kME_own >= KME_THRESHOLD ~ "Hub", abs(GS_Protection) >= GS_THRESHOLD ~ "GS fort",
      TRUE ~ "Autre"))
  ann_row_heat  <- data.frame(GeneType=factor(gene_info$GeneType[match(genes_heat,gene_info$Gene)],
                                              levels=c("Hub + GS fort","Hub","GS fort","Autre")), row.names=genes_heat)
  ann_colors_row <- list(GeneType=c("Hub + GS fort"="#C0392B","Hub"="#E8820C","GS fort"="#2471A3","Autre"="grey80"))
  col_order <- meta_aligned |> dplyr::arrange(Protection,Genotype) |> dplyr::pull(Sample)
  mat_z     <- mat_z[, col_order, drop=FALSE]
  title_str <- paste0("Module ", colr, " (", me, ") — top ", length(genes_heat), "/", n_total,
                      " genes (kME)\nFDR=", formatC(fdr_v,format="e",digits=2),
                      "  β=", round(beta_v,3), "  VST z-score  clamp±3")
  h_heat <- max(8, min(length(genes_heat) * 0.13 + 4, 30))
  base_path <- file.path(dir_modheat, paste0(colr,"_",me,"_gene_heatmap"))
  args_ph <- list(mat=mat_z, color=HEATMAP_DIVERGING,
                  annotation_col=ann_col_heat[col_order,,drop=FALSE],
                  annotation_row=ann_row_heat,
                  annotation_colors=c(ann_colors_heat, ann_colors_row),
                  cluster_rows=TRUE, cluster_cols=FALSE, clustering_method="ward.D2",
                  show_rownames=length(genes_heat)<=60, show_colnames=FALSE,
                  fontsize_row=max(4, 8-length(genes_heat)*0.04), fontsize=8,
                  border_color=NA, treeheight_row=25, treeheight_col=0,
                  main=title_str, silent=TRUE)
  for (ext in c("pdf","png","svg")) {
    fpath <- paste0(base_path,".",ext)
    if (ext == "pdf") {
      grDevices::pdf(fpath, width=14, height=h_heat, useDingbats=FALSE)
      ph <- do.call(pheatmap::pheatmap, args_ph); grDevices::dev.off()
    } else if (ext == "png") {
      args_ph_png <- args_ph; args_ph_png$filename <- fpath
      args_ph_png$width <- 14; args_ph_png$height <- h_heat; args_ph_png$silent <- FALSE
      do.call(pheatmap::pheatmap, args_ph_png)
    } else {
      ph <- do.call(pheatmap::pheatmap, args_ph)
      svglite::svglite(fpath, width=14, height=h_heat)
      grid::grid.newpage(); grid::grid.draw(ph$gtable); grDevices::dev.off()
    }
  }
  cat("✅ Heatmap:", colr, "[", length(genes_heat), "genes ×", ncol(mat_z), "samples]\n")
}
cat("\n====================================================\n")
cat("MODULE GENE HEATMAPS TERMINEES\n====================================================\n\n")

############################################################
## PARTIE 17) Enrichissement GO [C16] FDR cutoff 0.20
############################################################
cat("==== PARTIE 17 : GO enrichment ====\n")
load_pkgs("topGO", install_missing=INSTALL_MISSING_PKGS, bioc=TRUE)

dir_go_lists  <- file.path(dir_go, "01_gene_lists")
dir_go_tables <- file.path(dir_go, "02_tables")
dir_go_plots  <- file.path(dir_go, "03_plots")
for (d in c(dir_go, dir_go_lists, dir_go_tables, dir_go_plots))
  dir.create(d, recursive=TRUE, showWarnings=FALSE)

GO_PVALUE_CUTOFF     <- 0.05
GO_MIN_GENES         <- 4
TOP_TERMS_PLOT       <- 20
## [C16] GO_FDR_MODULE_CUTOFF 0.50 → 0.20
GO_FDR_MODULE_CUTOFF <- 0.50
IF_NONE_TAKE_TOP_N   <- 20

stopifnot(file.exists(go_annot_file))
go_annot <- readr::read_csv(go_annot_file, show_col_types=FALSE) |>
  dplyr::rename(gene=Gene_ID, go=Annotation) |> dplyr::distinct() |>
  dplyr::mutate(gene=as.character(gene), go=as.character(go))
go_df <- go_annot |>
  tidyr::separate_rows(go, sep="[,;\\s]+", convert=FALSE) |>
  dplyr::mutate(gene=trimws(gene), go=trimws(go)) |>
  dplyr::filter(!is.na(gene), gene!="", !is.na(go), go!="", grepl("^GO:\\d+$", go)) |>
  dplyr::distinct(gene, go)

genes_network <- colnames(datExpr)
gene_universe <- intersect(unique(go_df$gene), unique(genes_network))
if (length(gene_universe) == 0) stop("Aucun gène commun entre GO et réseau WGCNA.")
geneID2GO <- split(go_df$go[go_df$gene %in% gene_universe], go_df$gene[go_df$gene %in% gene_universe])
cat("Nb gènes univers GO:", length(gene_universe), "\n\n")

modules_from_mixed <- res_mixed |> dplyr::filter(!is.na(FDR)) |>
  dplyr::arrange(FDR, dplyr::desc(abs(betaProtection)))
modules_sig <- modules_from_mixed |> dplyr::filter(FDR < GO_FDR_MODULE_CUTOFF) |> dplyr::pull(ME)
if (length(modules_sig) == 0) {
  modules_sig <- head(modules_from_mixed$ME, IF_NONE_TAKE_TOP_N)
  cat("Aucun module FDR <", GO_FDR_MODULE_CUTOFF, "→ top", IF_NONE_TAKE_TOP_N, "\n")
}
modules_labels_v <- suppressWarnings(as.integer(gsub("^ME","", modules_sig)))
modules_ok       <- !is.na(modules_labels_v)
modules_final    <- data.frame(ModuleID=modules_sig[modules_ok],
                               ModuleLabel=modules_labels_v[modules_ok],
                               ModuleColor=WGCNA::labels2colors(modules_labels_v[modules_ok]))
write.csv(modules_final, file.path(dir_go, "Selected_modules_for_GO.csv"), row.names=FALSE)
cat("Modules retenus pour GO:\n"); print(modules_final); cat("\n")

do_GO_for_module <- function(module_color, module_label, out_dir_lists, out_dir_tables, out_dir_plots,
                             p_cutoff=0.05, min_genes=4, top_terms_plot=15) {
  genes_mod_all <- colnames(datExpr)[as.integer(net$colors) == module_label]
  genes_in      <- intersect(genes_mod_all, gene_universe)
  module_tag    <- paste0(module_color, "_label", module_label)
  for (d in c(file.path(out_dir_lists,module_tag), file.path(out_dir_tables,module_tag),
              file.path(out_dir_plots,module_tag)))
    dir.create(d, recursive=TRUE, showWarnings=FALSE)
  writeLines(sort(genes_mod_all), file.path(out_dir_lists, module_tag, "Module_all_genes.txt"))
  writeLines(sort(genes_in),      file.path(out_dir_lists, module_tag, "Module_genes_in_GO_universe.txt"))
  if (length(genes_in) < min_genes) { warning("Module ", module_color, ": trop peu de gènes annotés"); return(invisible(FALSE)) }
  geneList <- factor(as.integer(gene_universe %in% genes_in), levels=c(0,1))
  names(geneList) <- gene_universe
  collector <- list()
  for (ont in c("BP","MF","CC")) {
    cat("   -> Module", module_color, "| Ontologie", ont, "\n")
    GOdata      <- new("topGOdata", ontology=ont, allGenes=geneList,
                       annot=annFUN.gene2GO, gene2GO=geneID2GO)
    res_w01     <- runTest(GOdata, algorithm="weight01", statistic="fisher")
    res_classic <- runTest(GOdata, algorithm="classic",  statistic="fisher")
    allRes      <- GenTable(GOdata, weight01Fisher=res_w01, classicFisher=res_classic,
                            topNodes=length(score(res_w01)))
    allRes$weight01Fisher <- suppressWarnings(as.numeric(as.character(allRes$weight01Fisher)))
    allRes$classicFisher  <- suppressWarnings(as.numeric(as.character(allRes$classicFisher)))
    allRes$weight01Fisher[is.na(allRes$weight01Fisher)|allRes$weight01Fisher==0] <- 1e-300
    allRes$classicFisher[is.na(allRes$classicFisher)|allRes$classicFisher==0]   <- 1e-300
    term2genes    <- lapply(allRes$GO.ID, function(go_id)
      sort(unique(intersect(genesInTerm(GOdata, go_id)[[1]], genes_in))))
    gene_count    <- vapply(term2genes, length, integer(1))
    genes_in_term <- vapply(term2genes, function(v) paste(v, collapse=";"), character(1))
    keep <- which(allRes$weight01Fisher < p_cutoff & gene_count >= min_genes)
    if (!length(keep)) { message("Module ", module_color, ": aucun terme GO pour ", ont); next }
    out_ont <- allRes[keep,] |>
      dplyr::select(GO.ID,Term,Annotated,Significant,Expected,weight01Fisher,classicFisher) |>
      dplyr::mutate(genes_in_term=genes_in_term[match(GO.ID,allRes$GO.ID)],
                    gene_count=gene_count[match(GO.ID,allRes$GO.ID)],
                    Ontology=ont, ModuleColor=module_color, ModuleLabel=module_label,
                    log10_weight01=-log10(weight01Fisher),
                    Term_short=stringr::str_wrap(Term, width=42)) |>
      dplyr::arrange(weight01Fisher, dplyr::desc(gene_count))
    readr::write_csv(out_ont, file.path(out_dir_tables,module_tag,
                                        paste0("GO_",ont,"_",module_color,"_with_genes.csv")))
    collector[[ont]] <- out_ont
  }
  if (length(collector) == 0) { message("Module ", module_color, ": aucun enrichissement GO."); return(invisible(FALSE)) }
  merged_mod <- dplyr::bind_rows(collector)
  readr::write_csv(merged_mod, file.path(out_dir_tables,module_tag,paste0("GO_ALL_",module_color,"_with_genes.csv")))
  merged_plot <- merged_mod |> dplyr::group_by(Ontology) |>
    dplyr::arrange(weight01Fisher, dplyr::desc(gene_count), .by_group=TRUE) |>
    dplyr::slice_head(n=top_terms_plot) |> dplyr::ungroup() |>
    dplyr::mutate(Term_short=forcats::fct_reorder(Term_short,log10_weight01),
                  Ontology=factor(Ontology, levels=c("BP","MF","CC")))
  if (nrow(merged_plot) == 0) return(invisible(FALSE))
  p_dot <- ggplot2::ggplot(merged_plot, ggplot2::aes(x=log10_weight01,y=Term_short,size=gene_count,colour=log10_weight01)) +
    ggplot2::geom_segment(ggplot2::aes(x=0,xend=log10_weight01,y=Term_short,yend=Term_short),
                          linewidth=0.3, colour="grey82") +
    ggplot2::geom_point(alpha=0.92) +
    ggplot2::facet_wrap(~Ontology, scales="free_y", ncol=1) +
    ggplot2::scale_colour_gradient(low="#6BAED6", high="#CB181D",
                                   name=expression(-log[10](italic(p)))) +
    ggplot2::scale_size_continuous(range=c(3,9), name="Gene count", breaks=scales::breaks_pretty(n=3)) +
    ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(0,0.15))) +
    ggplot2::labs(title=paste0("GO enrichment — module ", module_color),
                  subtitle=paste0("topGO weight01 Fisher  |  top ", top_terms_plot, " terms/ontology"),
                  x=expression(-log[10](italic(p))~"weight01"), y=NULL) +
    theme_pub(base_size=9) +
    ggplot2::theme(axis.text.y=ggplot2::element_text(size=7.5),
                   strip.text=ggplot2::element_text(face="bold",size=9),
                   panel.spacing=ggplot2::unit(2,"mm"), legend.position="right")
  n_terms <- nrow(merged_plot)
  h_plot  <- max(5, min(n_terms * 0.22 + 2, 22))
  save_figure(p_dot, file.path(out_dir_plots,module_tag,paste0("GO_ALL_",module_color,"_dotplot")),
              width=8, height=h_plot)
  message("\u2705 GO terminé module ", module_color)
  invisible(TRUE)
}

for (i in seq_len(nrow(modules_final))) {
  try(do_GO_for_module(module_color=modules_final$ModuleColor[i],
                       module_label=modules_final$ModuleLabel[i],
                       out_dir_lists=dir_go_lists, out_dir_tables=dir_go_tables, out_dir_plots=dir_go_plots,
                       p_cutoff=GO_PVALUE_CUTOFF, min_genes=GO_MIN_GENES, top_terms_plot=TOP_TERMS_PLOT),
      silent=TRUE)
}
cat("====================================================\n")
cat("GO SIMPLE TERMINE\n====================================================\n")

############################################################
## PARTIE 18) Export Cytoscape
############################################################
cat("==== PARTIE 18 : Export Cytoscape ====\n")
dir.create(dir_cyto, recursive=TRUE, showWarnings=FALSE)
CYTO_TOP_GENES  <- 100
CYTO_TOM_CUTOFF <- 0.02

cyto_modules <- res_mixed |> dplyr::filter(FDR < 0.10, ModuleColor != "grey") |> dplyr::pull(ME)
if (length(cyto_modules) == 0) cyto_modules <- head(res_mixed$ME[order(res_mixed$FDR)], 4)

node_attrs <- hub_full_df |>
  dplyr::select(Gene,ModuleLabel,ModuleColor,kME_own,GS_Protection,GS_Pvalue) |>
  dplyr::mutate(IsHub=kME_own >= KME_THRESHOLD, IsHighGS=abs(GS_Protection) >= GS_THRESHOLD,
                IsPriority=IsHub & IsHighGS,
                GS_direction=dplyr::case_when(GS_Protection > GS_THRESHOLD ~ "Protected",
                                              GS_Protection < -GS_THRESHOLD ~ "NonProtected",
                                              TRUE ~ "Neutral"),
                NodeColor=dplyr::case_when(IsPriority~"#C0392B",IsHub~"#E8820C",IsHighGS~"#2471A3",TRUE~"grey80"),
                NodeSize=dplyr::case_when(IsPriority~50,IsHub~35,TRUE~20))
write.csv(node_attrs, file.path(dir_cyto, "Cytoscape_node_attributes_ALL.csv"), row.names=FALSE)

if (length(tom_files) > 0) {
  for (me in cyto_modules) {
    colr      <- map_me_color$ModuleColor[match(me, map_me_color$ME)]
    genes_mod <- hub_full_df |> dplyr::filter(ModuleColor==colr) |>
      dplyr::arrange(dplyr::desc(kME_own)) |> dplyr::slice_head(n=CYTO_TOP_GENES) |> dplyr::pull(Gene)
    if (length(genes_mod) < 2) next
    block_id <- 1
    if (length(net$blockGenes) > 1) {
      for (b in seq_along(net$blockGenes)) {
        if (any(colnames(datExpr)[net$blockGenes[[b]]] %in% genes_mod)) { block_id <- b; break }
      }
    }
    tom_file_b <- tom_files[grep(paste0("block.",block_id), tom_files)]
    if (length(tom_file_b) == 0) tom_file_b <- tom_files[1]
    TOM_env <- new.env(); load(tom_file_b, envir=TOM_env)
    TOM_mat <- as.matrix(get(ls(TOM_env)[1], envir=TOM_env))
    block_genes <- colnames(datExpr)[net$blockGenes[[block_id]]]
    idx_genes   <- match(genes_mod, block_genes); idx_genes <- idx_genes[!is.na(idx_genes)]
    genes_ok    <- block_genes[idx_genes]
    if (length(genes_ok) < 2) { next }
    TOM_sub <- TOM_mat[idx_genes, idx_genes, drop=FALSE]
    rownames(TOM_sub) <- genes_ok; colnames(TOM_sub) <- genes_ok
    edges <- which(upper.tri(TOM_sub) & TOM_sub >= CYTO_TOM_CUTOFF, arr.ind=TRUE)
    if (nrow(edges) == 0) next
    edge_df <- data.frame(fromNode=rownames(TOM_sub)[edges[,1]],
                          toNode=colnames(TOM_sub)[edges[,2]],
                          weight=TOM_sub[edges], direction="undirected") |>
      dplyr::arrange(dplyr::desc(weight))
    readr::write_csv(edge_df, file.path(dir_cyto, paste0(colr,"_",me,"_edges_TOM.csv")))
    readr::write_csv(node_attrs |> dplyr::filter(Gene %in% genes_ok),
                     file.path(dir_cyto, paste0(colr,"_",me,"_nodes.csv")))
    cat("✅ Cytoscape:", colr, "(", nrow(edge_df), "edges)\n")
    rm(TOM_mat, TOM_sub, TOM_env); gc()
  }
}

############################################################
## PARTIE 20) [C19] Connectivité différentielle (DC)
## Identifier les gènes qui changent de connectivité entre
## Protected et NonProtected — révèle les réseaux conditionnels.
## Référence : McKenzie et al. 2016, BMC Bioinformatics (DGCA)
############################################################
cat("==== PARTIE 20 : [C19] Connectivité différentielle ====\n")
dir.create(dir_dc, recursive=TRUE, showWarnings=FALSE)

idx_prot    <- which(meta_aligned$Protection == "Protected")
idx_nonprot <- which(meta_aligned$Protection == "NonProtected")

if (length(idx_prot) >= 5 && length(idx_nonprot) >= 5) {
  
  cat("Calcul softConnectivity Protected...\n")
  k_prot <- WGCNA::softConnectivity(
    datExpr[idx_prot, , drop=FALSE], power=softPower, type="signed", corFnc="bicor"
  )
  cat("Calcul softConnectivity NonProtected...\n")
  k_nonprot <- WGCNA::softConnectivity(
    datExpr[idx_nonprot, , drop=FALSE], power=softPower, type="signed", corFnc="bicor"
  )
  names(k_prot) <- names(k_nonprot) <- colnames(datExpr)
  
  DC_df <- data.frame(
    Gene           = colnames(datExpr),
    k_Protected    = k_prot,
    k_NonProtected = k_nonprot,
    DC             = k_prot - k_nonprot,
    ModuleColor    = colors,
    ModuleLabel    = as.integer(net$colors),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(GS_df, by="Gene") |>
    dplyr::left_join(gene_module_df[,c("Gene","kME_own")], by="Gene") |>
    dplyr::mutate(
      DC_zscore    = scale(DC)[,1],
      DC_direction = dplyr::case_when(DC >  quantile(DC, 0.90) ~ "Gain_Protected",
                                      DC <  quantile(DC, 0.10) ~ "Loss_Protected",
                                      TRUE                      ~ "Stable")
    ) |>
    dplyr::arrange(dplyr::desc(abs(DC)))
  
  write.csv(DC_df, file.path(dir_dc, "Differential_connectivity_all_genes.csv"), row.names=FALSE)
  
  ## Top DC per module
  DC_top <- DC_df |> dplyr::filter(ModuleColor != "grey") |>
    dplyr::group_by(ModuleColor) |>
    dplyr::arrange(dplyr::desc(abs(DC)), .by_group=TRUE) |>
    dplyr::slice_head(n=20) |> dplyr::ungroup()
  write.csv(DC_top, file.path(dir_dc, "Differential_connectivity_top20_per_module.csv"), row.names=FALSE)
  
  ## Figure globale DC distribution
  p_dc_hist <- ggplot2::ggplot(DC_df |> dplyr::filter(ModuleColor != "grey"),
                               ggplot2::aes(x=DC, fill=DC_direction)) +
    ggplot2::geom_histogram(bins=80, colour="white", linewidth=0.1, alpha=0.9) +
    ggplot2::geom_vline(xintercept=0, colour="black", linewidth=0.5) +
    ggplot2::scale_fill_manual(values=c("Gain_Protected"="#C0392B","Loss_Protected"="#2471A3","Stable"="grey75"),
                               name="DC direction") +
    ggplot2::labs(
      title    = "Distribution de la connectivité différentielle",
      subtitle = paste0("DC = k_Protected − k_NonProtected  |  signed bicor  |  β=", softPower),
      x = "DC (connectivity gain in Protected)", y = "Nombre de gènes",
      caption  = paste0("n Protected = ", length(idx_prot), "  |  n NonProtected = ", length(idx_nonprot))
    ) + theme_pub(base_size=11)
  save_figure(p_dc_hist, file.path(dir_dc, "DC_distribution"), width=7, height=4.5)
  
  ## Figure DC vs GS par module significatif
  sig_me <- head(res_mixed$ME[order(res_mixed$FDR)], 6)
  for (me in sig_me) {
    colr   <- map_me_color$ModuleColor[match(me, map_me_color$ME)]
    dc_mod <- DC_df |> dplyr::filter(ModuleColor == colr)
    if (nrow(dc_mod) < 5) next
    top_dc <- dc_mod |> dplyr::arrange(dplyr::desc(abs(DC))) |> dplyr::slice_head(n=20)
    p_dc_gs <- ggplot2::ggplot(dc_mod, ggplot2::aes(x=GS_Protection, y=DC)) +
      ggplot2::geom_point(colour=colr, size=1.5, alpha=0.5) +
      ggrepel::geom_text_repel(data=top_dc, ggplot2::aes(label=Gene),
                               size=2.4, max.overlaps=15, colour="#C0392B",
                               segment.colour="grey60", box.padding=0.4) +
      ggplot2::geom_hline(yintercept=0, colour="grey30", linewidth=0.4) +
      ggplot2::geom_vline(xintercept=0, colour="grey30", linewidth=0.4) +
      ggplot2::labs(
        title    = paste0("DC vs GS — module ", colr, " (", me, ")"),
        subtitle = "DC = connectivity gain in Protected | GS = bicor(gene, Protection)",
        x = "Gene Significance ~ Protection (GS)", y = "Differential Connectivity (DC)"
      ) + theme_pub(base_size=11)
    save_figure(p_dc_gs, file.path(dir_dc, paste0(colr,"_",me,"_DC_vs_GS")), width=6, height=5)
  }
  cat("✅ Connectivité différentielle terminée\n\n")
} else {
  cat("⚠️  Trop peu d'échantillons par groupe pour DC (min 5)\n\n")
}
############################################################
## SUMMARY & FINALISATION
############################################################
cat("\n====================================================\n")
cat("SCRIPT WGCNA v3 COMPLET TERMINE\n")
cat("  out_dir :", out_dir, "\n")
cat("====================================================\n\n")

cat("RESUME DES SORTIES :\n")
cat("  00_inputs/                   → données alignées\n")
cat("  01_normalization/            → DESeq2, VST, batch-corrected expr_mat\n")
cat("  02_QC_samples/plots/         → PCA (raw+BC), clustering, distance heatmap\n")
cat("  03_network/                  → TOM, network, MAD distrib, scale-free check,\n")
cat("                                  sensibilité paramètres\n")
cat("  04_module_trait/             → heatmaps, boxplots, lmerTest (primaire),\n")
cat("                                  permut p-values, méta-modules eigengene\n")
cat("  05_hub_genes/                → kME, GS, hub genes, scatter plots, heatmaps\n")
cat("  06_WGCNA_plots/              → soft-threshold, dendrogramme, TOM plot,\n")
cat("                                  bar/donut chart, scale-free distrib\n")
cat("  07_GO_simple/                → enrichissement topGO par module (FDR ≤ 0.20)\n")
cat("  08_Cytoscape/                → edges TOM + node attributes\n")
cat("  09_modulePreservation/       → Zsummary, MedianRank, figures [C18]\n")
cat("  10_diff_connectivity/        → DC par gène Protected vs NonProt [C19]\n")
cat("  11_network_stats/            → betweenness, transitivity, igraph [C20]\n")
cat("  12_bootstrap/                → stabilité hub genes bootstrap [C21]\n\n")

writeLines(c("R session information","=====================", capture.output(sessionInfo())),
           file.path(out_dir, "session_info.txt"))

message("\u2705 WGCNA v3 FINAL complet")
message("   Corrections critiques  : removeBatchEffect [C9], MAD top50% [C10],")
message("   scaleFreePlot [C11], sensibilité [C12], TOM plot [C13],")
message("   permutation p-val [C14], lmerTest primaire [C15], GO FDR 0.20 [C16]")
message("   Analyses nouvelles     : eigengene network [C17], modulePreservation [C18],")
message("   diff. connectivity [C19], igraph stats [C20], bootstrap [C21]")