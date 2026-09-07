############################################################
## SCRIPT V2 — INTEGRATION GENES + ANNOTATIONS BIOLOGIQUES
## ---------------------------------------------------------
## Objectif :
## 1. Lire un fichier d'annotation fonctionnelle
## 2. Lire automatiquement des sorties LOO / WGCNA
## 3. Ajouter les fonctions biologiques
## 4. Ajouter un nom lisible de gène/protéine
## 5. Exporter un tableau final intégré
############################################################

##############################
## PARTIE 1) PARAMETRES
##############################
## USER CONFIGURATION
## Functional annotation table for the DH13M14 reference
annotation_file <- "/path/to/annotation/DH13M14_functional_annotation.txt"

## Outputs produced by the upstream WGCNA/LOO integration and edgeR scripts
hub_nodes_file <- "/path/to/results/WGCNA/integration/Hub_network_nodes.csv"
tier1_both_file <- "/path/to/results/edgeR_LOO_topGO/02_DEG/DEG_Tier1_Both_Sides.csv"
tier1_both_up_file <- "/path/to/results/edgeR_LOO_topGO/02_DEG/Tier_Split_By_Direction/Tier1_both_Up.csv"

out_dir <- "/path/to/results/functional_annotation_integration"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##############################
## PARTIE 2) PACKAGES
##############################
load_pkgs <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) stop("Package manquant : ", p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

load_pkgs(c("dplyr", "readr", "stringr", "tidyr", "tibble", "purrr"))

##############################
## PARTIE 3) FONCTIONS UTILES
##############################

clean_empty <- function(x) {
  x <- as.character(x)
  x[x %in% c("-", "", "NA", "NaN")] <- NA
  x
}

extract_go_text <- function(x) {
  x <- clean_empty(x)
  ifelse(
    is.na(x),
    NA_character_,
    x |>
      stringr::str_replace_all("GO:\\d+\\(", "") |>
      stringr::str_replace_all("\\)", "") |>
      stringr::str_squish()
  )
}

clean_annotation_text <- function(x) {
  x <- clean_empty(x)
  ifelse(
    is.na(x),
    NA_character_,
    x |>
      stringr::str_replace_all("\\s+", " ") |>
      stringr::str_squish()
  )
}

detect_gene_col <- function(df) {
  candidates <- c("gene_id", "Gene", "gene", "name", "Gene_ID", "ID")
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) stop("Aucune colonne identifiant gène trouvée dans un fichier.")
  hit[1]
}

read_generic_table <- function(filepath) {
  ext <- tools::file_ext(filepath)
  
  if (!file.exists(filepath)) return(NULL)
  
  if (ext == "csv") {
    readr::read_csv(filepath, show_col_types = FALSE)
  } else if (ext %in% c("txt", "tsv")) {
    readr::read_tsv(filepath, show_col_types = FALSE)
  } else {
    stop("Format non supporté : ", filepath)
  }
}

## Essaie d'extraire un nom de protéine lisible
extract_gene_name <- function(swissprot_annot, nr_annot, kog_annot, kegg_annot) {
  src <- dplyr::coalesce(swissprot_annot, nr_annot, kog_annot, kegg_annot)
  src <- clean_annotation_text(src)
  
  out <- src
  
  ## enlever les préfixes d'accession type "XP_... PREDICTED:"
  out <- stringr::str_replace(out, "^[^ ]+\\s+", "")
  out <- stringr::str_replace(out, "^PREDICTED:\\s*", "")
  
  ## enlever les morceaux après OS=, OX=, GN=, PE=, SV=
  out <- stringr::str_replace(out, "\\s+OS=.*$", "")
  out <- stringr::str_replace(out, "\\s+OX=.*$", "")
  out <- stringr::str_replace(out, "\\s+GN=.*$", "")
  out <- stringr::str_replace(out, "\\s+PE=.*$", "")
  out <- stringr::str_replace(out, "\\s+SV=.*$", "")
  
  ## pour KEGG : garder surtout la partie après ";"
  out <- ifelse(
    !is.na(kegg_annot) & grepl(";", kegg_annot),
    stringr::str_trim(stringr::str_replace(kegg_annot, ".*;\\s*", "")),
    out
  )
  
  out <- stringr::str_squish(out)
  out[out %in% c("-", "", "NA")] <- NA
  
  out
}

##############################
## PARTIE 4) IMPORT ANNOTATION
##############################
cat("---- IMPORT ANNOTATION ----\n")
stopifnot(file.exists(annotation_file))

annot_raw <- readr::read_tsv(annotation_file, show_col_types = FALSE)
cat("Dimensions annotation :", nrow(annot_raw), "x", ncol(annot_raw), "\n")

required_cols <- c(
  "gene_id",
  "KOG_Annotation",
  "GO:BiologicalProcess",
  "GO:CellularComponent",
  "GO:MolecularFunction",
  "KEGG_Ortholog",
  "NR_Annotation",
  "Swissprot_Annotation"
)

missing_cols <- setdiff(required_cols, colnames(annot_raw))
if (length(missing_cols) > 0) {
  warning("Colonnes absentes dans le fichier d'annotation : ",
          paste(missing_cols, collapse = ", "))
}

annot_clean <- annot_raw |>
  dplyr::mutate(
    gene_id                 = as.character(gene_id),
    KOG_Annotation          = clean_annotation_text(KOG_Annotation),
    `GO:BiologicalProcess`  = extract_go_text(`GO:BiologicalProcess`),
    `GO:CellularComponent`  = extract_go_text(`GO:CellularComponent`),
    `GO:MolecularFunction`  = extract_go_text(`GO:MolecularFunction`),
    KEGG_Ortholog           = clean_annotation_text(KEGG_Ortholog),
    NR_Annotation           = clean_annotation_text(NR_Annotation),
    Swissprot_Annotation    = clean_annotation_text(Swissprot_Annotation)
  ) |>
  dplyr::distinct(gene_id, .keep_all = TRUE) |>
  dplyr::mutate(
    Gene_Name = extract_gene_name(
      Swissprot_Annotation,
      NR_Annotation,
      KOG_Annotation,
      KEGG_Ortholog
    ),
    Functional_Summary = dplyr::coalesce(
      Swissprot_Annotation,
      NR_Annotation,
      KOG_Annotation,
      KEGG_Ortholog
    )
  )

readr::write_csv(
  annot_clean,
  file.path(out_dir, "Annotation_table_cleaned.csv")
)

##############################
## PARTIE 5) IMPORT DES FICHIERS D'INTERET
##############################
cat("---- IMPORT DES FICHIERS D'INTERET ----\n")

files_to_import <- list(
  Hub_nodes      = hub_nodes_file,
  Tier1_both     = tier1_both_file,
  Tier1_both_Up  = tier1_both_up_file
)

input_tables <- purrr::imap(files_to_import, function(fp, nm) {
  if (is.na(fp) || !file.exists(fp)) return(NULL)
  
  df <- read_generic_table(fp)
  gene_col <- detect_gene_col(df)
  
  df |>
    dplyr::mutate(
      gene_id = as.character(.data[[gene_col]]),
      Source  = nm
    )
})

input_tables <- input_tables[!sapply(input_tables, is.null)]

if (length(input_tables) == 0) {
  stop("Aucun fichier d'entrée valide trouvé parmi hub_nodes_file / tier1_both_file / tier1_both_up_file.")
}

##############################
## PARTIE 6) STANDARDISATION DES TABLES
##############################
cat("---- STANDARDISATION ----\n")

standardize_table <- function(df) {
  cn <- colnames(df)
  
  ## colonnes optionnelles
  out <- df |>
    dplyr::mutate(
      ModuleColor = dplyr::if_else("ModuleColor" %in% cn, as.character(.data$ModuleColor), NA_character_),
      kME         = dplyr::if_else("kME_own" %in% cn, as.numeric(.data$kME_own),
                                   dplyr::if_else("kME" %in% cn, as.numeric(.data$kME), NA_real_)),
      log2FoldChange = dplyr::if_else("log2FoldChange" %in% cn, as.numeric(.data$log2FoldChange), NA_real_),
      padj        = dplyr::if_else("padj" %in% cn, as.numeric(.data$padj), NA_real_),
      Status      = dplyr::if_else("Status" %in% cn, as.character(.data$Status), NA_character_),
      Tier_Global = dplyr::if_else("Tier_Global" %in% cn, as.character(.data$Tier_Global), NA_character_),
      Tier_P      = dplyr::if_else("Tier_P" %in% cn, as.character(.data$Tier_P), NA_character_),
      Tier_N      = dplyr::if_else("Tier_N" %in% cn, as.character(.data$Tier_N), NA_character_),
      LOO_Score_P = dplyr::if_else("LOO_Score_P" %in% cn, as.numeric(.data$LOO_Score_P), NA_real_),
      LOO_Score_N = dplyr::if_else("LOO_Score_N" %in% cn, as.numeric(.data$LOO_Score_N), NA_real_)
    ) |>
    dplyr::select(
      gene_id, Source,
      ModuleColor, kME,
      log2FoldChange, padj, Status,
      Tier_Global, Tier_P, Tier_N,
      LOO_Score_P, LOO_Score_N,
      dplyr::everything()
    )
  
  out
}

input_std <- purrr::map(input_tables, standardize_table)
combined_input <- dplyr::bind_rows(input_std)

##############################
## PARTIE 7) FUSION DES SOURCES PAR GENE
##############################
cat("---- FUSION DES SOURCES ----\n")

## résumer les sources disponibles par gène
source_summary <- combined_input |>
  dplyr::group_by(gene_id) |>
  dplyr::summarise(
    Sources = paste(sort(unique(Source)), collapse = "; "),
    ModuleColor = dplyr::first(stats::na.omit(ModuleColor)),
    kME = dplyr::first(stats::na.omit(kME)),
    log2FoldChange = dplyr::first(stats::na.omit(log2FoldChange)),
    padj = dplyr::first(stats::na.omit(padj)),
    Status = dplyr::first(stats::na.omit(Status)),
    Tier_Global = dplyr::first(stats::na.omit(Tier_Global)),
    Tier_P = dplyr::first(stats::na.omit(Tier_P)),
    Tier_N = dplyr::first(stats::na.omit(Tier_N)),
    LOO_Score_P = dplyr::first(stats::na.omit(LOO_Score_P)),
    LOO_Score_N = dplyr::first(stats::na.omit(LOO_Score_N)),
    .groups = "drop"
  )

## petite fonction pour éviter les vecteurs vides
fix_empty <- function(x) ifelse(length(x) == 0, NA, x)

source_summary <- source_summary |>
  dplyr::rowwise() |>
  dplyr::mutate(
    ModuleColor = fix_empty(ModuleColor),
    kME = fix_empty(kME),
    log2FoldChange = fix_empty(log2FoldChange),
    padj = fix_empty(padj),
    Status = fix_empty(Status),
    Tier_Global = fix_empty(Tier_Global),
    Tier_P = fix_empty(Tier_P),
    Tier_N = fix_empty(Tier_N),
    LOO_Score_P = fix_empty(LOO_Score_P),
    LOO_Score_N = fix_empty(LOO_Score_N)
  ) |>
  dplyr::ungroup()

##############################
## PARTIE 8) AJOUT ANNOTATIONS BIOLOGIQUES
##############################
cat("---- AJOUT ANNOTATIONS ----\n")

final_df <- source_summary |>
  dplyr::left_join(
    annot_clean |>
      dplyr::select(
        gene_id,
        Gene_Name,
        Functional_Summary,
        `GO:BiologicalProcess`,
        `GO:MolecularFunction`,
        `GO:CellularComponent`,
        KEGG_Ortholog,
        NR_Annotation,
        Swissprot_Annotation,
        KOG_Annotation
      ),
    by = "gene_id"
  ) |>
  dplyr::arrange(
    dplyr::desc(!is.na(kME)),
    dplyr::desc(kME),
    padj
  )

##############################
## PARTIE 9) EXPORTS
##############################
cat("---- EXPORTS ----\n")

## table complète
readr::write_csv(
  final_df,
  file.path(out_dir, "Integrated_genes_with_biological_functions_FULL.csv")
)

## table compacte utile pour lecture rapide
final_short <- final_df |>
  dplyr::select(
    gene_id,
    Gene_Name,
    Sources,
    ModuleColor,
    kME,
    log2FoldChange,
    padj,
    Status,
    Tier_Global,
    Tier_P,
    Tier_N,
    LOO_Score_P,
    LOO_Score_N,
    Functional_Summary,
    `GO:BiologicalProcess`,
    `GO:MolecularFunction`,
    `GO:CellularComponent`,
    KEGG_Ortholog,
    NR_Annotation,
    Swissprot_Annotation
  )

readr::write_csv(
  final_short,
  file.path(out_dir, "Integrated_genes_with_biological_functions_SHORT.csv")
)

## table encore plus courte pour slide / lecture rapide
final_slide <- final_df |>
  dplyr::transmute(
    Gene = gene_id,
    Gene_Name = Gene_Name,
    Module = ModuleColor,
    kME = kME,
    log2FC = log2FoldChange,
    Tier = Tier_Global,
    Source = Sources,
    Biological_Process = `GO:BiologicalProcess`,
    Molecular_Function = `GO:MolecularFunction`,
    KEGG = KEGG_Ortholog,
    SwissProt = Swissprot_Annotation
  )

readr::write_csv(
  final_slide,
  file.path(out_dir, "Integrated_genes_for_presentation.csv")
)

##############################
## PARTIE 10) RESUME
##############################
cat("---- RESUME ----\n")

summary_df <- data.frame(
  Metric = c(
    "Genes_total",
    "With_Gene_Name",
    "With_GO_BP",
    "With_GO_MF",
    "With_GO_CC",
    "With_KEGG",
    "With_NR",
    "With_SwissProt",
    "With_kME",
    "With_log2FC",
    "With_Tier_Global"
  ),
  Value = c(
    nrow(final_df),
    sum(!is.na(final_df$Gene_Name)),
    sum(!is.na(final_df$`GO:BiologicalProcess`)),
    sum(!is.na(final_df$`GO:MolecularFunction`)),
    sum(!is.na(final_df$`GO:CellularComponent`)),
    sum(!is.na(final_df$KEGG_Ortholog)),
    sum(!is.na(final_df$NR_Annotation)),
    sum(!is.na(final_df$Swissprot_Annotation)),
    sum(!is.na(final_df$kME)),
    sum(!is.na(final_df$log2FoldChange)),
    sum(!is.na(final_df$Tier_Global))
  )
)

readr::write_csv(
  summary_df,
  file.path(out_dir, "Integration_summary.csv")
)

print(summary_df)

cat("\n✅ Terminé\n")
cat("Fichiers principaux :\n")
cat(" - Integrated_genes_with_biological_functions_FULL.csv\n")
cat(" - Integrated_genes_with_biological_functions_SHORT.csv\n")
cat(" - Integrated_genes_for_presentation.csv\n")