# =============================================================
# AML scRNA-seq pipeline - main script
# Author: Dominik Souillac
# Supervisor: Dr. Dominik Beck 
# =============================================================
# Usage: Rscript pipeline.R
# Edit the parameters in the section below, then run.
# =============================================================


# =============================================================
# PARAMETERS 
# =============================================================

# Input file: .dem.txt / .tsv / .csv (genes in rows, cells in columns)
# or .rds (a pre-built Seurat object)
INPUT_FILE  <- "data/my_sample.dem.txt"

# Identifier used in output filenames and figures
SAMPLE_ID   <- "MySample"

# Where to write the outputs
OUTPUT_DIR  <- "output/"

# Quality control thresholds
MIN_NFEATURE   <- 200    # min genes per cell
MAX_NFEATURE   <- 6000   # max genes per cell
MAX_PERCENT_MT <- 10     # max mitochondrial %

# Doublet detection (TRUE / FALSE)
REMOVE_DOUBLETS <- TRUE

# BoneMarrowMap mapping error threshold (Pass/Fail)
MAD_THRESHOLD  <- 2.5

# Paths to the BoneMarrowMap reference files
# These will be downloaded automatically if missing
REFERENCE_RDS  <- "reference/BoneMarrowMap_SymphonyReference.rds"
REFERENCE_UWOT <- "reference/BoneMarrowMap_uwot_model.uwot"


# =============================================================
# PIPELINE - no need to edit below
# =============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BoneMarrowMap)
  library(viridis)
})

set.seed(42)

cat("\n==========================================\n")
cat("  AML scRNA-seq pipeline\n")
cat("  Sample:", SAMPLE_ID, "\n")
cat("==========================================\n\n")


# --- 0. Create output directories ---
fig_dir <- file.path(OUTPUT_DIR, "figures")
csv_dir <- file.path(OUTPUT_DIR, "csv")
rds_dir <- file.path(OUTPUT_DIR, "rds")
for (d in c(fig_dir, csv_dir, rds_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# --- 1. Load the input data ---
cat(">>> Loading input...\n")
ext <- tools::file_ext(INPUT_FILE)

if (ext == "rds") {
  seu <- readRDS(INPUT_FILE)
  if (!inherits(seu, "Seurat")) stop("RDS file is not a Seurat object")
} else {
  sep <- ifelse(ext == "csv", ",", "\t")
  mat <- as.matrix(read.delim(INPUT_FILE, row.names = 1,
                              sep = sep, check.names = FALSE))
  storage.mode(mat) <- "numeric"
  seu <- CreateSeuratObject(counts = mat, project = SAMPLE_ID,
                            min.cells = 3, min.features = 0)
}
seu$sample_id <- SAMPLE_ID
cat("  Loaded:", ncol(seu), "cells x", nrow(seu), "genes\n")


# --- 2. Quality control ---
cat("\n>>> Quality control...\n")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Figure: nFeature_RNA histogram with QC thresholds
qc_data <- seu@meta.data
n_in  <- sum(qc_data$nFeature_RNA >= MIN_NFEATURE & qc_data$nFeature_RNA <= MAX_NFEATURE)
n_out <- sum(qc_data$nFeature_RNA <  MIN_NFEATURE | qc_data$nFeature_RNA >  MAX_NFEATURE)

p_qc <- ggplot(qc_data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 50, fill = "#3498DB", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = MIN_NFEATURE, color = "#E74C3C", linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = MAX_NFEATURE, color = "#E74C3C", linewidth = 1, linetype = "dashed") +
  annotate("text", x = MIN_NFEATURE, y = Inf, label = paste0(" min = ", MIN_NFEATURE),
           vjust = 1.5, hjust = 0, color = "#E74C3C", fontface = "bold", size = 4) +
  annotate("text", x = MAX_NFEATURE, y = Inf, label = paste0("max = ", MAX_NFEATURE, " "),
           vjust = 1.5, hjust = 1, color = "#E74C3C", fontface = "bold", size = 4) +
  labs(title = paste0(SAMPLE_ID, " - QC on nFeature_RNA before filtering"),
       subtitle = paste0("Kept: ", n_in, " cells   |   Filtered out: ", n_out, " cells"),
       x = "Number of genes detected per cell",
       y = "Number of cells") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(color = "gray30", hjust = 0.5))

ggsave(file.path(fig_dir, "01_QC_nFeature.png"),
       p_qc, width = 9, height = 5, dpi = 200, bg = "white")

# Filtering
n_before <- ncol(seu)
seu <- subset(seu, subset = nFeature_RNA > MIN_NFEATURE &
                            nFeature_RNA < MAX_NFEATURE &
                            percent.mt < MAX_PERCENT_MT)
cat("  Kept", ncol(seu), "/", n_before, "cells after QC\n")


# --- 3. Doublet detection ---
cat("\n>>> Doublet detection (scDblFinder)...\n")
sce <- as.SingleCellExperiment(seu)
sce <- scDblFinder(sce)
seu$doublet_score <- sce$scDblFinder.score
seu$is_doublet   <- sce$scDblFinder.class == "doublet"
n_dbl <- sum(seu$is_doublet)
cat("  Detected:", n_dbl, "doublets (",
    round(100 * n_dbl / ncol(seu), 1), "%)\n")

if (REMOVE_DOUBLETS) {
  seu <- subset(seu, cells = colnames(seu)[!seu$is_doublet])
  cat("  After removal:", ncol(seu), "cells\n")
}


# --- 4. Normalisation ---
cat("\n>>> Normalisation...\n")
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)


# --- 5. BoneMarrowMap reference mapping ---
cat("\n>>> BoneMarrowMap reference mapping...\n")

# Download reference if missing
if (!dir.exists("reference")) dir.create("reference")
if (!file.exists(REFERENCE_RDS)) {
  cat("  Downloading BoneMarrowMap reference (~5 GB, may take a few minutes)...\n")
  download.file(
    "https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_SymphonyReference.rds",
    REFERENCE_RDS, mode = "wb")
}
if (!file.exists(REFERENCE_UWOT)) {
  cat("  Downloading UWOT model...\n")
  download.file(
    "https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrow_RefMap_uwot_model.uwot",
    REFERENCE_UWOT, mode = "wb")
}

ref <- readRDS(REFERENCE_RDS)
ref$save_uwot_path <- REFERENCE_UWOT

query_mapped <- map_Query(
  exp_query      = GetAssayData(seu, assay = "RNA", layer = "counts"),
  metadata_query = seu@meta.data,
  ref_obj        = ref,
  vars           = "sample_id")

query_mapped <- query_mapped %>%
  calculate_MappingError(reference          = ref,
                         MAD_threshold      = MAD_THRESHOLD,
                         threshold_by_donor = TRUE,
                         donor_key          = "sample_id")

query_mapped <- predict_CellTypes(query_obj   = query_mapped,
                                   ref_obj     = ref,
                                   final_label = "predicted_CellType")
query_mapped <- predict_Pseudotime(query_obj   = query_mapped,
                                    ref_obj     = ref,
                                    final_label = "predicted_Pseudotime")

# Pass/Fail
n_pass <- sum(query_mapped$mapping_error_QC == "Pass")
n_fail <- sum(query_mapped$mapping_error_QC == "Fail")
cat("  Pass:", n_pass, "(", round(100 * n_pass / ncol(query_mapped), 1), "%)\n")
cat("  Fail:", n_fail, "\n")

# Keep only Pass cells for downstream figures
pass_cells <- colnames(query_mapped)[query_mapped$mapping_error_QC == "Pass"]
query_filtered <- subset(query_mapped, cells = pass_cells)


# --- 6. Save RDS objects ---
saveRDS(query_mapped,
        file.path(rds_dir, paste0("query_mapped_", SAMPLE_ID, ".rds")))
saveRDS(query_filtered,
        file.path(rds_dir, paste0("query_filtered_", SAMPLE_ID, ".rds")))


# --- 7. Figures ---
cat("\n>>> Generating figures...\n")

# Pass/Fail summary
pf_df <- data.frame(status = c("Pass", "Fail"), n = c(n_pass, n_fail))
p_pf <- ggplot(pf_df, aes(x = "", y = n, fill = status)) +
  geom_bar(stat = "identity", width = 0.5, color = "white") +
  geom_text(aes(label = paste0(n, " cells")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold") +
  scale_fill_manual(values = c("Pass" = "#27AE60", "Fail" = "#E74C3C")) +
  theme_classic(base_size = 12) +
  labs(title = paste0(SAMPLE_ID, " - mapping QC"),
       y = "Number of cells", x = "")
ggsave(file.path(fig_dir, "02_passfail.png"),
       p_pf, width = 6, height = 6, dpi = 200, bg = "white")

# UMAP by cell type
p_umap_ct <- DimPlot(query_filtered, reduction = "umap_projected",
                     group.by = "predicted_CellType",
                     label = TRUE, repel = TRUE, label.size = 2.5) +
  NoAxes() + NoLegend() +
  ggtitle(paste0(SAMPLE_ID, " - predicted cell type"))
ggsave(file.path(fig_dir, "03_UMAP_celltype.png"),
       p_umap_ct, width = 11, height = 8, dpi = 300, bg = "white")

# UMAP by pseudotime
p_umap_pt <- FeaturePlot(query_filtered, features = "predicted_Pseudotime",
                         reduction = "umap_projected") +
  scale_color_viridis_c(option = "plasma", name = "Pseudotime") +
  NoAxes() +
  ggtitle(paste0(SAMPLE_ID, " - predicted pseudotime"))
ggsave(file.path(fig_dir, "04_UMAP_pseudotime.png"),
       p_umap_pt, width = 11, height = 8, dpi = 300, bg = "white")

# Composition barplot
comp <- as.data.frame(table(query_filtered$predicted_CellType))
colnames(comp) <- c("cell_type", "n")
comp$pct <- 100 * comp$n / sum(comp$n)
comp <- comp[order(-comp$n), ]

p_comp <- ggplot(comp, aes(x = reorder(cell_type, n), y = pct,
                            fill = cell_type)) +
  geom_col(color = "black", linewidth = 0.2) +
  geom_text(aes(label = paste0(round(pct, 1), "% (n=", n, ")")),
            hjust = -0.1, size = 2.5) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  labs(title = paste0(SAMPLE_ID, " - cell type composition (Pass cells only)"),
       x = "", y = "Frequency (%)")
ggsave(file.path(fig_dir, "05_composition.png"),
       p_comp, width = 9, height = 8, dpi = 300, bg = "white")

# Pseudotime violin per cell type
pt_data <- query_filtered@meta.data
pt_data <- pt_data[!is.na(pt_data$predicted_Pseudotime) &
                   is.finite(pt_data$predicted_Pseudotime), ]
ct_order <- pt_data %>%
  group_by(predicted_CellType) %>%
  summarise(med = median(predicted_Pseudotime), .groups = "drop") %>%
  arrange(med) %>% pull(predicted_CellType)
pt_data$predicted_CellType <- factor(pt_data$predicted_CellType,
                                      levels = ct_order)

p_pt <- ggplot(pt_data, aes(x = predicted_CellType, y = predicted_Pseudotime,
                             fill = predicted_CellType)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, alpha = 0.9, outlier.size = 0.3) +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = paste0(SAMPLE_ID, " - pseudotime by cell type"),
       y = "Predicted pseudotime", x = "")
ggsave(file.path(fig_dir, "06_pseudotime_violin.png"),
       p_pt, width = 13, height = 6, dpi = 300, bg = "white")

# Marker dotplot (canonical hematopoietic markers)
markers <- c("SPINK2", "MEIS1", "CRHBP", "CD34", "MPO", "ELANE", "PRTN3",
             "LYZ", "S100A9", "CD14", "VCAN", "CLEC9A", "LILRA4",
             "GATA1", "HBB", "HBA1", "MS4A1", "CD79A", "MZB1",
             "CD3D", "CD3E", "CD8A", "NKG7", "GNLY")
markers_present <- intersect(markers, rownames(query_filtered))

if (length(markers_present) >= 3) {
  Idents(query_filtered) <- "predicted_CellType"
  p_dot <- DotPlot(query_filtered, features = markers_present,
                   cols = c("lightgrey", "#B71C1C"), dot.scale = 5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8)) +
    labs(title = paste0(SAMPLE_ID, " - canonical marker expression"),
         x = "Marker", y = "Predicted cell type")
  ggsave(file.path(fig_dir, "07_marker_dotplot.png"),
         p_dot, width = 13, height = 9, dpi = 300, bg = "white")
}


# --- 8. CSV outputs ---
cat("\n>>> Writing CSV outputs...\n")

# Cell-by-cell predictions
cell_pred <- data.frame(
  cell_barcode        = colnames(query_mapped),
  predicted_CellType  = query_mapped$predicted_CellType,
  predicted_Pseudotime = query_mapped$predicted_Pseudotime,
  mapping_error_QC    = query_mapped$mapping_error_QC,
  mapping_error_score = query_mapped$mapping_error_score,
  nFeature_RNA        = query_mapped$nFeature_RNA,
  nCount_RNA          = query_mapped$nCount_RNA,
  percent_mt          = query_mapped$percent.mt,
  is_doublet          = query_mapped$is_doublet
)
write.csv(cell_pred, file.path(csv_dir, "cell_predictions.csv"),
          row.names = FALSE)

# Composition summary
write.csv(comp, file.path(csv_dir, "composition_summary.csv"),
          row.names = FALSE)

# QC summary
qc_summary <- data.frame(
  sample_id        = SAMPLE_ID,
  n_cells_total    = ncol(query_mapped),
  n_cells_pass     = n_pass,
  n_cells_fail     = n_fail,
  pct_pass         = round(100 * n_pass / ncol(query_mapped), 2),
  median_nFeature  = median(query_mapped$nFeature_RNA),
  median_percent_mt = round(median(query_mapped$percent.mt), 2),
  n_doublets       = n_dbl
)
write.csv(qc_summary, file.path(csv_dir, "qc_summary.csv"),
          row.names = FALSE)


# --- Done ---
cat("\n==========================================\n")
cat("  Pipeline completed successfully\n")
cat("  Outputs:", normalizePath(OUTPUT_DIR), "\n")
cat("==========================================\n\n")
