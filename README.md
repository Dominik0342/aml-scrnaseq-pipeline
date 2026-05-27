# aml-scrnaseq-pipeline

A simple R pipeline for cell type classification of bone marrow scRNA-seq data, based on [BoneMarrowMap](https://github.com/andygxzeng/BoneMarrowMap).

Developed for the analysis of acute myeloid leukemia (AML) datasets as part of an Engineering Project at UTS.

## What it does

Given a single-cell count matrix, the pipeline:

1. Filters low-quality cells based on gene counts and mitochondrial content
2. Detects and removes doublets using `scDblFinder`
3. Normalises the data with Seurat
4. Projects each cell onto the BoneMarrowMap healthy bone marrow reference
5. Predicts a cell type label and a pseudotime value for each cell
6. Generates figures (UMAP, pseudotime, composition, marker dotplot) and CSV outputs

## Requirements

- R version 4.5 or higher
- At least 16 GB of RAM (BoneMarrowMap reference is ~5 GB)
- About 10 GB of free disk space

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/aml-scrnaseq-pipeline.git
cd aml-scrnaseq-pipeline
Rscript install_packages.R
```

This installs Seurat, scDblFinder, BoneMarrowMap and the other required packages. It can take 15-30 minutes the first time.

## How to use

1. Open `pipeline.R` and edit the parameters at the top of the file:

   - `INPUT_FILE`: path to your count matrix (`.dem.txt`, `.tsv`, `.csv` or `.rds`)
   - `SAMPLE_ID`: an identifier for your sample
   - `OUTPUT_DIR`: where to save the results
   - QC thresholds and other settings if needed

2. Run the script:

   ```bash
   Rscript pipeline.R
   ```

3. Check the outputs in `OUTPUT_DIR`:

   - `figures/`: PNG files (QC, UMAP, pseudotime, composition, marker dotplot)
   - `csv/`: tabular outputs (cell-by-cell predictions, composition summary, QC stats)
   - `rds/`: Seurat objects for further analysis

The first run will automatically download the BoneMarrowMap reference (about 5 GB).

## Input format

The pipeline accepts:

- **Text matrices** (`.dem.txt`, `.tsv`, `.csv`) with genes in the first column and cells in the header row. Values must be raw integer counts.
- **Seurat objects** saved in `.rds` format.

The format is detected from the file extension.

## Outputs

| File | Description |
|------|-------------|
| `figures/01_QC_violins.png` | QC metrics before filtering |
| `figures/02_passfail.png` | Pass/Fail mapping summary |
| `figures/03_UMAP_celltype.png` | UMAP coloured by predicted cell type |
| `figures/04_UMAP_pseudotime.png` | UMAP coloured by pseudotime |
| `figures/05_composition.png` | Cell type composition |
| `figures/06_pseudotime_violin.png` | Pseudotime per cell type |
| `figures/07_marker_dotplot.png` | Canonical marker expression |
| `csv/cell_predictions.csv` | Per-cell predictions |
| `csv/composition_summary.csv` | Cell type counts and percentages |
| `csv/qc_summary.csv` | QC summary statistics |

## Author

Dominik Souillac, supervised by Dr. Dominik Beck.
University of Technology Sydney, 2026.

