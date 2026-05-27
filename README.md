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
