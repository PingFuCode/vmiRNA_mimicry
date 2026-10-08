# vmiRNA_mimicry

R and Python scripts, together with figure source data, for computational analyses of viral and human microRNA (miRNA) mimicry.

The repository covers sequence similarity and permutation analyses, miRNA conservation, tissue expression, gene expression and pathway analyses, and cancer-associated analyses. Files are organized by figure number.

## Repository contents

This repository contains **16 R scripts**, **4 Python scripts**, and **`Source_Data.xlsx`**. The scripts are figure-specific research scripts; they do not form a single automated pipeline.

| File | Analysis or figure content |
| --- | --- |
| `Figure1B.R` | Distribution of miRNAs across viruses and associated figure annotations. |
| `Figure1BC.py` | Sequence-comparison statistics and permutation analyses. |
| `Figure2A.R` | Bar plots summarizing unique and non-unique seed-region counts; plotting values are embedded in the script. |
| `Figure 2A_uniquevmiRNA.py` | Summary of unique and shared seed regions in viral miRNA sequence collections. |
| `Figure2B.R` | Cross-species seed-sequence similarity visualizations using Jaccard indices. |
| `Figure2C.R` | Family-level comparisons of seed-sequence conservation. |
| `Figure 2E.R` | Heatmap of representative biological-process enrichment terms. |
| `Figure 2F.R` | Chromosomal distribution and enrichment analyses. |
| `Figure3A.R` | Tissue-expression heatmap and tissue annotations. |
| `Figure3B.R` | Comparisons of miRNA expression across tissue groups. |
| `Figure3D.R` | Expression plots and export of paired miRNA plotting data. |
| `Figure 4A.R` | Visualization of gene and protein interaction networks. |
| `Figure 4B.R` | Differential-expression analyses of public transcriptomic datasets. |
| `Figure 4C.py` | Permutation analysis of gene-expression direction. |
| `Figure 4D.R` | Pathway enrichment analysis and enrichment-network visualization. |
| `Figure 5A.R` | Disease-enrichment summary heatmap. |
| `Figure 5B.R` | Analysis of cancer-related gene categories. |
| `Figure 5C.R` | Genomic overlap analyses involving deletion regions, fragile regions, and gene annotations. |
| `Figure 5D.py` | Exploratory classification of cancer-stage groups using miRNA expression, with ROC/AUC summaries. |
| `Figure 5EF.R` | Differential-expression analyses of miRNAs and genes in gastric-cancer datasets. |
| `Source_Data.xlsx` | Figure-level source-data workbook with 23 worksheets. |

File names retain their original capitalization and spaces. A script can contain several analysis or plotting sections, so script names and workbook sheets do not always correspond one-to-one.

## Source data

[`Source_Data.xlsx`](Source_Data.xlsx) contains the following worksheets:

| Figure group | Worksheets |
| --- | --- |
| Figure 1 | `Fig1b`, `Fig1c`, `Fig1d` |
| Figure 2 | `Fig2a`, `Fig2b`, `Fig2c`, `Fig2e`, `Fig2f` |
| Figure 3 | `Fig3a`, `Fig3b`, `Fig3c`, `Fig3d`, `Fig3e` |
| Figure 4 | `Fig4a`, `Fig4b`, `Fig4c`, `Fig4d` |
| Figure 5 | `Fig5a`, `Fig5b`, `Fig5c`, `Fig5d`, `Fig5e`, `Fig5f` |

Worksheets include descriptive text and figure data. Inspect each sheet before importing it: descriptive rows and multiple table blocks may require selecting a specific range or skipping rows.

The workbook is a figure-data resource. Scripts generally read separate CSV, TSV, FASTA, Excel, expression-matrix, or genomic-annotation files from the original analysis directories; they do not automatically read their inputs from this workbook. Those external input files and intermediate objects are not all included here.

## Software dependencies

### Python

The Python scripts use Python 3 and the following third-party packages:

- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `scikit-learn`

A basic installation command is:

```sh
python -m pip install numpy pandas scipy matplotlib scikit-learn
```

Python standard-library modules used by the scripts require no separate installation.

### R

R dependencies vary by script. Package references in the scripts include:

| Purpose | Packages |
| --- | --- |
| Data handling | `tidyverse`, `data.table`, `readxl`, `reshape2` |
| Plotting and networks | `ggpubr`, `ggrepel`, `ggsci`, `gridExtra`, `patchwork`, `RColorBrewer`, `circlize`, `ComplexHeatmap`, `pheatmap`, `viridisLite`, `igraph`, `tidygraph`, `ggraph`, `ggVennDiagram`, `RIdeogram` |
| Statistical tests | `rstatix`, `lawstat` |
| Expression and enrichment analysis | `limma`, `edgeR`, `DESeq2`, `affy`, `clusterProfiler`, `fgsea`, `rrvgo` |
| Genomic annotations and data access | `AnnotationDbi`, `org.Hs.eg.db`, `GenomicRanges`, `GenomeInfoDb`, `ChIPseeker`, `karyoploteR`, `GEOquery`, `TCGAbiolinks`, `TxDb.Hsapiens.UCSC.hg19.knownGene`, `TxDb.Hsapiens.UCSC.hg38.knownGene` |

The scripts also reference individual tidyverse packages such as `dplyr`, `tidyr`, `readr`, `stringr`, `tibble`, `purrr`, and `ggplot2`, as well as `scales`. Base R packages such as `grid` do not need separate installation.

Install CRAN packages with `install.packages()` and Bioconductor packages with `BiocManager::install()`. Install the packages used by the particular script you intend to inspect or run. Some scripts use functions from packages loaded earlier in an interactive R session, so their initial `library()` statements may not be sufficient on their own.

Exact software versions and a dependency lockfile are not included in this repository.

## Working with the scripts

1. Clone or download the repository:

   ```sh
   git clone https://github.com/PingFuCode/vmiRNA_mimicry.git
   cd vmiRNA_mimicry
   ```

2. Choose the script for the figure or analysis of interest and inspect its data-loading and output sections.
3. Prepare the referenced input files and any intermediate objects. Check the expected column names, identifiers, worksheet names, and genome assembly where applicable.
4. Replace machine-specific paths, including `setwd()` calls and absolute `D:/...` or `E:/...` paths, with paths on your machine. Check output locations and create required directories.
5. Install and load the relevant dependencies, then run the required sections in R/RStudio or Python. Quote filenames that contain spaces.

