# Thac_ECM_Recount3
Exon–Gene Correlation Analysis (TCGA Recount3)

1) Overview

- This R script performs exon–gene expression correlation analysis across multiple TCGA cancer types using the recount3.
It identifies genes whose expression levels are correlated with a selected target exon.

- The workflow automatically retrieves exon-level and gene-level expression data, computes RPKM values, and calculates Spearman correlations between the target exon and all protein-coding genes for each TCGA project.

2) Features

- Automatic data retrieval: Uses recount3::create_rse_manual() to access TCGA expression data.

- RPKM normalization: Computes exon and gene expression as RPKM for cross-sample comparability.

- Correlation analysis: Calculates Spearman correlation between one selected exon and all genes.

- Pan-cancer coverage: Runs across multiple TCGA cancer types in a single execution.

- Structured outputs: Produces one correlation CSV file per TCGA project.

- User-friendly configuration: All key parameters (target gene, output directory, projects) are set in one section.

3) Requirements
- R Version: R ≥ 4.0.0 

- R Packages:

|Package| Purpose|
| ------------------------ | ---------------------------------------- |
| **recount**              | RPKM normalization and count processing  |
| **recount3**             | Retrieval of TCGA RSE data               |
| **GenomicRanges**        | Genomic coordinate handling              |
| **SummarizedExperiment** | Expression matrix structure              |
| **pheatmap**             | (Optional) visualization of correlations |
| **dplyr**                | Data manipulation                        |
| **biomaRt**              | Query Ensembl for protein-coding genes   |

4) Installation

5) Inputs

- TCGA Project List
  - Example:
```r
prj <- c("KIRC", "BRCA", "LUAD", "LIHC", "PAAD", "SKCM", "ACC", "UCEC")
```

- Exon List
  - A list of exon coordinates (gene, chromosome, start, end, strand).
```r
exon_list <- c(
  "TNC_chr9_115048260_115048532_-", 
  "FN1_chr2_215392931_215393203_-",
  "POSTN_chr13_37574572_37574652_-",
  "COL6A3_chr2_237378636_237379235_-"
)
```
- Target Exon Gene
  - The gene prefix (in lower or upper case) of the exon to correlate.
```r
target_exon_gene <- "tnc"
```
- File Paths
  - Set paths for TCGA data and working directory.
```r
prj_home     <- "data_sources/tcga"
prj_organism <- "human"
work_dir     <- "/share/Lab_Shaw/ThacProject/ECM_proj_4/"
```
- Output Settings
  - Setting your output version and date for tracking:
```r
out_ver  <- "v1"
out_date <- "20251111"
```

6) Outputs
- For each TCGA project, the script generates one CSV file containing all exon–gene correlation results:
```r
<work_dir>/<target_exon_gene>_correlation_matrix/
```

- Example:
```r
  /share/Lab_Shaw/ThacProject/ECM_proj_4/tnc_correlation_matrix/KIRC_TNC_correlation_matrix_20251111_v1.csv
```

