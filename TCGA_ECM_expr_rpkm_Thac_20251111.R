################################################################################
#                                USER SETTINGS                                 #
################################################################################

# Choose which exon-gene to correlate (must match prefix in exon_list)
# Examples: "tnc", "fn1", "postn", "vcan"
target_exon_gene <- "tnc"

# Set output version and date
out_ver  <- "v1"
out_date <- "20251111"

# Set project list (TCGA cancer types to process)
prj <- c(
  "KIRC", "BLCA", "OV", "PRAD", "BRCA", "LUAD", "UCEC", "THCA",
  "LUSC", "HNSC", "COAD", "LGG", "SKCM", "STAD", "LIHC", "KIRP",
  "CESC", "SARC", "ESCA", "PCPG", "PAAD", "LAML", "READ", "GBM",
  "TGCT", "THYM", "KICH", "MESO", "UVM", "ACC", "UCS", "DLBC", "CHOL"
)

# Set paths and organism info
prj_home     <- "data_sources/tcga"
prj_organism <- "human"
work_dir     <- "/share/Lab_Shaw/ThacProject/ECM_proj_4/"

# Define the exon list (must include your target gene’s exon)
exon_list <- c(
  "TNC_chr9_115048260_115048532_-", 
  "FN1_chr2_215392931_215393203_-",
  "POSTN_chr13_37574572_37574652_-",
  "COL6A3_chr2_237378636_237379235_-",
  "VCAN_chr5_83519349_83522309_+",
  "CLSTN1_chr1_9756481_9756510_-",
  "NRCAM_chr7_108191254_108191283_-",
  "FYN_chr6_111699515_111699670_-",
  "PICALM_chr11_85990250_85990378_-"
)

################################################################################
#                              LOAD PACKAGES                                   #
################################################################################
library(recount)
library(recount3)
library(GenomicRanges)
library(SummarizedExperiment)
library(pheatmap)
library(dplyr)
library(biomaRt)

################################################################################
#                          FUNCTION DEFINITIONS                                #
################################################################################

# ---------------------- Exon-level RPKM --------------------------------------
exon_level_RPKM <- function(prj, prj_home, prj_organism, exon_list) {
  rse_exon <- create_rse_manual(
    project = prj, 
    project_home = prj_home, 
    organism = prj_organism, 
    annotation = "gencode_v26", 
    type = "exon"
  )
  
  get_exon_index <- function(exon_str, rse) {
    exon_info <- unlist(strsplit(exon_str, "_"))
    chr <- exon_info[2]; start_pos <- as.numeric(exon_info[3])
    end_pos <- as.numeric(exon_info[4]); strand <- exon_info[5]
    which(
      seqnames(rowRanges(rse)) == chr &
        start(rowRanges(rse)) == start_pos &
        end(rowRanges(rse)) == end_pos &
        strand(rowRanges(rse)) == strand
    )[1]
  }
  
  exon_indices <- sapply(exon_list, get_exon_index, rse = rse_exon)
  exon_indices <- exon_indices[!is.na(exon_indices)]
  rse_selected <- rse_exon[exon_indices, ]
  
  assay(rse_selected, "counts_2") <- transform_counts(
    rse = rse_selected, 
    by = c("mapped_reads"),
    L = 100, round = FALSE
  )
  
  exon_lengths <- width(rowRanges(rse_selected))
  counts_2 <- assay(rse_selected, "counts_2")
  mapped_reads <- rse_selected$recount_qc.star.all_mapped_reads
  names(mapped_reads) <- colnames(rse_selected)
  
  exon_rpkm_matrix <- sweep(counts_2, 1, exon_lengths, "/")
  exon_rpkm_matrix <- sweep(exon_rpkm_matrix, 2, mapped_reads, "/")
  exon_rpkm_matrix <- exon_rpkm_matrix * 1e9
  
  colnames(exon_rpkm_matrix) <- rse_selected$tcga.tcga_barcode
  rownames(exon_rpkm_matrix) <- paste0(exon_list[!is.na(exon_indices)], "_RPKM")
  exon_rpkm_matrix <- exon_rpkm_matrix[rowSums(exon_rpkm_matrix) != 0, ]
  return(exon_rpkm_matrix)
}

# ---------------------- Gene-level RPKM --------------------------------------
gene_level_RPKM <- function(prj, prj_home, prj_organism, protein_coding_symbols) {
  all_genes <- unique(protein_coding_symbols)
  
  rse_gene <- create_rse_manual(
    project = prj, project_home = prj_home, organism = prj_organism,
    annotation = "gencode_v26", type = "gene"
  )
  
  rse_gene <- rse_gene[rowData(rse_gene)$gene_type == "protein_coding", ]
  keep_genes <- rowData(rse_gene)$gene_name %in% all_genes
  rse_filtered <- rse_gene[keep_genes, ]
  
  assay(rse_filtered, "counts") <- transform_counts(rse_filtered)
  gene_rpkm_matrix <- recount::getRPKM(rse_filtered)
  
  rownames(gene_rpkm_matrix) <- rowData(rse_filtered)$gene_name
  colnames(gene_rpkm_matrix) <- colData(rse_filtered)$tcga.tcga_barcode
  gene_rpkm_matrix <- gene_rpkm_matrix[rowSums(gene_rpkm_matrix) != 0, ]
  return(gene_rpkm_matrix)
}

# ---------------------- Correlation Function ---------------------------------
run_exon_gene_correlation <- function(merged_matrix, target_exon_gene, work_dir,
                                      out_date, out_ver, project_name) {
  gene_tag_upper <- toupper(target_exon_gene)
  gene_tag_lower <- tolower(target_exon_gene)
  
  prefix_pattern <- paste0("^", gene_tag_upper, "_chr")
  exon_row_candidates <- rownames(merged_matrix)[grepl(prefix_pattern, rownames(merged_matrix))]
  
  if (length(exon_row_candidates) == 0) {
    stop(paste0("No exon row found for gene tag: ", gene_tag_upper))
  }
  if (length(exon_row_candidates) > 1) {
    message("Multiple exon rows found for ", gene_tag_upper,
            ". Using the first one: ", exon_row_candidates[1])
  }
  
  exon_row_name <- exon_row_candidates[1]
  exon_expr <- as.numeric(merged_matrix[exon_row_name, ])
  
  results <- data.frame(gene = character(),
                        spearman_correlation = numeric(),
                        p_value = numeric(),
                        cancer_type = character(),
                        stringsAsFactors = FALSE)
  
  for (gene_name in rownames(merged_matrix)) {
    if (gene_name == exon_row_name) next
    if (grepl("chr", gene_name, ignore.case = TRUE)) next
    gene_expr <- as.numeric(merged_matrix[gene_name, ])
    cor_test <- suppressWarnings(cor.test(exon_expr, gene_expr, method = "spearman"))
    
    results <- rbind(results, data.frame(
      gene = gene_name,
      spearman_correlation = as.numeric(cor_test$estimate),
      p_value = as.numeric(cor_test$p.value),
      cancer_type = project_name
    ))
  }
  
  results <- results %>% arrange(p_value)
  out_dir <- paste0(work_dir, gene_tag_lower, "_correlation_matrix/")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_path <- paste0(
    out_dir, project_name, "_", gene_tag_upper,
    "_correlation_matrix_", out_date, "_", out_ver, ".csv"
  )
  write.csv(results, out_path, row.names = FALSE)
  message("Finished ", project_name, " - ", gene_tag_upper, " correlation.")
}

################################################################################
#                                MAIN SCRIPT                                   #
################################################################################

# Retrieve all protein coding genes from Ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
               filters = "biotype", values = "protein_coding", mart = ensembl)
protein_coding_symbols <- genes$hgnc_symbol[genes$hgnc_symbol != ""]

# Loop through TCGA projects
for (p in prj) {
  message("Processing project: ", p)
  
  exon_rpkm_matrix <- exon_level_RPKM(p, prj_home, prj_organism, exon_list)
  gene_rpkm_matrix <- gene_level_RPKM(p, prj_home, prj_organism, protein_coding_symbols)
  
  common_cols <- intersect(colnames(exon_rpkm_matrix), colnames(gene_rpkm_matrix))
  exon_rpkm_matrix <- exon_rpkm_matrix[, common_cols, drop = FALSE]
  gene_rpkm_matrix <- gene_rpkm_matrix[, common_cols, drop = FALSE]
  
  merged_matrix <- rbind(exon_rpkm_matrix, gene_rpkm_matrix)
  
  run_exon_gene_correlation(
    merged_matrix, target_exon_gene, work_dir, out_date, out_ver, project_name = p
  )
}

