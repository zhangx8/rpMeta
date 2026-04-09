#' Perform correlation analysis using RNA count matrix 
#' 
#' `rna_cor()` conducts correlation analysis between RNA expression and a continuous variable
#'
#' @param samples a data frame of sample information, including "Sample_ID", "Group" and covariates
#' @param rna_count_data  a data frame of RNASeq raw count matrix with rows as genes and columns as samples
#' @param result_dir absolute path where the results should be stored
#' @param continuous_var continuous variable (this must be included in the sample_file)
#' @param method correlation analysis method, "spearman" or "pearson", Default="spearman"
#' @param normalization normalization method, Default="default", if the value is "as_is", the data is already normalized 
#' @param mean_rna_count mean RNA count cutoff (per row/gene), Default=0
#' 
#' @return a data frame of correlation anlaysis result 
#' @examples
#' sample_file <- system.file("extdata", "apollo1.sample.info.tsv", package = "rpMeta")
#' rna_count_file <- system.file("extdata", "apollo1.rna.count.top6k.tsv", package = "rpMeta")
#' protein_quant_file <- system.file("extdata", "apollo1.protein.quant.top3k.tsv", package = "rpMeta")
#' inputs <- loadRpInput(sample_file, rna_count_file, protein_quant_file, format="tsv")
#' samples <- inputs$samples
#' rna_count_data <- inputs$rna
#' prot_quant_data <- inputs$protein
#' continuous_var <- "PkYrs"
#' method <- "spearman"
#' rna_norm <- "default"
#' result_dir <- file.path(getwd(), "rpMeta_example1")
#' if (!dir.exists(result_dir)) { dir.create(result_dir) }
#' # rna_cor(samples, rna_count_data, result_dir, continuous_var=continuous_var, method=method, rna_norm)
#' @export
rna_cor <- function(samples, rna_count_data, result_dir, continuous_var, 
                    method = "spearman", normalization = "default", 
                    mean_rna_count = 0) {
  cat("Running RNA correlation analysis using", method, "correlation\n")
 
  # Step 1: filter rna count matrix using "mean_rna_count" argument
  # "samples" is the data frame for sample information
  if (normalization == "as_is") {
    count_data <- rna_count_data
  } else {
    print(paste0("mean_rna_count cutoff: ", as.character(mean_rna_count)))
    if (mean_rna_count == 0) {
        count_data <- rna_count_data
    } else {
        rowmeans <- rowMeans(rna_count_data, na.rm = TRUE)
        inds <- which(rowmeans < mean_rna_count)
        count_data <- rna_count_data[-inds, ]
    }
  }
  
  # Step 2: Check if the specified continuous variable column exists in the samples table
  if (!(continuous_var %in% colnames(samples))) {
    stop(paste("The continuous variable", continuous_var, "was not found in the sample file."))
  }
  
  # Extract the continuous variable, and ensure it is numeric
  cont_var <- samples[[continuous_var]]
  if (!is.numeric(cont_var)) {
    cont_var <- as.numeric(as.character(cont_var))
  }
  
  # Step 3: Normalize the count data
  # If normalization is "local", use rna_normalize function (median centered log2(cpm+1))
  # Otherwise, perform a simple log2 transformation (plus one)
  if (normalization == "local") {
    count_data_norm <- rna_normalize(count_data)
  } else if (normalization == "as_is") {
    count_data_norm <- count_data
  } else {
    count_data_norm <- log2(count_data + 1)
  }
  
  # Initialize the results data frame
  n_genes <- nrow(count_data_norm)
  results <- data.frame(
    Gene      = rownames(count_data_norm),
    RNA_Cor  = numeric(n_genes),
    RNA_Pval  = numeric(n_genes),
    stringsAsFactors = FALSE
  )
  
  # compute correlation 
  for (i in seq_len(n_genes)) {
    gene_expr <- as.numeric(count_data_norm[i, ])
    # keep only the samples where both are finite
    ok <- which(is.finite(gene_expr) & is.finite(cont_var))
    if (length(ok) >= 3 && sd(gene_expr[ok]) > 0 && sd(cont_var[ok]) > 0) {
      test <- cor.test(gene_expr[ok], cont_var[ok], method = method, exact = FALSE)
      results$RNA_Cor[i]  <- unname(test$estimate)
      results$RNA_Pval[i] <- test$p.value
    } else {
      # not enough data or no variance, return NA
      results$RNA_Cor[i]  <- NA
      results$RNA_Pval[i] <- NA
    }
  }
  
  # Step 4: Adjust the p-values for multiple testing (Benjamini-Hochberg)
  results$RNA_FDR <- p.adjust(results$RNA_Pval, method = "BH")
  
  # Step 5: Write the results to a file and return the data frame.
  res_file <- paste0(result_dir, "/rpMeta.DEG.RNA.cor.", method, ".txt")
  write.table(results, file = res_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(results)
}

# CPM normalization 
rna_normalize <- function(count_data) {
    # Step 1: Calculate CPM (Counts Per Million)
    # Total counts per sample (column sums)
    total_counts_per_sample <- colSums(count_data)

    # Calculate CPM: count / total_counts * 1e6
    cpm_counts <- sweep(count_data, 2, total_counts_per_sample, "/") * 1e6

    # Step 2: Log2 transformation (adding a pseudocount to avoid log(0))
    log2_cpm_counts <- log2(cpm_counts + 1)

    # Step 3: Median centering for each gene
    median_centered_counts <- sweep(log2_cpm_counts, 1, apply(log2_cpm_counts, 1, median), "-")

    # Return the normalized matrix
    return(median_centered_counts)
}
