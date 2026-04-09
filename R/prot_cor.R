#' Perform correlation analysis using proteomics data file
#' 
#' `prot_cor()` conducts correlation analysis using spearman, pearson or kendall
#'
#' @param samples a data frame of sample information, including "Sample_ID", "Group" and covariates
#' @param prot_quant_data  a data frame of proteomics quantification matrix with rows as genes and columns as samples
#' @param result_dir absolute path where the results should be stored
#' @param continuous_var continuous variable (this must be included in the sample_file)
#' @param method correlation analysis method, "spearman" or "pearson", Default="spearman"
#' 
#' @return a data frame of correlation analysis result 
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
#' # prot_cor(samples, prot_quant_data, result_dir, continuous_var=continuous_var, method=method)
#' @export
prot_cor <- function(samples, prot_quant_data, result_dir, continuous_var, 
                       method = "spearman") {
  cat("Running protein correlation analysis using", method, "correlation\n")

  # samples is the sample table
  # prot_quant_data is the protein quantification table 
  
  # Step 1: Extract the continuous variable from the sample file and ensure it is numeric
  if (!(continuous_var %in% colnames(samples))) {
    stop(paste("The continuous variable", continuous_var, "was not found in the sample file."))
  }
  cont_var <- samples[[continuous_var]]
  if (!is.numeric(cont_var)) {
    cont_var <- as.numeric(as.character(cont_var))
  }
  
  # Step 2: assume quantification data is already normalized 
  quant_data_norm <- prot_quant_data
  
  # Step 3: Initialize results data frame
  n_proteins <- nrow(quant_data_norm)
  results <- data.frame(
    Gene  = rownames(quant_data_norm),
    Prot_Cor = numeric(n_proteins),
    Prot_Pval  = numeric(n_proteins),
    stringsAsFactors = FALSE
  )
  
  # Step 4: Compute the correlation for each protein
  for (i in seq_len(n_proteins)) {
    expr_values <- as.numeric(quant_data_norm[i, ])
    # only keep samples where both expr_values and cont_var are finite
    ok <- which(is.finite(expr_values) & is.finite(cont_var))
    # require at least 3 points and non-zero variance in both vectors
    if (length(ok) >= 3 &&
        sd(expr_values[ok]) > 0 &&
        sd(cont_var[ok]) > 0) {
      test <- cor.test(expr_values[ok], cont_var[ok], method = method, exact = FALSE)
      results$Prot_Cor[i]  <- unname(test$estimate)
      results$Prot_Pval[i] <- test$p.value
    } else {
      # not enough data or no variance → mark as NA
      results$Prot_Cor[i]  <- NA
      results$Prot_Pval[i] <- NA
    }
  }
  
  # Step 5: Adjust the p-values for multiple testing using Benjamini-Hochberg
  results$Prot_FDR <- p.adjust(results$Prot_Pval, method = "BH")
  
  # Step 6: Write the results to a file
  res_file <- paste0(result_dir, "/rpMeta.DEG.protein.cor.", method, ".txt")
  write.table(results, file = res_file, sep = "\t", row.names = FALSE, quote = FALSE)
  
  return(results)
}

