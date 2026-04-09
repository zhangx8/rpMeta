# Example user-defined DE function: wilcoxon_rank_sum_test
wilcoxon_rank_sum_test <- function(samples, count_data, group1, group2, covariates=NULL, normalization="default") {
    # Assume the 'samples' data.frame has a column 'Group' containing group labels
    group_labels <- factor(samples$Group)

    if (normalization == "default") {
        count_data_norm <- count_data
    }

    n_genes <- nrow(count_data_norm)
    results <- data.frame(gene_name = rownames(count_data_norm),
                          log_fc = numeric(n_genes),
                          p_value = numeric(n_genes),
                          FDR = numeric(n_genes),
                          stringsAsFactors = FALSE)

    # Loop over each gene (each row in count_data_norm)
    for (i in seq_len(n_genes)) {
        gene_expr <- as.numeric(count_data_norm[i, ])
        # Subset expression values for the two groups
        x <- gene_expr[group_labels == group1]
        y <- gene_expr[group_labels == group2]

        # Perform Wilcoxon rank sum test (Mann-Whitney U test)
        wt <- wilcox.test(x, y, exact = FALSE)
        p_val <- wt$p.value

        # Compute log fold change using medians (adding 1 to avoid division by zero)
        median1 <- median(x, na.rm = TRUE)
        median2 <- median(y, na.rm = TRUE)
        log_fc_value <- log2((median2 + 1) / (median1 + 1))

        results$log_fc[i] <- log_fc_value
        results$p_value[i] <- p_val
    }

    # Adjust p-values using FDR method
    results$FDR <- p.adjust(results$p_value, method = "fdr")

    return(results)
}

