#' Perform differrential expression analysis using RNA count matrix 
#' 
#' `rna_de()` conducts differential expression analysis using limma or DESeq2.
#' Covariates can be added to the model 
#'
#' @param samples a data frame of sample information, including "Sample_ID", "Group" and covariates
#' @param rna_count_data  a data frame of RNASeq raw count matrix with rows as genes and columns as samples
#' @param result_dir absolute path where the results should be stored
#' @param group1 name of the first group used in DE analysis, a.k.a condition1 or factor1
#' @param group2 name of the second group used in DE analysis, a.k.a condition2 or factor2
#' fold change is calculated as group2/group1
#' @param method DE analysis method, limma, DESeq2 or user defined function, Default="limma"
#' @param covariates covariates used in the model (this must be included in the sample_file)
#' @param normalization normalization method, Default="default" 
#' @param X11 output X11 or PDF for graphics, c("No", "Yes"), (Default="No")
#' @param ... Additional arguments. If a named argument `mean_rna_count` (mean RNA count cutoff per row/gene) is supplied (e.g. 5)
#' 
#' @return a data frame of DE anlaysis result 
#' @examples
#' sample_file <- system.file("extdata", "apollo1.sample.info.tsv", package = "rpMeta")
#' rna_count_file <- system.file("extdata", "apollo1.rna.count.top6k.tsv", package = "rpMeta")
#' protein_quant_file <- system.file("extdata", "apollo1.protein.quant.top3k.tsv", package = "rpMeta")
#' inputs <- loadRpInput(sample_file, rna_count_file, protein_quant_file, format="tsv")
#' samples <- inputs$samples
#' rna_count_data <- inputs$rna
#' prot_quant_data <- inputs$protein
#' group1 <- "TRU"
#' group2 <- "ProxProlif"
#' covariates <- c("Sex", "grade")
#' rna_method <- "limma"
#' rna_norm <- "default"
#' result_dir <- file.path(getwd(), "rpMeta_example")
#' mean_rna_count <- 5 
#' DESeq2_test <- "Wald" 
#' if (!dir.exists(result_dir)) { dir.create(result_dir) }
#' rna_de(samples, rna_count_data, result_dir, group1, group2, method=rna_method, covariates, rna_norm, X11="Yes", mean_rna_count=mean_rna_count, DESeq2_test = DESeq2_test)
#' @export
rna_de <- function(samples, rna_count_data, result_dir, group1, group2, 
                   method="limma", covariates=NULL, normalization="default", X11="No", ...) {
    dotArgs <- list(...)
    print(dotArgs)
    if (X11 != "No" && X11 != "Yes") {
        print("X11 argument musct by No or Yes")
        return("X11_Error")
    }
    
    # Step 1: filter rna count matrix using "mean_rna_count" argument
    # "samples" is the data frame for sample information 
    if (normalization == "as_is") {
	count_data <- rna_count_data
        print("rna_normalization method is *as_is*, only limma can be used.")
        method="limma"
    } else {
        mean_rna_count = 0
        if ("mean_rna_count" %in% names(dotArgs)) {
            mean_rna_count = dotArgs[["mean_rna_count"]]
        }
        print(paste0("mean_rna_count cutoff: ", as.character(mean_rna_count)))
        if (mean_rna_count == 0) {
            count_data <- rna_count_data
        } else {
            rowmeans <- rowMeans(rna_count_data, na.rm = TRUE)
            inds <- which(rowmeans < mean_rna_count)
            count_data <- rna_count_data[-inds, ]
        }
    }
    
    # Step 2: Define the group variable (for example, treatment vs control)
    group <- factor(samples$Group)
    
    # Step 3: Create design matrix (used for limma)
    design <- model.matrix(~0 + group)
    colnames(design)[1:2] <- c(group1, group2)
    
    # Step 4: Add covariates to design matrix if provided (for limma)
    if (!is.null(covariates)) {
        for (covariate_col in covariates) {
            covariate <- samples[[covariate_col]]
            if (any(is.na(covariate))) {
                stop(paste("Covariate column", covariate_col, "contains NAs. Please handle missing values before running the analysis."))
            }
            if (!is.numeric(covariate)) {
                covariate <- as.numeric(factor(covariate))
            }
            design <- cbind(design, covariate)
            colnames(design)[ncol(design)] <- covariate_col
        }
    }
    
    # Initialize variable to store DE results
    results_df <- NULL
    
    # Step 5: Fit the model using the chosen method
    if (is.character(method)) {
        if (method == "limma") {
            library(limma)
            cat("Using limma...\n")
            
            if (normalization == "local") {
                count_data_norm <- rna_normalize(count_data)
            } else {
                pdf_file <- paste0(result_dir, "/DEG.RNA.limma.pdf")
                if (X11 == "No") {
                    pdf(pdf_file)
                }
                fit <- voom(count_data, design, plot = TRUE)
                count_data_norm <- fit
                dev.off()
            }
            
            fit_lm <- lmFit(count_data_norm, design)
            contrasts_str <- paste0(group1, '-', group2)
            fit_lm_cm <- eval(parse(text=paste0("makeContrasts(", contrasts_str, ", levels = design)")))
            fit_lm_c <- contrasts.fit(fit_lm, fit_lm_cm)
            ebayes_fit_lm_c <- eBayes(fit_lm_c)

            # limma diagnosis plot 
            pdf_file <- file.path(result_dir, "rpMeta.rna.limma.MA.MDS.pdf")
            if (X11 == "No") {
                pdf(pdf_file)
            }
            plotMA(ebayes_fit_lm_c, main="RNA MA plot")
            plotMDS(count_data, labels=samples$Sample_ID, col=as.numeric(samples$Group))
            dev.off()
            
            p_values <- ebayes_fit_lm_c$p.value[,1]
            log_fc <- ebayes_fit_lm_c$coefficients[,1]
            gene_names <- rownames(count_data)
            results_df <- data.frame(gene_name = gene_names, log_fc = log_fc, p_value = p_values, stringsAsFactors = FALSE)
            results_df$FDR <- p.adjust(results_df$p_value, method = "fdr")
            
        } else if (method == "deseq2") {
            cat("Using DESeq2...\n")
            library(DESeq2)
            
            # Build DESeq2 design formula
            if (is.null(covariates)) {
                design_formula <- ~ Group
            } else {
                covariate_formula <- paste(covariates, collapse = " + ")
                design_formula <- as.formula(paste("~", covariate_formula, "+ Group"))
            }
            rownames(samples) <- samples$Sample_ID
            samples$Group <- factor(samples$Group)

            # get DESeq2 parameters from ... arguments 
            DESeq2_test <- ifelse("DESeq2_test" %in% names(dotArgs), dotArgs[["DESeq2_test"]], "Wald")
            DESeq2_fitType <- ifelse("DESeq2_fitType" %in% names(dotArgs), dotArgs[["DESeq2_fitType"]], "parametric")

            dds <- DESeqDataSetFromMatrix(countData = count_data, colData = samples, design = design_formula)
            dds <- DESeq(dds, test = DESeq2_test, fitType = DESeq2_fitType)
            res <- results(dds, contrast = c("Group", group1, group2))
            
            p_values <- res$pvalue
            log_fc <- res$log2FoldChange
            gene_names <- rownames(count_data)
            results_df <- data.frame(gene_name = gene_names, log_fc = log_fc, p_value = p_values, stringsAsFactors = FALSE)
            results_df$FDR <- p.adjust(results_df$p_value, method = "fdr")
            
        } else {
            stop("Invalid method. Choose either 'limma', 'deseq2' or provide a user-defined function.")
        }
    } else if (is.function(method)) {
        # User-supplied method for DE analysis.
        cat("Using user-defined DE function...\n")
        results_df <- method(samples, count_data, group1, group2, covariates, normalization)
        # Check that the returned data frame contains the required columns.
        required_cols <- c("gene_name", "log_fc", "p_value", "FDR")
        if (!all(required_cols %in% colnames(results_df))) {
            stop("The user-defined DE function must return a data frame with columns: gene_name, log_fc, p_value, FDR.")
        }
    } else {
        stop("Invalid method. Choose either 'limma', 'deseq2' or provide a user-defined function.")
    }
    
    # Step 6: Write the results to a file
    if (is.function(method)) {
        res_file <- paste0(result_dir, "/rpMeta.DEG.RNA.", "user_defined",  ".txt")
    } else {
        res_file <- paste0(result_dir, "/rpMeta.DEG.RNA.", as.character(method), ".txt")
    }
    write.table(results_df, file = res_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    return(results_df)
}

#' Perform RNASeq count matrix normalization using log2(CPM) and median center 
#' 
#' `rna_normalize()` conducts RNASeq count matrix normalization
#' @param count_data RNASeq count data frame, rows are genes and columns are samples. 
#' 
#' @return a data frame of normalized RNA quantification matrix 
#' @export
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

