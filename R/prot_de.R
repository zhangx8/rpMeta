#' Calculate differrential expression from proteomics data file
#' 
#' `prot_de()` conducts differential expression analysis using limma. 
#' Covariates can be added to the model. 
#'
#' @param samples a data frame of sample information, including "Sample_ID", "Group" and covariates
#' @param prot_quant_data  a data frame of proteomics quantification matrix with rows as genes and columns as samples
#' @param result_dir absolute path where the results should be stored
#' @param group1 name of the first group used in DE analysis, a.k.a condition1 or factor1
#' @param group2 name of the second group used in DE analysis, a.k.a condition2 or factor2
#' fold change is calculated as group2/group1
#' @param covariates covariates used in the model (this must be included in the sample_file)
#' @param method Differential expression analysis method, Default="limma"
#' @param X11 output X11 or PDF for graphics, c("No", "Yes"), (Default="No")
#' @param ... Additional arguments 
#' `limma_proportion` numeric value between 0 and 1, assumed proportion of genes which are differentially expressed 
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
#' prot_method <- "limma"
#' result_dir <- file.path(getwd(), "rpMeta_APOLLO1")
#' limma_proportion <- 0.02
#' if (!dir.exists(result_dir)) { dir.create(result_dir) }
#' prot_de(samples, prot_quant_data, result_dir, group1, group2, covariates, method=prot_method, X11="Yes", limma_proportion=limma_proportion)
#' @export
prot_de <- function(samples, prot_quant_data, result_dir, group1="TRU", group2="ProxProlif",
                    covariates=NULL, method="limma", X11="No", ...) {
    cat("Running protein DE analysis...\n")
    dotArgs <- list(...)
    if (X11 != "No" && X11 != "Yes") {
        print("X11 argument musct by No or Yes")
        return("X11_Error")
    }

    # samples is the sample table
    # prot_quant_data is the protein quantification table 
    
    # Step 1: Create group and design matrix (for limma)
    group <- factor(samples$Group)
    design <- model.matrix(~0 + group)
    colnames(design)[1:2] <- c(group1, group2)
    
    # Add covariates to design matrix, if provided
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
    
    cat("Design matrix:\n")
    print(head(design))
    
    if (nrow(design) != ncol(prot_quant_data)) {
        stop("Row dimension of design matrix doesn't match the column dimension of the quantification matrix.")
    }
    
    # Step 2: Differential Expression Analysis based on chosen method
    results_df <- NULL
    if (is.character(method)) {
        if (method == "limma") {
            library(limma)
            cat("Using limma...\n")

            # get limma parameters from ... arguments 
            limma_proportion <- ifelse("limma_proportion" %in% names(dotArgs), dotArgs[["limma_proportion"]], 0.01)
            
            # Fit the linear model and perform contrast analysis
            fit <- lmFit(prot_quant_data, design)
            contrasts_str <- paste0(group1, "-", group2)
            contrast.matrix <- eval(parse(text = paste0("makeContrasts(", contrasts_str, ", levels = design)")))
            fit2 <- contrasts.fit(fit, contrast.matrix)
            fit2 <- eBayes(fit2, proportion=limma_proportion)
            
            # limma diagnosis plot 
            pdf_file <- file.path(result_dir, "rpMeta.prot.limma.MA.MDS.pdf")
            if (X11 == "No") {
                pdf(pdf_file)
            }
            plotMA(fit2, main="Protein MA plot")
            plotMDS(prot_quant_data, labels=samples$Sample_ID, col=as.numeric(samples$Group))
            dev.off()
            
            p_values <- fit2$p.value[,1]
            log_fc <- fit2$coefficients[,1]
            gene_names <- rownames(prot_quant_data)
            results_df <- data.frame(gene_name = gene_names, log_fc = log_fc, p_value = p_values,
                                     stringsAsFactors = FALSE)
            results_df$FDR <- p.adjust(results_df$p_value, method = "fdr")
            
        } else {
            stop("Invalid method. Only 'limma' or a user-defined DE function are allowed.")
        }
    } else if (is.function(method)) {
        cat("Using user-defined DE function...\n")
        # Call the user-supplied function with the necessary parameters.
        # For consistency, the function is expected to accept:
        # (samples, prot_quant_data, group1, group2, covariates, result_dir)
        results_df <- method(samples, prot_quant_data, group1, group2, covariates)
        
        # Check that the returned data frame contains the required columns
        required_cols <- c("gene_name", "log_fc", "p_value", "FDR")
        if (!all(required_cols %in% colnames(results_df))) {
            stop("The user-defined DE function must return a data frame with columns: gene_name, log_fc, p_value, FDR.")
        }
    } else {
        stop("Invalid method parameter. Must be 'limma' or a user-defined function.")
    }
    
    # Step 5: Write results to a file and return the results data frame
    file_suffix <- ifelse(is.character(method), method, "user_defined")
    res_file <- paste0(result_dir, "/rpMeta.DEG.protein.", file_suffix, ".txt")
    write.table(results_df, file = res_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    return(results_df)
}

