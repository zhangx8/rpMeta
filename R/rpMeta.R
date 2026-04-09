#' Perform combined RNA protein differential expression analysis
#' 
#' `rpMeta()` Perform combined RNA protein differential expression analysis.
#'
#' @param sample_file a tab-delimited sample description file , in which rows are samples and features are columns (sample group and covariates)
#' @param rna_count_file  an RNA count .csv file with samples as columns and gene symbols as rows
#' @param protein_quant_file a proteomics quantification .csv file with samples as columns and gene symbols as rows
#' @param group1 name of the group1 from sample file used for DEG analysis, a.k.a condition1 or factor1
#' @param group2 name of the group2 from sample file used for DEG analysis, a.k.a condition2 or factor2
#' fold change is calculated as group2/group1
#' @param covariates character array of covariate names used in the model (this must be included in the sample_file, Default = NULL)
#' @param continuous_var name of the continuous variable to be analyzed 
#' if continuous_var is not NULL, rpMeta() will perform correlation analysis. 
#' If both group1 and group2 are not NULL, rpMeta() will perform differential expression analysis. 
#' if both group1 and group2 are not NULL and continuous_var is not NULL, rpMeta() will perform differential expression and correlation analysis.  
#' @param result_dir absolute path where the results should be stored, Default=cwd
#' @param rna_method RNA differential expression method, limma, DESeq2 or custom function, Default="limma"
#' @param rna_cor_method RNA correlation method (continuous variable), "spearman" or "pearson" , Default="spearman" 
#' @param prot_method protein differential expression method, limma or custom function, Default="limma"
#' @param prot_cor_method protein correlation method (continuous variable), "spearman", "pearson" or "kendall", Default="spearman"
#' @param combine_method P value combine method, e.g. Fisher, weighted Fisher, Pearson, Stouffer, Tippet(min), Wilkinson(max), Default = "Fisher" 
#' @param rna_norm normalization method for limma, 
#' e.g. default(voom) or local (Median centered log2(CPM+1)). Default="default" (voom)
#' if rna_norm == "as_is", it means RNA data is already normalized 
#' @param X11 output X11 or PDF for graphics, c("No", "Yes"), (Default="No")
#' @param max_na_frac maximum fraction of NAs in each row(gene), default=0.5
#' @param pathways Pathway database (Hallmark, KEGG, Cell_Type_Signature, IMMUNESIGDB or GO_BP), default=c("Hallmark")
#' @param path_plot_fdr Pathway plot FDR cutoff (default=0.05)
#' @param ... Additional arguments 
#' `heatmap_genesToDisplay` (e.g. c("BRCA1", "ATM", "CASP8")),
#' `heatmap_scale` character indicating if the values should be Z score normalized by row
#'                 values are "Yes" or "No"
#'                 This is passed in to rpHeatmap() function as rpscale argument 
#' `GenePlot_genes` (e.g. c("MET", "CDK12", "CASP1")) 
#' `DESeq2_test` (e.g. c("Wald", "LRT"))
#' `DESeq2_fitType` (e.g. c("parametric", "local", "mean", "glmGamPoi"))
#' `limma_proportion` numeric value between 0 and 1, assumed proportion of genes which are differentially expressed 
#' 
#' @return workflow execution status 
#' @examples
#' sample_file <- system.file("extdata", "apollo1.sample.info.tsv", package = "rpMeta")
#' rna_count_file <- system.file("extdata", "apollo1.rna.count.top6k.tsv", package = "rpMeta")
#' protein_quant_file <- system.file("extdata", "apollo1.protein.quant.top3k.tsv", package = "rpMeta")
#' input_format <- "tsv"
#' combine_method <- "Fisher"
#' group1 <- "TRU"
#' group2 <- "ProxProlif"
#' covariates <- c("Sex", "grade")
#' rna_method <- "limma"
#' prot_method <- "limma"
#' result_dir <- file.path(getwd(), "rpMeta_example1")
#' if (!dir.exists(result_dir)) { dir.create(result_dir) }
#' heatmap_genesToDisplay <- c("CASP1", "CDK1", "ELN", "DARS", "CYP7B1", "CKS1B")
#' heatmap_scale <- "Yes"
#' mean_rna_count <- 5
#' DESeq2_test <- "Wald"
#' limma_proportion <- 0.02
#' continuous_var <- NULL
#' X11 <- "Yes"
#' rpMeta(sample_file, rna_count_file, protein_quant_file, input_format=input_format, group1=group1, group2=group2, covariates=covariates, continuous_var=continuous_var, result_dir, rna_method=rna_method, rna_cor_method="spearman", combine_method="Fisher", rna_norm="local", X11=X11, heatmap_genesToDisplay=heatmap_genesToDisplay, heatmap_scale=heatmap_scale, mean_rna_count=mean_rna_count, DESeq2_test=DESeq2_test, limma_proportion=limma_proportion)
#' @export
rpMeta <- function(sample_file, rna_count_file, protein_quant_file, input_format=input_format, group1=NULL, group2=NULL, covariates=NULL, continuous_var=NULL, result_dir=getwd(), rna_method="limma", rna_cor_method="spearman",prot_method="limma", prot_cor_method="spearman", combine_method="Fisher", rna_norm="default", heatmap_scale="Yes", X11="No", max_na_frac=0.5, pathways=c("Hallmark"), path_plot_fdr=0.05, ...) {
    library(limma)
    # Step 1: parse ... arguments 
    print("Step 1: parsing dot arguments ...")
    dotArgs <- list(...)
    heatmap_genesToDisplay = ifelse("heatmap_genesToDisplay" %in% (dotArgs), dotArgs[["heatmap_genesToDisplay"]], 10)
    mean_rna_count       = ifelse("mean_rna_count" %in% (dotArgs), dotArgs[["mean_rna_count"]], 0)
    DESeq2_test <- ifelse("DESeq2_test" %in% names(dotArgs), dotArgs[["DESeq2_test"]], "Wald")
    DESeq2_fitType <- ifelse("DESeq2_fitType" %in% names(dotArgs), dotArgs[["DESeq2_fitType"]], "parametric")
    limma_proportion <- ifelse("limma_proportion" %in% names(dotArgs), dotArgs[["limma_proportion"]], 0.01)
    if (X11 != "No" && X11 != "Yes") {
        print("X11 argument must by No or Yes")
        return("X11_Error")
    }
    if (!dir.exists(result_dir)) { dir.create(result_dir) }

    # Step 2: load and validate inputs
    print("Step 2: load and validate input files (sample, rna, protein) ...")
    inputs <- loadRpInput(sample_file, rna_count_file, protein_quant_file, input_format, max_na_frac=max_na_frac)
    sampleDF <- inputs$samples
    rna_count_data <- inputs$rna
    prot_quant_data <- inputs$protein
    common_samples <- colnames(rna_count_data)
    print("common_samples")
    print(common_samples)
    # get normalized RNA quantitifcation table (for heatmap and gene plot)
    if (rna_norm == "as_is") {
        rna_norm_counts <- rna_count_data
    } else {
        rna_norm_counts <- rna_normalize(rna_count_data)   
    }
    # protein quantifiction matrix should be already normalized 
    prot_norm_counts <- prot_quant_data
    # set gene list for heatmap 
    if ("heatmap_genesToDisplay" %in% names(dotArgs)) {
        heatmap_genesToDisplay <- dotArgs[["heatmap_genesToDisplay"]]
    } else {
        heatmap_genesToDisplay <- 20
    }
    # set gene list for rpGenePlot
    if ("GenePlot_genes" %in% names(dotArgs)) {
        GenePlot_genes <- dotArgs[["GenePlot_genes"]]
    } else {
        GenePlot_genes <- 20
    }

    # Step 3a: run RNA and protein correlation test, merge results, calculate meta P value
    if (!is.null(continuous_var) ) {
        print(paste0("Step 3a: RNA and protein correlation test on ", continuous_var))
        rna_cor_res <- rna_cor(sampleDF, rna_count_data, result_dir, continuous_var=continuous_var, method=rna_cor_method, normalization=rna_norm)
	prot_cor_res <- prot_cor(sampleDF, prot_quant_data, result_dir, continuous_var=continuous_var, method=prot_cor_method)
        # merge common genes 
        m_cor <- merge(rna_cor_res, prot_cor_res, by="Gene", all.x=F, all.y=F)
        # calculate meta P values and FDR
        m_cor$Meta_PVal <- apply(m_cor, 1, combinep, combine_method)
        m_cor$Meta_FDR <- p.adjust(m_cor$Meta_PVal, method="fdr")
        cor_res_file <- paste0(result_dir, "/rpMeta.correlation.combined.txt")
        write.table(m_cor, file = cor_res_file, sep = "\t", row.names = FALSE, quote = FALSE)
        # Merge RNA and protein (all genes)
        m_cor_all <- merge(rna_cor_res, prot_cor_res, by="Gene", all.x=T, all.y=T)
        m_cor_all$Meta_PVal <- apply(m_cor_all, 1, combinep, combine_method)
        m_cor_all$Meta_FDR <- p.adjust(m_cor_all$Meta_PVal, method="fdr")
        cor_res_file_all <- paste0(result_dir, "/rpMeta.correlation.combined.all.txt")
        write.table(m_cor_all, file = cor_res_file_all, sep = "\t", row.names = FALSE, quote = FALSE)
        # plot heatmap 
        rpCorHeatmap(m_cor, rna_norm_counts, prot_norm_counts, sampleDF, continuous_var, result_dir, X11=X11, genesToDisplay=heatmap_genesToDisplay, rpscale=heatmap_scale)
        # rpCorGenePlot
        rpCorGenePlot(df = m_cor, rna_norm_counts = rna_norm_counts, prot_norm_counts = prot_norm_counts, samples = sampleDF, continuous_var, result_dir, GenePlot_genes, X11=X11)
	# rpCorPathway
        for (i in 1:length(pathways)) {
	    pathway <- pathways[i]
    	    pathway_res <- rpCorPathway(m_cor, result_dir, pathways = pathway, FDR_cutoff = 0.25, species = "Homo sapiens")
    	    rpCorPathwayPlot(pathway_res, result_dir, continuous_var, pathways = pathway, FDR_cutoff = path_plot_fdr, max_pathways = 30, X11=X11)
        }
        # print_gene_stat(rna_cor_res, prot_cor_res, m_cor)
        print_log(rna_cor_res, prot_cor_res, m_cor, sample_file, rna_count_file, protein_quant_file,
              group1, group2, covariates, continuous_var, result_dir,
              rna_method, rna_cor_method, prot_method, prot_cor_method,
              combine_method, rna_norm, max_na_frac, pathways, common_samples, 
              heatmap_genesToDisplay = dotArgs[["heatmap_genesToDisplay"]],
              heatmap_scale = heatmap_scale,
              mean_rna_count       = dotArgs[["mean_rna_count"]],
              DESeq2_test = dotArgs[["DESeq2_test"]],
              DESeq2_fitType = dotArgs[["DESeq2_fitType"]],
              limma_proportion = dotArgs[["limma_proportion"]])
        if (is.null(group1) || is.null(group2)) {
            print("Continuous variable only analysis is done.")
            return("Done")
        }
    }

    if (is.null(group1) || is.null(group2)) {
        print("group1 or group2 is NULL. DEG analysis can not be performed")
        return("Done")
    } 

    # Step 3b: perform DEG analysis on RNA and protein data, merge results, calculate meta P value
    print("Step 3b: perform DEG analysis on RNA and protein data, merge results, calculate meta P value.")
    # run RNA limma or deseq2
    rna_res <- rna_de(sampleDF, rna_count_data, result_dir, group1, group2, method=rna_method, covariates, rna_norm, X11=X11, mean_rna_count=mean_rna_count, DESeq2_test=DESeq2_test)
    # run protein limma 
    prot_res <- prot_de(sampleDF, prot_quant_data, result_dir, group1, group2, covariates, prot_method, X11=X11, limma_proportion=limma_proportion)

    # merge rna and protein (only common genes) 
    m <- merge(rna_res, prot_res, by="gene_name", all.x=F, all.y=F)
    header <- c("Gene", "RNA_log2FC", "RNA_Pval", "RNA_FDR", "Prot_log2FC", "Prot_Pval", "Prot_FDR")
    colnames(m) <- header
    print_log(rna_res, prot_res, m, sample_file, rna_count_file, protein_quant_file,
              group1, group2, covariates, continuous_var, result_dir,
              rna_method, rna_cor_method, prot_method, prot_cor_method,
              combine_method, rna_norm, max_na_frac, pathways, common_samples,
              heatmap_genesToDisplay = dotArgs[["heatmap_genesToDisplay"]],
              heatmap_scale = heatmap_scale,
              mean_rna_count       = dotArgs[["mean_rna_count"]],
              DESeq2_test = dotArgs[["DESeq2_test"]],
              DESeq2_fitType = dotArgs[["DESeq2_fitType"]],
              limma_proportion = dotArgs[["limma_proportion"]])
    # calculate meta P values and FDR
    m$Meta_PVal <- apply(m, 1, combinep, combine_method)
    m$Meta_FDR <- p.adjust(m$Meta_PVal, method="fdr")
    res_file <- paste0(result_dir, "/rpMeta.DEG.meta.txt")
    write.table(m, file = res_file, sep = "\t", row.names = FALSE, quote = FALSE)

    # Merge RNA and protein (all genes)
    m_all <- merge(rna_res, prot_res, by="gene_name", all.x=T, all.y=T)
    colnames(m_all) <- header
    m_all$Meta_PVal <- apply(m_all, 1, combinep, combine_method)
    m_all$Meta_FDR <- p.adjust(m_all$Meta_PVal, method="fdr")
    res_file_all <- paste0(result_dir, "/rpMeta.DEG.meta.all.txt")
    write.table(m_all, file = res_file_all, sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Step 4: rpvolcano plot (including scater plot of RNA vs Protein log2FC
    print("Step 4: rpVolcano plot")
    rpVolcano(m, result_dir, group1, group2, genesToDisplay=6, X11=X11)

    # Step 5: heatmap
    # Set up annotation colors for Group if provided via dotArgs
    print("Step 5: rpHeatmap plot")
    rpHeatmap(m, rna_norm_counts, prot_norm_counts, sampleDF, covariates=covariates, result_dir, X11=X11, genesToDisplay=heatmap_genesToDisplay, rpscale=heatmap_scale)

    # Step 6: plot scatter plot per gene for top 20 DEGs
    print("Step 6: rpGenePlot")
    rpGenePlot(df = m, rna_norm_counts = rna_norm_counts,prot_norm_counts = prot_norm_counts, samples = sampleDF, group1, group2, result_dir, GenePlot_genes, X11=X11)

    # Step 7: Venn diagrams for m and m_all
    print("Step 7: rpVenn plot")
    venn_data <- rpVenn(m, fdr_cutoff=0.05, result_dir, filename_suffix="", X11=X11)  # For common genes
    print(venn_data)

    # Step 8: pathway enrichment analysis
    print("Step 8: rpPathway analysis and plot")
    for (i in 1:length(pathways)) {
        pathway <- pathways[i]
        pathway_res <- rpPathway(m, result_dir, pathways = pathway, FDR_cutoff = 0.25, species = "Homo sapiens")
        rpPathwayPlot(pathway_res, result_dir, group1, group2, pathways = pathway, FDR_cutoff = path_plot_fdr, max_pathways = 30, X11=X11)
    }

    sessionInfo()
    return("Done.")
}

print_log <- function(
  rna,
  prot,
  merged, 
  sample_file,
  rna_count_file,
  protein_quant_file,
  group1,
  group2,
  covariates,
  continuous_var,
  result_dir,
  rna_method,
  rna_cor_method,
  prot_method,
  prot_cor_method,
  combine_method,
  rna_norm,
  max_na_frac,
  pathways,
  common_samples,
  heatmap_genesToDisplay = NULL,
  heatmap_scale = NULL,
  mean_rna_count       = NULL,
  DESeq2_test          = NULL,
  DESeq2_fitType       = NULL,
  limma_proportion     = NULL,
  writeTable           = TRUE
) {
  title <- basename(result_dir)
  num_rna_genes <- nrow(rna)
  num_prot_genes <- nrow(prot)
  num_common_genes <- nrow(merged)
  cat("-----------------------------\n")
  cat("Gene number stats: \n")
  cat(paste0("Number of RNA genes: ", as.character(num_rna_genes), "\n"))
  cat(paste0("Number of protein genes: ", as.character(num_prot_genes), "\n"))
  cat(paste0("Number of common genes: ", as.character(num_common_genes), "\n"))
  cat("-----------------------------\n")
  cat("\n")
  
  log <- matrix(
    ncol    = 2,
    byrow   = TRUE,
    c(
      "sample_file",           sample_file,
      "rna_count_file",        rna_count_file,
      "protein_quant_file",    protein_quant_file,
      "group1",                group1,
      "group2",                group2,
      "covariates",            paste(covariates, collapse=";"),
      "continuous_var",        ifelse(is.null(continuous_var), "", continuous_var),
      "result_dir",            result_dir,
      "rna_method",            deparse(substitute(rna_method)),
      "rna_cor_method",        rna_cor_method,
      "prot_method",           deparse(substitute(prot_method)),
      "prot_cor_method",       prot_cor_method,
      "combine_method",        combine_method,
      "rna_norm",              rna_norm,
      "max_na_frac",	       as.character(max_na_frac),
      "pathways", 	       paste(pathways, collapse=";"),
      "heatmap_genesToDisplay",if (is.null(heatmap_genesToDisplay)) "" else paste(heatmap_genesToDisplay, collapse=";"),
      "heatmap_scale",         heatmap_scale,
      "mean_rna_count",        if (is.null(mean_rna_count)) "" else as.character(mean_rna_count),
      "DESeq2_test",           if (is.null(DESeq2_test)) "" else DESeq2_test,
      "DESeq2_fitType",        if (is.null(DESeq2_fitType)) "" else DESeq2_fitType,
      "limma_proportion",      if (is.null(limma_proportion)) "" else as.character(limma_proportion),
      "num_rna_genes",         as.character(num_rna_genes),
      "num_prot_genes",        as.character(num_prot_genes),
      "num_common_genes",      as.character(num_common_genes),
      "common_samples",        paste(common_samples, collapse=",")
    )
  )

  colnames(log) <- c("Argument", "Value")
  
  # Print to console
  print(as.data.frame(log), row.names = FALSE)
  
  # Write CSV if desired
  if (writeTable) {
    write.csv(
      as.data.frame(log),
      file      = file.path(result_dir, paste0(title, ".log.csv")),
      row.names = FALSE,
      quote     = FALSE
    )
  }
}

