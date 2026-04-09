#' Perform pathway enrichment analysis
#'
#' `rpPathway()` performs KEGG, Cell_Type_Signature, IMMUNESIGDB, GO:BP, 
#' or Hallmark pathway enrichment analysis
#' Analysis is based on GSEA method using fgsea R package
#' Molecular Signatures Database (MSigDB) gene sets are used through msigdbr R package 
#' Subramanian A et al. Gene set enrichment analysis: a knowledge-based approach for 
#' interpreting genome-wide expression profiles. Proc Natl Acad Sci U S A. 2005 Oct 25;102(43)
#'
#' @param df Combined RNA-protein differential analysis result data frame.
#' @param result_dir Absolute path where the results should be stored.
#' @param FDR_cutoff false discovery rate (FDR) cutoff, default=0.25
#' @param pathways Pathway database (Hallmark, KEGG, Cell_Type_Signature, IMMUNESIGDB or GO_BP), default="Hallmark"
#' @param species Species, e.g., "Homo sapiens".
#'
#' @return a data frame of pathway enrichment results
#' @export
rpPathway <- function(df, result_dir, pathways = "Hallmark", FDR_cutoff = 0.25, species = "Homo sapiens") {
    library(fgsea)
    library(msigdbr)
    library(data.table)
    # library(ggplot2)
    library(dplyr)
    # library(ggpubr)

    # Ensure result directory exists
    if (!dir.exists(result_dir)) {
        dir.create(result_dir, recursive = TRUE)
    }

    # Fetch Gene Sets
    message("Fetching gene sets...")
    if (pathways == "Hallmark") {
        gene_sets <- msigdbr(species = species, category = "H")
    } else if (pathways == "Cell_Type_Signature") {
	gene_sets <- msigdbr(species = species, category = "C8")
    } else if (pathways == "KEGG") {
        gene_sets <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG_LEGACY")
    } else if (pathways == "IMMUNESIGDB") {
        gene_sets <- msigdbr(species = species, category = "C7", subcategory = "IMMUNESIGDB")
    } else if (pathways == "GO_BP") {
        gene_sets <- msigdbr(species = species, category = "C5", subcategory = "BP")
    } else {
        stop("Invalid pathway selection. Choose from 'Hallmark', 'KEGG', 'Cell_Type_Signature', 'IMMUNESIGDB' or 'GO_BP'.")
    }
    pathways_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)

    # Compute per‐category statistics
    rna_stat <- with(df, sign(RNA_log2FC) * -log10(RNA_Pval))
    names(rna_stat) <- df$Gene
    rna_stat <- sort(rna_stat[!is.na(rna_stat) & is.finite(rna_stat)], decreasing = TRUE)

    prot_stat <- with(df, sign(Prot_log2FC) * -log10(Prot_Pval))
    names(prot_stat) <- df$Gene
    prot_stat <- sort(prot_stat[!is.na(prot_stat) & is.finite(prot_stat)], decreasing = TRUE)

    meta_stat <- with(df, sign(RNA_log2FC) * -log10(Meta_PVal))
    names(meta_stat) <- df$Gene
    meta_stat <- sort(meta_stat[!is.na(meta_stat) & is.finite(meta_stat)], decreasing = TRUE)

    # Run fgsea for each category
    combined_results <- list()
    for (type in c("RNA","Protein","Meta")) {
        stat <- switch(type, RNA=rna_stat, Protein=prot_stat, Meta=meta_stat)
        if (length(stat)==0) next
        message("Running pathway analysis for ", type, "...")
        fgseaRes <- fgsea(pathways=pathways_list, stats=stat, minSize=15, maxSize=500)
        fgseaRes <- fgseaRes[fgseaRes$padj < FDR_cutoff, ]
        if (nrow(fgseaRes)==0) next
        fgseaRes <- fgseaRes[order(-fgseaRes$NES), ]
        fgseaRes$Category <- type
        combined_results[[type]] <- fgseaRes
    }

    if (length(combined_results)==0) {
        message("No significant pathways detected in any category.")
        return(NULL)
    }
    combined_df <- bind_rows(combined_results)

    # Save full results
    res_file <- file.path(result_dir, paste0("rpMeta.rpPathway.", tolower(pathways), ".txt"))
    fwrite(combined_df, file = res_file, sep = "\t")
    message("rpPathway results saved to: ", res_file)
    invisible(combined_df)
}

