#' Perform pathway enrichment analysis on correlation analysis result
#'
#' `rpCorPathway()` performs KEGG, Cell_Type_Signature, IMMUNESIGDB, GO:BP, 
#' or Hallmark pathway enrichment analysis on continuous variable correlation
#' results. Analysis is based on the GSEA method using the fgsea R package.
#' Molecular Signatures Database (MSigDB) gene sets are accessed via the 
#' msigdbr R package.
#'
#' @param df Combined RNA-protein correlation analysis result data frame.
#' @param result_dir Absolute path where the results should be stored.
#' @param FDR_cutoff False discovery rate (FDR) cutoff. Default is 0.25.
#' @param pathways Pathway database (Hallmark, KEGG, Cell_Type_Signature, 
#'   IMMUNESIGDB or GO_BP). Default is "Hallmark".
#' @param species Species name, e.g., "Homo sapiens".
#'
#' @return A data frame of pathway enrichment results.
#' @references
#' Subramanian A et al. (2005) "Gene set enrichment analysis: a 
#' knowledge-based approach for interpreting genome-wide expression profiles." 
#' Proc Natl Acad Sci U S A. 102(43):15545-15550.
#'
#' @export
rpCorPathway <- function(df, result_dir, pathways = "Hallmark", FDR_cutoff = 0.25, species = "Homo sapiens") {
    library(fgsea)
    library(msigdbr)
    library(data.table)
    library(dplyr)

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
        gene_sets <- msigdbr(species = species, category = "C5", subcategory = "GO:BP")
    } else {
        stop("Invalid pathway selection. Choose from 'Hallmark', 'KEGG', 'Cell_Type_Signature', 'IMMUNESIGDB' or 'GO_BP'.")
    }
    pathways_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)

    # Compute per‐category statistics with duplicate resolution and epsilon-padding
    
    # 1. RNA Statistic
    df_rna <- df[!is.na(df$RNA_Pval) & !is.na(df$Gene), ]
    if (nrow(df_rna) > 0) {
        # Sort by p-value (most significant first) and deduplicate
        df_rna <- df_rna[order(df_rna$RNA_Pval), ]
        df_rna_uniq <- df_rna[!duplicated(df_rna$Gene), ]
        
        p_vals <- df_rna_uniq$RNA_Pval
        p_vals[p_vals == 0] <- .Machine$double.xmin
        
        cor_val <- df_rna_uniq$RNA_Cor
        cor_val[is.na(cor_val)] <- 0
        
        rna_stat <- sign(cor_val) * -log10(p_vals)
        names(rna_stat) <- df_rna_uniq$Gene
        rna_stat <- sort(rna_stat[is.finite(rna_stat)], decreasing = TRUE)
    } else {
        rna_stat <- numeric()
    }

    # 2. Protein Statistic
    df_prot <- df[!is.na(df$Prot_Pval) & !is.na(df$Gene), ]
    if (nrow(df_prot) > 0) {
        # Sort by p-value (most significant first) and deduplicate
        df_prot <- df_prot[order(df_prot$Prot_Pval), ]
        df_prot_uniq <- df_prot[!duplicated(df_prot$Gene), ]
        
        p_vals <- df_prot_uniq$Prot_Pval
        p_vals[p_vals == 0] <- .Machine$double.xmin
        
        cor_val <- df_prot_uniq$Prot_Cor
        cor_val[is.na(cor_val)] <- 0
        
        prot_stat <- sign(cor_val) * -log10(p_vals)
        names(prot_stat) <- df_prot_uniq$Gene
        prot_stat <- sort(prot_stat[is.finite(prot_stat)], decreasing = TRUE)
    } else {
        prot_stat <- numeric()
    }

    # 3. Meta Statistic (using RNA correlation for directionality)
    df_meta <- df[!is.na(df$Meta_PVal) & !is.na(df$Gene), ]
    if (nrow(df_meta) > 0) {
        # Sort by p-value (most significant first) and deduplicate
        df_meta <- df_meta[order(df_meta$Meta_PVal), ]
        df_meta_uniq <- df_meta[!duplicated(df_meta$Gene), ]
        
        p_vals <- df_meta_uniq$Meta_PVal
        p_vals[p_vals == 0] <- .Machine$double.xmin
        
        cor_val <- df_meta_uniq$RNA_Cor
        cor_val[is.na(cor_val)] <- 0
        
        meta_stat <- sign(cor_val) * -log10(p_vals)
        names(meta_stat) <- df_meta_uniq$Gene
        meta_stat <- sort(meta_stat[is.finite(meta_stat)], decreasing = TRUE)
    } else {
        meta_stat <- numeric()
    }

    # Run fgsea for each category
    combined_results <- list()
    for (type in c("RNA","Protein","Meta")) {
        stat <- switch(type, RNA=rna_stat, Protein=prot_stat, Meta=meta_stat)
        if (length(stat) == 0) next
        message("Running pathway analysis for ", type, "...")
        fgseaRes <- fgsea(pathways=pathways_list, stats=stat, minSize=15, maxSize=500)
        fgseaRes <- fgseaRes[fgseaRes$padj < FDR_cutoff, ]
        if (nrow(fgseaRes) == 0) next
        fgseaRes <- fgseaRes[order(-fgseaRes$NES), ]
        fgseaRes$Category <- type
        combined_results[[type]] <- fgseaRes
    }

    if (length(combined_results) == 0) {
        message("No significant pathways detected in any category.")
        return(NULL)
    }
    combined_df <- bind_rows(combined_results)

    # Save full results
    res_file <- file.path(result_dir, paste0("rpMeta.rpCorPathway.", tolower(pathways), ".txt"))
    fwrite(combined_df, file = res_file, sep = "\t")
    message("rpCorPathway results saved to: ", res_file)
    invisible(combined_df)
}
