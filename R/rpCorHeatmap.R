# rpCorHeatmap.R
#' Plot RNA-protein heatmap for continuous variable correlation analysis 
#' 
#' `rpCorHeatmap()` plots an integrated heatmap for joint RNA and protein expression visualization (for continuous variable).
#' rpHeatmap.cor.meta.pdf displays RNA and protein expression side by side 
#' heatmap.cor.pdf displays RNA and protein expression separately 
#'
#' @param rpCor combined RNA protein correlation analysis result data frame 
#' @param rna_norm_counts an RNA count data frame with samples as columns and gene symbols as rows
#' @param prot_norm_counts a proteomics quantification data frame with samples as columns and gene symbols as rows
#' @param samples sample description data frame with group and covariates columns
#' @param continuous_var continuous variable used in the analysis
#' @param result_dir absolute path where the results should be stored, default="./"
#' @param rpscale Z score normalization or not (default="Yes")
#' @param boundN minimum and maximum value of the heatmap color scale 
#' @param genesToDisplay Either a numeric value indicating the top n meta DEGs (by Meta_PVal) to display or a vector/list of gene symbols.
#' @param X11 "No" (default) to write PDF, or "Yes" to plot to X11 device
#' @param ... Additional arguments. 
#' `legend` True or False, legend to be drawn or not  
#' `scale` character indicating if the values should be centered and
#'          scaled in either the row direction or the column direction,
#'          or none. Corresponding values are ‘"row"’, ‘"column"’ and ‘"none"’
#' 
#' @return None
#' @examples
#' # rpCorHeatmap(m, rna_norm_counts, prot_norm_counts, samples, continuous_var=continuous_var, result_dir, genesToDisplay=20)
#' @export
rpCorHeatmap <- function(rpCor, rna_norm_counts, prot_norm_counts, samples, continuous_var, result_dir = "./", rpscale="Yes", boundN = 2.0, genesToDisplay = 20, X11 = "No", ...) {
    library(pheatmap)
    library(gridExtra)
    dotArgs <- list(...)

    # validate X11
    if (!X11 %in% c("No","Yes")) {
        stop("X11 must be 'No' or 'Yes'")
    }

    # validate rpscale
    if (!rpscale %in% c("No","Yes")) {
        stop("rpscale must be 'No' or 'Yes'")
    }

    # Helper function to normalize counts
    minmax <- function(x, limit = boundN) {
        return(pmax(pmin(x, limit), -limit))
    }

    # Filter DEGs based on log2 fold change and FDR for RNA and Protein
    # cor_rna <- rpCor[rpCor$RNA_FDR < 0.05, ]
    # cor_prot <- rpCor[rpCor$Prot_FDR < 0.05, ]

    # Select top 10 up and downregulated DEGs for RNA and Protein
    top10_rna_up <- head(rpCor[order(-rpCor$RNA_Cor), ], 10)
    top10_rna_down <- head(rpCor[order(rpCor$RNA_Cor), ], 10)
    top10_prot_up <- head(rpCor[order(-rpCor$Prot_Cor), ], 10)
    top10_prot_down <- head(rpCor[order(rpCor$Prot_Cor), ], 10)

    # Filter meta DEGs based on adjusted meta P-values and select genes to display
    cor_meta <- rpCor
    if (is.numeric(genesToDisplay)) {
        displayGenes <- head(cor_meta[order(cor_meta$Meta_PVal), ], genesToDisplay)$Gene
    } else {
        displayGenes <- genesToDisplay
    }

    # Annotation Data (Group and Covariates)
    rownames(samples) <- samples$Sample_ID
    # cols <- c(continuous_var, "Sample_ID")
    cols <- c(continuous_var)
    annotation <- samples[, cols, drop = FALSE]

    # set annotation colors
    rp_colors = c(RNA = "lightblue", Protein = "darkgreen")
    annotation_colors_RP <- list(Type = rp_colors)

    # Define color palette for heatmap
    mycolors <- colorRampPalette(c("Navyblue", "White", "firebrick"))(400)

    # Open PDF to save RNA & Protein heatmaps
    if (X11=="No") pdf(file = paste0(result_dir, "/rpMeta.heatmap.cor.pdf"))

    print(top10_rna_up)
    print(top10_rna_down)
    rna_top10 <- rbind(rna_norm_counts[rownames(rna_norm_counts) %in% top10_rna_up$Gene, ],
                      rna_norm_counts[rownames(rna_norm_counts) %in% top10_rna_down$Gene, ])
    # if rpscale is Yes, Z score normalization by row
    if (rpscale == "Yes") rna_top10 <- t( scale(t(rna_top10), center=TRUE, scale=TRUE) )
    rna_top10_plot <- apply(rna_top10, c(1, 2), minmax, boundN)

    print(head(rna_top10_plot))
    print(head(annotation))

    pheatmap(rna_top10_plot, 
             annotation_col = annotation, 
             color = mycolors,
             # annotation_colors = annotation_colors,
             show_rownames = TRUE, 
             fontsize_row = 9, 
             fontsize_col = 7,
	     cluster_rows = FALSE,
	     cluster_cols = FALSE, 
             main = "RNA Top 10 Correlated Genes (Up/Down)", 
             colnames = NULL,
             ...)
    print("rna heatmap done")

    prot_top10 <- rbind(prot_norm_counts[rownames(prot_norm_counts) %in% top10_prot_up$Gene, ],
                        prot_norm_counts[rownames(prot_norm_counts) %in% top10_prot_down$Gene, ])
    # if rpscale is Yes, Z score normalization by row
    if (rpscale == "Yes") prot_top10 <- t( scale(t(prot_top10), center=TRUE, scale=TRUE) )
    prot_top10_plot <- apply(prot_top10, c(1, 2), minmax, boundN)
    pheatmap(prot_top10_plot, 
             annotation_col = annotation, 
             color = mycolors,
             # annotation_colors = annotation_colors, 
             show_rownames = TRUE, 
             fontsize_row = 9, 
             fontsize_col = 7,
             cluster_rows = FALSE,
             cluster_cols = FALSE, 
             main = "Protein Top 10 Correlated Genes (Up/Down)", 
             colnames = NULL,
             ...)
    print("protein heatmap done")
    if (X11=="No") dev.off()

    # Prepare Combined Heatmap Data for Meta DEGs (RNA and Protein)
    prepare_combined_heatmap <- function(rna_data, prot_data) {
        combined_matrix <- NULL
        genes = intersect(rownames(rna_data), rownames(prot_data))
        for (gene in genes) {
            gene_rna <- rna_data[rownames(rna_data) == gene, , drop = FALSE]
            gene_prot <- prot_data[rownames(prot_data) == gene, , drop = FALSE]
            if (nrow(gene_rna) == 0) {
                gene_rna <- matrix(NA, nrow = 1, ncol = ncol(rna_data))
                rownames(gene_rna) <- paste0(gene, "_R")
            } else {
                rownames(gene_rna) <- paste0(gene, "_R")
            }
            if (nrow(gene_prot) == 0) {
                gene_prot <- matrix(NA, nrow = 1, ncol = ncol(prot_data))
                rownames(gene_prot) <- paste0(gene, "_P")
            } else {
                rownames(gene_prot) <- paste0(gene, "_P")
            }
            combined_matrix <- rbind(combined_matrix, gene_rna, gene_prot)
        }
        return(combined_matrix)
    }

    # RNA quantification matrix 
    meta_norm_counts <- rna_norm_counts[rownames(rna_norm_counts) %in% displayGenes, ]
    if (rpscale == "Yes") meta_norm_counts <- t( scale(t(meta_norm_counts), center=TRUE, scale=TRUE) )
    meta_norm_counts_plot <- apply(meta_norm_counts, c(1, 2), minmax, boundN)
    # protein matrix 
    meta_norm_counts_prot <- prot_norm_counts[rownames(prot_norm_counts) %in% displayGenes, ]
    if (rpscale == "Yes") meta_norm_counts_prot <- t( scale(t(meta_norm_counts_prot), center=TRUE, scale=TRUE) )
    meta_norm_counts_prot_plot <- apply(meta_norm_counts_prot, c(1, 2), minmax, boundN)
    # Combined heatmap for selected meta DEGs (displayGenes)
    combined_display_meta <- prepare_combined_heatmap(meta_norm_counts_plot, meta_norm_counts_prot_plot)

    # Add row annotation for RNA/Protein
    row_type <- ifelse(grepl("_R$", rownames(combined_display_meta)), "RNA", "Protein")
    row_annotation <- data.frame(Type = row_type)
    rownames(row_annotation) <- rownames(combined_display_meta)

    # Reorder combined_display_meta based on RNA clustering
    rna_rows <- grep("_R$", rownames(combined_display_meta))
    rna_data_subset <- combined_display_meta[rna_rows, , drop = FALSE]
    hc <- hclust(dist(rna_data_subset))
    ordered_genes <- sub("_R$", "", rownames(rna_data_subset)[hc$order])
    all_genes <- sub("_[RP]$", "", rownames(combined_display_meta))
    gene_order_factor <- factor(all_genes, levels = ordered_genes)
    combined_display_meta_ordered <- combined_display_meta[order(gene_order_factor), ]

    # sort columns by continuous_var
    annotation <- annotation[colnames(combined_display_meta_ordered), , drop=FALSE]
    cont_order <- order(annotation[[continuous_var]], na.last=TRUE)
    combined_display_meta_ordered <- combined_display_meta_ordered[, cont_order, drop=FALSE]
    annotation               <- annotation[cont_order, , drop=FALSE]

    if (X11=="No") pdf(file = paste0(result_dir, "/rpMeta.rpHeatmap.cor.meta.pdf"))
    pheatmap(combined_display_meta_ordered, 
             annotation_col = annotation, 
             annotation_colors = annotation_colors_RP,
             annotation_row = row_annotation, 
             color = mycolors,
             show_rownames = TRUE, 
             fontsize_row = 9, 
             fontsize_col = 6, 
             cluster_rows = FALSE, 
             cluster_cols = FALSE, 
             main = "rpCorHeatmap", 
             colnames = NULL,
             ...)
    if (X11=="No") dev.off()
}

