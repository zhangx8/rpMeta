# rpHeatmap.R
#' Plot RNA-protein heatmap plot
#' 
#' `rpHeatmap()` plots an integrated heatmap for joint RNA and protein expression visualization.
#' rpMeta.rpHeatmap.meta.[timestamp].pdf displays RNA and protein expression side by side
#' samples are clustered within each group by RNA expression 
#' rpMeta.heatmap.DEG.[timestamp].pdf displays RNA and protein expression separately 
#'
#' @param rpDEG combined RNA protein differential analysis result data frame 
#' @param rna_norm_counts an RNA count data frame with samples as columns and gene symbols as rows
#' @param prot_norm_counts a proteomics quantification data frame with samples as columns and gene symbols as rows
#' @param samples sample description data frame with group and covariates columns
#' @param covariates covariates used in the model, Default=NULL
#' @param result_dir absolute path where the results should be stored, default="./"
#' @param rpscale Z score normalization or not (default="Yes")
#' @param boundN minimum and maximum value of the heatmap color scale 
#' @param genesToDisplay Either a numeric value indicating the top n meta DEGs (by Meta_PVal) to display or a vector/list of gene symbols.
#' @param X11 "No" (default) to write PDF, or "Yes" to plot to X11 device
#' @param rpColCluster option to cluster colums(samples), "group" or "all", default="group"
#' @param ... Additional arguments. 
#' `legend` True or False, legend to be drawn or not  
#' `scale` character indicating if the values should be centered and
#'          scaled in either the row direction or the column direction,
#'          or none. Corresponding values are ‘"row"’, ‘"column"’ and ‘"none"’
#' 
#' @return None
#' @examples
#' # rpHeatmap(m, rna_norm_counts, prot_norm_counts, samples, covariates=covariates, result_dir, genesToDisplay=20)
#' @export
rpHeatmap <- function(rpDEG, rna_norm_counts, prot_norm_counts, samples, covariates = NULL, result_dir = "./", rpscale="Yes", boundN = 2.0, genesToDisplay = 20, X11 = "No", rpColCluster="group", ...) {
    library(pheatmap)
    library(gridExtra)
    df <- rpDEG
    dotArgs <- list(...)

    # print out user's dotArgs that uses pheatmap() arguments 
    pheat_args <- intersect(names(dotArgs), names(formals(pheatmap)))
    if (length(pheat_args)) {
        message(
          "You have overridden pheatmap() defaults for: ",
          paste(pheat_args, collapse = ", ")
        )
    }

    # validate X11
    if (!X11 %in% c("No","Yes")) {
        stop("X11 must be 'No' or 'Yes'")
    }

    # validate rpscale
    if (!rpscale %in% c("No","Yes")) {
        stop("rpscale must be 'No' or 'Yes'")
    }

    # validate rpscale
    if (!rpColCluster %in% c("group","all")) {
        stop("rpColCluster must be 'group' or 'all'")
    }

    # Helper function to normalize counts
    minmax <- function(x, limit = boundN) {
        return(pmax(pmin(x, limit), -limit))
    }

    # Filter DEGs based on log2 fold change and FDR for RNA and Protein
    # deg_rna <- df[df$RNA_FDR < 0.05 & abs(df$RNA_log2FC) > 0.263, ]
    # deg_prot <- df[df$Prot_FDR < 0.05 & abs(df$Prot_log2FC) > 0.263, ]
    deg_rna <- df[df$RNA_FDR < 0.05, ]
    deg_prot <- df[df$Prot_FDR < 0.05, ]

    # Select top 10 up and downregulated DEGs for RNA and Protein
    top10_rna_up <- head(deg_rna[order(-deg_rna$RNA_log2FC), ], 10)
    top10_rna_down <- head(deg_rna[order(deg_rna$RNA_log2FC), ], 10)
    top10_prot_up <- head(deg_prot[order(-deg_prot$Prot_log2FC), ], 10)
    top10_prot_down <- head(deg_prot[order(deg_prot$Prot_log2FC), ], 10)

    # Filter meta DEGs based on adjusted meta P-values and select genes to display
    deg_meta <- df[df$Meta_FDR < 0.05, ]
    if (is.numeric(genesToDisplay)) {
        displayGenes <- head(deg_meta[order(deg_meta$Meta_PVal), ], genesToDisplay)$Gene
    } else {
        displayGenes <- genesToDisplay
    }

    # Annotation Data (Group and Covariates)
    rownames(samples) <- samples$Sample_ID
    if (is.null(covariates)) {
        cols <- c("Group", "Sample_ID")
        print("covariates is null")
        # print(cols)
    } else {
        cols <- c("Group", covariates)
    }
    annotation <- samples[, cols, drop = FALSE]

    # set annotation colors
    grp_colors <- c("Black", "Red")
    unique_groups <- unique(annotation$Group)
    if (length(grp_colors) == length(unique_groups)) {
        names(grp_colors) <- unique_groups
    }
    rp_annotation_colors <- list(Group = grp_colors)
    rp_colors = c(RNA = "lightblue", Protein = "darkgreen")
    rp_annotation_colors <- list(Group = grp_colors, Type = rp_colors)
    if (!"annotation_colors" %in% names(dotArgs)) {
        dotArgs[["annotation_colors"]] <- rp_annotation_colors
    }

    # Define color palette for heatmap
    mycolors <- colorRampPalette(c("Navyblue", "White", "firebrick"))(400)

    # Open PDF to save RNA & Protein heatmaps
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    if (X11=="No") pdf(file = paste0(result_dir, "/rpMeta.heatmap.DEG.by.set.", timestamp, ".pdf"))

    rna_top10 <- rbind(rna_norm_counts[rownames(rna_norm_counts) %in% top10_rna_up$Gene, ],
                      rna_norm_counts[rownames(rna_norm_counts) %in% top10_rna_down$Gene, ])
    # if rpscale is Yes, Z score normalization by row
    if (rpscale == "Yes") rna_top10 <- t( scale(t(rna_top10), center=TRUE, scale=TRUE) )
    rna_top10_plot <- apply(rna_top10, c(1, 2), minmax, boundN)

    pheatmap(rna_top10_plot, 
             annotation_col = annotation,
             color = mycolors,
             show_rownames = TRUE, 
             fontsize_row = 9, 
             fontsize_col = 7, 
             main = "RNA Top 10 DEGs (Up/Down)", 
             colnames = NULL,
             ...)

    prot_top10 <- rbind(prot_norm_counts[rownames(prot_norm_counts) %in% top10_prot_up$Gene, ],
                        prot_norm_counts[rownames(prot_norm_counts) %in% top10_prot_down$Gene, ])
    # if rpscale is Yes, Z score normalization by row
    if (rpscale == "Yes") prot_top10 <- t( scale(t(prot_top10), center=TRUE, scale=TRUE) )
    prot_top10_plot <- apply(prot_top10, c(1, 2), minmax, boundN)
    pheatmap(prot_top10_plot, 
             annotation_col = annotation, 
             color = mycolors,
             show_rownames = TRUE, 
             fontsize_row = 9, 
             fontsize_col = 7, 
             main = "Protein Top 10 DEGs (Up/Down)", 
             colnames = NULL,
             ...)
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

    if (rpColCluster == "group") {
        # Reorder columns: Cluster within each group and then concatenate
        sample_groups <- annotation[colnames(combined_display_meta_ordered), "Group"]
        names(sample_groups) <- rownames(annotation[colnames(combined_display_meta_ordered), ])
        unique_groups <- unique(sample_groups)
        ordered_cols <- c()
        for (grp in unique_groups) {
            grp_cols <- names(sample_groups)[sample_groups == grp]
            if (length(grp_cols) > 1) {
                d <- dist(t(combined_display_meta_ordered[, grp_cols, drop = FALSE]))
                hc_cols <- hclust(d)
                ordered_cols <- c(ordered_cols, grp_cols[hc_cols$order])
            } else {
               ordered_cols <- c(ordered_cols, grp_cols)
            }
        }
        combined_display_meta_ordered <- combined_display_meta_ordered[, ordered_cols, drop = FALSE]
        annotation_ordered <- annotation[ordered_cols, , drop = FALSE]
        ann_col_reord <- c(sort(covariates), "Group")
        annotation_ordered <- annotation_ordered[, ann_col_reord, drop = FALSE]
        rp_cluster_cols <- FALSE
    } else {
        combined_display_meta_ordered <- combined_display_meta_ordered
        annotation_ordered <- annotation
        ann_col_reord <- c(sort(covariates), "Group")
        annotation_ordered <- annotation_ordered[, ann_col_reord, drop = FALSE]
        rp_cluster_cols <- TRUE
    }

    if (X11=="No") pdf(file = paste0(result_dir, "/rpMeta.rpHeatmap.meta.", timestamp, ".pdf"))
    pheatmap(combined_display_meta_ordered, 
             annotation_col = annotation_ordered, 
             annotation_row = row_annotation, 
             color = mycolors,
             show_rownames = TRUE, 
             fontsize_row = 9, 
             fontsize_col = 6, 
             cluster_rows = FALSE, 
             cluster_cols = rp_cluster_cols, 
             main = "rpHeatmap", 
             colnames = NULL,
             ...)
    if (X11=="No") dev.off()
}

