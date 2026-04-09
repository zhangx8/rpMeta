# rpCorGenePlot.R
#' `rpCorGenePlot()` plot RNA–protein scatterplot and boxplot for top 20 DEGs
#'
#' @param df                combined RNA–protein differential analysis result data frame 
#' @param rna_norm_counts   an RNA count data frame with samples as columns and gene symbols as rows
#' @param prot_norm_counts  an proteomics count data frame with samples as columns and gene symbols as rows
#' @param samples           sample description data frame with group and covariates columns
#' @param continuous_var    name of the continuous variable
#' @param result_dir        absolute path where the results should be stored, default="./"
#' @param genes             a vector of gene symbols, e.g. c("MET","BRCA2") or top n by Meta_FDR (default=20)
#' @param X11               "No" (default) to write PDF, or "Yes" to plot to X11 device
#'
#' @return A data.frame with columns \code{gene_name}, \code{RNA_pvalue}, \code{Prot_pvalue}, \code{Meta_pvalue}, \code{Spearman_rho}, \code{Spearman_pvalue}
#' @examples
#' # rpCorGenePlot(df=m_cor, rna_norm_counts=rna_norm_counts, prot_norm_counts=prot_norm_counts,
#' #            samples=samples, continuous_var, result_dir)
#' @export
rpCorGenePlot <- function(df, rna_norm_counts, prot_norm_counts, samples, continuous_var, result_dir="./", genes = 20, X11 = "No") {
    library(ggplot2)
    library(ggpubr)
    library(cowplot)

    # validate X11
    if (!X11 %in% c("No","Yes")) {
        stop("X11 must be 'No' or 'Yes'")
    }

    # select genes
    if (is.numeric(genes)) {
        top_meta <- head(df[order(df$Meta_FDR),], genes)
        plot_genes <- top_meta$Gene
    } else if (is.vector(genes)) {
        plot_genes <- genes
        top_meta <- df[df$Gene %in% plot_genes,]
    } else {
        top_meta <- head(df[order(df$Meta_FDR),], 20)
        plot_genes <- top_meta$Gene
    }

    # ensure sample order
    if (!identical(samples$Sample_ID, colnames(rna_norm_counts))) {
        stop("Sample order mismatch between samples table and RNA matrix")
    }
    if (!identical(samples$Sample_ID, colnames(prot_norm_counts))) {
        stop("Sample order mismatch between samples table and protein matrix")
    }

    # check gene availability
    avail   <- intersect(plot_genes, intersect(rownames(rna_norm_counts), rownames(prot_norm_counts)))
    missing <- setdiff(plot_genes, avail)
    if (length(missing) > 0) {
        message("These genes not found in both matrices and will be skipped: ", paste(missing, collapse=", "))
    }
    if (length(avail) == 0) {
        stop("No valid genes remain for plotting")
    }
    plot_genes <- avail
    top_meta   <- top_meta[top_meta$Gene %in% plot_genes,]

    # prepare storage for p-values and correlations
    rna_p_vals       <- setNames(numeric(length(plot_genes)), plot_genes)
    prot_p_vals      <- setNames(numeric(length(plot_genes)), plot_genes)
    meta_p_vals      <- setNames(numeric(length(plot_genes)), plot_genes)
    cor_vals         <- setNames(numeric(length(plot_genes)), plot_genes)
    spearman_p_vals  <- setNames(numeric(length(plot_genes)), plot_genes)
    gene_plots       <- list()

    for (g in plot_genes) {
        rvals <- as.numeric(rna_norm_counts[g,])
        pvals <- as.numeric(prot_norm_counts[g,])
        cont_vals <- samples[[continuous_var]]

        # pull p-values from the merged df
        rna_p_vals[g]  <- df$RNA_Pval[df$Gene == g]
        prot_p_vals[g] <- df$Prot_Pval[df$Gene == g]
        meta_p_vals[g] <- df$Meta_PVal[df$Gene == g]

        # compute genewise Spearman correlation and p-value (on multiple samples, internally drop missing values)
        test <- cor.test(rvals, pvals, method="spearman", exact=FALSE)
        cor_vals[g]        <- test$estimate
        spearman_p_vals[g] <- test$p.value

        pdf_df <- data.frame(Sample=colnames(rna_norm_counts), RNA=rvals, Protein=pvals, Continuous = cont_vals)

        p <- ggplot(pdf_df, aes(x = RNA, y = Protein, color = Continuous)) +
          geom_point(size = 3, alpha = 0.8) +
          scale_color_gradient(
            low  = "magenta",
            high = "navyblue",
            name = continuous_var
          ) +
        labs(
          title    = paste("Gene:", g),
          subtitle = paste("Spearman Corr:", round(cor_vals[g], 3)),
          x        = "RNA expression",
          y        = "Protein expression"
        ) +
        theme_bw() +
        theme(
          plot.title     = element_text(face = "bold", size = 14),
          plot.subtitle  = element_text(size = 12),
          legend.position = "right"
        )

        # sp <- ggscatter(pdf_df, x="Protein", y="RNA",
        #                 color="Group", palette=c("firebrick","dodgerblue"),
        #                 size=3, alpha=0.6, ggtheme=theme_bw()) +
        #       ggtitle(paste("Gene:",g)) +
        #       theme(plot.title=element_text(face="bold", size=14)) +
        #       annotate("text", x=-Inf, y=Inf,
        #                label=paste("Spearman:", round(cor_vals[g],3)),
        #                hjust=-0.1, vjust=1.5, size=5)

        gene_plots[[g]] <- p
    }

    # write or display PDF with timestamp
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    pdf_file <- file.path(result_dir, paste0("rpMeta.rpCorGenePlot_", timestamp, ".pdf"))
    if (X11=="No") pdf(pdf_file, width=10, height=8)
    for (p in gene_plots) print(p)
    if (X11=="No") dev.off()

    # return the summary table
    result_df <- data.frame(
        gene_name         = plot_genes,
        RNA_pvalue        = unname(rna_p_vals[plot_genes]),
        Prot_pvalue       = unname(prot_p_vals[plot_genes]),
        Meta_pvalue       = unname(meta_p_vals[plot_genes]),
        Spearman_rho      = unname(cor_vals[plot_genes]),
        Spearman_pvalue   = unname(spearman_p_vals[plot_genes]),
        stringsAsFactors = FALSE
    )
    return(result_df)
}

