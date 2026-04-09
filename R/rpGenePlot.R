# rpGenePlot.R
#' `rpGenePlot()` plot RNA–protein scatterplot and boxplot for top 20 DEGs
#'
#' @param df                combined RNA–protein differential analysis result data frame 
#' @param rna_norm_counts   an RNA count data frame with samples as columns and gene symbols as rows
#' @param prot_norm_counts  an proteomics count data frame with samples as columns and gene symbols as rows
#' @param samples           sample description data frame with group and covariates columns
#' @param group1            name of the first group used in DE analysis
#' @param group2            name of the second group used in DE analysis
#' @param result_dir        absolute path where the results should be stored, default="./"
#' @param genes             a vector of gene symbols, e.g. c("MET","BRCA2") or top n by Meta_FDR (default=20)
#' @param X11               "No" (default) to write PDF, or "Yes" to plot to X11 device
#'
#' @return A data.frame with columns \code{gene_name}, \code{RNA_pvalue}, \code{Prot_pvalue}, \code{Meta_pvalue}, \code{Spearman_rho}, \code{Spearman_pvalue}
#' @examples
#' deg_file <- system.file("extdata","DEG.combined.txt",package="rpMeta")
#' m <- read.table(deg_file,header=TRUE,stringsAsFactors=FALSE)
#' # rpGenePlot(df=m, rna_norm_counts=rna_norm_counts, prot_norm_counts=prot_norm_counts,
#' #            samples=samples, group1, group2, result_dir)
#' @export
rpGenePlot <- function(df, rna_norm_counts, prot_norm_counts, samples, group1, group2, result_dir="./", genes = 20, X11 = "No") {
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

        # pull p-values from the merged df
        rna_p_vals[g]  <- df$RNA_Pval[df$Gene == g]
        prot_p_vals[g] <- df$Prot_Pval[df$Gene == g]
        meta_p_vals[g] <- df$Meta_PVal[df$Gene == g]

        # compute Spearman correlation and p-value
        test <- cor.test(rvals, pvals, method="spearman", exact=FALSE)
        cor_vals[g]        <- test$estimate
        spearman_p_vals[g] <- test$p.value

        pdf_df <- data.frame(Sample=colnames(rna_norm_counts), RNA=rvals, Protein=pvals, Group=samples$Group)

        sp <- ggscatter(pdf_df, x="RNA", y="Protein",
                        color="Group", palette=c("firebrick","dodgerblue"),
                        size=3, alpha=0.6, ggtheme=theme_bw()) +
              ggtitle(paste("Gene:",g)) +
              theme(plot.title=element_text(face="bold", size=14)) +
              annotate("text", x=-Inf, y=Inf,
                       label=paste("Spearman:", round(cor_vals[g],3)),
                       hjust=-0.1, vjust=1.5, size=5)

        xpl <- ggboxplot(pdf_df, x="Group", y="RNA",
                         xlab="", ylab="Protein",
                         color = "Group", fill = "Group",
                         palette=c("dodgerblue","firebrick"),
                         alpha=0.5, ggtheme=theme_bw()) +
               theme_minimal() +
               rremove("legend") +
               theme(axis.title.x=element_blank(),
                     axis.text.x=element_blank(),
                     axis.ticks.x=element_blank()) +
               coord_flip()

        ypl <- ggboxplot(pdf_df, x="Group", y="Protein",
                         xlab="", ylab="RNA",
                         color = "Group", fill = "Group",
                         palette=c("dodgerblue","firebrick"),
                         alpha=0.5, ggtheme=theme_bw()) +
               theme_minimal() +
               rremove("legend") +
               theme(axis.title.y=element_blank(),
                     axis.text.y=element_blank(),
                     axis.ticks.y=element_blank())

        # assemble a small panel with all three p‐values in top right
        pp_text <- paste0(
            "RNA p=",   signif(rna_p_vals[g],3), "\n",
            "Prot p=",  signif(prot_p_vals[g],3), "\n",
            "rpMeta p=",signif(meta_p_vals[g],3)
        )
        blank <- cowplot::ggdraw() +
                 cowplot::draw_label(pp_text,
                                     x=1, y=1, hjust=1, vjust=1,
                                     fontface="italic", size=12)

        gene_plots[[g]] <- plot_grid(xpl, blank, sp, ypl,
                                     ncol=2, align="hv",
                                     rel_widths=c(2,1),
                                     rel_heights=c(1,2))
    }

    # write or display PDF with timestamp
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    pdf_file <- file.path(result_dir, paste0("rpMeta.rpGenePlot_", timestamp, ".pdf"))
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

