#' Create pathway enrichment plots
#'
#' `rpPathwayPlot()` create a combined RNA/protein/meta pathway enrichment analysis bubble plot.
#'
#' @param df Combined RNA-protein differential analysis result data frame.
#' @param result_dir Absolute path where the results should be stored.
#' @param group1 name of the first group used in DE analysis
#' @param group2 name of the second group used in DE analysis
#' @param pathways Pathway database (Hallmark, KEGG, Cell_Type_Signature, IMMUNESIGDB or GO_BP). Default="Hallmark"
#' @param FDR_cutoff false discovery rate (FDR) cutoff for plot, default=0.05
#' @param max_pathways maximum number of pathways to be displayed in the plot, default=30
#' @param X11 "No" (default) to write PDF, or "Yes" to plot to X11 device
#'
#' @return None
#' @export
rpPathwayPlot <- function(df, result_dir, group1, group2, pathways = "Hallmark", FDR_cutoff = 0.05, max_pathways = 30, X11 = "No") {
    library(data.table)
    library(ggplot2)
    library(dplyr)
    library(ggpubr)

    # Ensure result directory exists
    if (!dir.exists(result_dir)) {
        dir.create(result_dir, recursive = TRUE)
    }

    # validate X11
    if (!X11 %in% c("No","Yes")) {
        stop("X11 must be 'No' or 'Yes'")
    }

    # filter results using FDR cutoff 
    combined_df <- df[df$padj < FDR_cutoff, ]

    # limit & reorder by Meta padj to top n, then order by NES
    meta_df <- combined_df %>%
      filter(Category == "Meta") %>%
      arrange(pval)
    if (nrow(meta_df) > max_pathways) meta_df <- head(meta_df, max_pathways)
    meta_df <- meta_df %>% arrange(NES)
    pathway_levels <- meta_df$pathway

    combined_df <- combined_df %>%
      filter(pathway %in% pathway_levels) %>%
      mutate(pathway = factor(pathway, levels = pathway_levels))

    # Generate Bubble Plot
    pdf_file <- file.path(result_dir, paste0("rpMeta.rpPathway.", tolower(pathways), ".pdf"))
    if (X11=="No") pdf(pdf_file, width = 10, height = 12)
    title <- paste0('rpPathway plot (', group2, '/', group1, ')')

    bubble_plot <- ggplot(combined_df, aes(x=Category, y=pathway, size=-log10(pval), color=NES)) +
      geom_point(alpha=0.8) +
      scale_color_gradient2(low="dodgerblue", mid="white", high="firebrick", midpoint=0) +
      scale_size_continuous(range=c(2,10)) +
      labs(title=title, x="", y="Pathway", color="NES", size="-log10(P-value)") +
      theme_minimal() +
      theme(
        axis.ticks.x.top  = element_line(),
        axis.text.x.top   = element_text(face="bold", size=14, vjust=0),
        axis.text.x       = element_blank(),    # remove default bottom x-axis text
        axis.title.x      = element_blank(),
        plot.title        = element_text(face="bold", size=16),
        axis.title.y      = element_text(face="bold", size=14)
      ) +
      scale_x_discrete(position="top", labels=c("RNA","Protein","Joint"))

    print(bubble_plot)
    if (X11=="No") dev.off()

    message("rpPathway plot saved to: ", pdf_file)
}

