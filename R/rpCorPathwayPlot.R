#' Create pathway enrichment plots for correlation analysis 
#'
#' `rpCorPathwayPlot()` creates a combined RNA/protein/meta pathway enrichment
#' analysis bubble plot for correlation analysis.
#'
#' @param df Pathway enrichment result data frame.
#' @param result_dir Absolute path where the results should be stored.
#' @param continuous_var Continuous variable used in the analysis.
#' @param pathways Pathway database (Hallmark, KEGG, Cell_Type_Signature,
#'   IMMUNESIGDB or GO_BP). Default = "Hallmark".
#' @param FDR_cutoff False discovery rate (FDR) cutoff for plot. Default is 0.05.
#' @param max_pathways Maximum number of pathways to be displayed in the plot.
#'   Default is 30.
#' @param X11 "No" (default) to write PDF, or "Yes" to plot to X11 device.
#'
#' @return None
#' @export
rpCorPathwayPlot <- function(df, result_dir, continuous_var, pathways = "Hallmark", 
                             FDR_cutoff = 0.05, max_pathways = 30, X11 = "No") {
    library(data.table)
    library(ggplot2)
    library(dplyr)

    # Ensure dataframe is not NULL or empty
    if (is.null(df) || nrow(df) == 0) {
        warning("The input pathway result is NULL or empty. Skipping plot generation.")
        return(NULL)
    }

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
    meta_df_nes_sorted <- meta_df %>% arrange(NES)
    pathway_levels <- meta_df_nes_sorted$pathway

    # Filter to keep only the pathways that are significant in meta analysis.
    # All pathways in pathway_levels will be displayed (even if RNA or Protein
    # FDR is > FDR_cutoff) to maintain alignment.
    combined_df_filt <- df[df$pathway %in% pathway_levels, ]

    # Explicitly set the levels to ensure the y-axis order matches Meta results
    combined_df_filt$pathway <- factor(combined_df_filt$pathway, levels = pathway_levels)
    # set Category order 
    combined_df_filt$Category <- factor(combined_df_filt$Category, levels = c("RNA", "Protein", "Meta"))
    # rename "Meta" to "Joint"
    combined_df_filt$Category <- gsub("Meta", "Joint", combined_df_filt$Category)
    combined_df_filt$Category <- factor(combined_df_filt$Category, levels = c("RNA", "Protein", "Joint")) 

    # Generate Bubble Plot
    pdf_file <- file.path(result_dir, paste0("rpMeta.rpCorPathway.", tolower(pathways), ".pdf"))
    if (X11=="No") pdf(pdf_file, width = 10, height = 12)
    title <- paste0('rpCorPathway plot (continuous variable: ', continuous_var, ')')

    bubble_plot <- ggplot(combined_df_filt, aes(x=Category, y=pathway, size=-log10(pval), color=NES)) +
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
      scale_x_discrete(position="top")

    print(bubble_plot)
    if (X11=="No") dev.off()

    message("rpCorPathway plot saved to: ", pdf_file)
}
