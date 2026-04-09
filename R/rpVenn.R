#' Plot Venn diagram of RNA, protein and meta DEGs
#' 
#' `rpVenn` Plot Venn diagram of RNA, protein and meta DEGs
#'
#' @param rpDEG combined RNA protein differential analysis result data frame 
#' @param fdr_cutoff The FDR threshold RNA, protein and meta DEGs in rpVenn plot 
#' @param result_dir absolute path where the results should be stored, default="./"
#' @param filename_suffix file name suffix of the output pdf file  
#' @param X11 "No" (default) to write PDF, or "Yes" to plot to X11 device
#' 
#' @return A list of counts, RNA, protein, meta, two way overlaps, three way overlaps
#' @examples
#' deg_file <- system.file("extdata", "DEG.combined.txt", package = "rpMeta")
#' rpDEG <- read.table(deg_file, header=TRUE, stringsAsFactors=FALSE)
#' result_dir <- file.path(getwd(), "rpMeta_APOLLO1")
#' if (!dir.exists(result_dir)) { dir.create(result_dir) }
#' rpVenn(rpDEG, fdr_cutoff=0.05, result_dir)
#' @export
rpVenn <- function(rpDEG, fdr_cutoff=0.05, result_dir="./", filename_suffix="", X11="No") {
    library(VennDiagram)
    df <- rpDEG

    # validate X11
    if (!X11 %in% c("No","Yes")) {
        stop("X11 must be 'No' or 'Yes'")
    }

    # Filter RNA, protein, and meta DEGs based on FDR < 0.05, remove NAs
    deg_rna <- na.omit(df$Gene[df$RNA_FDR < fdr_cutoff])
    deg_prot <- na.omit(df$Gene[df$Prot_FDR < fdr_cutoff])
    deg_meta <- na.omit(df$Gene[df$Meta_FDR < fdr_cutoff])
    venn_values <- list(
      RNA     = deg_rna,
      Protein = deg_prot,
      Meta    = deg_meta
    )

    # Check if there are any DEGs in each set
    if (length(deg_rna) == 0 || length(deg_prot) == 0 || length(deg_meta) == 0) {
        stop("One or more DEG sets are empty. Venn diagram cannot be created.")
    }

    # Create a Venn diagram of the three sets
    venn_plot <- venn.diagram(
        x = venn_values, 
        category.names = c("RNA", "Protein", "Meta"),
        filename = NULL,
        output = TRUE,
        fill = c("cornflowerblue", "green", "red"),
        alpha = 0.5,
        cex = 2,
        cat.cex = 2,
        cat.fontface = "bold",
        main = "rpVenn",
        main.cex = 2
    )

    # Save the plot to a PDF file
    if (filename_suffix == "") {
        venn_pdf <- paste0(result_dir, "/rpMeta.rpVenn.pdf")
    } else {
        venn_pdf <- paste0(result_dir, "/rpMeta.rpVenn.", filename_suffix, ".pdf")
    }
    if (X11=="No") pdf(venn_pdf)
    grid.draw(venn_plot)
    dev.off()
    cat("rpVenn plot saved to:", venn_pdf, "\n")

    # basic sizes
    n_rna  <- length(deg_rna)
    n_prot <- length(deg_prot)
    n_meta <- length(deg_meta)

    # pairwise overlaps
    rna_prot  <- intersect(deg_rna, deg_prot)
    rna_meta  <- intersect(deg_rna, deg_meta)
    prot_meta <- intersect(deg_prot, deg_meta)

    n_rna_prot  <- length(rna_prot)
    n_rna_meta  <- length(rna_meta)
    n_prot_meta <- length(prot_meta)

    # three‐way overlap
    all_three   <- Reduce(intersect, list(deg_rna, deg_prot, deg_meta))
    n_all_three <- length(all_three)

    # bundle into a named vector (or data.frame)
    venn_counts <- c(
      RNA           = n_rna,
      Protein       = n_prot,
      Meta          = n_meta,
      RNA_Protein   = n_rna_prot,
      RNA_Meta      = n_rna_meta,
      Protein_Meta  = n_prot_meta,
      All_Three     = n_all_three
    )

    # return the counts
    venn_counts
}

