# rpVolcano.R
#' Plot RNA-protein volcano plot 
#' 
#' `rpVolcano()` plots RNA-protein volcano plots and fold change.
#'
#' @param rpDEG Combined RNA protein differential analysis result data frame.
#' @param result_dir Absolute path where the results should be stored.
#' @param group1 name of the first group used in DE analysis
#' @param group2 name of the second group used in DE analysis
#' @param FDR_cutoff FDR cutoff in volcano plot, Default = 0.05.
#' @param log2FC_cutoff log2 Fold Change cutoff in volcano plot, Default = 1 (FoldChange = 2).
#' @param genesToDisplay  a numeric value (default 5) indicating the top n up‐ and down‐regulated genes (based on Meta_FDR) to annotate.
#' @param X11 "No" (default) to write PDF, or "Yes" to plot to X11 device
#'
#' @return Null
#' @examples
#' deg_file <- system.file("extdata", "DEG.combined.txt", package = "rpMeta")
#' m <- read.table(deg_file, header = TRUE, stringsAsFactors = FALSE)
#' result_dir <- file.path(getwd(), "rpMeta_APOLLO1")
#' group1 <- "TRU"
#' group2 <- "ProxProlif"
#' if (!dir.exists(result_dir)) { dir.create(result_dir) }
#' rpVolcano(m, result_dir, group1, group2, genesToDisplay = 5, X11="Yes")
#' @export
rpVolcano <- function(rpDEG, result_dir, group1, group2, FDR_cutoff = 0.05, log2FC_cutoff = 1, genesToDisplay = 5, X11 = "No") {
  
  # Ensure output directory exists
  if (!dir.exists(result_dir)) {
    dir.create(result_dir, recursive = TRUE)
  }

  # validate X11
  if (!X11 %in% c("No","Yes")) {
      stop("X11 must be 'No' or 'Yes'")
  }
  
  df <- rpDEG
  
  pdf_file <- file.path(result_dir, "rpMeta.rpVolcano.pdf")
  if (X11=="No") pdf(pdf_file, width = 10, height = 12)

  # leave room in the top outer margin for a main title:
  par(oma = c(0, 0, 2, 0))
  
  # Set layout: panel 1 spans the top row; panels 2 and 3 are in the bottom row
  layout(matrix(c(1, 1, 2, 3), nrow = 2, byrow = TRUE))

  ## Panel 1: Scatter Plot of RNA vs Protein log2 Fold Change
  valid_indices <- complete.cases(df$RNA_log2FC, df$Prot_log2FC, df$Meta_FDR)
  rna_log2FC <- df$RNA_log2FC[valid_indices]
  prot_log2FC <- df$Prot_log2FC[valid_indices]
  meta_fdr <- df$Meta_FDR[valid_indices]
  
  colors <- ifelse(meta_fdr < 0.05, "dodgerblue", "grey")
  title1 <- paste0('Joint Fold Change (', group2, '/', group1, ')')
  plot(prot_log2FC, rna_log2FC,
       xlab = "Protein log2 Fold Change",
       ylab = "RNA log2 Fold Change",
       main = title1,
       pch = 16, cex = 0.5, col = colors)
  legend("bottomright", 
         legend = c("DEG (Meta FDR < 0.05)", "None DEG"),
         col = c("dodgerblue", "grey"), pch = 16, pt.cex = 0.5)

  # draw the main title in the outer margin
  mtext("rpVolcano", outer = TRUE, cex = 1.5, font = 2)
  
  model <- lm(rna_log2FC ~ prot_log2FC)
  abline(model, col = "firebrick", lwd = 2)
  abline(h = 0, v = 0, col = "black", lty = 2)
  
  spearman_corr <- cor(rna_log2FC, prot_log2FC, method = "spearman")
  legend("topleft", legend = paste("Spearman Corr:", round(spearman_corr, 2)),
         bty = "n", col = "black", text.col = "black")
  
  quad1 <- sum(rna_log2FC > 0 & prot_log2FC > 0)
  quad2 <- sum(rna_log2FC < 0 & prot_log2FC > 0)
  quad3 <- sum(rna_log2FC < 0 & prot_log2FC < 0)
  quad4 <- sum(rna_log2FC > 0 & prot_log2FC < 0)
  
  text(x = max(prot_log2FC) * 0.75, y = max(rna_log2FC) * 0.75, labels = quad1, col = "black", cex = 1.2)
  text(x = min(prot_log2FC) * 0.75, y = max(rna_log2FC) * 0.75, labels = quad2, col = "black", cex = 1.2)
  text(x = min(prot_log2FC) * 0.75, y = min(rna_log2FC) * 0.75, labels = quad3, col = "black", cex = 1.2)
  text(x = max(prot_log2FC) * 0.75, y = min(rna_log2FC) * 0.75, labels = quad4, col = "black", cex = 1.2)

  ## Panel 2: RNA Volcano Plot
  # Color rules for RNA volcano plot:
  # Downregulated RNA DEGs: if Meta_FDR < 0.05 -> "navyblue", else "dodgerblue"
  # Upregulated RNA DEGs: if Meta_FDR < 0.05 -> "firebrick", else "tomato"
  plot(df$RNA_log2FC, -log10(df$RNA_Pval), xlim = c(-5, 5),
       col = ifelse(df$RNA_FDR < FDR_cutoff & df$RNA_log2FC < -log2FC_cutoff & df$Meta_FDR < FDR_cutoff, "navyblue",
                    ifelse(df$RNA_FDR < FDR_cutoff & df$RNA_log2FC < -log2FC_cutoff, "dodgerblue",
                           ifelse(df$RNA_FDR < FDR_cutoff & df$RNA_log2FC > log2FC_cutoff & df$Meta_FDR < FDR_cutoff, "firebrick",
                                  ifelse(df$RNA_FDR < FDR_cutoff & df$RNA_log2FC > log2FC_cutoff, "tomato", "grey")))),
       pch = 16, cex = 0.6,
       xlab = "RNA log2 Fold Change", ylab = "-log10(RNA P-value)",
       main = "RNA")
  
  abline(h = -log10(0.01), col = "black", lty = 2)
  abline(v = c(-log2FC_cutoff, log2FC_cutoff), col = "black", lty = 2)
  
  # Gene Annotation for RNA Volcano Plot
  upGenes <- subset(df, RNA_FDR < FDR_cutoff & RNA_log2FC > 0)
  upGenes <- upGenes[order(upGenes$Meta_FDR), ]
  upGenes <- head(upGenes, genesToDisplay)
  
  downGenes <- subset(df, RNA_FDR < FDR_cutoff & RNA_log2FC < 0)
  downGenes <- downGenes[order(downGenes$Meta_FDR), ]
  downGenes <- head(downGenes, genesToDisplay)
  
  selected_genes <- rbind(upGenes, downGenes)
  
  if (nrow(selected_genes) > 0) {
    for (i in 1:nrow(selected_genes)) {
      xg <- selected_genes$RNA_log2FC[i]
      yg <- -log10(selected_genes$RNA_Pval[i])
      # Highlight the gene's point with a circle
      points(xg, yg, pch = 21, col = "black", bg = adjustcolor("orange", alpha.f = 0.2), cex = 1.8, lwd = 2)
      # points(xg, yg, pch = 21, col = "orange", cex = 1.8, lwd = 2)
      offset_x <- if (xg >= 0) 0.8 else -0.8
      offset_y <- 0.5
      arrows(x0 = xg + offset_x, y0 = yg + offset_y, x1 = xg, y1 = yg, length = 0.1, col = "black", lwd = 1)
      text(x = xg + offset_x, y = yg + offset_y, labels = as.character(selected_genes$Gene[i]),
           pos = if (xg >= 0) 4 else 2, cex = 0.8, col = "black")
    }
  }
 
  ## Panel 3: Protein Volcano Plot
  plot(df$Prot_log2FC, -log10(df$Prot_Pval), xlim = c(-5, 5),
       col = ifelse(df$Prot_FDR < FDR_cutoff & df$Prot_log2FC < -log2FC_cutoff & df$Meta_FDR < FDR_cutoff, "navyblue",
                    ifelse(df$Prot_FDR < FDR_cutoff & df$Prot_log2FC < -log2FC_cutoff, "dodgerblue",
                           ifelse(df$Prot_FDR < FDR_cutoff & df$Prot_log2FC > log2FC_cutoff & df$Meta_FDR < FDR_cutoff, "firebrick",
                                  ifelse(df$Prot_FDR < FDR_cutoff & df$Prot_log2FC > log2FC_cutoff, "tomato", "grey")))),
       pch = 16, cex = 0.6,
       xlab = "Protein log2 Fold Change", ylab = "-log10(Protein P-value)",
       main = "Protein")
 
  abline(h = -log10(0.01), col = "black", lty = 2)
  abline(v = c(-log2FC_cutoff, log2FC_cutoff), col = "black", lty = 2)
  
  # Gene Annotation for Protein Volcano Plot
  protUp <- subset(df, Prot_FDR < FDR_cutoff & Prot_log2FC > 0)
  protUp <- protUp[order(protUp$Meta_FDR), ]
  protUp <- head(protUp, genesToDisplay)
  
  protDown <- subset(df, Prot_FDR < FDR_cutoff & Prot_log2FC < 0)
  protDown <- protDown[order(protDown$Meta_FDR), ]
  protDown <- head(protDown, genesToDisplay)
  
  selected_prot <- rbind(protUp, protDown)
  
  if (nrow(selected_prot) > 0) {
    for (i in 1:nrow(selected_prot)) {
      xp <- selected_prot$Prot_log2FC[i]
      yp <- -log10(selected_prot$Prot_Pval[i])
      # Highlight the gene's point with a circle
      points(xp, yp, pch = 21, col = "black", bg = adjustcolor("orange", alpha.f = 0.2), cex = 1.8, lwd = 2)
      # points(xp, yp, pch = 21, col = "black", bg = "yellow", cex = 1.8, lwd = 2)
      offset_x <- if (xp >= 0) 0.8 else -0.8
      offset_y <- 0.5
      arrows(x0 = xp + offset_x, y0 = yp + offset_y, x1 = xp, y1 = yp, length = 0.1, col = "black", lwd = 1)
      text(x = xp + offset_x, y = yp + offset_y, labels = as.character(selected_prot$Gene[i]),
           pos = if (xp >= 0) 4 else 2, cex = 0.8, col = "black")
    }
  }
  
  if (X11=="No") dev.off()
  cat("rpVolcano plots saved to:", pdf_file, "\n")
}

