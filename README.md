# rpMeta: Integrated RNA and Protein Meta-Analysis

rpMeta is an R package designed for the joint differential expression analysis and visualization of integrated RNA-Seq (transcriptomics) and mass spectrometry-based proteomics data.

# Overview

While various statistical tools exist to tackle differential expression for a single platform (e.g., edgeR, limma, DESeq2), rpMeta provides a unified framework to leverage both transcription and translation information.

By performing joint analysis of RNA-Seq and MS-based proteomics datasets, rpMeta helps researchers identify biologically relevant results with increased statistical power, particularly in studies with small sample sizes.

# Key Features

- Integrated DEG Analysis: Jointly test for differential expression between groups or associations with continuous variables.

- Flexible Meta-Analysis: Combine P-values from RNA and protein datasets using multiple methods (Fisher, Stouffer, etc.).

- Rich Visualizations: Generate publication-ready figures including:

  - rpHeatmaps: Adjacent RNA/Protein expression views.

  - rpVolcano & rpVenn Plots: Visualizing single-platform and joint significance.

  - rpGenePlot: Correlation scatter plots and boxplots for specific genes.

  - Pathway Analysis: Integrated GSEA bubble plots (rpPathwayPlot).

# Installation

rpMeta is tested on R version 4.3 or higher.

# 1. Install Dependencies

Before installing rpMeta, ensure the following Bioconductor and CRAN packages are installed:

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor dependencies
BiocManager::install(c("limma", "DESeq2", "pheatmap", "msigdbr", "fgsea"))

# Install CRAN dependencies
install.packages(c("ggplot2", "gridExtra", "VennDiagram", "ggpubr", "cowplot", "data.table", "dplyr"))


# 2. Install rpMeta

You can install the package directly from the source tarball available in the Releases section.

**From the R console:**

install.packages("rpMeta_1.0.0.tar.gz", repos = NULL, type = "source")


**From the Terminal:**

R CMD INSTALL rpMeta_1.0.0.tar.gz


# Documentation

A comprehensive tutorial including example data analysis is available as a Word document in our [Latest Release](https://github.com/zhangx8/rpMeta/releases/latest).

**License:** MIT License

**Authors:** Xijun Zhang and Matthew D. Wilkerson
