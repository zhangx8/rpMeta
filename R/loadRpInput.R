#' Load and validate input files for rpMeta
#'
#' `loadRpInput()` reads a sample table, an RNA count matrix and a protein quant
#' matrix (all in CSV or TSV), validates them, reports summaries, finds the
#' overlapping samples, subsets & re‑orders each table, and returns a named list
#' of data.frames.
#'
#' @param sample_file        Path to the sample info file (CSV or TSV)
#' @param rna_count_file     Path to the RNA counts file (CSV or TSV; row names = genes)
#' @param protein_quant_file Path to the protein quant file (CSV or TSV; row names = genes)
#' @param format             `"csv"` or `"tsv"` (all three files use this format)
#' @param max_na_frac        maximum fraction of NAs in each row(gene)
#'
#' @return A named list of three data.frames: samples, rna, protein
#' @examples
#' # inputs <- loadRpInput("samples.tsv", "counts.tsv", "prot.tsv", format="tsv")
#' @export
loadRpInput <- function(
  sample_file,
  rna_count_file,
  protein_quant_file,
  format = c("tsv", "csv"),
  max_na_frac = 0.5
) {
  format <- match.arg(format)
  delim  <- if (format=="csv") "," else "\t"
  read_df <- function(f, row.names=FALSE) {
    if (format=="csv") {
      read.csv(f, header=TRUE, stringsAsFactors=FALSE,
               check.names=FALSE,
               row.names = if (row.names) 1 else NULL)
    } else {
      read.table(f, header=TRUE, sep="\t", stringsAsFactors=FALSE,
                 check.names=FALSE,
                 row.names = if (row.names) 1 else NULL)
    }
  }

  ## Step 1: read in sample, rna and protein files 
  samples  <- read_df(sample_file,      row.names=FALSE)
  rna      <- read_df(rna_count_file,   row.names=TRUE)
  protein  <- read_df(protein_quant_file, row.names=TRUE)

  ## Step 1.5: remove rows with NA fraction >= max_na_frac
  num_col_rna <- ncol(rna)
  max_num_na_rna <- floor(num_col_rna*max_na_frac)
  num_na_rna <- rowSums(is.na(rna))
  rna <- rna[ num_na_rna <= max_num_na_rna, ] 
  num_col_prot <- ncol(protein)
  max_num_na_prot <- floor(num_col_prot*max_na_frac)
  num_na_prot <- rowSums(is.na(protein))
  protein <- protein[ num_na_prot <= max_num_na_prot, ]

  ## Step 2: Validate 
  # 2a) sample file
  # reqCols <- c("Sample_ID","Group")
  reqCols <- c("Sample_ID")
  if (!all(reqCols %in% colnames(samples))) {
    stop("`sample_file` must contain columns: ", paste(reqCols,collapse=", "))
  }
  # 2b) count matrices: must have rownames and colnames
  if (is.null(rownames(rna))   || ncol(rna)==0)   stop("RNA file must have genes as rownames and samples as columns")
  if (is.null(rownames(protein)) || ncol(protein)==0) stop("Protein file must have genes as rownames and samples as columns")
  # 2c) numeric entries
  if (!all(vapply(rna, is.numeric, TRUE)))
    stop("All RNA count columns must be numeric")
  if (!all(vapply(protein, is.numeric, TRUE)))
    stop("All protein quant columns must be numeric")

  ## Step 3:  Summaries before filtering
  message(">>> Raw input summary:")
  message("  RNA:      genes = ", nrow(rna),   ", samples = ", ncol(rna))
  message("  Protein:  genes = ", nrow(protein),", samples = ", ncol(protein))
  message("  Samples:  total = ", nrow(samples),
          ", groups = ", paste0(names(table(samples$Group)), "(", table(samples$Group),")", collapse=", "))

  ## Step 4: Find overlapping samples
  sampIDs <- samples$Sample_ID
  common  <- Reduce(intersect, list(sampIDs, colnames(rna), colnames(protein)))
  message(">>> Overlapping samples: ", length(common),
          " (out of ", length(sampIDs),")\n  ", paste(common, collapse=", "))

  if (length(common)==0) stop("No samples in common!")

  ## Step 5:  Subset & reorder
  samples_f <- samples[match(common, samples$Sample_ID), , drop=FALSE]
  rna_f     <- rna   [    , common, drop=FALSE]
  protein_f <- protein[, common, drop=FALSE]

  ## Step 6: Final summaries
  message(">>> Final filtered summary:")
  message("  RNA:      genes = ", nrow(rna_f),   ", samples = ", ncol(rna_f))
  message("  Protein:  genes = ", nrow(protein_f),", samples = ", ncol(protein_f))
  message("  Samples:  total = ", nrow(samples_f),
          ", groups = ", paste0(names(table(samples_f$Group)), "(", table(samples_f$Group),")", collapse=", "))

  ## Step 7: Return data frames 
  invisible(list(
    samples = samples_f,
    rna     = rna_f,
    protein = protein_f
  ))
}

