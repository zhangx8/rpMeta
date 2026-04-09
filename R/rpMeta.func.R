#' Combine RNA and proteomics DE analysis P values
#'
#' This function combines two p-values (RNA and protein) using one of several
#' meta-analysis methods:
#'
#' 1. **Fisher** (default): 
#'    \deqn{X^2 = -2 \sum_{i=1}^{k} \ln(p_i) \sim \chi^2_{2k}}
#'
#' 2. **Pearson (sensitive to large P values)**: 
#'    \deqn{X^2 = -2 \sum_{i=1}^{k} \ln(1 - p_i) \sim \chi^2_{2k}}
#'
#' 3. **Stouffer**: 
#'    \deqn{Z = \frac{1}{\sqrt{k}} \sum_{i=1}^{k} Z_i, \quad Z_i = \Phi^{-1}(1 - p_i)}
#'    The resulting \eqn{Z} follows approximately \eqn{N(0,1)}.
#'
#' 4. **Weighted_Stouffer**: 
#'    \deqn{Z = \frac{\sum_{i=1}^{k} w_i Z_i}{\sqrt{\sum_{i=1}^{k} w_i^2}}, 
#'          \quad Z_i = \Phi^{-1}(1 - p_i)}
#'    where weights \eqn{w_i} are provided as \code{rna_weight} and \code{prot_weight}.
#'    The resulting \eqn{Z} follows approximately \eqn{N(0,1)}.
#'
#' 5. **Tippet** (minimum p-value): 
#'    \deqn{\min(p_1, \dots, p_k) \sim \mathrm{Beta}(1, k) \quad 
#'          \Rightarrow \quad p_{\text{combined}} = 1 - [1 - \min(p_i)]^k}
#'
#' 6. **Wilkinson** (maximum p-value): 
#'    \deqn{\max(p_1, \dots, p_k) \sim \mathrm{Beta}(k, 1) \quad
#'          \Rightarrow \quad p_{\text{combined}} = [\max(p_i)]^k}
#'
#' @param pvals A named numeric vector of length 2, with names "RNA_Pval" and "Prot_Pval".
#' @param method A character string specifying the method to combine the p-values.
#'   One of \code{c("Fisher", "Pearson", "Stouffer", "Weighted_Stouffer", "Tippet", "Wilkinson")}.
#'   Defaults to \code{"Fisher"}.
#' @param rna_weight Numeric weight for the RNA p-value when using \code{method = "Weighted_Stouffer"}.
#'   Default is 1.
#' @param prot_weight Numeric weight for the protein p-value when using \code{method = "Weighted_Stouffer"}.
#'   Default is 1.
#'
#' @return A single numeric value representing the combined p-value.
#'
#' @examples
#' example_data <- c(RNA_Pval = 0.05, Prot_Pval = 0.01)
#' combinep(example_data, method = "Fisher")
#' combinep(example_data, method = "Pearson")
#' combinep(example_data, method = "Stouffer")
#' combinep(example_data, method = "Weighted_Stouffer", rna_weight = sqrt(30), prot_weight = sqrt(20))
#' combinep(example_data, method = "Tippet")
#' combinep(example_data, method = "Wilkinson")
#'
#' @references
#' Yoon, S., Baik, B., Park, T. et al. (2021) "Powerful p-value combination methods to detect 
#' incomplete association." Sci Rep 11, 6980. \doi{10.1038/s41598-021-86465-y}
#'
#' Toro-Domínguez D et al. (2021) "A survey of gene expression meta-analysis: methods and applications." 
#' Brief Bioinform. 22(2):1694-1705.
#'
#' @export
combinep <- function(row,
                     method = c("Fisher", "Pearson", "Stouffer", 
                                "Weighted_Stouffer", "Tippet", "Wilkinson"),
                     rna_weight = 1,
                     prot_weight = 1) {
  
    method <- match.arg(method)
    
    # Extract p-values
    p_rna <- as.numeric(row["RNA_Pval"])
    p_prot <- as.numeric(row["Prot_Pval"])
    if (is.na(p_rna)) {
        return(p_prot)
    }
    if (is.na(p_prot)) {
        return(p_rna)
    }
    
    p_vec <- c(p_rna, p_prot)
    k <- length(p_vec)  # should be 2 in this case
    
    # Combine p-values based on the selected method
    p_combined <- switch(
      method,
      
      "Fisher" = {
	# Fisher's method: -2 * sum(log(p_i)) ~ chi-squared with 2k degrees of freedom
	stat <- -2 * sum(log(p_vec))
	df <- 2 * k
	pchisq(stat, df = df, lower.tail = FALSE)
      },
      
      "Pearson" = {
	# Pearson's method: -2 * sum(log(1 - p_i)) ~ chi-squared with 2k degrees of freedom
	stat <- -2 * sum(log(1 - p_vec))
	df <- 2 * k
	pchisq(stat, df = df, lower.tail = FALSE)
      },
      
      "Stouffer" = {
	# Stouffer's method: convert p-values to Z-scores then combine
	z_vals <- stats::qnorm(1 - p_vec)
	z_sum <- sum(z_vals) / sqrt(k)
	stats::pnorm(z_sum, lower.tail = FALSE)
      },
      
      "Weighted_Stouffer" = {
	# Weighted Stouffer's method: weight the Z-scores by rna_weight and prot_weight
	w <- c(rna_weight, prot_weight)
	z_vals <- stats::qnorm(1 - p_vec)
	z_comb <- sum(w * z_vals) / sqrt(sum(w^2))
	stats::pnorm(z_comb, lower.tail = FALSE)
      },
      
      "Tippet" = {
	# Tippet's method: use the minimum p-value adjusted for the number of tests
	p_min <- min(p_vec)
	1 - (1 - p_min)^k
      },
      
      "Wilkinson" = {
	# Wilkinson's method: use the maximum p-value adjusted for the number of tests
	p_max <- max(p_vec)
	(p_max)^k
      }
    )
    
    return(p_combined)
}

# Set min and max value for a normalized expression matrix for heatmap visualization 
minmax <- function(x, limit=2.0) {
    min = -abs(limit)
    max = abs(limit)
    if( x<=min ) {
        return(min)
    } else if ( x<max && x>min ) {
        return(x)
    } else {
        return(max);
    }
}

