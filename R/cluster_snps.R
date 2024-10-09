#' Clustering SNPs based on p-value and proximity with a PLINK C+T file.
#'
#' This function clusters SNPs within a given window size based on their P-value and proximity. It iterates through
#' each chromosome, finds the SNP with the smallest P-value, and groups all SNPs within the specified window size
#' around this SNP into a cluster.
#'
#' @param df A data.frame containing SNP data with columns for SNP (SNP ID), CHR (chromosome), BP (base pair position), and P (p-value).
#' @param window_size An integer specifying the window size around each SNP (in base pairs) within which other SNPs are considered for clustering. Default to 1e6.
#' @return A data.frame containing the clustered SNPs with an additional column 'ClusterSize' indicating the number of SNPs in each cluster.
#'
#' @details The function processes each chromosome independently. It orders the SNPs by their base pair positions, identifies
#'          the SNP with the smallest P-value, and clusters all SNPs within the specified window size around this SNP. The process
#'          is repeated until all SNPs are assigned to a cluster.
#'
#' @examples
#' df <- data.frame(SNP=c("rs1", "rs2", "rs3", "rs4", "rs5"),
#'                  CHR=c(1, 1, 1, 1, 2),
#'                  BP=c(100000, 150000, 200000, 250000, 300000),
#'                  P=c(0.01, 0.02, 0.03, 0.04, 0.05))
#' window_size <- 50000
#' clustered_snps <- cluster_snps(df, window_size)
#'
#' @export
cluster_snps <- function(df, window_size=1e6) {
  df <- df[order(df$CHR, df$BP), ]
  result <- data.frame(SNP = character(), CHR = integer(), BP = integer(), P = numeric(), ClusterSize = integer())
  for (chr in unique(df$CHR)) {
    df_chr <- df[df$CHR == chr, ]
    while (nrow(df_chr) > 0) {
      min_p_snp <- df_chr[which.min(df_chr$P), ]
      window_snps <- df_chr[abs(df_chr$BP - min_p_snp$BP) <= window_size, ]
      df_chr <- df_chr[!df_chr$SNP %in% window_snps$SNP, ]
      result <- rbind(result, c(min_p_snp, ClusterSize = nrow(window_snps)))
    }
  }
  return(result)
}
