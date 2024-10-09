#' Clustering second data frame based on closest SNP centers from first data frame
#'
#' This function performs clustering of SNPs in a second data.frame based on the closest SNP centers defined in a first data.frame.
#' Both data.frames should include SNP, BP, and CHR columns. This function scales CHR and BP to ensure distinctiveness across chromosomes
#' and employs Euclidean distance to find the nearest cluster centers from the first data.frame for each SNP in the second data.frame.
#'
#' @param df1 A data.frame representing the output of a plink clump with parameters r2=0.01.
#'           It contains columns for SNP, BP (base pair position), CHR (chromosome), and P (p-value).
#' @param df2 A data.frame similar to df1, representing a plink output with a less stringent r2 value, typically r2=0.5,
#'           including columns for SNP, BP, CHR, and P.
#' @return A modified version of df2 where each SNP is annotated with a 'cluster' index corresponding to the closest
#'         SNP center from df1 based on scaled CHR and BP values.
#'
#' @details The function first standardizes the CHR and BP columns by multiplying CHR by 10000 and dividing BP by 1e6.
#'          This standardization helps to manage the scale differences between chromosome numbers and base pair positions.
#'          After standardization, it calculates the Euclidean distances between each SNP in df2 to all SNP centers in df1,
#'          assigns each SNP in df2 to the nearest center from df1, and adds a new column 'cluster' to df2 to reflect this assignment.
#'
#' @examples
#' df1 <- data.frame(SNP=c("rs1", "rs2"), CHR=c(1, 1), BP=c(150000, 250000), P=c(0.001, 0.002))
#' df2 <- data.frame(SNP=c("rs1", "rs3", "rs2", "rs4"), CHR=c(1,1,1,1),
#'                   BP=c(150000,160000,250000,260000),
#'                   P=c(0.001,0.003,0.002, 0.004))
#' clustered_df2 <- clump_cluster(df1, df2)
#'
#' @importFrom stats dist
#' @export
clump_cluster <- function(df1, df2) {

df1$CHR_scaled <- df1$CHR * 10000
df1$BP_scaled <- df1$BP / 1e6

df2$CHR_scaled <- df2$CHR * 10000
df2$BP_scaled <- df2$BP / 1e6

coords1 <- df1[, c("CHR_scaled", "BP_scaled")]
coords2 <- df2[, c("CHR_scaled", "BP_scaled")]

distances <- as.matrix(dist(rbind(coords1, coords2)))[-(1:nrow(coords1)), 1:nrow(coords1)]

cluster_indices <- apply(distances, 1, which.min)

df2$cluster <- cluster_indices

return(df2)
}
