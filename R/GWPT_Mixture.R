#' @title Genome-Wide Pleiotropy Test for mixture model
#' @description
#' This function performs a genome-wide pleiotropy test (GWPT) after Mendelian randomization. It offers an option for a two-mixture model, where the residual is chosen as the smaller one resulting from the two causal effect estimates from two mixtures.
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bX A matrix of effect estimates from the exposure GWAS.
#' @param bXse A matrix of standard errors of effect estimates from the exposure GWAS.
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last column corresponds to the outcome.
#' @param theta1 The causal effect estimate of the first mixture.
#' @param theta.cov1 The covariance matrix of the causal effect estimate of the first mixture.
#' @param theta2 The causal effect estimate of the second mixture.
#' @param theta.cov2 The covariance matrix of the causal effect estimate of the second mixture.
#' @param LD.blcok A vector of indices of LD blocks.
#' @import data.table
#' @return A list with two components:
#' \item{BETA}{The estimated residual values.}
#' \item{SE}{The standard errors of the residual estimates.}
#'
#' @export

GWPT_Mixture <- function(by, byse, bX, bXse, Rxy, theta1, theta.cov1, theta2, theta.cov2, LD.blcok) {
bZ <- cbind(by, bX)
bZse <- cbind(byse, bXse)
vartheta1 <- c(1, -theta1)
vartheta2 <- c(1, -theta2)
residual1 <- c(bZ %*% vartheta1)
residual2 <- c(bZ %*% vartheta2)
var_residual1 <- var_residual2 <- residual1 * 0
Covtheta1 <- diag(length(vartheta1) - 1) * 0
Covtheta1[which(theta1 != 0), which(theta1 != 0)] <- theta.cov1[which(theta1 != 0), which(theta1 != 0)]
Covtheta1 <- as.matrix(Matrix::bdiag(0, Covtheta1))
Covtheta2 <- diag(length(vartheta2) - 1) * 0
Covtheta2[which(theta2 != 0), which(theta2 != 0)] <- theta.cov2[which(theta2 != 0), which(theta2 != 0)]
Covtheta2 <- as.matrix(Matrix::bdiag(0, Covtheta2))

cat("Calculating the residuals of two mixtures:\n")
pb <- txtProgressBar(min = 0, max = length(var_residual1), style = 3)
for (i in 1:nrow(bZ)) {
setTxtProgressBar(pb, i)
beta <- as.vector(bZ[i, ])
se <- as.vector(bZse[i, ])
G <- t(Rxy * se) * se
var_residual1[i] <- sum(vartheta1 * c(G %*% vartheta1)) + sum(beta * c(Covtheta1 %*% beta))
var_residual2[i] <- sum(vartheta2 * c(G %*% vartheta2)) + sum(beta * c(Covtheta2 %*% beta))
}
close(pb)

cat("Classifying the whole genome into two mixtures:\n")
dt <- data.table(
idx = seq_len(length(residual1)),
residual1 = residual1,
residual2 = residual2,
var_residual1 = var_residual1,
var_residual2 = var_residual2,
LD.blcok = LD.blcok
)
dt_block <- dt[, .(
count1 = sum(abs(residual1) < abs(residual2)),
count2 = sum(abs(residual1) >= abs(residual2))
), by = LD.blcok]
dt_block[, cluster := ifelse(count1 >= count2, 1, 2)]
dt <- merge(dt, dt_block[, .(LD.blcok, cluster)], by = "LD.blcok", all.x = TRUE)
dt[cluster == 1, `:=`(
residual = residual1,
var_residual = var_residual1
)]
dt[cluster == 2, `:=`(
residual = residual2,
var_residual = var_residual2
)]
setorder(dt, idx)

A=data.frame(SNP=rownames(bX),BETA=dt$residual,SE=sqrt(dt$var_residual),Cluster=dt$cluster)
return(A)
}
