#' @title Univariable Mendelian Randomization using Bias-corrected Estimating Equations
#' @description
#' This function removes weak instrument bias using Bias-corrected Estimating Equations and identifies uncorrelated horizontal pleiotropy (UHP) and correlated horizontal pleiotropy (CHP) through two distinct methods. UHP is detected using the IPOD algorithm, where outliers are interpreted as UHP. For CHP, a two-mixture regression model is applied. Both the IPOD algorithm and the Mixture method support the inclusion of correlated instrumental variables using an LD matrix and provide advanced options for exposure selection and horizontal pleiotropy correction.
#'
#' @param by A vector of effect estimates from the outcome GWAS.
#' @param bX A vector of effect estimates from the exposure GWAS.
#' @param byse A vector of standard errors of effect estimates from the outcome GWAS.
#' @param bXse A vector of standard errors of effect estimates from the exposure GWAS.
#' @param LD The linkage disequilibrium (LD) matrix. Default is the identity matrix, assuming independent instrumental variables (IVs).
#' @param Rxy The correlation matrix of estimation errors of exposures and outcome GWAS. The last element corresponds to the outcome.
#' @param cluster.index A vector indicating the LD block indices each IV belongs to. The length is equal to the number of IVs, and values are the LD block indices.
#' @param reliability.thres A threshold for the minimum value of the reliability ratio. If the original reliability ratio is less than this threshold, only part of the estimation error is removed so that the working reliability ratio equals this threshold.
#' @param Method Method for handling horizontal pleiotropy. Options are \code{"IPOD"} and \code{"Mixture"}.
#' @param tauvec When choosing \code{"IPOD"}, the candidate vector of tuning parameters for the MCP penalty function. Default is \code{seq(3, 30, by=3)}.
#' @param rho When choosing \code{"IPOD"}, the tuning parameter in the nested ADMM algorithm. Default is \code{2}.
#' @param main.cluster.thres When choosing \code{"Mixture"}, a threshold for weights belonging to the first category. To prevent instability caused by small-effect IVs falling into both categories, we slightly lower the voting threshold for the first category to below 0.5, ensuring it remains dominant. Default is \code{0.4}.
#' @param min.cluster.size When choosing \code{"Mixture"}, a threshold for the minimum number of IVs in the second cluster. If our initial check reveals that the number is below this threshold, the IPOD algorithm will be applied. Default is \code{10}.
#' @param robust.se When choosing \code{"Mixture"}, an indicator of whether the robust covariance estimate is applied to calculate the empirical covariance matrix of causal effect estimates from subsampling results. Default is \code{T}.
#' @param max.iter Maximum number of iterations for causal effect estimation. Defaults to \code{100}.
#' @param max.eps Tolerance for stopping criteria. Defaults to \code{0.001}.
#' @param maxdiff The maximum difference between the MRBEE causal estimate and the initial estimator. Defaults to \code{1.5}.
#' @param ebic.theta EBIC factor on causal effect. Default is \code{0}.
#' @param ebic.gamma EBIC factor on horizontal pleiotropy Default is \code{2}.
#' @param sampling.time Number of resampling times. Default is \code{100}.
#' @param sampling.iter Number of iterations per resampling. Default is \code{5}.
#' @param theta.ini Initial value of theta. If \code{FALSE}, the default method is used to estimate. Default is \code{FALSE}.
#' @param gamma.ini Initial value of gamma. Default is \code{FALSE}.
#'
#' @importFrom MASS rlm
#' @importFrom CppMatrix matrixInverse matrixMultiply matrixVectorMultiply matrixEigen matrixListProduct
#' @importFrom Matrix Matrix solve chol bdiag
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom mixtools regmixEM
#' @importFrom FDRestimation p.fdr
#'
#' @return A list containing the results of the MRBEEX.UV analysis:
#' \describe{
#'   \item{\code{theta}}{Causal effect estimate.}
#'   \item{\code{theta.se}}{Standard error of the causal effect estimate.}
#'   \item{\code{theta.cov}}{Covariance matrix of the causal effect estimate.}
#'   \item{\code{theta.pratt}}{Pratt index estimate of exposure.}
#'   \item{\code{gamma}}{Estimate of horizontal pleiotropy.}
#'   \item{\code{gamma.pratt}}{Pratt index estimate of horizontal pleiotropy.}
#'   \item{\code{Bic}}{A vector or matrix recording the Bayesian Information Criterion (BIC) values.}
#'   \item{\code{theta.ini}}{Initial value of \code{theta} used in the estimation procedure.}
#'   \item{\code{gamma.ini}}{Initial value of \code{gamma} used in the estimation procedure.}
#'   \item{\code{reliability.adjust}}{Estimated reliability-adjusted values.}
#'   \item{\code{thetalist}}{List of \code{theta} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{gammalist}}{List of \code{gamma} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{theta1}}{Causal effect estimate for the first mixture component (when \code{Method="Mixture"}).}
#'   \item{\code{theta2}}{Causal effect estimate for the second mixture component (when \code{Method="Mixture"}).}
#'   \item{\code{theta.se1}}{Standard error of \code{theta1}.}
#'   \item{\code{theta.se2}}{Standard error of \code{theta2}.}
#'   \item{\code{theta.cov1}}{Covariance matrix of \code{theta1}.}
#'   \item{\code{theta.cov2}}{Covariance matrix of \code{theta2}.}
#'   \item{\code{theta.pratt1}}{Pratt index estimates of exposures in the first mixture.}
#'   \item{\code{theta.pratt2}}{Pratt index estimates of exposures in the second mixture.}
#'   \item{\code{thetalist1}}{List of \code{theta1} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{thetalist2}}{List of \code{theta2} estimates recorded during each iteration in the subsampling procedure.}
#'   \item{\code{cluster1}}{Indices of individual IVs in the first mixture component.}
#'   \item{\code{cluster2}}{Indices of individual IVs in the second mixture component.}
#'   \item{\code{fit.mixture}}{Return of fit.mixture in the lastest iteration.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage of MRBEEX.UV function
#' result <- MRBEEX_UV(by, bX, byse, bXse)
#' }
MRBEEX_UV=function(by,bX,byse,bXse,LD="identity",Rxy,cluster.index=c(1:length(by)),
        reliability.thres=0.8,
        Method=c("IPOD","Mixture"),ebic.theta=0,
        tauvec=seq(3,30,by=3),rho=2,ebic.gamma=2,
        main.cluster.thres=0.48,min.cluster.size=10,robust.se=T,
        max.iter=100,max.eps=0.001,maxdiff=3,
        sampling.time=100,sampling.iter=10,
        theta.ini=F,gamma.ini=F){
if(LD[1]!="identity"){
##########################################################################
if(Method[1]=="IPOD"){
A=MRBEE_IPOD_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,cluster.index=cluster.index,ebic.gamma=ebic.gamma,max.iter=max.iter,max.eps=max.eps,tauvec=tauvec,reliability.thres=reliability.thres,rho=rho,maxdiff=maxdiff,sampling.time=sampling.time,sampling.iter=sampling.iter,theta.ini=theta.ini,gamma.ini=gamma.ini)
}
##########################################################################
if(Method[1]=="Mixture"){
A=MRBEE_Mixture_UV(by=by,bX=bX,byse=byse,bXse=bXse,LD=LD,Rxy=Rxy,main.cluster.thres=main.cluster.thres,cluster.index=cluster.index,reliability.thres=reliability.thres,sampling.time=sampling.time,min.cluster.size=min.cluster.size,robust.se=robust.se,ebic.theta=ebic.theta)
if(A$IsIPOD==T){
A=list()
A$IsMixture="Initial check suggests only one pathway"
sink()
print(A$IsMixture)
}
}
}

if(LD[1]=="identity"){
if(Method[1]=="IPOD"){
A=MRBEE_IPOD_UV_Independent(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,rho=rho,tauvec=tauvec,max.iter=max.iter,max.eps=max.eps,ebic.gamma=ebic.gamma,maxdiff=maxdiff,sampling.time=sampling.time,theta.ini=theta.ini,gamma.ini=gamma.ini,sampling.iter=sampling.time,reliability.thres=reliability.thres)
}
##########################################################################
if(Method[1]=="Mixture"){
A=MRBEE_Mixture_Independent_UV(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,reliability.thres=reliability.thres,sampling.time=sampling.time,main.cluster.thres=main.cluster.thres,min.cluster.size=min.cluster.size,robust.se=robust.se,ebic.theta=ebic.theta)
if(A$IsIPOD==T){
A=list()
A$IsMixture="Initial check suggests only one pathway"
sink()
print(A$IsMixture)
}
}
##########################################################################
}
return(A)
}


