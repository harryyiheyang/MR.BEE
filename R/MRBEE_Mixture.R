MRBEE_Mixture=function(by,bX,byse,bXse,LD,Rxy,cluster.index=c(1:length(by)),main.cluster.thres=0.48,min.cluster.size=5,reliability.thres=0.8,sampling.time=100,robust.se=T,ebic.theta=1){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=nrow(bX)
p=ncol(bX)
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
tilde.y=as.vector(TC%*%by)
tilde.X=as.matrix(TC%*%bX)
tilde.R=as.matrix(TC%*%LD)
r=reliability.adj(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
############################ Initial Estimate #######################
fit.init=fit.mixture=regmixEM(y=tilde.y,x=tilde.X,k=2,addintercept=F,maxit=m)
max.cluster=ifelse(sum(fit.init$posterior[,1]>main.cluster.thres)>(m/2),1,2)
cluster.1=which(fit.init$posterior[,max.cluster]>main.cluster.thres)
cluster.2=setdiff(1:m,cluster.1)
########################### Check if applied mixture model #######################
if(length(cluster.2)>=min.cluster.size){
cluster.ratio=c(length(cluster.1),length(cluster.2))/m
indvalid1=cluster.1
indvalid2=cluster.2
Rxysum1=biasterm(RxyList=RxyList,indvalid1)
XtX1=matrixMultiply(t(tilde.X[cluster.1,]),tilde.X[cluster.1,])
Xty1=matrixVectorMultiply(t(tilde.X[cluster.1,]),tilde.y[cluster.1])
theta1=c(solve(XtX1-Rxysum1[1:p,1:p])%*%(Xty1-Rxysum1[1:p,p+1]))
Rxysum2=biasterm(RxyList=RxyList,indvalid2)
XtX2=matrixMultiply(t(tilde.X[cluster.2,]),tilde.X[cluster.2,])
Xty2=matrixVectorMultiply(t(tilde.X[cluster.2,]),tilde.y[cluster.2])
theta2=c(solve(XtX2-Rxysum2[1:p,1:p])%*%(Xty2-Rxysum2[1:p,p+1]))
Bic=MRFit(by=tilde.y,bX=tilde.X,theta1=theta1,theta2=theta2,cluster.1=cluster.1,cluster.2=cluster.2,df1=p,df2=p)+2*p*(log(p*2)*ebic.theta+log(m))
Bic=Bic/m
############################### inference #########################
names(theta1)=names(theta2)=colnames(bX)
ThetaList1=ThetaList2=matrix(0,sampling.time,p)
colnames(ThetaList1)=colnames(ThetaList2)=colnames(bX)
cat("Bootstrapping process ...\n")
sink(tempfile())
for(j in 1:sampling.time) {
indicator <- FALSE
tryCatch({
cluster.sampling <- sample(1:max(cluster.index), 0.5*max(cluster.index), replace = F)
indj=which(cluster.index%in%cluster.sampling)
indj=sort(indj)
LDj=LD[indj,indj]
Thetaj <- solve(LDj)
TCj=chol(Thetaj)
tilde.yj=as.vector(TCj%*%by[indj])
tilde.Xj=as.matrix(TCj%*%bX[indj,])
tilde.Rj=as.matrix(TCj%*%LDj)
theta1j=theta1*runif(1,0.95,1.05)*0.95
theta2j=theta2*runif(1,0.95,1.05)*0.95
fit.mixturej=regmixEM(y=tilde.yj,x=tilde.Xj,k=2,addintercept=F,lambda=cluster.ratio,beta=cbind(theta1,theta2)*0.9,maxit=100,epsilon=1e-3)
max.cluster=ifelse(sum(fit.mixturej$posterior[,1]>main.cluster.thres)>(length(indj)/2),1,2)
cluster.1j=which(fit.mixturej$posterior[,max.cluster]>main.cluster.thres)
cluster.2j=which(fit.mixturej$posterior[,max.cluster]<=main.cluster.thres)
Rxysum1j=biasterm(RxyList=RxyList,indj[cluster.1j])
XtX1j=matrixMultiply(t(tilde.Xj[cluster.1j,]),tilde.Xj[cluster.1j,])
Xty1j=matrixVectorMultiply(t(tilde.Xj[cluster.1j,]),tilde.yj[cluster.1j])
theta1j=c(solve(XtX1j-Rxysum1j[1:p,1:p])%*%(Xty1j-Rxysum1j[1:p,p+1]))
Rxysum2j=biasterm(RxyList=RxyList,indj[cluster.2j])
XtX2j=matrixMultiply(t(tilde.Xj[cluster.2j,]),tilde.Xj[cluster.2j,])
Xty2j=matrixVectorMultiply(t(tilde.Xj[cluster.2j,]),tilde.yj[cluster.2j])
theta2j=c(solve(XtX2j-Rxysum2j[1:p,1:p])%*%(Xty2j-Rxysum2j[1:p,p+1]))
ThetaVecj=cbind(theta1j,theta2j)
ThetaNormj=colMeans((ThetaVecj-cbind(theta1,theta1))^2)
theta1j=ThetaVecj[,which.min(ThetaNormj)]
theta2j=ThetaVecj[,which.max(ThetaNormj)]

ThetaList1[j, ] <- theta1j
ThetaList2[j, ] <- theta2j
}, error = function(e) {
# Error handling block
cat("Error occurred: ", e$message, "\n")
indicator <<- TRUE  # Set indicator to TRUE if an error occurs
j <<- j - 1  # Decrement the iteration counter to retry
})
if (indicator) {
next  # Retry the current iteration
}
}
indtheta1=which(theta1!=0)
indtheta2=which(theta2!=0)
theta.se1=colSD(ThetaList1)
theta.cov1=cov(ThetaList1)
theta.se2=colSD(ThetaList2)
theta.cov2=cov(ThetaList2)
if(robust.se==T){
theta.se1=colSDMAD(ThetaList1)
theta.cov1=covmad(ThetaList1)
theta.se2=colSDMAD(ThetaList2)
theta.cov2=covmad(ThetaList2)
}
colnames(theta.cov2)=rownames(theta.cov2)=names(theta.se2)=colnames(bX)
colnames(theta.cov1)=rownames(theta.cov1)=names(theta.se1)=colnames(bX)
A=list()
A$theta1=theta1
A$theta2=theta2
A$theta.se1=theta.se1
A$theta.cov1=theta.cov1
A$theta.se2=theta.se2
A$theta.cov2=theta.cov2
A$reliability.adjust=r
A$thetalist1=ThetaList1
A$thetalist2=ThetaList2
A$cluster1=cluster.1
A$cluster2=cluster.2
A$fit.mixture=fit.mixture
A$theta.pratt1=getPratt(bX=bX[cluster.1,],by=by[cluster.1],bXse=bXse[cluster.1,],byse=byse[cluster.1],Theta=solve(LD[cluster.1,cluster.1]),theta=theta1,Rxy=Rxy)
A$theta.pratt2=getPratt(bX=bX[cluster.2,],by=by[cluster.2],bXse=bXse[cluster.2,],byse=byse[cluster.2],Theta=solve(LD[cluster.2,cluster.2]),theta=theta2,Rxy=Rxy)
A$IsIPOD=F
A$Bic=Bic
sink()
return(A)
}else{
A=list()
A$IsIPOD=T
sink()
return(A)
}
}
