MRBEE_IPOD_SuSiE_Independent=function(by,bX,byse,bXse,Rxy,Lvec=c(1:min(10,nrow(bX))),pip.thres=0.5,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,susie.iter=100,ebic.theta=0,ebic.gamma=1,reliability.thres=0.5,rho=2,maxdiff=1.5,sampling.time=100,sampling.iter=5,theta.ini=F,gamma.ini=F){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=nrow(bX)
p=ncol(bX)
cluster.index=c(1:m)
BtB=matrixMultiply(t(bX),bX)
BtB=t(BtB)/2+BtB/2
r=reliability.adj(bX,bXse,Theta="identity",thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
RxyList=IVweight(byse,bXse,Rxy)
Rxyall=biasterm(RxyList=RxyList,c(1:m))
############################ Initial Estimate #######################
if(theta.ini[1]==F){
fit0=MRBEE_IMRP(by=by,bX=bX,byse=byse,bXse=bXse,Rxy=Rxy,var.est="variance",FDR="Sidak")
gamma.ini=gamma.ini1=fit0$gamma
theta.ini=theta.ini1=fit0$theta
}else{
gamma.ini=gamma.ini1=gamma.ini/byse1
theta.ini=theta.ini1=theta.ini
}
############################## Tuning Parameter ######################
w=length(tauvec)
q=length(Lvec)
Btheta=array(0,c(p,w,q))
Bgamma=array(0,c(m,w,q))
Bbic=matrix(0,w,q)
for(v in length(Lvec):1){
fit.theta=NULL
for(j in length(tauvec):1){
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1*0
error=1
iter=1
empirical.variance=5
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-gamma
XtX=BtB
Xty=matrixVectorMultiply(t(bX),res.theta)
yty=sum(res.theta^2)
fit.theta=susie_suff_stat(XtX=BtB,Xty=Xty,yty=yty,n=m,L=Lvec[v],residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=T,max_iter=susie.iter,s_init=fit.theta)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.thres)
indtheta=which(theta!=0)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
if((norm(theta,"2")/norm(theta.ini1,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
}
gamma=(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1)/(1+rho)
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
Btheta[,j,v]=theta
Bgamma[,j,v]=gamma1
df1=sum(gamma1!=0)
df2=min(Lvec[v],sum(theta!=0))
res=c(by-matrixVectorMultiply(bX,theta)-gamma)
rss=sum(res^2)/(m-df1-df2)
Bbic[j,v]=log(rss)*m+df2*(log(m)+log(p)*ebic.theta)+(log(m)+ebic.gamma*log(m))*df1
}
}
Bbic=Bbic/m
jstar=bimin(Bbic)[1]
vstar=bimin(Bbic)[2]
theta.ini=Btheta[,jstar,vstar]
gamma.ini=Bgamma[,jstar,vstar]
indvalid=which(gamma.ini==0)
res=by-matrixVectorMultiply(bX,theta.ini)-gamma.ini
rss=sum(res^2)
######################## Final Estimate #################################
error=1
iter=1
gamma1=gamma
delta=gamma1
empirical.variance=5
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
if(length(indvalid)==m){
Rxysum=Rxyall
}else{
Rxysum=Rxyall-biasterm(RxyList=RxyList,setdiff(1:m,indvalid))
}
res.theta=by-gamma
XtX=BtB
Xty=matrixVectorMultiply(t(bX),res.theta)
yty=sum(res.theta^2)
fit.theta=susie_suff_stat(XtX=BtB,Xty=Xty,yty=yty,n=m,L=Lvec[vstar],residual_variance=empirical.variance,estimate_prior_method="EM",intercept=F,estimate_residual_variance=T,max_iter=susie.iter)
theta=coef.susie(fit.theta)[-1]*(fit.theta$pip>pip.thres)
indtheta=which(theta!=0)
if(length(indtheta)==1){
xtx=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=xty/xtx
}
if(length(indtheta)>1){
XtX=XtX[indtheta,indtheta]-Rxysum[indtheta,indtheta]
Xty=Xty[indtheta]-Rxysum[indtheta,p+1]
theta[indtheta]=c(solve(XtX)%*%Xty)
}
if((norm(theta,"2")/norm(theta.ini1,"2"))>maxdiff){
theta=theta/norm(theta,"2")*maxdiff*norm(theta.ini1,"2")
}
########################### update gamma ############################
gamma=(by-matrixVectorMultiply(bX,theta)-delta+rho*gamma1)/(1+rho)
gamma1=mcp(gamma+delta/rho,tauvec[jstar]/rho)
delta=delta+rho*(gamma-gamma1)
res=c(by-matrixVectorMultiply(bX,theta)-gamma1)
rss=sum(res^2)
iter=iter+1
if(iter>3){
error=max(abs(theta-theta1))
}
}
############################### inference #########################
gamma=gamma1
theta=theta
names(theta)=colnames(bX)
names(gamma)=rownames(bX)
indtheta=which(theta!=0)
indgamma=which(gamma1!=0)
indvalid=which(gamma1==0)
res=by-matrixVectorMultiply(bX,theta)-gamma1

ThetaList=matrix(0,sampling.time,p)
GammaList=matrix(0,sampling.time,m)
cat("Bootstrapping process:\n")
pb <- txtProgressBar(min = 0, max = sampling.time, style = 3)
for(j in 1:sampling.time) {
setTxtProgressBar(pb, j)
indj=sample(m,0.5*m,replace=F)
indj <- sort(indj)
BtBj <- matrixMultiply(t(bX[indj,]),bX[indj,])
BtBj=t(BtBj)/2+BtBj/2
thetaj=theta
gammaj=gamma1j=gamma
deltaj=gammaj*0
for(jiter in 1:sampling.iter){
indvalidj <- which(gamma1j==0)
indvalidj <- intersect(indvalidj, indj)
Rxysumj <- biasterm(RxyList = RxyList, indvalidj)
res.thetaj=by[indj]-gammaj[indj]
XtXj=BtBj
Xtyj=matrixVectorMultiply(t(bX[indj,]),res.thetaj)
ytyj=sum(res.thetaj^2)
fit.thetaj=susie_suff_stat(XtX=BtBj,Xty=Xtyj,yty=ytyj,n=length(indvalidj),L=Lvec[vstar],estimate_prior_method="EM",intercept=F,estimate_residual_variance=T,max_iter=sampling.iter,s_init=fit.theta)
thetaj=coef.susie(fit.thetaj)[-1]*(fit.thetaj$pip>pip.thres)
indthetaj=which(thetaj!=0)
if(length(indthetaj)==1){
xtxj=XtXj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,p+1]
thetaj[indthetaj]=xtyj/xtxj
}
if(length(indthetaj)>1){
XtXj=XtXj[indthetaj,indthetaj]-Rxysumj[indthetaj,indthetaj]
Xtyj=Xtyj[indthetaj]-Rxysumj[indthetaj,p+1]
thetaj[indthetaj]=c(solve(XtXj)%*%Xtyj)
}
if((norm(thetaj, "2") / norm(theta.ini1, "2")) > maxdiff) {
thetaj <- thetaj / norm(thetaj, "2") * maxdiff * norm(theta.ini1, "2")
}
gammaj[indj]=(by[indj]-matrixVectorMultiply(bX[indj, ],thetaj)-deltaj[indj]+rho*gamma1j[indj])/(1+rho)
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}
ThetaList[j, ] <- thetaj
GammaList[j, ] = gamma1j
}
close(pb)
theta.se=colSD(ThetaList)*sqrt((m-length(indtheta))/(m-length(indtheta)-length(indgamma)))
theta.cov=cov(ThetaList)*(m-length(indtheta))/(m-length(indtheta)-length(indgamma))
colnames(theta.cov)=rownames(theta.cov)=names(theta.se)=colnames(bX)

A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$theta.cov=theta.cov
A$Bic=Bbic
A$reliability.adjust=r
A$susie.theta=fit.theta
A$thetalist=ThetaList
A$gammalist=GammaList
A$theta.pip=colMeans(ThetaList!=0)
A$theta.pratt=getPratt(bX=bX,by=by,bXse=bXse,byse=byse,theta=theta,Rxy=Rxy)
A$gamma.pratt=pleiotropyPratt(by=by,pleiotropy=gamma)
return(A)
}
