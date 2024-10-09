MRBEE_IPOD_UV=function(by,bX,byse,bXse,LD=LD,Rxy,cluster.index,tauvec=seq(3,50,by=2),max.iter=100,max.eps=0.001,ebic.gamma=1,rho=2,maxdiff=1.5,sampling.time=100,sampling.iter=5,theta.ini=F,gamma.ini=F,reliability.thres=0.8){
########################### Basic information #######################
by=by/byse
byseinv=1/byse
bX=bX*byseinv
bXse=bXse*byseinv
byse1=byse
byse=byse/byse
m=length(bX)
LD=Matrix(LD,sparse=T)
Theta=solve(LD)
TC=chol(Theta)
RC=as.matrix(TC%*%LD)
byinv=as.vector(Theta%*%by)
bXinv=as.vector(Theta%*%bX)
tilde.y=as.vector(TC%*%by)
tilde.X=as.vector(TC%*%bX)
Bt=t(bXinv)
BtB=sum(Bt*bX)
Thetarho=solve(LD+rho*diag(m))
r=reliability.adj.uv(bX,bXse,Theta=Theta,thres=reliability.thres)
r=c(r,1)
Rxy=t(t(Rxy)*r)*r
############################ Initial Estimate #######################
if(theta.ini[1]==F){
fit0=varbvs(X=RC,Z=tilde.X,y=tilde.y,verbose=F,maxiter=100)
gamma.ini=fit0$beta*(fit0$pip>0.5)
theta.ini=fit0$beta.cov[-1]
}
############################## Tuning Parameter ######################
w=length(tauvec)
Btheta=c(1:w)
Bgamma=matrix(0,m,w)
Bbic=tauvec
for(j in length(tauvec):1){
error=1
iter=1
theta=theta.ini
gamma=gamma.ini
gamma1=gamma
delta=gamma1
while(error>max.eps&iter<max.iter){
theta1=theta
indvalid=which(gamma1==0)
Hinv=1/(BtB-sum(bXse[indvalid]^2)*Rxy[1,1])
g=sum(Bt*(by-as.vector(LD%*%gamma)))-sum(bXse[indvalid])*Rxy[2,1]
theta=g*Hinv
########################### update gamma ############################
gamma=as.vector(Thetarho%*%c(by-bX*theta-delta+rho*gamma1))
gamma1=mcp(gamma+delta/rho,tauvec[j]/rho)
delta=delta+rho*(gamma-gamma1)

iter=iter+1
if(iter>3){
error=abs(theta-theta1)
}
}
Btheta[j]=theta
Bgamma[,j]=gamma1
df1=sum(gamma1!=0)
res=by-bX*theta-as.vector(LD%*%gamma)
rss=sum(res*(Theta%*%res))/(m-df1-1)
Bbic[j]=log(rss)*m+(log(m)+ebic.gamma*log(m))*df1+log(m)
}
Bbic=Bbic/m
######################## Inference #################################
jstar=which.min(Bbic)
theta=Btheta[jstar]
gamma=Bgamma[,jstar]
error=1
iter=1
theta=theta
gamma=gamma
names(gamma)=rownames(bX)
indgamma=which(gamma!=0)
indvalid=which(gamma==0)
res=by-bX*theta-as.vector(LD%*%gamma)
ThetaList=c(1:sampling.time)
for(j in 1:sampling.time){
cluster.sampling <- sample(1:max(cluster.index), 0.5*max(cluster.index), replace = F)
indj=which(cluster.index%in%cluster.sampling)
indj=sort(indj)
LDj=LD[indj,indj]
Thetaj=solve(LDj)
Thetarhoj=solve(LDj+rho*diag(length(indj)))
BtB=sum(bX[indj]*(Thetaj%*%bX[indj]))
gammaj=gamma1j=gamma
deltaj=0*gamma1j
for(iterj in 1:sampling.iter){
indvalidj=which(gamma1j==0)
indvalidj=intersect(indvalidj,indj)
Hinv=1/(BtB-sum(bXse[indvalidj]^2)*Rxy[1,1])
g=sum(bX[indj]*(by[indj]-as.vector(LD[indj,]%*%gammaj)))-sum(bXse[indvalidj])*Rxy[2,1]
thetaj=g*Hinv
resgammaj=as.vector(Thetarhoj%*%(by[indj]-bX[indj]*thetaj-deltaj[indj]+rho*gamma1j[indj]))
gammaj[indj]=resgammaj
gamma1j=mcp(gammaj+deltaj/rho,tauvec[jstar]/rho)
deltaj=deltaj+rho*(gammaj-gamma1j)
}
ThetaList[j]=thetaj
}
theta.se=sd(ThetaList)*sqrt((m-length(theta))/(m-length(theta)-length(indgamma)))

A=list()
A$theta=theta
A$gamma=gamma*byse1
A$theta.se=theta.se
A$Bic=Bbic
A$theta.ini=theta.ini
A$gamma.ini=gamma.ini
A$reliability.adjust=r
A$theta.pratt=getPratt.uv(bX=bX,by=by,bXse=bXse,byse=byse,Theta=Theta,theta=theta,Rxy=Rxy)
A$gamma.pratt=pleiotropyPratt(by=by,pleiotropy=gamma,Theta=Theta,LD=LD)
return(A)
}
