MRBEE_IMRP_UV=function(by,bX,byse,bXse,Rxy,max.iter=30,max.eps=1e-4,pv.thres=0.05,var.est="variance",FDR=T,adjust.method="Sidak"){
by=by/byse
byseinv=1/byse
bx=bX*byseinv
bxse=bXse*byseinv
byse1=byse
byse=byse/byse
n=length(by)
RxyList=IVweight(byse,bxse,Rxy)
########## Initial Estimation ############
fit=MASS::rlm(by~bx-1)
theta=fit$coefficient
theta1=10000
e=c(by-bx*theta)
indvalid=which(abs(e)<=3*stats::mad(e))
indvalid=validadj(abs(e),indvalid,0.5)
########## Iteration ###################
error=1
iter=0
while(error>max.eps&iter<max.iter){
theta1=theta
e=c(by-bx*theta)
pv=imrpdetect(x=e,theta=theta,RxyList=RxyList,var.est=var.est,FDR=FDR,adjust.method=adjust.method,indvalid=indvalid)
indvalid=which(pv>pv.thres)
if (length(indvalid) < length(pv) * 0.5) {
indvalid.cut = which(pv > stats::quantile(pv, 0.5))
indvalid = union(indvalid, indvalid.cut)
}
h=sum(bx[indvalid]^2)-sum(bxse[indvalid]^2*Rxy[1,1])
g=sum(bx[indvalid]*by[indvalid])-Rxy[1,2]*sum(bxse[indvalid]*byse[indvalid])
theta=g/h
iter=iter+1
if(iter>5) error=sqrt(sum((theta-theta1)^2))
}
adjf=n/(length(indvalid)-1)
Hat=outer(bx[indvalid],bx[indvalid])/h
Hat=1-diag(Hat)
Hat[Hat<0.5]=0.5
e[indvalid]=e[indvalid]/Hat
E=-bx[indvalid]*e[indvalid]+bxse[indvalid]*byse[indvalid]-bxse[indvalid]^2*theta
vartheta=sum(E^2)/h^2*adjf
A=list()
A$theta=theta
A$vartheta=vartheta
r=c(by-bx*theta)*byse1
r[indvalid]=0
names(r)=rownames(bx)
A$gamma=r
A$theta.se=sqrt(vartheta)
A$theta.pratt=getPratt.uv(bX=bX,by=by,bXse=bXse,byse=byse,theta=theta,Rxy=Rxy)
A$gamma.pratt=pleiotropyPratt(by=by,pleiotropy=r)
return(A)
}
