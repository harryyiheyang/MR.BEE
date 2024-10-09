#' Generating simulated data for Mendelian randomization simulation
#'
#' This function generates simulated data for Mendelian Randomization (MR) analysis,
#' considering genetic effects, estimation errors, and horizontal pleiotropy.
#' It allows for different distributions of genetic effects and pleiotropy,
#' and accommodates both independent and correlated instrumental variables (IVs).
#'
#' @param theta An (px1) vector of causal effects.
#' @param m The number of instrumental variables (IVs).
#' @param Rbb An (pxp) correlation matrix of genetic effects.
#' @param Ruv An ((p+1)x(p+1)) correlation matrix of residuals in outcome and exposures; the outcome is the first one.
#' @param Rnn An ((p+1)x(p+1)) correlation matrix of sample overlap; the outcome is the first one.
#' @param Nxy An ((p+1)x1) vector of GWAS sample sizes; the outcome is the first one.
#' @param Hxy An ((p+1)x1) vector of heritabilities; the outcome is the first one.
#' @param LD An (mxm) correlation matrix of the IVs or "identity" indicating independent IVs.
#' @param non.zero.frac An (px1) vector with all entries in (0,1]; each entry is the probability of deltaj such that betaj=betaj'*deltaj.
#' @param UHP.frac A number indicating the fraction of IVs affected by UHP.
#' @param UHP.var A number indicating the variance attributed to UHP.
#' @param UHP.dis Distribution of pleiotropy effects: "normal" (default), "uniform",  "t" distribution (with degree of freedom 5).
#' @param CHP.frac A number indicating the fraction of IVs affected by CHP.
#' @param CHP.effect A vector of effects corresponding to the variables correlated with the correlated horizontal pleiotropy.
#' @param effect.dis Distribution of genetic effects: "normal" (default), "uniform",  "t" distribution (with degree of freedom 5).
#'
#' @return A list containing simulated GWAS effect sizes for exposures (bX), their standard errors (bXse),
#'         the GWAS effect size for the outcome (by), its standard error (byse), the pleiotropy effects (pleiotropy), and the true effects.
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
#' @importFrom stats rt sd var runif rnorm
#'
#' @export

summary_generation=function(theta,m,Rbb,Ruv,Rnn,Nxy,Hxy,LD="identity",non.zero.frac,UHP.frac=0,CHP.frac=0,UHP.var=0.5,UHP.dis="uniform",CHP.effect=c(1,rep(0,length(theta)-1)),effect.dis="normal"){

if(CHP.frac==0){
p=length(theta)
hy=Hxy[p+1];hx=Hxy[1:p];Ny=Nxy[p+1];Nx=Nxy[1:p];

################ generate genetic effect ##########################
if(effect.dis=="normal"){
b0=MASS::mvrnorm(m,rep(0,p),Rbb)
b0=b0/sqrt(m);b0=t(t(b0)*sqrt(hx/non.zero.frac))
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
### Ensure that there must be 0.9*non.zero.frac nonzero elements
b0=b00*b0
b0=subjecth2(b0,hx)
}
if(effect.dis=="uniform"){
b0=MASS::mvrnorm(m,rep(0,p),Rbb)
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
### Ensure that there must be 0.9*non.zero.frac nonzero elements
b0=b00*b0
b0=subjecth2(b0,hx)
}
if(effect.dis=="t"){
b0=mvtnorm::rmvt(m,Rbb,df=5)
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
### Ensure that there must be 0.9*non.zero.frac nonzero elements
b0=b00*b0
b0=subjecth2(b0,hx)
}

################### Generate Estimation Error #######################
Syy=verrorvar(Rbb=Rbb,Suu=Ruv[-1,-1],Suv=Ruv[1,-1],hx=diag(hx),theta=theta,hy=hy,pleiotropy.var=UHP.var)
Duv=c(rep(1,p),sqrt(Syy))
Sigmauv=diag(Duv)%*%Ruv%*%diag(Duv)
Vxy=Sigmauv*Rnn

E=MASS::mvrnorm(m,rep(0,p+1),Vxy)
for(i in 1:(p+1)){
e=E[,i]
e=(e-mean(e))/sd(e)*sqrt(Vxy[i,i])
E[,i]=e
}
E=E%*%diag(1/sqrt(Nxy))

if(LD[1]=="identity"){
E1=E
bX=b0+E1[,-(p+1)]
by=b0%*%theta+E1[,p+1]
bX=as.matrix(bX)
by=as.vector(by)
by0=c(b0%*%theta)
}

if(LD[1]!="identity"){
C=matrixsqrt(LD)$w
E1=C%*%E
bX=LD%*%b0+E1[,-(p+1)]
by=LD%*%b0%*%theta+E1[,p+1]
bX=as.matrix(bX)
by=as.vector(by)
by0=c(b0%*%theta)
}

############################    generate pleiotropy    #####################
s=0*by;ra=1;pleiotropy0=0*by
if(UHP.frac!=0){
indpleio=sample(1:m,UHP.frac*m,replace=F)
bu=b0[,1]*0
if(length(indpleio)==0){
indpleio=round(m/2)
}
if(UHP.dis=="normal"){
bu1=rnorm(length(indpleio),0,1)
bu[indpleio]=bu1
}
if(UHP.dis=="uniform"){
bu1=runif(length(indpleio),0,1)
bu[indpleio]=bu1
}
if(UHP.dis=="t"){
bu1=stats::rt(length(indpleio),5,ncp=0)
bu[indpleio]=bu1
}

s=as.vector(bu)
s1=as.vector(b0%*%theta)
ra=var(s)/var(s1)/UHP.var
if(LD[1]=="identity"){
by=by+s/sqrt(ra)
pleiotropy0=s/sqrt(ra)
by0=s1+s/sqrt(ra)
}
if(LD[1]!="identity"){
by=by+LD%*%s/sqrt(ra)
pleiotropy0=s/sqrt(ra)
by0=s1+s/sqrt(ra)
}
}
by=c(by)
byse=rep(sd(E1[,p+1]),m)
bXse=matrix(1,m,p)
for(jj in 1:p){bXse[,jj]=sd(E1[,jj])}
}


if(CHP.frac>0){
m1=round(m*CHP.frac)
m=m-m1
################ generate the first cluster ##########################
p=length(theta)
hy=Hxy[p+1];hx=Hxy[1:p];Ny=Nxy[p+1];Nx=Nxy[1:p];
if(effect.dis=="normal"){
b0=MASS::mvrnorm(m,rep(0,p),Rbb)
b0=b0/sqrt(m);b0=t(t(b0)*sqrt(hx/non.zero.frac))
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
b0=b00*b0
b0=subjecth2(b0,hx*(1-CHP.frac))
}
if(effect.dis=="uniform"){
b0=MASS::mvrnorm(m,rep(0,p),Rbb)
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
b0=b00*b0
b0=subjecth2(b0,hx*(1-CHP.frac))
}
if(effect.dis=="t"){
b0=mvtnorm::rmvt(m,Rbb,df=5)
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
b0=b00*b0
b0=subjecth2(b0,hx*(1-CHP.frac))
}
Syy=verrorvar(Rbb=Rbb,Suu=Ruv[-1,-1],Suv=Ruv[1,-1],hx=diag(hx),theta=theta,hy=hy,pleiotropy.var=UHP.var*(1-CHP.frac))
Duv=c(rep(1,p),sqrt(Syy))
Sigmauv=diag(Duv)%*%Ruv%*%diag(Duv)
Vxy=Sigmauv*Rnn
E=MASS::mvrnorm(m,rep(0,p+1),Vxy)
for(i in 1:(p+1)){
e=E[,i]
e=(e-mean(e))/sd(e)*sqrt(Vxy[i,i])
E[,i]=e
}
E=E%*%diag(1/sqrt(Nxy))
pleiotropy0=c(1:m)*0
if(UHP.frac!=0){
indpleio=sample(1:m,UHP.frac*m,replace=F)
if(length(indpleio)==0){
indpleio=round(m/2)
}
bu=b0[,1]*0
if(UHP.dis=="normal"){
bu1=rnorm(length(indpleio),0,1)
bu[indpleio]=bu1
}
if(UHP.dis=="uniform"){
bu1=runif(length(indpleio),0,1)
bu[indpleio]=bu1
}
if(UHP.dis=="t"){
bu1=stats::rt(length(indpleio),5,ncp=0)
bu[indpleio]=bu1
}
s=as.vector(bu)
s1=as.vector(b0%*%theta)
ra=var(s)/var(s1)/UHP.var
pleiotropy0=s/sqrt(ra)
}
EMixture1=E
pleiotropy0Mixture1=pleiotropy0
b0Mixture1=b0
eta0Mixture1=c(b0%*%theta)
################ generate the second cluster ##########################
theta=CHP.effect
m=m1
hy=Hxy[p+1];hx=Hxy[1:p];Ny=Nxy[p+1];Nx=Nxy[1:p];
if(effect.dis=="normal"){
b0=MASS::mvrnorm(m,rep(0,p),Rbb)
b0=b0/sqrt(m);b0=t(t(b0)*sqrt(hx/non.zero.frac))
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
b0=b00*b0
b0=subjecth2(b0,hx*CHP.frac)
}
if(effect.dis=="uniform"){
b0=MASS::mvrnorm(m,rep(0,p),Rbb)
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
b0=b00*b0
b0=subjecth2(b0,hx*CHP.frac)
}
if(effect.dis=="t"){
b0=mvtnorm::rmvt(m,Rbb,df=5)
cutoff=qnorm(1-non.zero.frac)
b00=MASS::mvrnorm(m,rep(0,p),Rbb)
for(i in 1:p){
b00[,i]=b00[,i]>cutoff[i]
}
b00=checkzero(b00,nonzero.frac=non.zero.frac*0.9)
b0=b00*b0
b0=subjecth2(b0,hx*CHP.frac)
}
Syy=verrorvar(Rbb=Rbb,Suu=Ruv[-1,-1],Suv=Ruv[1,-1],hx=diag(hx),theta=theta,hy=hy,pleiotropy.var=UHP.var*CHP.frac)
Duv=c(rep(1,p),sqrt(Syy))
Sigmauv=diag(Duv)%*%Ruv%*%diag(Duv)
Vxy=Sigmauv*Rnn
E=MASS::mvrnorm(m,rep(0,p+1),Vxy)
for(i in 1:(p+1)){
e=E[,i]
e=(e-mean(e))/sd(e)*sqrt(Vxy[i,i])
E[,i]=e
}
E=E%*%diag(1/sqrt(Nxy))
pleiotropy0=c(1:m)*0
if(UHP.frac!=0){
indpleio=sample(1:m,UHP.frac*m,replace=F)
if(length(indpleio)==0){
indpleio=round(m/2)
}
bu=b0[,1]*0
if(UHP.dis=="normal"){
bu1=rnorm(length(indpleio),0,1)
bu[indpleio]=bu1
}
if(UHP.dis=="uniform"){
bu1=runif(length(indpleio),0,1)
bu[indpleio]=bu1
}
if(UHP.dis=="t"){
bu1=stats::rt(length(indpleio),5,ncp=0)
bu[indpleio]=bu1
}
s=as.vector(bu)
s1=as.vector(b0%*%theta)
ra=var(s)/var(s1)/UHP.var
pleiotropy0=s/sqrt(ra)
}
EMixture2=E
pleiotropy0Mixture2=pleiotropy0
b0Mixture2=b0
eta0Mixture2=c(b0%*%theta)

E=rbind(EMixture1,EMixture2)
pleiotropy0=c(pleiotropy0Mixture1,pleiotropy0Mixture2)
b0=rbind(b0Mixture1,b0Mixture2)
eta0=c(eta0Mixture1,eta0Mixture2)

if(LD[1]=="identity"){
E1=E
bX=b0+E[,1:p]
by=eta0+E[,p+1]+pleiotropy0
}else{
RC=matrixsqrt(LD)$w
E1=RC%*%E
bX=LD%*%b0+E1[,1:p]
by=LD%*%(eta0+pleiotropy0)+E1[,p+1]
}
by0=eta0
by=c(by)
byse=rep(sd(E1[,p+1]),length(by))
bXse=matrix(1,length(by),p)
for(jj in 1:p){bXse[,jj]=sd(E1[,jj])}
}

A=list(bX=bX,bXse=bXse,by=by,byse=byse,pleiotropy=pleiotropy0,bX0=b0,by0=by0)
return(A)
}
