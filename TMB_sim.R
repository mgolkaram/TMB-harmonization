require(BB)
require(MASS)
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1.300000
b<-35.600000
m<-a/b

tissue = 'LUSC'

d<-WES$WES[which(WES$Cancer.Type==tissue)]*b
res <- as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #Fn<-ecdf(d)
N<-2000
rand.samples.x = rep(NA,N)
rand.samples.y = rep(NA,N)
  #res$par[1]<-265.22 #LUSC
  #res$par[1]<-202.53 #LUAD
  #res$par[1]<-1.5336370/m #test
  #res$par[2]<-84.372 #LUSC
  #res$par[2]<-123.57 #LUAD
  #res$par[2]<-0.8825977/m #test
  
constant_noise= 0.05
S<-sample(d,replace = TRUE,size = N)
for(i in 1:N){
  POISSON_NOISE=rpois(1,1)  
  #POISSON_NOISE=0
  
  #rand.samples.x[i] = (rnorm(1,mean  = (res[1]),sd = (res[2])))
  #rand.samples.x[i] =   as.numeric(quantile(Fn, runif(1)))/m
  rand.samples.x[i]=S[i]
  rand.samples.y[i] = rand.samples.x[i]*m+rnorm(1,0,sd = constant_noise+sqrt((rand.samples.x[i])/b*a))+POISSON_NOISE
  rand.samples.x[i]=S[i]+POISSON_NOISE
  
}
#rand.samples.x[i]=rand.samples.x[i]+P_noise/m}
idx<-which((rand.samples.y>=0) & rand.samples.x>=0)
X<-as.data.frame(cbind(round(rand.samples.x[idx]),round(rand.samples.y[idx])))
  
## some pretty colors
library(RColorBrewer)
k <- 30
my.cols <- rev(brewer.pal(k, "RdYlBu"))

## compute 2D kernel density, see MASS book, pp. 130-131
z <- kde2d(X[,1]/b, X[,2]/a, n=50)

plot(WES$WES[which(WES$Cancer.Type==tissue)],WES$TST500[which(WES$Cancer.Type==tissue)], pch=19, cex=.4, xlim = c(0,25), ylim = c(0,25), xlab = 'WES', ylab = 'TSO500')
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
#abline(h=mean(X[,2]), v=mean(X[,1]), lwd=2)
#legend("topleft", paste("R=", round(cor(X)[1,2],2)), bty="n")
legend("topleft", c('TCGA','model'), lty = c(NA,1),pch = c(19, NA), bty = 'n')

qqplot(X[,2],WES$TST500[which(WES$Cancer.Type==tissue)]*a,xlab = "model", ylab = "TCGA")
lines(c(0,250),c(0,250))

qqplot(X[,1],WES$WES[which(WES$Cancer.Type==tissue)]*b,xlab = "model", ylab = "TCGA")
lines(c(0,1200),c(0,1200))
 
perc.exp<-as.vector(quantile(WES$TST500[which(WES$Cancer.Type==tissue)]*a,probs = seq(0.01,1,length.out = 20)))
f<-ecdf(X[,2])
plot(seq(1,100,length.out = 20),100*f(perc.exp), xlab = 'model', ylab = 'TCGA', pch=19,ylim = c(0,100))
lines(c(0,100),c(0,100))

##### tissue dependant slope panel vs panel
require(BB)
require(MASS)
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-2.00000
b<-35.600000
c<-1.0
m<-a/b
M<-c/b
TMB.sim <- function(tissue){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]*b
  #res <- as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #Fn<-ecdf(d)
  N<-500
  rand.samples.x = rep(NA,N)
  rand.samples.y = rep(NA,N)
  rand.samples.z = rep(NA,N)
  constant_noise= 0.5
  S=sample(d,size = N,replace = TRUE)
  for(i in 1:N){
    POISSON_NOISE=rpois(1,1)
    #rand.samples.x[i] = (rnorm(1,mean  = (res[1]),sd = (res[2])))
    #rand.samples.x[i] =   as.numeric(quantile(Fn, runif(1)))/m
    rand.samples.x[i] = S[i]
    rand.samples.y[i] = rand.samples.x[i]*m+rnorm(1,0,sd = constant_noise+sqrt((rand.samples.x[i])/b*a)) + POISSON_NOISE
    rand.samples.z[i] = rand.samples.x[i]*M+rnorm(1,0,sd = constant_noise+sqrt((rand.samples.x[i])/b*c)) + POISSON_NOISE
    rand.samples.x[i] = rand.samples.x[i] + POISSON_NOISE
  }
  #rand.samples.x[i]=rand.samples.x[i]+P_noise/m}
  idx<-which((rand.samples.y>=0) & rand.samples.z>=0)
  X<-as.data.frame(cbind(round(rand.samples.z[idx]),round(rand.samples.y[idx])))
  return(round(lm(X[,2]~X[,1])$coefficients[2],3))}



tissue.list<-as.character(unique(WES$Cancer.Type))
slopes<-sapply(X = tissue.list,function(t) replicate(50,TMB.sim(t)))
#colVar <- function(x)(apply(x, 2, FUN = var))
#boxplot(slopes[,order(colMeans(slopes))]/m,las = 2,col = 'red',ylab = 'slopes',ylim = c(0,0.04)/m)
boxplot(slopes[,order(colSums(slopes))]/(a/c),las = 2,col = 'red',ylab = 'slopes',ylim = c(0,1.2))




####
TMB.sim <- function(tissue){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]*b
  #res <- as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #Fn<-ecdf(d)
  N<-500
  rand.samples.x = rep(NA,N)
  rand.samples.y = rep(NA,N)
  constant_noise= 0.05
  S=sample(d,size = N,replace = TRUE)
  for(i in 1:N){
    POISSON_NOISE=rpois(1,1)
    #rand.samples.x[i] = (rnorm(1,mean  = (res[1]),sd = (res[2])))
    #rand.samples.x[i] =   as.numeric(quantile(Fn, runif(1)))/m
    rand.samples.x[i] = S[i]
    rand.samples.y[i] = rand.samples.x[i]*m+rnorm(1,0,sd = constant_noise+sqrt((rand.samples.x[i])*a/b)) + POISSON_NOISE
    rand.samples.x[i] = rand.samples.x[i] + POISSON_NOISE
  }
  #rand.samples.x[i]=rand.samples.x[i]+P_noise/m}
  idx<-which((rand.samples.y>=0) & rand.samples.x>=0)
  X<-as.data.frame(cbind(round(rand.samples.x[idx]),round(rand.samples.y[idx])))
  return(round(lm(X[,2]~X[,1])$coefficients[2],3))}

tissue.list<-as.character(unique(WES$Cancer.Type))
slopes<-sapply(X = tissue.list,function(t) replicate(50,TMB.sim(t)))
colVar <- function(x)(apply(x, 2, FUN = var))
#boxplot(slopes[,order(colMeans(slopes))]/m,las = 2,col = 'red',ylab = 'slopes',ylim = c(0,0.04)/m)
boxplot(slopes[,order(colVar(slopes))]/m,las = 2,col = 'red',ylab = 'slopes',ylim = c(0,2))

true.slopes<-sapply(tissue.list,FUN = function(x) {
  round(lm((WES$TST500[which(WES$Cancer.Type==x)])~(WES$WES[which(WES$Cancer.Type==x)]))$coefficients[2],3)})     

true.slopes.bs=list()
for(i in seq(500)){
  true.slopes.bs[[i]]<-sapply(tissue.list,FUN = function(x) {
    bs.size = 30
    idx=sample(1:length(which(WES$Cancer.Type==x)),size = bs.size)
    var.x<-(WES$WES[which(WES$Cancer.Type==x)])[idx]
    var.y<-(WES$TST500[which(WES$Cancer.Type==x)])[idx]
    round(lm(var.y~var.x)$coefficients[2],3)})     
}
true.slopes.bs=do.call(rbind, true.slopes.bs)

require(MALDIquant)
bs <- function(x,u,v){
  v<-v[order(u,decreasing = FALSE)]
  u<-u[order(u,decreasing = FALSE)]
  num = 20
  rand.index<-runif(num,min(u),max(u))
  u.sampled<-as.vector(scale(u[match.closest(rand.index,u)]))
  v.sampled<-as.vector(scale(v[match.closest(rand.index,u)]))
  return(round(lm(v.sampled ~ u.sampled)$coefficients[2],3))
}

bootstrap.slopes<-sapply(tissue.list,FUN = function(x) {
  U<-WES$WES[which(WES$Cancer.Type==x)]
  V<-WES$TST500[which(WES$Cancer.Type==x)]
  mean(sapply(X = 1:1000,FUN = bs,u=U, v=V))})   


#true.mean<-sapply(tissue.list,FUN = function(x) {
#  mean(WES$WES[which(WES$Cancer.Type==x)])}   )


#text(x=seq(dim(slopes)[2]),y = true.slopes[order(colMeans(slopes))]*m,'_',col = 'blue',cex = 2)
tmp<-order(colVar(slopes))
#tmp<-order(true.slopes)
text(x=seq(dim(slopes)[2]),y = true.slopes[tmp],'_',col = 'blue',cex = 2)

text(x=seq(dim(slopes)[2]),y = bootstrap.slopes[tmp],'_',col = 'purple',cex = 2)

abline(h=1,lty=2)
#lines(lowess(true.slopes[tmp]),lwd = 2, col = 'blue')
#lines(lowess(bootstrap.slopes[tmp]),lwd = 2, col = 'purple')
##lines(lowess(bootstrap.slopes[tmp]/2+true.slopes[tmp]*m/2),lwd = 2, col = 'green')
##lines(lowess(true.mean[tmp]),lwd = 2, col = 'blue')


#lines(lowess(colMeans(slopes[,order(colVar(slopes))])/m),lwd = 2)

legend("topleft",c('model prediction','observations','bootstrap regression\n (after scaling)'),col = c('black','blue','purple'),lwd = c(2,2,2),lty = c(1,1,1))

boxplot(true.slopes.bs[,order(colVar(slopes))],las = 2,col = 'blue',ylab = 'slopes',ylim = c(0,2), names = tissue.list[tmp])
abline(h=1,lty=2)

barplot(unlist(lapply(tissue.list,function(x){mean(WES$WES[WES$Cancer.Type==x])}))[tmp],names.arg = tissue.list[tmp],las=2,ylab = 'average TMB (WES)')

##### regression bias ####

TMB.sim <- function(tissue){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]*b
  #res <- as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #Fn<-ecdf(d)
  N<-50
  rand.samples.x = rep(NA,N)
  rand.samples.y = rep(NA,N)
  constant_noise= 0.5
  S=sample(d,size = N,replace = TRUE)
  for(i in 1:N){
    POISSON_NOISE=rpois(1,1)
    #rand.samples.x[i] = (rnorm(1,mean  = (res[1]),sd = (res[2])))
    #rand.samples.x[i] =   as.numeric(quantile(Fn, runif(1)))/m
    rand.samples.x[i] = S[i]
    rand.samples.y[i] = rand.samples.x[i]*m+rnorm(1,0,sd = constant_noise+sqrt((rand.samples.x[i])*a/b)) + POISSON_NOISE
    rand.samples.x[i] = rand.samples.x[i] + POISSON_NOISE
  }
  #rand.samples.x[i]=rand.samples.x[i]+P_noise/m}
  idx<-which((rand.samples.y>=-0) & rand.samples.x>=-0)
  X<-as.data.frame(cbind(round(rand.samples.x[idx])/b,round(rand.samples.y[idx])/a))
  return(round(lm(X[,2]~X[,1])$coefficients[1],3))}

tissue.list<-as.character(unique(WES$Cancer.Type))
offset<-sapply(X = tissue.list,function(t) replicate(50,TMB.sim(t)))


true.offset<-sapply(tissue.list,FUN = function(x) {
  A<-(WES$TST500[which(WES$Cancer.Type==x)])
  B<-(WES$WES[which(WES$Cancer.Type==x)])
  round(lm(A~B)$coefficients[1],3)})     

boxplot(true.offset,offset,col = 'red',names = c('TCGA','model'),ylab = 'offset',outline = FALSE)
boxplot(true.offset,offset,offset.no.noise,offset.no.driver, offset.no.driver.no.noise,col = 'red',names = c('TCGA','model','w/o noise','w/o driver\nmutations','w/o noise or\ndriver mutations'),ylab = 'offset',outline = FALSE)

abline(h = 1,lty = 2)
abline(h = 0,lty = 2)

#### part II

TMB.sim <- function(x){
  #d<-WES$WES[which(WES$Cancer.Type==tissue)]*b
  d<-WES$WES*b
  #res <- as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #Fn<-ecdf(d)
  N<-50
  rand.samples.x = rep(NA,N)
  rand.samples.y = rep(NA,N)
  constant_noise= 0.5
  S=sample(d,size = N,replace = TRUE)
  for(i in 1:N){
    POISSON_NOISE=rpois(n = 1,lambda = 0)
    #rand.samples.x[i] = (rnorm(1,mean  = (res[1]),sd = (res[2])))
    #rand.samples.x[i] =   as.numeric(quantile(Fn, runif(1)))/m
    rand.samples.x[i] = S[i]
    rand.samples.y[i] = rand.samples.x[i]*m+0*rnorm(1,0,sd = constant_noise+sqrt((rand.samples.x[i])*x/b)) + POISSON_NOISE
    rand.samples.x[i] = rand.samples.x[i] + POISSON_NOISE
  }
  #rand.samples.x[i]=rand.samples.x[i]+P_noise/m}
  idx<-which((rand.samples.y>=0) & rand.samples.x>=0)
  X<-as.data.frame(cbind((rand.samples.x[idx])/b,(rand.samples.y[idx])/x))
  return(round(lm(X[,2]~X[,1])$coefficients[1],3))}


plot(seq(from = 0.1,by = 0.1,to = 4),rowMeans(matrix(ncol = 500, byrow = TRUE, unlist(lapply(X = seq(from = 0.1,by = 0.1,to = 4),function(t) replicate(500,TMB.sim(t)))))),col = 'red', type = 'l', lwd = 2, xlab = 'Panel size (L) [Mbp]', ylab = 'offset',ylim = c(-1,10))
lines(seq(from = 0.1,by = 0.1,to = 4),rowMeans(matrix(ncol = 500, byrow = TRUE, unlist(lapply(X = seq(from = 0.1,by = 0.1,to = 4),function(t) replicate(500,TMB.sim(t)))))),col = 'black', type = 'l', lwd = 2, xlab = 'Panel size (L) [Mbp]', ylab = 'offset')
lines(seq(from = 0.1,by = 0.1,to = 4),rowMeans(matrix(ncol = 500, byrow = TRUE, unlist(lapply(X = seq(from = 0.1,by = 0.1,to = 4),function(t) replicate(500,TMB.sim(t)))))),col = 'green', type = 'l', lwd = 2, xlab = 'Panel size (L) [Mbp]', ylab = 'offset')
lines(seq(from = 0.1,by = 0.1,to = 4),rowMeans(matrix(ncol = 500, byrow = TRUE, unlist(lapply(X = seq(from = 0.1,by = 0.1,to = 4),function(t) replicate(500,TMB.sim(t)))))),col = 'blue', type = 'l', lwd = 2, xlab = 'Panel size (L) [Mbp]', ylab = 'offset')

legend(x = 1, y = 6, legend = c('full model','w/o noise','w/o driver mutations','w/o noise or driver mutations'), col=c('red','black','green', 'blue'),lty = rep(1,4))

### Survival analysis

OS<-read.csv('~/Desktop/Rizvi_tmb_res.csv')#read OS data
rsp<-ifelse(OS$response=='PR',1,0)
OS$response<-rsp
OS$tmb<-OS$tmb/b
plot(OS$tmb,rsp,xlab = "Somatic mutations [1/Mbp]", ylab = "P(Response to PD-L1 inhibitor)",pch = 3,ylim = c(0,1.1))

f = function(rate, shape = 10) {
  sum((OS$response-pgamma(OS$tmb, shape = shape, rate = rate))^2)
}

# arbitrary start (probability parameter is important to remain as close as possible)
start0 <- list("rate"=0.1) 

fcn <- function(x) f(x[1])   
res<-0
res <- spg(par=as.numeric(start0), fn=fcn, lower=c(0), upper=c(500), control=list( maxit=1000 ) ,quiet = TRUE )
#res

x<-seq(min(OS$tmb),max(OS$tmb))
lines(x,pgamma(x, shape = shape, rate = res$par),type = 'l',lwd = 2)

wait_time_iter<-function(shape){# test error for different wait time > n for iter over n

  # arbitrary start (probability parameter is important to remain as close as possible)
  start0 <- list("rate"=0.1) 
  
  fcn <- function(x) f(x[1], shape)   
  res<-0
  res <- spg(par=as.numeric(start0), fn=fcn, lower=c(0), upper=c(500), control=list( maxit=1000 ) ,quiet = TRUE )
  return(res$value)
}
plot(1:10,sapply(wait_time_iter,X = 1:10),xlab = "number of epitopes (n)", ylab = "OLS cost function value", type = 'b', lwd = 2)

f.sampling = function(rate, shape = 1) {
  sampling.rate = 0.1
  idx<-sample(1:length(OS$response), ceiling(length(OS$response)*sampling.rate), replace = FALSE)
  sum((OS$response[-idx]-pgamma(OS$tmb[-idx], shape = shape, rate = rate))^2)
}

fit.sampling <- function(shape){
  # arbitrary start (probability parameter is important to remain as close as possible)
  start0 <- list("rate"=0.1) 
  
  fcn <- function(x) f.sampling(x[1], shape)   
  res<-0
  res <- spg(par=as.numeric(start0), fn=fcn, lower=c(0), upper=c(500), control=list( maxit=1000 ) ,quiet = TRUE )
  return(res$value)
}

OLS.sampling<-matrix(data = 0,nrow = 3,ncol = 10)
for(i in 1:3) OLS.sampling[i,]<-sapply(1:10,function(x) fit.sampling(x))
boxplot(OLS.sampling,xlab = "number of epitopes (n)", ylab = "OLS cost function value")
##
shape=1

Psi_rate<-res$par
Psi<-function(t) pgamma(t,shape,Psi_rate)
tmb.range = seq(0,100,1)
plot(tmb.range,Psi(tmb.range),type = 'l',lwd =3,xlab = 'TMB',ylab='Response',xlim =c(0,40),pch = '+')
lines(OS$tmb,OS$response,type = 'p',cex = 1, pch='+')



tissue='LUSC'
d<-WES$WES[which(WES$Cancer.Type==tissue)]
#res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
myDensity<-density(d)
Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]

Theta<-function(t,L,tau){ 
  WES<-t*b
  constant_noise= 0.5
  tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L
  KEEP<-tmp>0
  mean(ifelse(tmp[KEEP]>tau,1,0))
}
  
  
P<-function(L,tau) sum(sapply(seq(min(OS$tmb),max(OS$tmb),1/b), function(t) Phi(t)*Psi(t)*Theta(t,L,tau)))
Z<-function(L,tau) sum(sapply(seq(min(OS$tmb),max(OS$tmb),1/b), function(t) Phi(t)*Theta(t,L,tau)))  
i<-0

tau.range<-seq(1,20,0.2)
L.range<-seq(0.1,2,0.01)

H<-matrix(0,nrow = length(L.range),ncol = length(tau.range))
i<-0
for(L in L.range){
  i<-i+1
  j<-0
  for(tau in tau.range){
    j<-j+1
    H[i,j]<-P(L,tau)/Z(L,tau)
    
  }
}
row.names(H)<-L.range
require(gplots)
#my_palette<-colorRampPalette(c('brown', 'red', 'orange', 'yellow', 'green', 'blue', 'violet'))
#my_palette<-colorRampPalette(c('brown', brewer.pal(11, "Spectral"),'violet'))
my_palette<-colorRampPalette(c('brown', brewer.pal(9, "Set1"),'violet'))
par(cex.main=0.8)
heatmap.2(H,trace = 'none',col = my_palette,dendrogram='none',Rowv=FALSE,Colv=FALSE)

WES.rsp<-sapply(X = tau.range, FUN = function(t) P(35.6,t)/Z(35.6,t)) # WES response
heatmap.2(as.matrix(rbind(WES.rsp,WES.rsp)), trace = 'none',col = my_palette,dendrogram='none',Rowv=FALSE,Colv=FALSE)

heatmap.2((1-H)/(1-0.3),trace = 'none',col = my_palette,dendrogram='none',Rowv=FALSE,Colv=FALSE)
heatmap.2((1-as.matrix(rbind(WES.rsp,WES.rsp)))/(1-0.3), trace = 'none',col = my_palette,dendrogram='none',Rowv=FALSE,Colv=FALSE)


auc<-sapply(seq(0.1,4,0.1), function(l){
  var<-seq(0,1,0.01)
  for(i in var){mat>i
  as.numeirc(Phi(t)*Psi(t)*Theta(t,L,tau)>i)
  }
  FP<-
  TP<-
  AUC(FP,TP)
})
  

### sample size effect

OS<-read.csv('~/Desktop/Rizvi_tmb_res.csv')#read OS data
rsp<-ifelse(OS$response=='PR',1,0)
OS$response<-rsp
OS$tmb<-OS$tmb/b
plot(OS$tmb,OS$response,xlab = "Somatic mutations [1/Mbp]", ylab = "P(Response to PD-L1 inhibitor)",pch = 3,ylim = c(0,1.1))

f = function(rate, shape = 3) {
  sum((OS$response-pgamma(OS$tmb, shape = shape, rate = rate))^2)
}

# arbitrary start (probability parameter is important to remain as close as possible)
start0 <- list("rate"=0.1) 

fcn <- function(x) f(x[1])   
res<-0
res <- spg(par=as.numeric(start0), fn=fcn, lower=c(0), upper=c(500), control=list( maxit=1000 ) ,quiet = TRUE )
#res

x<-seq(min(OS$tmb),max(OS$tmb))
lines(x,pgamma(x, shape = shape, rate = res$par),type = 'l',lwd = 2)

wait_time_iter<-function(shape){# test error for different wait time > n for iter over n
  
  # arbitrary start (probability parameter is important to remain as close as possible)
  start0 <- list("rate"=0.1) 
  
  fcn <- function(x) f(x[1], shape)   
  res<-0
  res <- spg(par=as.numeric(start0), fn=fcn, lower=c(0), upper=c(500), control=list( maxit=1000 ) ,quiet = TRUE )
  return(res$value)
}
plot(1:10,sapply(wait_time_iter,X = 1:10),xlab = "number of epitopes (n)", ylab = "OLS cost function value", type = 'b', lwd = 2)


fit.sampling <- function(shape){
  # arbitrary start (probability parameter is important to remain as close as possible)
  start0 <- list("rate"=0.1) 
  
  fcn <- function(x) f.sampling(rate = x[1], shape)   
  res<-0
  res <- spg(par=as.numeric(start0), fn=fcn, lower=c(0), upper=c(500), control=list( maxit=1000 ) ,quiet = TRUE )
  return(res$value)
}

OLS.sampling<-matrix(0,10,10)
for(i in 1:10) OLS.sampling[i,]<-sapply(1:10,function(x) fit.sampling(x))
boxplot(OLS.sampling, ylab = 'cost function', xlab= 'shape', main = 'subsampling')



f.sampling = function(rate = rate, shape = 3, n) {
  X<-OS
  i<-1:n
  sum((X$response[i]-pgamma(X$tmb[i], shape = shape, rate = rate))^2)
}

patient.sampling <- function(sampling.rate){
  shape = 3
  start0 <- list("rate"=0.1) 
  fcn <- function(x) f.sampling(rate = x[1], n = sampling.rate, shape = shape)   
  res<-0
  res <- spg(par=as.numeric(start0), fn=fcn, lower=c(0), upper=c(500), control=list( maxit=1000 ) ,quiet = TRUE )
  print(res)
  return(pgamma(10, rate = res$par, shape = shape))
}

size.sampling<-matrix(0,30,34)

for(i in 1:30) {
  OS<-OS[sample(nrow(OS)),]
  size.sampling[i,]<-sapply(1:34,function(x) patient.sampling(x))
}
boxplot((1-size.sampling)/(1-median(size.sampling[,34])), ylab = 'HR', xlab= 'trial size',col = 'red')


#### binary model and PPV analysis
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1.300000
b<-35.600000
m<-a/b
threshold = 6
tissue.list=c('LUSC', 'LUAD','PRAD','BRCA')
#tissue.list=c('LUSC')

Psi<-function(t) ifelse(t>threshold,1,0)

Theta<-function(t,L,tau){ 
  WES<-t*b
  constant_noise= 0.5
  tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L
  KEEP<-tmp>0
  mean(ifelse(tmp[KEEP]>tau,1,0))
}

iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.1,4,0.01)
  tmb.range = seq(0,20,1)
  d<-WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*Theta(t,L,tau)))
  Pos<-function(L,tau) sum(sapply(seq(min(tmb.range),max(tmb.range),1/b), function(t) Phi(t)*Theta(t,L,tau)))
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/Pos(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'PPV', type = 'l', lwd = 3, ylim = c(0,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/Pos(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'PPV', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}

legend(2,0.4, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')

#### binary model and NPV analysis

iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.1,4,0.01)
  tmb.range = seq(0,20,1)
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))))
  Neg<-function(L,tau) sum(sapply(seq(min(tmb.range),max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))))
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/Neg(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'NPV', type = 'l', lwd = 3, ylim = c(0,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/Neg(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'NPV', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}

legend(3,0.4, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')

#### panel accuracy analysis
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1.300000
b<-35.600000
m<-a/b
threshold = 8
#threshold = 14

#tissue.list=c('LUSC','BRCA','PRAD')
tissue.list=c('LUSC','LUAD')
mix=TRUE
#tissue.list=c('LUSC')

ERROR = 0

Psi<-function(t) ifelse(t>threshold,1,0)

Theta<-function(t,L,tau){
  WES<-t*b
  constant_noise= 0.5
  tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L+ERROR
  KEEP<-tmp>0
  mean(ifelse(tmp[KEEP]>tau,1,0))
}

iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.1,4,0.01)
  tmb.range = seq(0,20,1)
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*Theta(t,L,tau)),na.rm = TRUE)
  Pos<-function(L,tau) sum(sapply(seq(min(tmb.range),max(tmb.range),1/b), function(t) Phi(t)*Theta(t,L,tau)),na.rm = TRUE)
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/Pos(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'PPV', type = 'l', lwd = 3, ylim = c(0,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/Pos(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'PPV', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}

legend(3,0.4, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')

#### binary model and NPV analysis

iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.1,4,0.01)
  tmb.range = seq(0,20,1)
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  Neg<-function(L,tau) sum(sapply(seq(min(tmb.range),max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/Neg(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'NPV', type = 'l', lwd = 3, ylim = c(0,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/Neg(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'NPV', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}

legend(3,0.4, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')


#### binary model and PPA analysis
threshold = 8
tissue.list=c('LUSC','LUAD','PRAD','BRCA')

iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  FN<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPA', type = 'l', lwd = 3, ylim = c(0.6,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPA', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}



if(mix==TRUE){
  ### checkmate227
  checkmate227.sqms.ratio <- c(floor(0.29*480),floor(0.71*480))
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  FN<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPA', type = 'l', lwd = 3, col = iter+1)
  iter<-iter+1
  tissue.list<-c(tissue.list, 'CheckMate227 Mixture')
  
  ### threshold setting study
  #checkmate227.sqms.ratio <- c(floor(0.5*480),floor(0.5*480))
  #panel.range = seq(0.2,4,0.05)
  #tmb.range = seq(0,100,1)
  ##d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
  #d<-TS$WES_NonSynonymousTMB
  ##res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #myDensity<-density(d)
  #Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  #FN<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  #lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPA', type = 'l', lwd = 3, col = iter+1)
  #tissue.list<-c(tissue.list, 'Threshold Study Mixture')
  #iter<-iter+1
  
  legend(1,0.8, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')
}


#### binary model and NPA analysis
threshold = 8
tissue.list=c('LUSC','LUAD','PRAD','BRCA')

iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  FP<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/(TN(x,threshold)+FP(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'NPA', type = 'l', lwd = 3, ylim = c(0.3,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/(TN(x,threshold)+FP(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'NPA', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}


if(mix==TRUE){
  ### checkmate227
  checkmate227.sqms.ratio <- c(floor(0.29*480),floor(0.71*480))
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  FP<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  lines(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/(TN(x,threshold)+FP(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'NPA', type = 'l', lwd = 3, col = iter+1)
  iter<-iter+1
  tissue.list<-c(tissue.list, 'CheckMate227 Mixture')
  
  ### threshold setting study
  #checkmate227.sqms.ratio <- c(floor(0.5*480),floor(0.5*480))
  #panel.range = seq(0.2,4,0.05)
  #tmb.range = seq(0,100,1)
  #d<-TS$WES_NonSynonymousTMB
  ##d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
  ##res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #m#yDensity<-density(d)
  #Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  #FP<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  #lines(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/(TN(x,threshold)+FP(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'NPA', type = 'l', lwd = 3, col = iter+1)
  #tissue.list<-c(tissue.list, 'Threshold Study Mixture')
  #iter<-iter+1
  
  legend(1,0.6, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')
}



#### binary model and OPA analysis
threshold = 8
tissue.list=c('LUSC','LUAD','PRAD','BRCA')
#threshold = 8.2
iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  Tot<-function(L,tau) sum(sapply(seq(min(tmb.range),max(tmb.range),1/b), function(t) Phi(t)),na.rm = TRUE)
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) (TN(x,threshold)+TP(x,threshold))/Tot(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'OPA', type = 'l', lwd = 3, ylim = c(0.6,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) (TN(x,threshold)+TP(x,threshold))/Tot(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'OPA', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}

legend(3,0.4, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')

if(mix==TRUE){
### checkmate227
    checkmate227.sqms.ratio <- c(floor(0.29*480),floor(0.71*480))
    panel.range = seq(0.2,4,0.05)
    tmb.range = seq(0,100,1)
    d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
    #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
    myDensity<-density(d)
    Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
    TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
    TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
    Tot<-function(L,tau) sum(sapply(seq(min(tmb.range),max(tmb.range),1/b), function(t) Phi(t)),na.rm = TRUE)
    lines(panel.range, sapply(panel.range,FUN = function(x) (TN(x,threshold)+TP(x,threshold))/Tot(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'OPA', type = 'l', lwd = 3, col = iter+1)
    iter<-iter+1
    tissue.list<-c(tissue.list, 'CheckMate227 Mixture')
    
### threshold setting study
    #checkmate227.sqms.ratio <- c(floor(0.5*480),floor(0.5*480))
    #panel.range = seq(0.2,4,0.05)
    #tmb.range = seq(0,100,1)
    #d<-TS$WES_NonSynonymousTMB
    ##d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
    ##res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
    #myDensity<-density(d)
    #Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
    #TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
    #TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
    #Tot<-function(L,tau) sum(sapply(seq(min(tmb.range),max(tmb.range),1/b), function(t) Phi(t)),na.rm = TRUE)
    #lines(panel.range, sapply(panel.range,FUN = function(x) (TN(x,threshold)+TP(x,threshold))/Tot(x,threshold)), xlab = 'Panel size (L) [Mbp]', ylab = 'OPA', type = 'l', lwd = 3, col = iter+1)
    #tissue.list<-c(tissue.list, 'Threshold Study Mixture')
    #iter<-iter+1
    
    legend(1,0.8, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')
}



#### binary model and PPV analysis ########
threshold = 2
tissue.list=c('LUSC','LUAD','PRAD','BRCA')
#tissue.list=c('PRAD')

iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  FP<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FP(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPV', type = 'l', lwd = 3, ylim = c(0,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FP(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPV', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}



if(mix==TRUE){
  ### checkmate227
  checkmate227.sqms.ratio <- c(floor(0.29*480),floor(0.71*480))
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  FP<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FP(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPV', type = 'l', lwd = 3, col = iter+1)
  iter<-iter+1
  tissue.list<-c(tissue.list, 'CheckMate227 Mixture')
  
  ### threshold setting study
  #checkmate227.sqms.ratio <- c(floor(0.5*480),floor(0.5*480))
  #panel.range = seq(0.2,4,0.05)
  #tmb.range = seq(0,100,1)
  ##d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
  #d<-TS$WES_NonSynonymousTMB
  ##res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #myDensity<-density(d)
  #Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  #FN<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  #lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPA', type = 'l', lwd = 3, col = iter+1)
  #tissue.list<-c(tissue.list, 'Threshold Study Mixture')
  #iter<-iter+1
  
  #legend(1,0.8, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')
}

###### binary model and NPV analysis ########
threshold = 12
tissue.list=c('LUSC','LUAD','PRAD','BRCA')

iter<-0
for(tissue in tissue.list){
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  FN<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  if(iter == 0){
    plot(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/(TN(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'NPV', type = 'l', lwd = 3, ylim = c(0.4,1), main = paste('threshold =', threshold))
  } else {
    lines(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/(TN(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'NPV', type = 'l', lwd = 3, col = iter+1)
  }
  iter <- iter + 1
}



if(mix==TRUE){
  ### checkmate227
  checkmate227.sqms.ratio <- c(floor(0.29*480),floor(0.71*480))
  panel.range = seq(0.2,4,0.05)
  tmb.range = seq(0,100,1)
  d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  TN<-function(L,tau) sum(sapply(seq(min(tmb.range),tau,1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  FN<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  lines(panel.range, sapply(panel.range,FUN = function(x) TN(x,threshold)/(TN(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'NPV', type = 'l', lwd = 3, col = iter+1)
  iter<-iter+1
  tissue.list<-c(tissue.list, 'CheckMate227 Mixture')
  
  ### threshold setting study
  #checkmate227.sqms.ratio <- c(floor(0.5*480),floor(0.5*480))
  #panel.range = seq(0.2,4,0.05)
  #tmb.range = seq(0,100,1)
  ##d<-c(sample(WES$WES[which(WES$Cancer.Type=='LUAD')],checkmate227.sqms.ratio[2]),sample(WES$WES[which(WES$Cancer.Type=='LUSC')],checkmate227.sqms.ratio[1]))
  #d<-TS$WES_NonSynonymousTMB
  ##res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #myDensity<-density(d)
  #Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #TP<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(Theta(t,L,tau))),na.rm = TRUE)
  #FN<-function(L,tau) sum(sapply(seq(tau,max(tmb.range),1/b), function(t) Phi(t)*(1-Theta(t,L,tau))),na.rm = TRUE)
  #lines(panel.range, sapply(panel.range,FUN = function(x) TP(x,threshold)/(TP(x,threshold)+FN(x,threshold))), xlab = 'Panel size (L) [Mbp]', ylab = 'PPA', type = 'l', lwd = 3, col = iter+1)
  #tissue.list<-c(tissue.list, 'Threshold Study Mixture')
  #iter<-iter+1
  
  #legend(1,0.8, lty = rep(1,iter-1), col = 1:iter, legend = tissue.list, bty = 'n')
}


























#####

#### Universal cutoff study
# Part 1.
require(BB)
require(MASS)
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1
b<-35.600000
m<-a/b
tissue.list=c('LUSC','LUAD','PRAD','BRCA')
i=1
surv.P.small = rep(0,length(tissue.list))
surv.P.large = rep(0,length(tissue.list))
surv.P.wes = rep(0,length(tissue.list))
surv.P.baseline = rep(0,length(tissue.list))
median.tmb = rep(0,length(tissue.list))
threshold = 2
for(tissue in tissue.list){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #Psi<-function(t) ifelse(t>threshold,1,0)
  #Psi<-function(t) ifelse(t>threshold,0.8,0.2)
  Psi<-function(t) pgamma(t,shape=1,rate = 0.04404708)
  #Psi<-function(t) pgamma(t,shape=2,rate = 0.1200846)
  
  #Psi<-function(t) ifelse(t<20,t/20,1)
  Theta<-function(t,L,tau){ 
    WES<-t*b
    constant_noise= 0.5
    tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L
    KEEP<-tmp>0
    mean(ifelse(tmp[KEEP]>tau,1,0))
  }
  
  
  P<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Psi(t)*Theta(t,L,tau)))
  Z<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Theta(t,L,tau)))  
  
  surv.P.small[i]<-P(1,threshold)/Z(1,threshold)
  surv.P.large[i]<-P(2,threshold)/Z(2,threshold)
  surv.P.wes[i]<-P(b,threshold)/Z(b,threshold)
  surv.P.baseline[i]<-P(b,0)/Z(b,0)
  #print(P(b,threshold))
  #print(Z(b,threshold))
  #small.panel.LOR[i]<-log(surv.P.small/(1-surv.P.small))
  #large.panel.LOR[i]<-log(surv.P.large/(1-surv.P.large))
  median.tmb[i]<-median(WES$WES[WES$Cancer.Type==tissue])
  i=i+1
}
idx<-order(median.tmb)
plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
#barplot(as.matrix(rbind(c(small.panel.LOR),c(large.panel.LOR))),beside = TRUE, n)
plot(median.tmb[idx], surv.P.small[idx],type = 'b', ylim = c(0,1), col='blue', ylab = 'survival probability (PPV)', xlab = 'median TMB', main = paste('threshold = ', threshold), pch = c(2,1,3,4))
lines(median.tmb[idx], surv.P.large[idx],type = 'b',col = 'red', pch = c(2,1,3,4))
lines(median.tmb[idx], surv.P.wes[idx],type = 'b',col = 'black', pch = c(2,1,3,4))
lines(median.tmb[idx], surv.P.baseline[idx],type = 'b',col = 'purple', lty = 2, pch = c(2,1,3,4))
#legend(x = 4, y = 0.35, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES', 'w/o TMB'), pch = c('*','*','*', '*'), col = c('blue', 'red', 'black','purple'), bty = 'n')
#legend(x = 4, y = 0.35, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES', 'w/o TMB'), col = c('blue', 'red', 'black','purple'), bty = 'n', lty = c(1,1,1,2))
#legend(x = 4, y = 0.35, tissue.list, bty = 'n', pch = 4:1)

#plot(median.tmb[idx], (1-surv.P.small[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',ylim = c(0,1), col='blue', ylab = 'HR', xlab = 'median TMB', main = paste('threshold = ', threshold))
#lines(median.tmb[idx], (1-surv.P.large[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',col = 'red')
#lines(median.tmb[idx], (1-surv.P.wes[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',col = 'black')
#legend(x = 4, y = 1, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES'), pch = c('*','*','*'), col = c('blue', 'red', 'black'), bty = 'n')


# part 2
require(BB)
require(MASS)
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1.300000
b<-35.600000
m<-a/b
tissue.list=c('LUSC','LUAD','PRAD','BRCA')
i=1
surv.P.small = rep(0,length(tissue.list))
surv.P.large = rep(0,length(tissue.list))
surv.P.wes = rep(0,length(tissue.list))
surv.P.baseline = rep(0,length(tissue.list))
median.tmb = rep(0,length(tissue.list))
threshold.list = rep(0,length(tissue.list))
for(tissue in tissue.list){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  threshold.list[i] = as.numeric(quantile(d, probs = seq(0, 1, by= 0.01))[81])
  threshold = threshold.list[i]
  print(threshold)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #Psi<-function(t) ifelse(t>threshold,1,0)
  #Psi<-function(t) ifelse(t>threshold,0.8,0.2)
  #Psi<-function(t) pgamma(t,shape=1,rate = threshold)
  #Psi<-function(t) pgamma(t,shape=1,rate = 0.04404708)
  RATE<-seq(0,2,by = 0.0001)[which.min(sapply(seq(0,2,by = 0.0001),FUN = function(x) abs(0.5-pgamma(threshold,shape = 1, rate = x))))]
  Psi<-function(t) pgamma(t,shape=1,rate = RATE)
  #Psi<-function(t) pgamma(t,shape=2,rate = 0.1200846)
  plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
  
  #Psi<-function(t) ifelse(t<20,t/20,1)
  Theta<-function(t,L,tau){ 
    WES<-t*b
    constant_noise= 0.5
    tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L
    KEEP<-tmp>0
    mean(ifelse(tmp[KEEP]>tau,1,0))
  }
  
  
  P<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Psi(t)*Theta(t,L,tau)))
  Z<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Theta(t,L,tau)))  
  
  surv.P.small[i]<-P(1,threshold)/Z(1,threshold)
  surv.P.large[i]<-P(2,threshold)/Z(2,threshold)
  surv.P.wes[i]<-P(b,threshold)/Z(b,threshold)
  surv.P.baseline[i]<-P(b,0)/Z(b,0)
  #print(P(b,threshold))
  #print(Z(b,threshold))
  #small.panel.LOR[i]<-log(surv.P.small/(1-surv.P.small))
  #large.panel.LOR[i]<-log(surv.P.large/(1-surv.P.large))
  median.tmb[i]<-median(WES$WES[WES$Cancer.Type==tissue])
  i=i+1
}
idx<-order(threshold.list)
#plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
#barplot(as.matrix(rbind(c(small.panel.LOR),c(large.panel.LOR))),beside = TRUE, n)
plot(threshold.list[idx], surv.P.small[idx],type = 'b', ylim = c(0,1), pch = c(2,1,4,3), col='blue', ylab = 'survival probability (PPV)', xlab = 'threshold', main = paste('threshold = top 20% for each histology'))
lines(threshold.list[idx], surv.P.large[idx],type = 'b',col = 'red', pch = c(2,1,4,3))
lines(threshold.list[idx], surv.P.wes[idx],type = 'b', col = 'black', pch = c(2,1,4,3))
lines(threshold.list[idx], surv.P.baseline[idx],type = 'b', col = 'purple', lty = 2, pch = c(2,1,4,3))

#legend(x = 4, y = 0.35, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES', 'w/o TMB'), pch = c('*','*','*', '*'), col = c('blue', 'red', 'black','purple'), bty = 'n')

# plot(median.tmb[idx], (1-surv.P.small[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',ylim = c(0,1), col='blue', ylab = 'HR', xlab = 'median TMB', main = paste('threshold = top 20% for each histology'))
# lines(median.tmb[idx], (1-surv.P.large[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',col = 'red')
# lines(median.tmb[idx], (1-surv.P.wes[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',col = 'black')
# legend(x = 4, y = 1, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES'), pch = c('*','*','*'), col = c('blue', 'red', 'black'), bty = 'n')


###### survival probability as a function panel size for different histologies
require(BB)
require(MASS)
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1.300000
b<-35.600000
m<-a/b
tissue.list=c('LUSC','LUAD','PRAD','BRCA')
#tissue.list=c('LUSC')

i=1
surv.P = list()
for(tissue in tissue.list){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  #threshold.list[i] = as.numeric(quantile(d, probs = seq(0, 1, by= 0.01))[81])
  threshold = 6
  #print(threshold)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #Psi<-function(t) ifelse(t>threshold,1,0)
  #Psi<-function(t) ifelse(t>threshold,0.8,0.2)
  #Psi<-function(t) pgamma(t,shape=1,rate = threshold)
  Psi<-function(t) pgamma(t,shape=1,rate = 0.04404708)
  #RATE<-seq(0,2,by = 0.0001)[which.min(sapply(seq(0,2,by = 0.0001),FUN = function(x) abs(0.5-pgamma(threshold,shape = 1, rate = x))))]
  #Psi<-function(t) pgamma(t,shape=1,rate = RATE)
  #Psi<-function(t) pgamma(t,shape=2,rate = 0.1200846)
  plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
  
  #Psi<-function(t) ifelse(t<20,t/20,1)
  Theta<-function(t,L,tau){ 
    WES<-t*b
    constant_noise= 0.5
    tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L
    KEEP<-tmp>0
    mean(ifelse(tmp[KEEP]>tau,1,0))
  }
  
  
  P<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Psi(t)*Theta(t,L,tau)))
  Z<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Theta(t,L,tau)))  
  
  surv.P[[i]]<-sapply(seq(0.1,4,by = 0.1), FUN = function(x) P(x,threshold)/Z(x,threshold))
  
  #print(Z(b,threshold))
  #small.panel.LOR[i]<-log(surv.P.small/(1-surv.P.small))
  #large.panel.LOR[i]<-log(surv.P.large/(1-surv.P.large))
  i=i+1
}

plot(seq(0.1,4,by = 0.1),surv.P[[1]],type = 'l', lwd = 3,ylab = 'survival probability', xlab = 'Panel size (L) [Mbp]', main = paste('threshold = 6') , ylim = c(0,0.6))
lines(seq(0.1,4,by = 0.1),surv.P[[2]],type = 'l', lwd = 3,ylab = 'survival probability', xlab = 'Panel size (L) [Mbp]', main = paste('threshold = 6'), col = 'red')
lines(seq(0.1,4,by = 0.1),surv.P[[3]],type = 'l', lwd = 3,ylab = 'survival probability', xlab = 'Panel size (L) [Mbp]', main = paste('threshold = 6'), col ='green')
lines(seq(0.1,4,by = 0.1),surv.P[[4]],type = 'l', lwd = 3,ylab = 'survival probability', xlab = 'Panel size (L) [Mbp]', main = paste('threshold = 6'), col = 'blue' )



### part 3 fraction of patients
require(BB)
require(MASS)
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1
b<-35.600000
m<-a/b
tissue.list=c('LUSC','LUAD','PRAD','BRCA')
i=1
frac.P.small = rep(0,length(tissue.list))
frac.P.large = rep(0,length(tissue.list))
frac.P.wes = rep(0,length(tissue.list))
median.tmb = rep(0,length(tissue.list))
threshold = 6
for(tissue in tissue.list){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  Psi<-function(t) ifelse(t>threshold,1,0)
  #Psi<-function(t) ifelse(t>threshold,0.8,0.2)
  #Psi<-function(t) pgamma(t,shape=1,rate = 0.04404708)
  #Psi<-function(t) pgamma(t,shape=2,rate = 0.1200846)
  
  #Psi<-function(t) ifelse(t<20,t/20,1)
  Theta<-function(t,L,tau){ 
    WES<-t*b
    constant_noise= 0.5
    tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L
    KEEP<-tmp>0
    mean(ifelse(tmp[KEEP]>tau,1,0))
  }
  
  
  P<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Psi(t)*Theta(t,L,tau)))
  Z<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Theta(t,L,tau)))  
  
  frac.P.small[i]<-Z(1,threshold)/sum(sapply(FUN = Phi,seq(0,100,1/b)))
  frac.P.large[i]<-Z(2,threshold)/sum(sapply(FUN = Phi,seq(0,100,1/b)))
  frac.P.wes[i]<-Z(b,threshold)/sum(sapply(FUN = Phi,seq(0,100,1/b)))
  median.tmb[i]<-median(WES$WES[WES$Cancer.Type==tissue])
  i=i+1
}
idx<-order(median.tmb)
plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
#barplot(as.matrix(rbind(c(small.panel.LOR),c(large.panel.LOR))),beside = TRUE, n)
plot(median.tmb[idx], frac.P.small[idx],type = 'b', ylim = c(0,1), col='blue', ylab = 'fraction of patients selected for treatment', xlab = 'median TMB', main = paste('threshold = ', threshold), pch = c(2,1,3,4))
lines(median.tmb[idx], frac.P.large[idx],type = 'b',col = 'red', pch = c(2,1,3,4))
lines(median.tmb[idx], frac.P.wes[idx],type = 'b',col = 'black', pch = c(2,1,3,4))


#### survival by keeping the treatment fraction constant
require(BB)
require(MASS)
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1.300000
b<-35.600000
m<-a/b

Theta<-function(t,L,tau){ 
  WES<-t*b
  constant_noise= 0.5
  tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L
  KEEP<-tmp>0
  mean(ifelse(tmp[KEEP]>tau,1,0))
}

threshold_finder=function(size,tissue){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  
  Z<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Theta(t,L,tau)))  
  
  fnToFindRoot = function(r) {
    return((Z(size,r)/Z(size,0.01)-0.2)^2)
    #return((sin(r)-0.5)^2)

  }
  # arbitrary starting values
  #r0 = 0
  # minimise the function to get the parameter estimates
  #rootSearch = optim(r0, fnToFindRoot,method = 'BFGS', hessian=TRUE)
  return(seq(0,15,by = 0.1)[which.min(sapply(fnToFindRoot,X = seq(0,15,by = 0.1)))])
  #return(rootSearch$par)
}


tissue.list=c('LUSC','LUAD','PRAD','BRCA')
i=1
threshold.list.small = rep(0,length(tissue.list))
threshold.list.large = rep(0,length(tissue.list))
threshold.list.wes = rep(0,length(tissue.list))

surv.P.small = rep(0,length(tissue.list))
surv.P.large = rep(0,length(tissue.list))
surv.P.wes = rep(0,length(tissue.list))
surv.P.baseline = rep(0,length(tissue.list))
median.tmb = rep(0,length(tissue.list))

for(tissue in tissue.list){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #Psi<-function(t) ifelse(t>threshold,1,0)
  #Psi<-function(t) ifelse(t>threshold,0.8,0.2)
  #Psi<-function(t) pgamma(t,shape=1,rate = threshold)
  #Psi<-function(t) pgamma(t,shape=1,rate = 0.04404708)
  Psi<-function(t) pgamma(t,shape=10,rate = 0.8696841)
  #RATE<-seq(0,2,by = 0.0001)[which.min(sapply(seq(0,2,by = 0.0001),FUN = function(x) abs(0.5-pgamma(threshold,shape = 1, rate = x))))]
  #Psi<-function(t) pgamma(t,shape=1,rate = RATE)
  #Psi<-function(t) pgamma(t,shape=2,rate = 0.1200846)
 # plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
  
  #Psi<-function(t) ifelse(t<20,t/20,1)

  
  P<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Psi(t)*Theta(t,L,tau)))
  Z<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Theta(t,L,tau)))  
  
  threshold.list.small[i]=threshold_finder(1,tissue)
  threshold.list.large[i]=threshold_finder(2,tissue)
  threshold.list.wes[i]=threshold_finder(b,tissue)
  #threshold.list.baseline[i]=
  
  
  surv.P.small[i]<-P(1,threshold.list.small[i])/Z(1,threshold.list.small[i])
  surv.P.large[i]<-P(2,threshold.list.large[i])/Z(2,threshold.list.large[i])
  surv.P.wes[i]<-P(b,threshold.list.wes[i])/Z(b,threshold.list.wes[i])
  surv.P.baseline[i]<-P(b,0)/Z(b,0)
  #print(P(b,threshold))
  #print(Z(b,threshold))
  #small.panel.LOR[i]<-log(surv.P.small/(1-surv.P.small))
  #large.panel.LOR[i]<-log(surv.P.large/(1-surv.P.large))
  median.tmb[i]<-median(WES$WES[WES$Cancer.Type==tissue])
  i=i+1
}
idx<-order(threshold.list.small)
plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
#barplot(as.matrix(rbind(c(small.panel.LOR),c(large.panel.LOR))),beside = TRUE, n)
plot(threshold.list.small[idx], surv.P.small[idx],type = 'b', ylim = c(0,0.8), xlim = c(0,15), pch = c(2,1,4,3), col='blue', ylab = 'survival probability (PPV)', xlab = 'threshold', main = paste('threshold set to keep treatment fraction constant at 20%'))
lines(threshold.list.large[idx], surv.P.large[idx],type = 'b',col = 'red', pch = c(2,1,4,3))
lines(threshold.list.wes[idx], surv.P.wes[idx],type = 'b', col = 'black', pch = c(2,1,4,3))
#lines(threshold.list.baseline[idx], surv.P.baseline[idx],type = 'b', col = 'purple', lty = 2, pch = c(2,1,4,3))

#legend(x = 4, y = 0.35, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES', 'w/o TMB'), pch = c('*','*','*', '*'), col = c('blue', 'red', 'black','purple'), bty = 'n')

# plot(median.tmb[idx], (1-surv.P.small[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',ylim = c(0,1), col='blue', ylab = 'HR', xlab = 'median TMB', main = paste('threshold = top 20% for each histology'))
# lines(median.tmb[idx], (1-surv.P.large[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',col = 'red')
# lines(median.tmb[idx], (1-surv.P.wes[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',col = 'black')
# legend(x = 4, y = 1, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES'), pch = c('*','*','*'), col = c('blue', 'red', 'black'), bty = 'n')



plot(c(1,2,4),c(surv.P.small[4],surv.P.large[4],surv.P.wes[4]),type = 'b', ylim = c(0,0.8), pch = c(4,4,4), col='blue', ylab = 'survival probability (PPV)', xlab = 'Panel size (L) [Mbp]', main = paste('threshold set to keep treatment fraction constant'))
lines(c(1,2,4),c(surv.P.small[3],surv.P.large[3],surv.P.wes[3]),type = 'b', pch = c(4,4,4), col='red', ylab = 'survival probability (PPV)', xlab = 'Panel size (L) [Mbp]', main = paste('threshold set to keep treatment fraction constant'))
lines(c(1,2,4),c(surv.P.small[2],surv.P.large[2],surv.P.wes[2]),type = 'b', pch = c(4,4,4), col='red', ylab = 'survival probability (PPV)', xlab = 'Panel size (L) [Mbp]', main = paste('threshold set to keep treatment fraction constant'))
lines(c(1,2,4),c(surv.P.small[1],surv.P.large[1],surv.P.wes[1]),type = 'b', pch = c(4,4,4), col='red', ylab = 'survival probability (PPV)', xlab = 'Panel size (L) [Mbp]', main = paste('threshold set to keep treatment fraction constant'))


plot(c(threshold.list.small[1],threshold.list.large[1],threshold.list.wes[1]), c(surv.P.small[1],surv.P.large[1],surv.P.wes[1]),type = 'b', ylim = c(0,1), xlim = c(0,15), pch = c(4,4,4), col = c('blue','red','black'), ylab = 'survival probability (PPV)', xlab = 'threshold', main = paste('threshold set to keep treatment fraction constant at 20%'))
lines(c(threshold.list.small[2],threshold.list.large[2],threshold.list.wes[2]), c(surv.P.small[2],surv.P.large[2],surv.P.wes[2]),type = 'b',col = c('blue','red','black'), pch = c(3,3,3))
lines(c(threshold.list.small[3],threshold.list.large[3],threshold.list.wes[3]), c(surv.P.small[3],surv.P.large[3],surv.P.wes[3]),type = 'b', col = c('blue','red','black'), pch = c(2,2,2))
lines(c(threshold.list.small[4],threshold.list.large[4],threshold.list.wes[4]), c(surv.P.small[4],surv.P.large[4],surv.P.wes[4]),type = 'b', col = c('blue','red','black'), pch = c(1,1,1))

#####
#### the treatment fraction by keeping survival constant
require(BB)
require(MASS)
WES<-read.csv('~/Box/Sven`s TMB project/FOCR.TMB.9104.csv')
a<-1.300000
b<-35.600000
m<-a/b

Theta<-function(t,L,tau){ 
  WES<-t*b
  constant_noise= 0.5
  tmp<-(WES*(L/b)+rnorm(100,0,sd = constant_noise+sqrt(WES*L/b))+rpois(100,1))/L
  KEEP<-tmp>0
  mean(ifelse(tmp[KEEP]>tau,1,0))
}

threshold_finder=function(size,tissue){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #Psi<-function(t) pgamma(t,shape=1,rate = 0.04404708)
  Psi<-function(t) pgamma(t,shape=10,rate = 0.8696841)
  Z<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Theta(t,L,tau)))  
  P<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Psi(t)*Theta(t,L,tau)))
  
  fnToFindRoot = function(r) {
    return((P(size,r)/Z(size,r)-0.5)^2)
    #return((sin(r)-0.5)^2)
    
  }
  # arbitrary starting values
  #r0 = 0
  # minimise the function to get the parameter estimates
  #rootSearch = optim(r0, fnToFindRoot,method = 'BFGS', hessian=TRUE)
  return(seq(0,15,by = 0.1)[which.min(sapply(fnToFindRoot,X = seq(0,15,by = 0.1)))])
  #return(rootSearch$par)
}


tissue.list=c('LUSC','LUAD','PRAD','BRCA')
i=1
threshold.list.small = rep(0,length(tissue.list))
threshold.list.large = rep(0,length(tissue.list))
threshold.list.wes = rep(0,length(tissue.list))

surv.P.small = rep(0,length(tissue.list))
surv.P.large = rep(0,length(tissue.list))
surv.P.wes = rep(0,length(tissue.list))
surv.P.baseline = rep(0,length(tissue.list))
median.tmb = rep(0,length(tissue.list))

for(tissue in tissue.list){
  d<-WES$WES[which(WES$Cancer.Type==tissue)]
  #res<-as.numeric(fitdistr(d,densfun = "normal")$estimate)
  myDensity<-density(d)
  Phi<-function(t) myDensity$y[which.min(abs(myDensity$x-t))]
  #Psi<-function(t) ifelse(t>threshold,1,0)
  #Psi<-function(t) ifelse(t>threshold,0.8,0.2)
  #Psi<-function(t) pgamma(t,shape=1,rate = threshold)
  #Psi<-function(t) pgamma(t,shape=1,rate = 0.04404708)
  Psi<-function(t) pgamma(t,shape=10,rate = 0.8696841)
  #RATE<-seq(0,2,by = 0.0001)[which.min(sapply(seq(0,2,by = 0.0001),FUN = function(x) abs(0.5-pgamma(threshold,shape = 1, rate = x))))]
  #Psi<-function(t) pgamma(t,shape=1,rate = RATE)
  #Psi<-function(t) pgamma(t,shape=2,rate = 0.1200846)
  # plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
  
  #Psi<-function(t) ifelse(t<20,t/20,1)
  
  
  P<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Psi(t)*Theta(t,L,tau)))
  Z<-function(L,tau) sum(sapply(seq(0,100,1/b), function(t) Phi(t)*Theta(t,L,tau)))  
  
  threshold.list.small[i]=threshold_finder(1,tissue)
  threshold.list.large[i]=threshold_finder(2,tissue)
  threshold.list.wes[i]=threshold_finder(b,tissue)
  #threshold.list.baseline[i]=
  
  
  surv.P.small[i]<-Z(1,threshold.list.small[i])/Z(1,0.01)
  surv.P.large[i]<-Z(2,threshold.list.large[i])/Z(2,0.01)
  surv.P.wes[i]<-Z(b,threshold.list.wes[i])/Z(b,0.01)
  #surv.P.baseline[i]<-P(b,0)/Z(b,0)
  #print(P(b,threshold))
  #print(Z(b,threshold))
  #small.panel.LOR[i]<-log(surv.P.small/(1-surv.P.small))
  #large.panel.LOR[i]<-log(surv.P.large/(1-surv.P.large))
  median.tmb[i]<-median(WES$WES[WES$Cancer.Type==tissue])
  i=i+1
}
idx<-order(threshold.list.small)
#plot(1:40,Psi(1:40),type = 'l', lwd = 3, xlab = 'TMB', ylab = 'Response',ylim = c(0,1))
#barplot(as.matrix(rbind(c(small.panel.LOR),c(large.panel.LOR))),beside = TRUE, n)
plot(threshold.list.small[idx], surv.P.small[idx],type = 'b', ylim = c(0,0.8), xlim = c(0,15), pch = c(2,1,4,3), col='blue', ylab = 'treatment fraction', xlab = 'threshold', main = paste('threshold set to keep survival constant at 50%'))
lines(threshold.list.large[idx], surv.P.large[idx],type = 'b',col = 'red', pch = c(2,1,4,3))
lines(threshold.list.wes[idx], surv.P.wes[idx],type = 'b', col = 'black', pch = c(2,1,4,3))
#lines(threshold.list.baseline[idx], surv.P.baseline[idx],type = 'b', col = 'purple', lty = 2, pch = c(2,1,4,3))

#legend(x = 4, y = 0.35, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES', 'w/o TMB'), pch = c('*','*','*', '*'), col = c('blue', 'red', 'black','purple'), bty = 'n')

# plot(median.tmb[idx], (1-surv.P.small[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',ylim = c(0,1), col='blue', ylab = 'HR', xlab = 'median TMB', main = paste('threshold = top 20% for each histology'))
# lines(median.tmb[idx], (1-surv.P.large[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',col = 'red')
# lines(median.tmb[idx], (1-surv.P.wes[idx])/(1-surv.P.baseline[idx]),type = 'b', pch = '*',col = 'black')
# legend(x = 4, y = 1, c('panel size = 1 Mbp','panel size = 2 Mbp', 'WES'), pch = c('*','*','*'), col = c('blue', 'red', 'black'), bty = 'n')

plot(c(threshold.list.small[1],threshold.list.large[1],threshold.list.wes[1]), c(surv.P.small[1],surv.P.large[1],surv.P.wes[1]),type = 'b', ylim = c(0,0.8), xlim = c(0,15), pch = c(4,4,4), col = c('blue','red','black'), ylab = 'treatment fraction', xlab = 'threshold', main = paste('threshold set to keep survival constant at 50%'))
lines(c(threshold.list.small[2],threshold.list.large[2],threshold.list.wes[2]), c(surv.P.small[2],surv.P.large[2],surv.P.wes[2]),type = 'b',col = c('blue','red','black'), pch = c(3,3,3))
lines(c(threshold.list.small[3],threshold.list.large[3],threshold.list.wes[3]), c(surv.P.small[3],surv.P.large[3],surv.P.wes[3]),type = 'b', col = c('blue','red','black'), pch = c(2,2,2))
lines(c(threshold.list.small[4],threshold.list.large[4],threshold.list.wes[4]), c(surv.P.small[4],surv.P.large[4],surv.P.wes[4]),type = 'b', col = c('blue','red','black'), pch = c(1,1,1))

#plot(c(1,2,4),c(surv.P.small[4],surv.P.large[4],surv.P.wes[4]),type = 'b', ylim = c(0,0.8), pch = c(4,4,4), col='blue', ylab = 'survival probability (PPV)', xlab = 'Panel size (L) [Mbp]', main = paste('threshold set to keep treatment fraction constant'))
#lines(c(1,2,4),c(surv.P.small[3],surv.P.large[3],surv.P.wes[3]),type = 'b', pch = c(4,4,4), col='red', ylab = 'survival probability (PPV)', xlab = 'Panel size (L) [Mbp]', main = paste('threshold set to keep treatment fraction constant'))
#lines(c(1,2,4),c(surv.P.small[2],surv.P.large[2],surv.P.wes[2]),type = 'b', pch = c(4,4,4), col='red', ylab = 'survival probability (PPV)', xlab = 'Panel size (L) [Mbp]', main = paste('threshold set to keep treatment fraction constant'))
#lines(c(1,2,4),c(surv.P.small[1],surv.P.large[1],surv.P.wes[1]),type = 'b', pch = c(4,4,4), col='red', ylab = 'survival probability (PPV)', xlab = 'Panel size (L) [Mbp]', main = paste('threshold set to keep treatment fraction constant'))

