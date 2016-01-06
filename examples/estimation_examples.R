#==================#
# estimation_example.R -- time-varying estimation for four methods
#  (1) the dynamic Kendall's tau 
#  (2) the kernel Pearson CLIME estimator 
#  (3) the kernel graphical Lasso estimator 
#  (4) the kernel neighborhood selection estimator 
#==================#


library("huge")
library("corpcor")
library("mvtnorm")
library("JGL")
source('funs_est.R')


n = 200;
d = 20;
num.data = 1;
num.test = 10;
lambda = c(0.15,0.12,0.1,0.085,0.07,0.05,0.02,0.01) # tuning parameter
h = 1.35/n^(1/5); #bandwidth

nsim=1
FPR=array(0,dim=c(nsim,length(lambda),4))
FNR=array(0,dim=c(nsim,length(lambda),4))
eps=array(0,dim=c(nsim,3,4))



for(sim in 1:nsim){
  cat(sim,"in",nsim,"\n")
  
  ptm=proc.time()

# generate data
Xn = matrix(0,d,n*num.data)
Zn = matrix(0,1,n*num.data)
data = data.gen(n = n,d = d, num.data)
Xn = data$X
Zn = data$Z

test.idx = sort(sample(1:n,num.test))
z.test = data$zn[test.idx];
Omega.true = data$Omega
Omega.true.test = vector('list',num.test)
for (itest in 1:num.test){
Omega.true.test[[itest]] = Omega.true[[test.idx[itest]]]
}

# estimate


Omega.kendall = vector("list", length(z.test))
Omega.Sample = vector("list", length(z.test))
Omega.ZWL = vector("list", length(z.test))
Omega.Glasso = vector("list", length(z.test))
adj.kendall = vector("list", length(z.test))
adj.Sample = vector("list", length(z.test))
adj.ZWL = vector("list", length(z.test))
adj.Glasso = vector("list", length(z.test))
for (i in 1:length(z.test)){
  cat('Estimate Omega ', i, 'in',length(z.test),"\n")
  out.Sigma = tvKT(Zn,Xn,z=z.test[i],h = h);
  
  # the dynamic Kendall's tau 
  out1 = caliClime.path(out.Sigma$Sigma,lambda = lambda)
  Omega.kendall[[i]] = out1$Omega
  adj.kendall[[i]] = out1$adj
  
  # the kernel Pearson CLIME estimator 
  out2 = caliClime.path(out.Sigma$Sample.est,lambda = lambda)
  Omega.Sample[[i]] = out2$Omega
  adj.Sample[[i]] = out2$adj
  
  # the kernel graphical Lasso estimator 
  out3 = ZWL(out.Sigma$Sample.est,lambda = lambda)
  Omega.ZWL[[i]] = out3$icov
  adj.ZWL[[i]] = out3$path
  
  # the kernel neighborhood selection estimator 
  out4 = Glasso.KX(Zn,Xn,z=z.test[i],h,lambda=c(0.3,0.2,0.15,0.12,0.1,0.085,0.05,0.02))
  Omega.Glasso[[i]] = out4$Omega
  adj.Glasso[[i]] = out4$adj
  
  
}



#error evaluate
err.kendall = errs.avg(length(lambda),length(z.test),Omega.kendall,adj.kendall,Omega.true.test)
err.Sample = errs.avg(length(lambda),length(z.test),Omega.Sample,adj.Sample,Omega.true.test)
FPR[sim,,1]=err.kendall$FPR
FNR[sim,,1]=err.kendall$FNR
eps[sim,,1]=err.kendall$eps
FPR[sim,,2]=err.Sample$FPR
FNR[sim,,2]=err.Sample$FNR
eps[sim,,2]=err.Sample$eps


ZWL.ers=errs.avg(length(lambda),length(z.test),Omega.ZWL,adj.ZWL,Omega.true.test)
FPR[sim,,3]=ZWL.ers$FPR
FNR[sim,,3]=ZWL.ers$FNR
eps[sim,,3]=ZWL.ers$eps

Glasso.ers=errs.avg(length(lambda),length(z.test),Omega.Glasso,adj.Glasso,Omega.true.test)
FPR[sim,,4]=Glasso.ers$FPR
FNR[sim,,4]=Glasso.ers$FNR
eps[sim,,4]=Glasso.ers$eps


print(proc.time()-ptm)
}


save(FPR,FNR,eps,file="ersS01N200.RData")

# plot ROC curve
a11=sort((FPR[,,1]));a12=sort(1-(FNR[,,1]));a11=c(0,a11,1);a12=c(0,a12,1)
a21=sort((FPR[,,2]));a22=sort(1-(FNR[,,2]));a21=c(0,a21,1);a22=c(0,a22,1)
a31=sort((FPR[,,3]));a32=sort(1-(FNR[,,3]));a31=c(0,a31,1);a32=c(0,a32,1)
a41=sort((FPR[,,4]));a42=sort(1-(FNR[,,4]));a41=c(0,a41,1);a42=c(0,a42,1)


pdf("ROCsimS01N200.pdf")
plot(a11, a12,type='l',lwd=3, xlab='FPR', ylab= 'TPR', main="")
lines(a21, a22,type="l",lwd=3,lty=2, col=2)
lines(a31, a32,type="l",lwd=3,lty=4, col='blue')
lines(a41, a42,type="b",lwd=3,lty=6, col='green')
legend("bottomright",c("Kendall","Pearson","ZLW","KSAX"),lty=c(1,2,4,6),col=c(1,2,'blue','green'),lwd=c(3,3,3,3))
dev.off() 
