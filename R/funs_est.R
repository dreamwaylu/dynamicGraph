# functions for times-varying graph estimation

library("huge")
library("corpcor")
library("mvtnorm")
library("igraph")
library('fastclime')
library('clime')
library("glmnet")
library('glasso')


####################################
# Functions for Data Generation
####################################


#=============================#
# graph.gen
# Generate the time-varying precision matrix Markov to a banded graph
#
# IN: num.itv -- num of time intervals, num.band -- size of bands, L,U -- Unif[L,U]
# OUT:  list of precision matrix and adj matrix
#=============================#
graph.gen = function(d, num.itv = 5,num.band = 2, L=0.5, U = 0.9){
S=matrix(0,d,d)
for (i in 1:d){
  for (j in 1:d){
    if(j > i & (abs(j-i)<=num.band | abs(j-i)>= d- num.band))
    S[j,i] =1;
  }
}
tridiag.idx=which(S == 1)


num.edge = 25; 
num.change = 10; 

graph.idx = matrix(0, num.edge, num.itv)
theta = matrix(0, length(tridiag.idx), 2*num.itv+1)

graph.idx[,1] = sort(sample(tridiag.idx, num.edge))
theta[which(tridiag.idx %in% graph.idx[,1]),1] = runif(num.edge,L,U)
theta[which(tridiag.idx %in% graph.idx[,1]),2] = runif(num.edge,L,U)
                     
                     
for (i in 2:num.itv){
  notin.idx = tridiag.idx[which(! tridiag.idx %in% graph.idx[,i-1])]
  add.idx = sample(notin.idx, num.change);
  remove.idx = sample(graph.idx[,i-1], num.change)
  theta[which(tridiag.idx %in% add.idx),2*i] = runif(num.change,L,U)
  theta[which( ((tridiag.idx %in% graph.idx[,i-1])&(! tridiag.idx %in% remove.idx)) ), 2*i-1] = runif(num.edge - num.change,L,U)
  theta[which( ((tridiag.idx %in% graph.idx[,i-1])&(! tridiag.idx %in% remove.idx)) ), 2*i] = runif(num.edge - num.change,L,U)
  graph.idx[,i] = tridiag.idx[which( ((tridiag.idx %in% graph.idx[,i-1])&(! tridiag.idx %in% remove.idx)) | (tridiag.idx %in% add.idx) )];
}
theta[which(tridiag.idx %in% add.idx),2*num.itv+1] = runif(num.change,L,U)

return(list(theta = theta, graph = graph.idx, tri = tridiag.idx))
}


#=============================#
# prec.gen
# Generate the precision matrix enforcing psd
#=============================#
prec.gen = function(z,theta,tridiag.idx,d,num.itv = 5){
   z.idx = floor(2*z*num.itv)+1
   theta.z = theta[,z.idx] + 2*num.itv*(theta[,z.idx+1] - theta[,z.idx])%*%diag(z-(z.idx-1)/(2*num.itv));
   theta.z.psd = theta.z 
   for(i in 1:length(z)){
   Prec.z = matrix(0,d,d)
   Prec.z[tridiag.idx] = theta.z[,i]
   Prec.z = Prec.z + t(Prec.z)
   eig.prec=eigen(Prec.z, symmetric = TRUE, only.values = TRUE)
   eigmin = abs(min(eig.prec$values)) + 1;
   theta.z.psd[,i] = theta.z[,i]/eigmin;
   }
   return(theta.z.psd)
}

#=============================#
# prec2.gen
# Generate the precision matrix enforcing psd
#=============================#
prec.gen2 = function(z,theta,tridiag.idx,d,num.itv = 5,edge,e.strength){
  j0 = edge[1]
  k0 = edge[2]
  z.idx = floor(2*z*num.itv)+1
  theta.z = theta[,z.idx] + 2*num.itv*(theta[,z.idx+1] - theta[,z.idx])%*%diag(z-(z.idx-1)/(2*num.itv));
  theta.z.psd = theta.z 
  for(i in 1:length(z)){
    Prec.z = matrix(0,d,d)
    Prec.z[tridiag.idx] = theta.z[,i]
    Prec.z = Prec.z + t(Prec.z)
    Prec.z[j0,k0] = e.strength
    Prec.z[k0,j0] = e.strength
    eig.prec=eigen(Prec.z, symmetric = TRUE, only.values = TRUE)
    eigmin = abs(min(eig.prec$values)) + 1;
    theta.z.psd[,i] = Prec.z[tridiag.idx]/eigmin;
  }
  return(theta.z.psd)
}

#=============================#
# prec3.gen
# Generate the precision matrix enforcing psd
#=============================#
prec.gen3 = function(z,theta,tridiag.idx,d,num.itv = 5,edge.list,e.strength){
  z.idx = floor(2*z*num.itv)+1
  theta.z = theta[,z.idx] + 2*num.itv*(theta[,z.idx+1] - theta[,z.idx])%*%diag(z-(z.idx-1)/(2*num.itv));
  theta.z.psd = theta.z 
  for(i in 1:length(z)){
    Prec.z = matrix(0,d,d)
    Prec.z[tridiag.idx] = theta.z[,i]
    Prec.z = Prec.z + t(Prec.z)
    for(iedge in 1:length(edge.list)){
      j.idx = edge.list[j]
      j0 = floor((j.idx-1)/d)+1; k0 = j.idx%%d
      Prec.z[k0,j0] = e.strength
      Prec.z[j0,k0] = e.strength
    }
    eig.prec=eigen(Prec.z, symmetric = TRUE, only.values = TRUE)
    eigmin = abs(min(eig.prec$values)) + 1;
    theta.z.psd[,i] = Prec.z[tridiag.idx]/eigmin;
  }
  return(theta.z.psd)
}

#=============================#
# covtrue
# Generate the covariance matrix
#=============================#
covtrue = function(z,graph){
  theta = prec.gen(z,graph$theta,graph$tri,d)
    SigmaInv = matrix(0,d,d);
    SigmaInv[graph.n$tri] = theta
    SigmaInv = SigmaInv + t(SigmaInv) + diag(d);  
    Cov = solve(SigmaInv)
    Cov.diag = 1/sqrt(diag(Cov));
    Cov = diag(Cov.diag)%*% Cov %*% diag(Cov.diag)
    Omega = diag(1/Cov.diag)%*% SigmaInv %*% diag(1/Cov.diag)
    return(list(Cov = Cov, Omega = Omega))
}

#=============================#
# data.gen
# Generate the iid samples for the time-varying graphical model
#
# IN: n--sample size, d--dim, num.data--samples from each fixed time, num.itv--nume of time intervals
# OUT:  list of data ana parameters
#=============================#
data.gen = function(n,d,num.data=1,graph = NULL,num.itv = 5){
  zn = sort(runif(n,0,1))
  if(is.null(graph)){
    graph.n = graph.gen(d)
  }else
  {graph.n = graph}
  theta = prec.gen(zn,graph.n$theta,graph.n$tri,d,num.itv = num.itv)
  Xn = matrix(0,d,n*num.data)
  Zn = matrix(0,1,n*num.data)
  Omega = vector('list',n)
  Sigma = vector('list',n)
  for(i in 1:n){
    SigmaInv = matrix(0,d,d);
    SigmaInv[graph.n$tri] = theta[,i]
    SigmaInv = SigmaInv + t(SigmaInv) + diag(d);  
    Cov = solve(SigmaInv)
    Cov.diag = 1/sqrt(diag(Cov));
    Cov = diag(Cov.diag)%*% Cov %*% diag(Cov.diag)
    Sigma[[i]] = Cov
    Omega[[i]] = diag(1/Cov.diag)%*% SigmaInv %*% diag(1/Cov.diag)
    Xn[,((i-1)*num.data+1):(i*num.data)] = rmvnorm(n = num.data,rep(0,d),Cov)
    Zn[,((i-1)*num.data+1):(i*num.data)] = zn[i]
  }
  return(list(X = Xn, Z = Zn, Omega = Omega,Cov = Sigma, zn = zn, graph = graph.n))
}


#=============================#
# data.gen2
# Generate the iid samples for the time-varying graphical model
# use prec.gen2
#=============================#
data.gen2 = function(n,d,edge,e.strength,graph = NULL,num.data=1){
  j0 = edge[1]
  k0 = edge[2]
  zn = sort(runif(n,0,1))
  if(is.null(graph)){
    graph.n = graph.gen(d)
  }else
  {graph.n = graph}
  
  
  theta = prec.gen2(zn,graph.n$theta,graph.n$tri,d,edge = edge,e.strength = e.strength)
  Xn = matrix(0,d,n*num.data)
  Zn = matrix(0,1,n*num.data)
  Omega = vector('list',n)
  Sigma = vector('list',n)
  for(i in 1:n){
    SigmaInv = matrix(0,d,d);
    SigmaInv[graph.n$tri] = theta[,i]
    SigmaInv = SigmaInv + t(SigmaInv) + diag(d);  
    Cov = solve(SigmaInv)
    Cov.diag = 1/sqrt(diag(Cov));
    Cov = diag(Cov.diag)%*% Cov %*% diag(Cov.diag)
    Sigma[[i]] = Cov
    Omega[[i]] = diag(1/Cov.diag)%*% SigmaInv %*% diag(1/Cov.diag)
    Xn[,((i-1)*num.data+1):(i*num.data)] = rmvnorm(n = num.data,rep(0,d),Cov)
    Zn[,((i-1)*num.data+1):(i*num.data)] = zn[i]
  }
  return(list(X = Xn, Z = Zn, Omega = Omega,Cov = Sigma, zn = zn, graph = graph.n))
}

#=============================#
# data.gen3
# Generate the iid samples for the time-varying graphical model
# use prec.gen3
#=============================#
data.gen3 = function(n,d,edge.list,e.strength,graph = NULL,num.data=1){
  zn = sort(runif(n,0,1))
  if(is.null(graph)){
  graph.n = graph.gen(d)
  }else
  {graph.n = graph}
  theta = prec.gen3(zn,graph.n$theta,graph.n$tri,d,edge.list = edge.list,e.strength = e.strength)
  Xn = matrix(0,d,n*num.data)
  Zn = matrix(0,1,n*num.data)
  Omega = vector('list',n)
  Sigma = vector('list',n)
  for(i in 1:n){
    SigmaInv = matrix(0,d,d);
    SigmaInv[graph.n$tri] = theta[,i]
    SigmaInv = SigmaInv + t(SigmaInv) + diag(d);  
    Cov = solve(SigmaInv)
    Cov.diag = 1/sqrt(diag(Cov));
    Cov = diag(Cov.diag)%*% Cov %*% diag(Cov.diag)
    Sigma[[i]] = Cov
    Omega[[i]] = diag(1/Cov.diag)%*% SigmaInv %*% diag(1/Cov.diag)
    Xn[,((i-1)*num.data+1):(i*num.data)] = rmvnorm(n = num.data,rep(0,d),Cov)
    Zn[,((i-1)*num.data+1):(i*num.data)] = zn[i]
  }
  return(list(X = Xn, Z = Zn, Omega = Omega,Cov = Sigma, zn = zn, graph = graph.n))
}

####################################
# Functions for Estimation
####################################

# kernel function
kern = function(x,h){
  return(15/(16*h)*(1-(x/h)^2)^2*ifelse(abs((x/h)) <1, 1, 0))
}



#=============================#
# tvKT
# Calculate the kernel smoohted Kendall's Tau 
#
# IN: Zn--index var, Xn--variables, z--time, h--bandwidth
# OUT:  list of Kendall's tau matrix
#=============================#
tvKT = function(Zn,Xn,z,h){
   n = dim(Xn)[2]
   d = dim(Xn)[1]
   nonzero.idx = which(abs(Zn-z)<h);
   if (length(nonzero.idx) == 0){
     cat('Warning: kernel is zero')
   }
   mat.KT = matrix(0,d,d)
   mat.Sample = matrix(0,d,d)
   mat.ker = matrix(0,n,n)
   for (i1 in nonzero.idx){
     for  (i2 in nonzero.idx){
       if(i2>i1){
       mat.ker[i1,i2] = kern((Zn[i1] - z),h) * kern((Zn[i2] - z),h);
       }
     }
   }
   U.ker = sum(mat.ker)
   for (j in 1:(d-1)){
     for (k in (j+1):d ){
   tmp.sign = 0
   tmp.sample = 0
   for (i1 in nonzero.idx){
     tmp.sample = tmp.sample + kern((Zn[i1] - z),h) * Xn[j,i1]*Xn[k,i1];
     for  (i2 in nonzero.idx){
       if(i2>i1){
       tmp.sign = tmp.sign + sign(Xn[j,i1] - Xn[j,i2])*sign(Xn[k,i1] - Xn[k,i2])*mat.ker[i1,i2]
       }
     }
   }
   mat.KT[j,k] = tmp.sign/U.ker
   mat.Sample[j,k] = tmp.sample/sum(kern((Zn - z),h))
     }
   }
   #mark
   Sig.est = sin(pi*mat.KT/2)

   Sig.est = (Sig.est + t(Sig.est));
   mat.Sample = (mat.Sample +t(mat.Sample));
   diag(Sig.est) = 1;
   diag(mat.Sample) =1;
   Ken.est = mat.KT + t(mat.KT)
   diag(Ken.est) = 1
   
   return(list(Sigma = Sig.est,kendall = Ken.est, Sample.est = mat.Sample, U.ker = U.ker))
}

#=============================#
# tvKT.big
# Calculate the kernel smoohted Kendall's Tau 
# function for larger size of n,d
#=============================#
tvKT.big = function(Zn,Xn,z,h){
  n = dim(Xn)[2]
  d = dim(Xn)[1]
  nonzero.idx = which(abs(Zn-z)<h);
  if (length(nonzero.idx) == 0){
    cat('Warning: kernel is zero')
  }
  mat.KT = matrix(0,d,d)
  mat.ker = matrix(0,1,n)
  for (i1 in nonzero.idx){
        mat.ker[i1] = kern((Zn[i1] - z),h) * kern((Zn[i2] - z),h);
      }
  for (j in 1:(d-1)){
    for (k in (j+1):d ){
      tmp.sign = 0
      tmp.ker = 0
      for (i1 in nonzero.idx){
        for  (i2 in nonzero.idx){
          if(i2>i1){
            tmp.ker = tmp.ker + mat.ker[i1]*mat.ker[i2]
            tmp.sign = tmp.sign + sign(Xn[j,i1] - Xn[j,i2])*sign(Xn[k,i1] - Xn[k,i2])*mat.ker[i1]*mat.ker[i2]
          }
        }
      }
      mat.KT[j,k] = tmp.sign/tmp.ker
    }
  }
  #mark
  Sig.est = sin(pi*mat.KT/2)
  
  Sig.est = (Sig.est + t(Sig.est));
  diag(Sig.est) = 1;
  Ken.est = mat.KT + t(mat.KT)
  diag(Ken.est) = 1
  
  return(list(Sigma = Sig.est,kendall = Ken.est))
}

#=============================#
# caliClime
# Calibrated CLIME
# estimate the inverse covariance matrix
#=============================#
caliClime = function(Sig,lambda,gamma = 0.5){
  d = dim(Sig)[1]
  Omega = matrix(0,d,d)
  for(k in 1:d){
  e = matrix(0,d,1);
  e[k] = 1;
  obj = -rbind(matrix(1,2*d,1),gamma)
  mat1 = cbind(Sig, -Sig, -lambda*matrix(1,d,1))
  mat2 = cbind(-Sig, Sig, -lambda*matrix(1,d,1))
  mat3 = cbind(matrix(1,1,2*d), -1)
  mat = rbind(mat1,mat2,mat3)
  rhs = rbind(e,-e,0)
  output = fastlp(obj,mat,rhs)
  Omega[,k] = output[1:d]- output[(d+1):(2*d)]
  }
  Omega = (Omega+t(Omega))/2;
  return(Omega)
} 

#=============================#
# caliClime
# Calibrated CLIME
# estimate the path inverse covariance matrix along lambda
#=============================#
caliClime.path = function(Sig,lambda,gamma = 0.5){
  d = dim(Sig)[1]
  nl = length(lambda)
  Omega.path = vector('list',nl)
  adj.path = vector('list',nl)
  for(i in 1:nl){
  Omega = matrix(0,d,d)
  for(k in 1:d){
    e = matrix(0,d,1);
    e[k] = 1;
    obj = -rbind(matrix(1,2*d,1),gamma)
    mat1 = cbind(Sig, -Sig, -lambda[i]*matrix(1,d,1))
    mat2 = cbind(-Sig, Sig, -lambda[i]*matrix(1,d,1))
    mat3 = cbind(matrix(1,1,2*d), -1)
    mat = rbind(mat1,mat2,mat3)
    rhs = rbind(e,-e,0)
    output = fastlp(obj,mat,rhs)
    Omega[,k] = output[1:d] - output[(d+1):(2*d)]
  }
  Omega.path[[i]] = (Omega+t(Omega))/2;
  adj.ls = ifelse(abs(Omega.path[[i]])>0.0001,1,0)
  diag(adj.ls)=0
  adj.path[[i]] = adj.ls
  }
  return(list(Omega = Omega.path, adj = adj.path))
}

#=============================#
# ZWL
# The kernel graphical Lasso estimator: Precision estimator using GLasso
# Method proposed in Zhou, Wasserman, Lafferty,2010
#=============================#
ZWL = function(Sigma,lambda){
  d = dim(Sigma)[1]
  eig = eigen(Sigma) 
  Sigma = eig$vectors%*%diag(ifelse(eig$values>0.000001,eig$values,0.000001))%*%t(eig$vectors) 
  #out3 = huge(Sigma,lambda=lambda,method="glasso")
#   out3 = glassopath(Sigma,rholist=sort(lambda))
#   nlambda = length(lambda)
#   icov = vector('list',nlambda)
#   path = vector('list',nlambda)
#   for(i in 1:nlambda){
#     icov[[i]] = out3$wi[,,i]
#     adj.ls = ifelse(abs(icov[[i]])>0.0001,1,0)
#     diag(adj.ls)=0
#     path[[i]] = adj.ls    
#   }    
#   return(list(icov = icov, path = path))
  
  out3 = huge(Sigma,lambda=lambda,method="glasso")
  return(out3)
}

#=============================#
# Precision estimator using CLIME
#=============================#
clime.naive = function(X,lambda){
  fit.clime = clime(X, lambda = lambda)
  clime.icov = fit.clime$Omega
  adj.path = vector('list',length(lambda))
  for(k in 1:length(lambda)){
    adj.ls = ifelse(abs(clime.icov[[k]])>0.0001,1,0)
    diag(adj.ls)=0
    adj.path[[k]] = adj.ls
  }
  
  return(list(Omega = clime.icov, adj = adj.path))
}

#=============================#
# Glasso.KX 
# The kernel Pearson CLIME estimator using node-wise regression
# Method proposed in Kolar, Ahmed, Xing, 2010
#=============================#
Glasso.KX = function(Zn,Xn,z,h,lambda){
  n = dim(Xn)[2]
  d = dim(Xn)[1]
  nonzero.idx = which(abs(Zn-z)<h);
  if (length(nonzero.idx) == 0){
    cat('Warning: kernel is zero')
  }
  
  Omega = vector('list',length(lambda))
  tmp.Omega = array(0,dim= c(d,d,length(lambda)))

  weight.kern = kern(Zn[nonzero.idx]-z,h)/mean(kern(Zn[nonzero.idx]-z,h))
  for(j in 1:d){
    rst = glmnet(x = t(Xn[-j,nonzero.idx]), y = t(Xn[j,nonzero.idx]),  standardize = FALSE, weights = weight.kern, lambda = lambda, intercept=FALSE)
    beta = coef(rst)
    for(i in 1:length(lambda))
      tmp.Omega[-j,j,i] = beta[-1,i]
  }
  for(i in 1:length(lambda))
    Omega[[i]] = -(tmp.Omega[,,i] + t(tmp.Omega[,,i]))/2 + diag(d)

  adj.path = vector('list',length(lambda))
  for(k in 1:length(lambda)){
    adj.ls = ifelse(abs(Omega[[k]])>0.0001,1,0)
    diag(adj.ls)=0
    adj.path[[k]] = adj.ls
  }
  return(list(Omega = Omega, adj = adj.path))
}

############## Error Evaluation #############
errs=function(icov.ls,adj.ls,icov.true){
  adj.true=ifelse(abs(icov.true)>0.0001,1,0)
  diag(adj.true)=0
  FP=sapply(adj.ls,function(mat) sum((mat-adj.true)>0)/2)
  FN=sapply(adj.ls,function(mat) sum((mat-adj.true)<0)/2)
  
  d=dim(icov.true)[1]
  TotalNum=d*(d-1)/2
  EdgeNum=sum(adj.true)/2
  FPR=FP/(TotalNum-EdgeNum)
  FNR=FN/(EdgeNum)
  
  min.idx=which.min(FPR+FNR)
  eps1=max(colSums(abs(icov.true-icov.ls[[min.idx]])))
  eps2=max(abs(eigen(icov.true-icov.ls[[min.idx]])$values))
  epsF=sqrt(sum((icov.true-icov.ls[[min.idx]])^2))
  
  return(list(FPR=FPR,FNR=FNR,eps1=eps1,eps2=eps2,epsF=epsF,min.idx=min.idx))
}

errs.avg = function(nlambda,nz,Omega.kendall,adj.kendall,Omega.true){
  FNR.mat.kendall = matrix(0,nlambda,nz)
  FPR.mat.kendall = matrix(0,nlambda,nz)
  eps.mat = matrix(0,3,nz)
  for(i in 1:nz){
    kendall.errs = errs(Omega.kendall[[i]],adj.kendall[[i]],Omega.true[[i]])
    FNR.mat.kendall[,i] =  kendall.errs$FNR
    FPR.mat.kendall[,i] =  kendall.errs$FPR
    eps.mat[,i] =  c(kendall.errs$eps1,kendall.errs$eps2,kendall.errs$epsF)
  }
  FNR = apply(FNR.mat.kendall,1,min)
  FPR = apply(FPR.mat.kendall,1,min)
  eps = apply(eps.mat,1,min)
  return(list(FNR = FNR, FPR = FPR, eps = eps))
}

errs.avg2 = function(nlambda,nz,Omega.kendall,adj.kendall,Omega.true){
  FNR.mat.kendall = matrix(0,nlambda,nz)
  FPR.mat.kendall = matrix(0,nlambda,nz)
  eps.mat = matrix(0,3,nz)
  for(i in 1:nz){
    kendall.errs = errs(Omega.kendall,adj.kendall,Omega.true[[i]])
    FNR.mat.kendall[,i] =  kendall.errs$FNR
    FPR.mat.kendall[,i] =  kendall.errs$FPR
    eps.mat[,i] =  c(kendall.errs$eps1,kendall.errs$eps2,kendall.errs$epsF)
  }
  FNR = apply(FNR.mat.kendall,1,min)
  FPR = apply(FPR.mat.kendall,1,min)
  eps = apply(eps.mat,1,min)
  return(list(FNR = FNR, FPR = FPR, eps = eps))
}

marg.tranf = function(X,type){
  switch(type,
    quad = sign(X)*abs(X)^3,
    Phi = pnorm(X),
    Cauchy = qt(pnorm(X),df=20)
  )
}

contam = function(X,rho){
  d = dim(X)[1]
  n = dim(X)[2]
  knoise=round(rho*d*n)
  X[sample(d*n,knoise)]=sample(c(-1,1),1)*rnorm(knoise,3,3)
  return(X)
}

