# functions for times-varying graph estimation

#=============================#
# UKern
# Generate the U-statistics for kernel, omega(z,z')
#=============================#
UKern = function(Zn,zlist,h){
  nz = length(zlist)
  n = length(Zn)
  K.Ustat = vector('list',nz)
  Ustat = matrix(0,nz,1)
  for(k in 1:nz){
    nonzero.idx = which(abs(Zn- zlist[k])<h);
  mat.ker = matrix(0,n,n)
  for (i1 in nonzero.idx){
    for  (i2 in nonzero.idx){
      if(i2>i1){
        mat.ker[i1,i2] = kern((Zn[i1] - zlist[k]),h) * kern((Zn[i2] - zlist[k]),h);
      }
    }
  }
  K.Ustat[[k]]  = mat.ker
  Ustat[k] = 2*sum(mat.ker)/(n*(n-1))
  }
  return(list(mat = K.Ustat, vec = Ustat))
}


#=============================#
# UKern.boot
# Generate the bootstrapped U-statistics for kernel
#=============================#
UKern.boot = function(Zn,zlist,h,xi){
  nz = length(zlist)
  n = length(Zn)
  K.Ustat = vector('list',nz)
  Ustat = matrix(0,nz,1)

  for(k in 1:nz){
    mat.ker = matrix(0,n,n)
    nonzero.idx = which(abs(Zn-zlist[k])<h);
    for (i1 in nonzero.idx){
      for  (i2 in nonzero.idx){
        if(i2>i1){
          mat.ker[i1,i2] = kern((Zn[i1] - zlist[k]),h) * kern((Zn[i2] - zlist[k]),h)*(xi[i1]+xi[i2]);
        }
      }
    }
    K.Ustat[[k]]  = mat.ker
    Ustat[k] = 2*sum(mat.ker)/(n*(n-1))
  }
  return(list(mat = K.Ustat, vec = Ustat))
}

#=============================#
# Score
# Score statistic Omega_j(Sigma*Omega_k-ek)
#=============================#
Score = function(Sigma,Omega,j,k){
  beta = Omega[,k]
  beta[j] = 0;
  e = matrix(0,d,1)
  e[k] = 1
  S = Omega[j,]%*%(Sigma%*%beta - e)
  return(S)
}

Var.est = function(Xn,Zn,Tau,Omega,z,j,k,K.Ustat= NULL){
  n = length(Zn);
  d = dim(Xn)[1]
  nonzero.idx = which(abs(Zn-z)<h);
  idx1 = 1:floor(n/2)
  idx2 = (floor(n/2)+1):n
  nonzero.idx = which(abs(Zn-z)<h);
  nonzero.idx1 = nonzero.idx[which(nonzero.idx<=max(idx1))]
  nonzero.idx2 = nonzero.idx[which(nonzero.idx>max(idx1))]
  Var.Theta = matrix(0,d,d)
  if(is.null(K.Ustat)){
  mat.ker = matrix(0,n,n)
  for (i1 in nonzero.idx){
    for  (i2 in nonzero.idx){
      if(i2>i1){
        mat.ker[i1,i2] = kern((Zn[i1] - z),h) * kern((Zn[i2] - z),h);
      }
    }
  }
  }else{
    mat.ker = K.Ustat
  }
  U.ker = 2*sum(mat.ker)/(n*(n-1))
  
  for(s in 1:d){
    for(t in s:d){
      tau1 = matrix(0,length(idx1),1)
      for(i1 in nonzero.idx1){
        tmp = 0
        for(i2 in nonzero.idx2){
          tmp = tmp + mat.ker[i1,i2]*sign(Xn[s,i1] - Xn[s,i2])*sign(Xn[t,i1] - Xn[t,i2])
        }
        tau1[i1] = sqrt(h)*(tmp/length(idx2) - Tau[s,t])
      }
      Var.Theta[s,t] = mean(tau1^2)*pi*cos((pi/2)*Tau[s,t])
    }
  }
  Var.Theta = (Var.Theta + t(Var.Theta) ) - diag(diag(Var.Theta))
  Var.est = Omega[j,]%*%Var.Theta%*%Omega[,k]/(U.ker)^2
  return(list(var=Var.est,Var.Theta = Var.Theta, U.ker = U.ker))
}


# Calculate Sign in the Kernel smoothed Kendall's tau
sign.comp = function(Xn){
  n = dim(Xn)[2]
  d = dim(Xn)[1]
  mat.sign = array(0,dim = c(d,d,n,n))
  for (j in 1:(d-1)){
    for (k in (j+1):d ){
      for (i1 in nonzero.idx){
        for  (i2 in nonzero.idx){
          if(i2>i1){
            mat.sign[j,k,i1,i2] =  sign(Xn[j,i1] - Xn[j,i2])*sign(Xn[k,i1] - Xn[k,i2])
          }
        }
      }
    }
  }
  return(mat.sign)
  
  
}

# Construct the d-by-d Kendall's tau matrix for all dimensions
tvKT.mat = function(Zn,Xn,z,h,K.Ustat = NULL,mat.sign){
  n = dim(Xn)[2]
  d = dim(Xn)[1]
  nonzero.idx = which(abs(Zn-z)<h);
  if (length(nonzero.idx) == 0){
    cat('Warning: kernel is zero')
  }
  mat.KT = matrix(0,d,d)
  #  mat.Sample = matrix(0,d,d)
  if(is.null(K.Ustat)){
    mat.ker = matrix(0,n,n)
    for (i1 in nonzero.idx){
      for  (i2 in nonzero.idx){
        if(i2>i1){
          mat.ker[i1,i2] = kern((Zn[i1] - z),h) * kern((Zn[i2] - z),h);
        }
      }
    }
  }else{
    mat.ker = K.Ustat
  }
  U.ker = sum(mat.ker)
  for (j in 1:(d-1)){
    for (k in (j+1):d ){
      tmp = 0
      for (i1 in nonzero.idx){
        #        tmp.sample = tmp.sample + kern((Zn[i1] - z),h) * Xn[j,i1]*Xn[k,i1];
        for  (i2 in nonzero.idx){
          if(i2>i1){
            tmp = tmp +  sum(mat.sign[j,k,,]*mat.ker)
          }
        }
      }
      mat.KT[j,k] = tmp/U.ker
      #mat.Sample[j,k] = tmp.sample/sum(kern((Zn - z),h))
    }
  }
  Sig.est = sin(pi*mat.KT/2)
  
  Sig.est = (Sig.est + t(Sig.est));
  diag(Sig.est) = 1;

  #mat.Sample = (mat.Sample +t(mat.Sample))/2;
  
  return(list(Sigma = Sig.est,kendall = mat.KT))
  
  
}

Xi.mat = function(xi){
  n = length(xi)
  Xi = matrix(0,n,n)
  for (i1 in 1:n){
    #        tmp.sample = tmp.sample + kern((Zn[i1] - z),h) * Xn[j,i1]*Xn[k,i1];
    for  (i2 in 1:n){
      if(i2>i1){
        Xi[i1,i2] =  xi[i1] + xi[i2]
      }
    }
  }
  return(Xi)
}


# Fast calculation of  the d-by-d bootstrap Kendall's tau matrix for all dimensions
tvKT.boot.list.fast =  function(Xi,sign.mat,K.Ustat){
  d = dim(sign.mat)[1]
  n = dim(sign.mat)[3]  
  Sigma.list = vector('list', length(K.Ustat))
  U.list = matrix(0,length(K.Ustat),1)
  for(i in 1:length(K.Ustat)){
    mat.KT = matrix(0,d,d)
    mat.ker = K.Ustat[[i]]
    tmp2 = mat.ker*Xi
    U.ker = sum(tmp2)
    for (j in 1:(d-1)){
      for (k in (j+1):d ){
        tmp1 = sign.mat[j,k,,]*mat.ker*tmp2
        mat.KT[j,k] = sum(tmp1)/U.ker
      }
    }
    Sig.est = sin(pi*mat.KT/2)
    
    Sig.est = (Sig.est + t(Sig.est));
    diag(Sig.est) = 1;
    U.list[i] = 2*U.ker/(n*(n-1))
    Sigma.list[[i]] = Sig.est
  }
  return(list(Sigma = Sigma.list, ker = U.list))
}

#Calculate the kernel smoothed Kendall's tau, another method
tvKT2 = function(Zn,Xn,z,h,K.Ustat = NULL){
  n = dim(Xn)[2]
  d = dim(Xn)[1]
  nonzero.idx = which(abs(Zn-z)<h);
  if (length(nonzero.idx) == 0){
    cat('Warning: kernel is zero')
  }
  mat.KT = matrix(0,d,d)
  #  mat.Sample = matrix(0,d,d)
  if(is.null(K.Ustat)){
    mat.ker = matrix(0,n,n)
    for (i1 in nonzero.idx){
      for  (i2 in nonzero.idx){
        if(i2>i1){
          mat.ker[i1,i2] = kern((Zn[i1] - z),h) * kern((Zn[i2] - z),h);
        }
      }
    }
  }else{
    mat.ker = K.Ustat
  }
  U.ker = sum(mat.ker)
  for (j in 1:(d-1)){
    for (k in (j+1):d ){
      tmp.sign = 0
      tmp.sample = 0
      for (i2 in nonzero.idx){
        #        tmp.sample = tmp.sample + kern((Zn[i1] - z),h) * Xn[j,i1]*Xn[k,i1];
        for  (i1 in nonzero.idx){
          if(i2>i1){
            tmp.sign = tmp.sign + sign(Xn[j,i1] - Xn[j,i2])*sign(Xn[k,i1] - Xn[k,i2])*mat.ker[i1,i2]
          }
        }
      }
      mat.KT[j,k] = tmp.sign/U.ker
      #mat.Sample[j,k] = tmp.sample/sum(kern((Zn - z),h))
    }
  }
  Sig.est = sin(pi*mat.KT/2)
  
  Sig.est = (Sig.est + t(Sig.est));
  diag(Sig.est) = 1;
  #mat.Sample = (mat.Sample +t(mat.Sample))/2;
  
  return(list(Sigma = Sig.est,kendall = mat.KT, U.ker.boot = 2*U.ker/(n*(n-1))))
}

#Calculate the kernel smoothed Kendall's tau
tvKT.boot = function(Zn,Xn,z,h,xi,K.Ustat = NULL){
  n = dim(Xn)[2]
  d = dim(Xn)[1]
  nonzero.idx = which(abs(Zn-z)<h);
  if (length(nonzero.idx) == 0){
    cat('Warning: kernel is zero')
  }
  mat.KT = matrix(0,d,d)
#  mat.Sample = matrix(0,d,d)
  if(is.null(K.Ustat)){
  mat.ker = matrix(0,n,n)
  for (i1 in nonzero.idx){
    for  (i2 in nonzero.idx){
      if(i2>i1){
        mat.ker[i1,i2] = kern((Zn[i1] - z),h) * kern((Zn[i2] - z),h)*(xi[i1]+xi[i2]);
      }
    }
  }
  }else{
    mat.ker = K.Ustat
  }
  U.ker = sum(mat.ker)
  for (j in 1:(d-1)){
    for (k in (j+1):d ){
      tmp.sign = 0
      tmp.sample = 0
      for (i1 in nonzero.idx){
#        tmp.sample = tmp.sample + kern((Zn[i1] - z),h) * Xn[j,i1]*Xn[k,i1];
        for  (i2 in nonzero.idx){
          if(i2>i1){
            tmp.sign = tmp.sign + sign(Xn[j,i1] - Xn[j,i2])*sign(Xn[k,i1] - Xn[k,i2])*mat.ker[i1,i2]
          }
        }
      }
      mat.KT[j,k] = tmp.sign/U.ker
      #mat.Sample[j,k] = tmp.sample/sum(kern((Zn - z),h))
    }
  }
  Sig.est = sin(pi*mat.KT/2)

  Sig.est = (Sig.est + t(Sig.est));
  diag(Sig.est) = 1;
  diag(mat.Sample) =1;
  #mat.Sample = (mat.Sample +t(mat.Sample))/2;

  return(list(Sigma = Sig.est,kendall = mat.KT, U.ker.boot = 2*U.ker/(n*(n-1))))
}


#Calculate the bootstrapped kernel smoothed Kendall's tau, another method
tvKT2.boot = function(Xn,xi){
  n = dim(Xn)[2]
  d = dim(Xn)[1]
  for (j in 1:(d-1)){
    for (k in (j+1):d ){
      tmp.sign = 0
      for (i1 in 1:(n-1)){
        #        tmp.sample = tmp.sample + kern((Zn[i1] - z),h) * Xn[j,i1]*Xn[k,i1];
        for  (i2 in i1:n){
          if(i2>i1){
            tmp.sign = tmp.sign + sign(Xn[j,i1] - Xn[j,i2])*sign(Xn[k,i1] - Xn[k,i2])*(xi[i1]+xi[i2])
          }
        }
      }
      mat.KT[j,k] = 2*tmp.sign/(n*(n-1))
      #mat.Sample[j,k] = tmp.sample/sum(kern((Zn - z),h))
    }
  }
  Sig.est = sin(pi*mat.KT/2)
  
  Sig.est = (Sig.est + t(Sig.est));
  diag(Sig.est) = 1;
  #mat.Sample = (mat.Sample +t(mat.Sample))/2;
  
  return(Sig.est)
}


tvKT.boot.list = function(Zn,Xn,zlist,h,xi,K.Ustat){
  nz = length(zlist)
  Sigma.list = vector('list',nz)
  for(i in 1:nz){
    Sigma.list[[i]]=tvKT.boot(Zn,Xn,zlist[i],h,xi,K.Ustat[[i]])$Sigma
  }
  return(Sigma.list)
}

#Calculate the bootstrapped score test statistics
Score.boot = function(Zn,Xn,h,xi,Omega,K.Ustat,j,k){
  n = length(Zn)
  rst = tvKT.boot(Zn,Xn,z,h,xi,K.Ustat)
  S.boot = sqrt(n*h)*K.Ustat*Score(rst$Sigma,Omega,j,k)
  return(S.boot)
}

#Calculate the maximum bootstrapped score test statistics
max.Score.boot = function(Sigma.list,Omega.list,K.Ustat.list,edge.list,z.list,h){
  num.edge = length(edge.list)
  nz = length(z.list)
  d= dim(Xn)[1]
  n= dim(Xn)[2]
  Test.super = matrix(0,nz,num.edge)
  for(i in 1:nz){
    for(j in  1:num.edge){
      j.idx = edge.list[j]
      j0 = floor((j.idx-1)/d)+1; k0 = j.idx%%d
      if(k0 == 0) k0 = d
      Test.super[i,j] = sqrt(n*h)*K.Ustat.list[i]*Score(Sigma.list[[i]],Omega.list[[i]],k0,j0)
    }
  }
  return(list(max = max(Test.super),stat = Test.super))
}

#Calculate the maximum score test statistics
max.Score = function(Sigma.list,Omega.list,K.Ustat.list,edge.list,z.list,h){
  num.edge = length(edge.list)
  nz = length(z.list)
  d= dim(Xn)[1]
  n= dim(Xn)[2]
  Test.super = matrix(0,nz,num.edge)
  for(i in 1:nz){
    for(j in 1:num.edge){
      j.idx = edge.list[j]
      j0 = floor((j.idx-1)/d)+1; k0 = j.idx%%d
      if(k0 == 0) k0 = d
      Uker = K.Ustat.list[i]
      Test.super[i,j] = sqrt(n*h)*Uker*Score(Sigma.list[[i]],Omega.list[[i]],k0,j0)
    }
  }
  return(list(max = max(Test.super),stat = Test.super))
}
