#==================#
# inference_example.R -- inferential analysis for the time-varying graphical model
#==================#

library("huge")
library("corpcor")
library("mvtnorm")
library("JGL")
source('funs_est.R')
source('funs_test.R')


n = 500;
d = 20;
num.data = 1;
num.test = 10;
num.band = 4;
h = 0.3/n^(1/5); #bandwidth
lambda = 0.2*(h^2 + sqrt(log(d/h)/(n*h)))


################# Super Graph Test #####################
omega.z0 = c(0,0.4,0.6,0.9)
z.alpha = qnorm(1-0.05/2)
alpha.level = 0.95
n.monte = 500
num.boot = 500
num.neighbor = 4;
test.stat = matrix(0,length(omega.z0),n.monte)
prob.stat = matrix(0,length(omega.z0),n.monte)
quant.stat = matrix(0,length(omega.z0),n.monte)
########## Super Graph ###########
Super=matrix(0,d,d)
for (i in 1:d){
  for (j in 1:d){
    if(j < i | (abs(j-i)<=num.neighbor | abs(j-i)>= d- num.neighbor))
      Super[j,i] =1;
  }
}
Super.null = which(Super == 0)
edge.list = Super.null



######## Super z #############
z.list =  seq(from = 0, to = 1, length.out = 10)
  





for(iz in 1:length(omega.z0)){
  ######## Graph Generation #####
  graph.true = graph.gen(d,num.itv = 3, num.band = 4,L = omega.z0[iz], U = min(0.9,2*omega.z0[iz]))
  
  for (monte in 1:n.monte){
    ptm=proc.time()   
    cat(monte,'in',n.monte,iz,'in',length(omega.z0),'\n')
    
    edge.list.now = edge.list;
    
    e.signal.z0 = omega.z0[iz]
    # generate data
    Xn = matrix(0,d,n*num.data)
    Zn = matrix(0,1,n*num.data)
    data = data.gen(n = n,d = d, graph = graph.true, num.itv =3)
    Xn = data$X
    Zn = data$Z
    
    # estimation
    Sigma.list = vector("list", length(z.list))
    Omega.list = vector("list", length(z.list))
    Sign.list = vector("list", length(z.list))


    K.Ustat.list = UKern(Zn,z.list,h)
    UK.vec = K.Ustat.list$vec
    UK.mat = K.Ustat.list$mat
    for (i in 1:length(z.list)){
      out.tmp = tvKT2(Zn,Xn,z=z.list[i], h = h, K.Ustat = UK.mat[[i]])
      Sigma.list[[i]] = out.tmp$Sigma;
      Omega.list[[i]] = caliClime(out.tmp$Sigma,lambda = lambda)
    }

    Test.super.mat = max.Score(Sigma.list,Omega.list,UK.vec,edge.list,z.list,h)
    Test.super = Test.super.mat$max
    
    cat('1st time',proc.time()-ptm,'\n')
    # bootstrap estimate quantile by monte carol
    super.stat.boot = matrix(0,num.boot,1)
    
    ptm = proc.time()
    for(iboot in 1:num.boot){
      cat('boot',iboot,'in',num.boot,'\n')
      xi = rnorm(n)
      K.boot.list = UKern.boot(Zn,z.list,h,xi)
      UKboot.vec = K.boot.list$vec
      UKboot.mat = K.boot.list$mat
      
      
      Sigma.boot.list = vector("list", length(z.list))
      for(izz in 1:length(z.list)){
        out.tmp = tvKT2(Zn,Xn,z=z.list[izz], h = h, K.Ustat = UKboot.mat[[izz]])
        Sigma.boot.list[[izz]] = out.tmp$Sigma
      }
      
      boot.super.mat = max.Score(Sigma.boot.list,Omega.list,UKboot.vec,edge.list,z.list,h)
      super.stat.boot[iboot] = boot.super.mat$max
    }
    super.qunatile = quantile(super.stat.boot,prob = alpha.level)
    cat('10boots time',proc.time()-ptm,'\n')
    test.stat[iz,monte] = Test.super
    quant.stat[iz,monte] = super.qunatile
    prob.stat[iz,monte] = ifelse(Test.super < super.qunatile, 1,0)
  }
}

save(test.stat,prob.stat,quant.stat,file="superT01N200.RData")
