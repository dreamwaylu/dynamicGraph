# dynamicGraph
This is a R implementation of the paper [Post-regularization Inference for Dynamic Nonparanormal Graphical Models](http://arxiv.org/abs/1512.08298) by Junwei Lu, Mladen Kolar and Han Liu.

The paper proposes a kernel smoothed Kendall’s tau correlation matrix based estimator for the dynamic nonparanormal graphical models, which allows us to model high dimensional heavy-tailed systems and the evolution of their latent network structures.

## Estimation of Graphical Models
We provide the R codes for four estimation methods for the graphical models.
* [The kernel smoothed Kendall's tau estimator](http://arxiv.org/abs/1512.08298)
 * A robust estimator for dynamic nonparanormal graphical models.
* [The kernel Pearson CLIME estimator](http://sites.stat.psu.edu/~rli/research/A20n120.pdf)
 * A graph estimator combining the kernel smoothed Pearson sample covariance with CLIME estimator.
* [The kernel graphical Lasso estimator](http://link.springer.com/article/10.1007/s10994-010-5180-0#/page-1)
  * A graph estimator combining the kernel smoothed Pearson sample covariance with  graphical Lasso maximum likelihood estimator.
* [The kernel neighborhood selection estimator](http://arxiv.org/abs/0812.5087)
  * A kernel smoothed nodewise regression method for graph estimation.
  
## Inference of Graphical Models
We provide the R codes for three hypothesis tests for the dynamic graphical models
* Edge presence test
 * H0: the edge (j,k) is zero at a fixed time
* Super-graph test
 * H0: the true graph is a subgraph of given graph G at a fixed time
* Uniform edge presence test
 * H0: the true graph is a subgraph of given graph G for all time period
 
## Usage
Estimate the correlation matrix:
```
out.Sigma = tvKT(Zn, Xn, z = z0, h = h)
```
The kernel smoothed Kendall's tau estimator:
```
out1 = caliClime.path(out.Sigma$Sigma, lambda = lambda)$Omega
```
The kernel Pearson CLIME estimator:
```
out2 = caliClime.path(out.Sigma$Sample.est, lambda = lambda)$Omega
```
The kernel graphical Lasso estimator:
```
out3 = ZWL(out.Sigma$Sample.est, lambda = lambda)$icov
```
The kernel neighborhood selection estimator:
```
out4 = Glasso.KX(Zn, Xn, z = z0, h = h,lambda = lambda$Omega)$Omega
```
Inference for the graphical model:

The edge set to test is ```edge.list``` and the time period to test is ```z.list```. We can set them to be single or multiple for three tests above. The covariance and inverse correlation estimators for multiple times are ```Sigma.list``` and ```Omega.list```.
The testing statistic
```
UK.vec = UKern(Zn,z.list,h)$vec
Test.super = max.Score(Sigma.list,Omega.list,UK.vec,edge.list,z.list,h)$max
```



