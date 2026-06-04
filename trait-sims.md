---
title: "trait-sims"
output:
  html_document:
    keep_md: yes
---



##1. Simulate neutral trait data for 2 traits that are correlated (so I decide Va and C)


``` r
#140 individuals, 2 traits
myM = matrix(0, nrow=140, ncol=2)

##pull in one of our actual kinship matrices
load('data/Norwich_qpc_inputs.rda')

#make up some Vas and Cs. 
myV = matrix(c(1,.5,.5,1), nrow=2)

set.seed(3)

myZ = matrix_normal(M = myM, U = myK, V = myV)

par(mfrow=c(1,3))
hist(myZ[,1], col = mycol[1], border="white", breaks=10, main = "", xlab = "trait 1")
hist(myZ[,2], col = mycol[2], border="white", breaks=10, main= "", xlab = "trait 2")
plot(myZ[,1], myZ[,2], xlab = "trait 1", ylab = "trait 2", bty="n", col = mycol[3], lwd=2)
```

![](trait-sims_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

``` r
cor.test(myZ[,1], myZ[,2])
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  myZ[, 1] and myZ[, 2]
## t = 3.93, df = 138, p-value = 0.0001338
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.1597654 0.4590034
## sample estimates:
##       cor 
## 0.3172595
```

##2. Simulate trait data for 2 correlated traits where trait 1 is neutral and trait 2  is under selection.


``` r
myM = matrix(0, nrow=140, ncol=2)
myV = matrix(c(1,.8,.8,1), nrow=2)

#neutral traits
myZ = matrix_normal(myM, myK, myV)

#selected traits -- could add separate selection by adding different PCs into the mix. Yes, this is hacky. But right now what I'm doing in is adding the selection effect is the position on PC2*the variance explained by PC3 on trait 2. For trait 1, I'm multiplying this by the correlation btw traits
selectionEffect = matrix(c(myEig$vectors[,3]*myEig$values[3]*.8, myEig$vectors[,3]*myEig$values[3]), nrow=140) 
myZs = myZ + selectionEffect

par(mfrow=c(1,3))
hist(myZs[,1], col = mycol[1], border="white", breaks=10, main = "", xlab = "trait 1")
hist(myZs[,2], col = mycol[2], border="white", breaks=10, main= "", xlab = "trait 2")
plot(myZs[,1], myZs[,2], xlab = "trait 1", ylab = "trait 2", bty="n", col = mycol[3], lwd=2)
```

![](trait-sims_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

``` r
cor.test(myZ[,1], myZ[,2])
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  myZ[, 1] and myZ[, 2]
## t = 16.963, df = 138, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.7598380 0.8694286
## sample estimates:
##       cor 
## 0.8221068
```

``` r
##does it work to create differences detectable by qpc
myQpc1 = calcQpc(myZ = myZ[,1], myLambdas = myEig$values, myU = myEig$vectors, myM = 1:10, myL = 70:140)
myQpc1s = calcQpc(myZ = myZs[,1], myLambdas = myEig$values, myU = myEig$vectors, myM = 1:10, myL = 70:140)
myQpc2 = calcQpc(myZ = myZ[,2], myLambdas = myEig$values, myU = myEig$vectors, myM = 1:10, myL = 70:140)
myQpc2s = calcQpc(myZ = myZs[,2], myLambdas = myEig$values, myU = myEig$vectors, myM = 1:10, myL = 70:140)

myQpc1$pvals[3]
```

```
## [1] 0.9887872
```

``` r
myQpc1s$pvals[3]
```

```
## [1] 0.09973492
```

``` r
myQpc2$pvals[3]
```

```
## [1] 0.9509037
```

``` r
myQpc2s$pvals[3]
```

```
## [1] 0.03481998
```



