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

#make up some Vas and Cs. Not sure how to do this super well.
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
