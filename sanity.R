library(scCompReg)
library(tictoc)
library(R.matlab)


file.mat <- readMat( '/Users/Sophia/Desktop/BioStats/RAd4_sc_data.mat')
X <- Matrix(file.mat[['X']], sparse=TRUE)
A <- Matrix(file.mat[['D']], sparse=TRUE)
peakO <- Matrix(log10(file.mat[['PeakO']] + 1), sparse=TRUE)

tic('cnmf')
output = cnmf(peakO,
         X,
         A,
         k=3,
         alpha=0.5,
         beta_max_scale=5,
         verbose=T)
toc()

rm(list=ls())

t2 = as.numeric(read.table('/Users/Sophia/Desktop/BioStats/scRNA-result.txt', header=F))
t1 = as.numeric(read.table('/Users/Sophia/Desktop/BioStats/scATAC-result.txt', header=F))

par(mfrow=c(1,2))
hist(t1)
hist(output$cluster1)
c1 = output$cluster1
c2 = output$cluster2

library(mclust)
adjustedRandIndex(t1, c1)
adjustedRandIndex(t2, c2)

hist(t2)
hist(c2)




