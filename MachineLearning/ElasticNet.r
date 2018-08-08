#!/usr/bin/Rscript
library(glmnet)
library(ROCR)
trainData=read.table("../data/train.dat")
testData=read.table("../data/test.dat")

fit <- cv.glmnet(as.matrix(trainData[-ncol(data)]), trainPhenotype[ncol(data)], 
  nfolds=10, alpha=1, keep=TRUE, foldid=groupid, parallel=FALSE, type.measure="auc", family="binomial")
print(max(fit$cvm))
pred=predict(fit,as.matrix(testData))
print(unlist(performance(prediction(pred,testPhenotype),measure="auc")@y.values))
