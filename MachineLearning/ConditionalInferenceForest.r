#!/usr/bin/Rscript
library(methods)
library(mlr)
trainData=read.table("../data/train.dat")
testData=read.table("../data/test.dat")

lrn = makeLearner("classif.cforest", predict.type = "prob")
rdesc = makeResampleDesc(method="CV", stratify=TRUE)
task=makeClassifTask(id="psapsc", data=trainData, target="PsAPsC")
r=resample(learner=lrn, task=task, resampling=rdesc, show.info=FALSE, measures=auc)
print(r$aggr)
pred=predict(train(lrn, task), newdata=testData)
print(performance(pred, measures=auc))
