rm(list=ls())

####################################
#
# chris-create package
#
####################################
#install.packages("roxygen2")
setwd("C:/Users/chris/OneDrive - Imperial College London/PHD/PHD Ensemble Package/New Package/EnsemblePCClass/")
#devtools::install("../EnsembleBaseClass")
require("EnsembleBaseClass")

####################################
#
# chris-test : 
#
####################################
# Required package for benchmarking dataset
require(mlbench)  # install.packages("mlbench")  

# require ML BENCH
data(PimaIndiansDiabetes)
# Data preparation
dataset <- PimaIndiansDiabetes
dataset$diabetes <- as.numeric(ifelse(dataset$diabetes=="pos",1,0))
myformula <- diabetes ~ insulin + age + triceps
perc.train <- 0.7
index.train <- sample(1:nrow(dataset), size = round(perc.train*nrow(dataset)))
data.train <- dataset[index.train,]
# balance dataset 
num.negative <- sum(data.train$diabetes==0)
num.positive <- sum(data.train$diabetes==1)
remove.n.negative <- num.negative-num.positive
data.train <- data.train[-which(data.train$diabetes==0)[1:remove.n.negative],]
data.predict <- dataset[-index.train,]

#devtools::install("EnsemblePCClass")
# chris-change: TESTINGS
sources <-paste("./R/",c("epcclass.R","integrator.R"),sep="")
for (src in sources){
  source(src)
}

control <- epcclass.baselearner.control(baselearners=c("svm","xgboost","rf",
                                                       "penreg","nnet","knn",
                                                       "bart"))
ensemble <- epcclass(formula =  myformula, data = data.train, 
                     baselearner.control= control, ncores=8)
plot(ensemble)
newpred <- predict(ensemble,newdata=data.predict)
mean((newpred-data.predict$diabetes)^2)
table(round(newpred), data.predict$diabetes)

batchPred <- predict(ensemble$est$baselearner.cv.batch, newdata= data.predict)
mean((batchPred[,2]-data.predict$diabetes)^2)                 
table(round(batchPred[,2]), data.predict$diabetes)
