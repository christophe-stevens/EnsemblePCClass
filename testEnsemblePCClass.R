rm(list=ls())

####################################
#
# chris-create package
#
####################################
#install.packages("roxygen2")
#devtools::create("EnsembleBaseClass")
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
dataset$diabetes <- as.factor(ifelse(dataset$diabetes=="pos",1,0))
myformula <- diabetes ~ insulin + age + triceps
perc.train <- 0.7
index.train <- sample(1:nrow(dataset), size = round(perc.train*nrow(dataset)))
data.train <- dataset[index.train,]
data.predict <- dataset[-index.train,]

devtools::install("EnsemblePCClass")
# chris-change: TESTINGS
sources <-paste("./R/",c("epcclass.R","integrator.R"),sep="")
for (src in sources){
  source(src)
}

control <- epcclass.baselearner.control(baselearners=c("penreg","svm","knn","nnet","gbm","rf"))
ensemble <- epcclass(formula =  myformula, data = data.train, baselearner.control= control)
plot(ensemble)
pred <- round(predict(ensemble,newdata=data.predict))
table(pred, data.predict$diabetes)
