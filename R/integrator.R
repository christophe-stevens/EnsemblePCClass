# defining base classes and generic methods to support sweep operation; TODO: move these base classes to EnsembleBase?

# base class for sweep configurations
setClass("Classification.Sweep.Config", slots=c(n="OptionalNumeric"), contains="VIRTUAL")
# base class for output of sweep training
setClass("Classification.Sweep.FitObj", slots = c(config="Classification.Sweep.Config", est="ANY", pred="matrix"), contains = "VIRTUAL")
# generic method for training sweep operations
setGeneric("Classification.Sweep.Fit", function(object, X, y, print.level=1) standardGeneric("Classification.Sweep.Fit"))
# class for output of sweep cv training
Classification.Sweep.CV.FitObj <- setClass("Classification.Sweep.CV.FitObj"
  , slots = c(sweep.list="list", pred="matrix", partition="OptionalInteger")
)

# sweep cv training function
Classification.Sweep.CV.Fit <- function(config, X, y, partition, print.level=1) {
  nfolds <- max(partition)
  if (length(partition)!=nrow(X)) stop("length of fold parameter does not match number of rows in data")
  pred <- array(NA, dim=dim(X))
  sweep.list <- list()
  n.min <- +Inf
  for (i in 1:nfolds) {
    if (print.level>=1) cat("processing fold", i, "of", nfolds, "\n")
    index_predict <- which(partition==i)
    X_train <- X[-index_predict,]
    X_predict <- X[index_predict,]
    y_train <- y[-index_predict]
    y_predict <- y[index_predict]
    regtmp <- Classification.Sweep.Fit(object=config, X=X_train, y=y_train, print.level=print.level-1)
    pred[index_predict,1:regtmp@config@n] <- predict(regtmp, Xnew=X_predict)
    sweep.list[[i]] <- regtmp
    n.min <- min(n.min, regtmp@config@n) # TODO: loop-carried, dependency, cannot be parallelized in current form
  }
  # this is needed specifically for pcr sweep since we are allowing for n in each fold to be determined using the try-error method
  # alternatively, we can keep this function general, and expect the PCR sweep to define its own CV method
  for (i in 1:nfolds) {
    sweep.list[[i]]@config@n <- n.min
  }
  pred <- pred[,1:n.min,drop=F]
  
  ret <- Classification.Sweep.CV.FitObj(sweep.list=sweep.list, pred=pred, partition=partition)
  return (ret)
}
# sweep cv predict
predict.Classification.Sweep.CV.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  nfolds <- length(object@sweep.list)
  ret <- 0
  for (n in 1:nfolds) {
    ret <- ret + predict(object@sweep.list[[n]], Xnew=Xnew)
  }
  ret <- ret/nfolds
  return (ret)
}

### PCR sweep methods and classes ###
Classification.Sweep.PCR.Config <- setClass("Classification.Sweep.PCR.Config", contains = "Classification.Sweep.Config")
Classification.Sweep.PCR.FitObj <- setClass("Classification.Sweep.PCR.FitObj", contains = "Classification.Sweep.FitObj")
# training
setMethod("Classification.Sweep.Fit", "Classification.Sweep.PCR.Config",
  function(object, X, y) {
    maxpc <- min(dim(X))
    X.trsf <- log(X / ( 1 - X))
    pcr <- prcomp(X.trsf)  # avoid having probability but use count using loggOdd here returns error due to +/-Inf
    Xpc <- predict(pcr, X.trsf)
    predmat <- array(NA, dim=c(nrow(X.trsf),maxpc))
    beta.list <- list()
    pca.reg.data <- cbind(data.frame(y=y), as.data.frame(Xpc))
    for (i in 1:maxpc) {
      # add a try wrapper, trap error and exit loop; this allows us to go to higher pc's and not arbitrarily choose a cutoff
      fml <- formula(paste("y ~ ",paste("PC",1:i,collapse="+",sep=""),sep=""))
      glmmdl <- try(glm( fml, data  = pca.reg.data ,family=binomial(link = "logit"),
                                               singular.ok = F))
      if (inherits(glmmdl, "try-error")) {
        cat("singularity encountered at i=", i, "\n")
        i <- i-1
        break
      }
      predmat[,i] <- as.vector(round(glmmdl$fitted))
      beta.list[[i]] <- as.vector(glmmdl$coeff)
    }
    object@n <- i
    est <- list(pcr=pcr, beta.list=beta.list, dimX=dim(X.trsf))
    ret <- Classification.Sweep.PCR.FitObj(config=object, est=est, pred=predmat)
    return (ret)
  }
)
# prediction
predict.Classification.Sweep.PCR.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  maxpc <- object@config@n
  Xnew.trsf <- log(Xnew / ( 1 - Xnew))
  Xpc.new <- cbind(1,predict(object@est$pcr, Xnew.trsf)[,1:maxpc,drop=F])
  newpred <- array(NA, dim=c(nrow(Xnew.trsf),maxpc))
  for (i in 1:maxpc) {
    xb <- Xpc.new[,1:(i+1)]%*%object@est$beta.list[[i]]
    newpred[,i] <- exp(xb)/(exp(xb)+1)
  }
  return (newpred)
}

### PCR integrator methods and classes ###
## Select.MinErr operation, TODO: this must be exported in EnsembleCV and imported in this package, or moved to EnsembleBase
Classification.Select.MinErr.Config <- setClass("Classification.Select.MinErr.Config", contains = "Classification.Select.Config")
Classification.Select.MinErr.FitObj <- setClass("Classification.Select.MinErr.FitObj", contains = "Classification.Select.FitObj")
setMethod("Classification.Select.Fit", "Classification.Select.MinErr.Config",
  function(object, X, y, print.level=1) {
    error <- apply(X, 2, FUN= function(x) {object@errfun(x, y)})
    error.min <- min(error)
    index.min <- which(error==error.min)[1] # in case we find multiple optimal indexes
    est <- list(index.min=index.min, error.min=error.min, error=error)
    ret <- Classification.Select.MinErr.FitObj(config=object, est=est, pred=X[,index.min])
    return (ret)
  }
)
predict.Classification.Select.MinErr.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  return (Xnew[,object@est$index.min])
}
## end of Select.MinErr operation

Classification.Integrator.PCR.SelMin.Config <- setClass("Classification.Integrator.PCR.SelMin.Config", slots = c(partition="integer"), contains = "Classification.Integrator.Config")
Classification.Integrator.PCR.SelMin.FitObj <- setClass("Classification.Integrator.PCR.SelMin.FitObj", contains = "Classification.Integrator.FitObj")
setMethod("Classification.Integrator.Fit", "Classification.Integrator.PCR.SelMin.Config",
  function(object, X, y, print.level=1) {
    # 1) cv sweep - pcr method
    my.pcr.config <- Classification.Sweep.PCR.Config(n=NULL)
    est.pcr <- Classification.Sweep.CV.Fit(my.pcr.config, X=X, y=y, partition=object@partition, print.level=print.level)
    # 2) select - min method
    my.select.config <- Classification.Select.MinErr.Config(errfun=object@errfun)
    est.select <- Classification.Select.Fit(my.select.config, X=est.pcr@pred, y=y)
    
    est <- list(pcr=est.pcr, select=est.select)
    ret <- Classification.Integrator.PCR.SelMin.FitObj(config=object, est=est, pred=est.select@pred)
    return (ret)
  }
)
predict.Classification.Integrator.PCR.SelMin.FitObj <- function(object, Xnew=NULL, ...) {
  if (is.null(Xnew)) return (object@pred)
  newpred.pcr <- predict(object@est$pcr, Xnew)
  newpred <- predict(object@est$select, newpred.pcr)
  attr(newpred, "newpred.pcr") <- newpred.pcr
  return (newpred)
}




