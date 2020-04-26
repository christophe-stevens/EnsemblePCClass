epcclass.baselearner.control <- function(baselearners=c("nnet","penreg","rf","svm","gbm","knn")
  , baselearner.configs=make.configs(baselearners, type="classification"), npart=1, nfold=5) {
    return (list(configs=baselearner.configs, npart=npart, nfold=nfold))
}

epcclass.integrator.control <- function(errfun=rmse.error, nfold=5, method=c("default")) {
  return (list(errfun=rmse.error, nfold=nfold, method=method))
}

epcclass.set.filemethod <- function(formula, data, instance.list, type="classification") FALSE # TODO: to be implemented

epcclass <- function(formula, data
  , baselearner.control=epcclass.baselearner.control()
  , integrator.control=epcclass.integrator.control()
  , ncores=1, filemethod=FALSE, print.level=1
  , preschedule = TRUE
  , schedule.method = c("random", "as.is", "task.length"), task.length) {
  if (integrator.control$method!="default") stop("invalid PCR integration method")
  ncores.max <- try(detectCores(),silent=T)
  mycall <- match.call()
  if (!inherits(ncores.max,"try-error")) ncores <- min(ncores,ncores.max)
  if (print.level>=1 && ncores>1) cat("running in parallel mode, using", ncores, "cores\n")

  # training a batch of baselearners
  partitions.bl <- generate.partitions(baselearner.control$npart, nrow(data), baselearner.control$nfold)
  my.instance.list <- make.instances(baselearner.control$configs, partitions.bl)

  # TODO: determining of base learner estimation objects must be saved to disk or not
  if (missing(filemethod)) filemethod <- epcclass.set.filemethod(formula, data, my.instance.list)
  
  if (print.level>=1) cat("CV training of base learners...\n")
  est.baselearner.cv.batch <- Classification.CV.Batch.Fit(my.instance.list, formula, data, ncores=ncores, filemethod=filemethod, print.level=print.level, preschedule = preschedule, schedule.method = schedule.method, task.length = task.length)
  if (print.level>=1) cat("finished CV training of base learners\n")
  Xcv <- est.baselearner.cv.batch@pred
  y <- as.character(data[,all.vars(formula)[1]]) # TODO: more robust way of extracting y (here and inside base learner functions)
  
  partition.int <- generate.partition(nrow(data), integrator.control$nfold)
  my.integrator.config <- Classification.Integrator.PCR.SelMin.Config(errfun=integrator.control$errfun, partition=partition.int)
  est.integrator <- Classification.Integrator.Fit(my.integrator.config, X=Xcv, y=y, print.level=print.level)
  pred <- est.integrator@pred
  
  ret <- list(call=mycall, formula=formula, instance.list=my.instance.list, integrator.config=my.integrator.config, method=integrator.control$method
    , est=list(baselearner.cv.batch=est.baselearner.cv.batch, integrator=est.integrator)
    , y=y, pred=pred, filemethod=filemethod)
  class(ret) <- "epcclass"
  if (filemethod) class(ret) <- c(class(ret), "epcclass.file")
  return (ret)
}

predict.epcclass <- function(object, newdata=NULL, ncores=1, preschedule = TRUE, ...) {
  if (is.null(newdata)) return (object$pred)
  if (object$method=="default") {
    newpred.baselearner.cv.batch <- predict(object$est$baselearner.cv.batch, newdata, ncores=ncores, preschedule = preschedule, ...)
    newpred <- predict(object$est$integrator, Xnew=newpred.baselearner.cv.batch, config.list=object$est$baselearner.batch@config.list, ...)
  } else {
    stop("invalid PCR integration method")
  }
  return (as.numeric(newpred))
}

plot.epcclass <- function(x, ...) {
  errfun <- x$integrator.config@errfun 
  error <- errfun(as.numeric(x$pred), as.numeric(x$y))
  oldpar <- par(mfrow=c(1,2))
  plot(x$est$baselearner.cv.batch, errfun=errfun, ylim.adj = error)
  abline(h=error, lty=2)
  pcr.errors <- x$est$integrator@est$select@est$error
  #pcr.pred <- x$est$integrator@est$pcr@pred
  #pcr.errors <- apply(pcr.pred, 2, errfun, x$y)
  plot(pcr.errors, type="l", xlab="Number of Principal Components", ylab="CV Error", main="Integrator Performance")
  par(oldpar)
}

epcclass.save <- function (obj, file) {
  if (!("epcclass" %in% class(obj))) 
    stop("invalid object class (must be epcclass & epcclass.file)")
  if (missing(file)) 
    stop("must provide file argument")
  tmpfiles <- obj$est$baselearner.cv.batch@tmpfiles
  
  if (is.null(tmpfiles)) { # ordinary save
    save(obj, file = file)
  } else {
    tmpfile.new <- tempfile()
    save(obj, file = tmpfile.new, compress = F)
    all.files <- c(tmpfile.new, tmpfiles)
    all.files.basename <- basename(all.files)
    tmpdir <- paste0("./.", basename(tempfile("dir")), "/")
    dir.create(tmpdir)
    all.files.new <- paste0(tmpdir, all.files.basename)
    file.copy(all.files, all.files.new)
    meta <- list(filename.mainobj = all.files.basename[1], filenames.batchobj = all.files.basename[1 + 1:length(tmpfiles)])
    save(meta, file = paste0(tmpdir, "meta"), compress = FALSE)
    tar(file, files = tmpdir, compression = "gzip")
    unlink(tmpdir, recursive = TRUE)
  }
}

epcclass.load <- function (file) {
  env <- new.env()
  #load(file, envir = env)
  loadret <- suppressWarnings(try(load(file, envir = env), silent = TRUE))
  
  if (class(loadret) == "try-error") { # filemethod load
    filepaths <- untar(file, list = T)
    basenames <- basename(filepaths)
    dirnames <- dirname(filepaths)
    if (length(unique(dirnames)) > 1) 
      stop("unexpected multiple directories in tar filepaths")
    metafile.index <- which(basenames == "meta")
    extdir <- dirnames[1]
    untar(file)
    meta <- NULL
    load(filepaths[metafile.index])
    mainfile.index <- which(basenames == meta$filename.mainobj)
    load(filepaths[mainfile.index])
    if (!identical(class(obj), c("epcclass", "epcclass.file"))) 
      stop("invalid object class (must be epcclass & epcclass.file)")
    basenames.ordered <- basename(obj$est$baselearner.cv.batch@tmpfiles)
    filepaths.ordered <- paste(extdir, basenames.ordered, sep = "/")
    tmpfiles.new <- tempfile(rep("file", length(filepaths.ordered)))
    file.copy(from = filepaths.ordered, to = tmpfiles.new)
    unlink(filepaths.ordered, recursive = TRUE)
    unlink(extdir, recursive = TRUE)
    obj$est$baselearner.cv.batch@tmpfiles <- tmpfiles.new
    n.instance <- length(obj$est$baselearner.cv.batch@instance.list@instances)
    for (i in 1:n.instance) {
      partid <- obj$est$baselearner.cv.batch@instance.list@instances[[1]]@partid
      nfold <- length(unique(obj$est$baselearner.cv.batch@instance.list@partitions[, partid]))
      for (j in 1:nfold) {
        obj$est$baselearner.cv.batch@fitobj.list[[i]]@fitobj.list[[j]]@est <- tmpfiles.new[obj$est$baselearner.cv.batch@tmpfiles.index.list$start[i] + j - 1]
      }
    }
    return(obj)
  } else { # ordinary load
    loadedObjects <- objects(env, all.names = TRUE)
    stopifnot(length(loadedObjects) == 1)
    return (env[[loadedObjects]])
  }
}

print.epcclass <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
}

summary.epcclass <- function(object, ...) {
  #summary.baselearner <- NULL
  #summary.integrator <- summary(object$est$integrator)
  #summary.integrator <- NULL
  #ret <- list(baselearner=summary.baselearner, integrator=summary.integrator)
  
  n.instance <- length(object$instance.list@instances)
  maxpc <- object$est$integrator@est$pcr@sweep.list[[1]]@config@n
  index.min <- object$est$integrator@est$select@est$index.min
  error.min <- object$est$integrator@est$select@est$error.min
  tvec <- object$est$baselearner.cv.batch@tvec
  
  ret <- list(n.instance=n.instance, maxpc=maxpc, index.min=index.min, error.min=error.min, tvec = tvec)
  
  class(ret) <- "summary.epcclass"
  return (ret)
}

print.summary.epcclass <- function(x, ...) {
  cat("number of base learner instances:", x$n.instance, "\n")
  cat("maximum number of PC's considered:", x$maxpc, "\n")
  cat("optimal number of PC's:", x$index.min, "\n")
  cat("minimum error:", x$error.min, "\n")
}


#  ret <- list(call=mycall, formula=formula, instance.list=my.instance.list, integrator.config=my.integrator.config, method=integrator.control$method
#    , est=list(baselearner.cv.batch=est.baselearner.cv.batch, integrator=est.integrator)
#    , y=y, pred=pred, filemethod=filemethod)



