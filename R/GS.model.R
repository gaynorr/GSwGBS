GS.model <-
function(phenoTrain, genoTrain, genoPredict, n.core="auto", ntree=1000, fineTune=T, methods=c("RRBLUP", "GAUSS", "PLSR", "ELNET", "RF")){
  
  genoTrain = as.matrix(genoTrain)
  genoPredict = as.matrix(genoPredict)
  methods = toupper(methods)
  if(length(methods)>1){
    methods = c(methods, "AVE")
  }
  
  #check data
  if(!is.numeric(phenoTrain)){
    stop("phenoTrain must be numeric")
  }
  if(!is.vector(phenoTrain)){
    stop("phenoTrain must be a vector")
  }
  if(any(methods=="GAUSS") & any(row.names(genoTrain)%in%row.names(genoPredict))){
    stop("unique row names for genoPredict and genoTrain must be used with method='GAUSS'")
  }
  if(!is.numeric(genoTrain)){
    stop("genoTrain must be numerically coded")
  }
  if(any(is.na(genoTrain))){
    stop("genoTrain cannot contain missing values")
  }
  if(!is.numeric(genoPredict)){
    stop("genoPredict must be numerically coded")
  }
  if(any(is.na(genoPredict))){
    stop("genoPredict cannot contain missing values")
  }
  if(ncol(genoTrain)!=ncol(genoPredict)){
    stop("marker number in genoTrain and genoPredict must be equal")
  }
  
  #load required packages
  if(n.core=="auto"){
     n.core=parallel::detectCores()
  }
  if(any(methods=="RRBLUP") | any(methods=="GAUSS")){
    require(rrBLUP)
  }
  if(any(methods=="PLSR")){
    require(pls)
  }
  if(any(methods=="ELNET")){
    require(glmnet)
  }
  if(any(methods=="RF")){
    require(randomForest)
	ntree = ceiling(ntree/n.core)
  }
  
  #prepare output
  output = matrix(nrow=nrow(genoPredict), ncol=length(methods))
  storage.mode(output) = "numeric"
  dimnames(output) = list(row.names(genoPredict), methods)
  output = as.data.frame(output)
  
  #build functions for GS methods
  RRBLUP.F = function(...){
    cat("Executing RR-BLUP...", "\n")
    model = rrBLUP::mixed.solve(phenoTrain, Z=genoTrain)
    tmp = genoPredict %*% model$u
    prediction = as.vector(tmp) + model$beta
    invisible(as.vector(prediction))
  }
  GAUSS.F = function(...){
    cat("Executing Gaussian Kernel...", "\n")
    pheno = c(phenoTrain, rep(NA, length(genoPredict[,1])))
    pheno.names = c(row.names(genoTrain), row.names(genoPredict))
    data = data.frame(pheno.names, pheno)
    geno = rbind(genoTrain, genoPredict)
    if(.Platform$OS.type=="windows") gk.core=1 else gk.core=n.core
    model = rrBLUP::kin.blup(data=data, geno=names(data)[1], pheno=names(data)[2],
							 GAUSS=TRUE, K=as.matrix(dist(geno)), n.core=gk.core)
    prediction = model[[4]][(length(phenoTrain)+1):length(pheno)]
    prediction = prediction + mean(phenoTrain)
    invisible(as.vector(prediction))
  }
  PLSR.F = function(...){
    cat("Executing PLS Regression...", "\n")
    pls::pls.options(parallel=cl)
    model = pls::plsr(phenoTrain~genoTrain, ncomp=25, validation="CV")
    error = pls::MSEP(model, "adjCV", intercept=FALSE)
    error.min = which.min(as.vector(error[[1]]))
    if(error.min>=24){
      cat("     Rebuilding Model Using ncomp=50...", "\n")
      model = pls::plsr(phenoTrain~genoTrain, ncomp=50, validation="CV")
      error = pls::MSEP(model, "adjCV", intercept=FALSE)
      error.min = which.min(as.vector(error[[1]]))
    }
    if(error.min>=49){
      cat("     Rebuilding Model Using ncomp=100...", "\n")
      model = pls::plsr(phenoTrain~genoTrain, ncomp=100, validation="CV")
      error = pls::MSEP(model, "adjCV", intercept=FALSE)
      error.min = which.min(as.vector(error[[1]]))
    }
    cat("     Components = ", error.min, "\n")
    prediction = predict(model, genoPredict, ncomp = error.min)
    invisible(as.vector(prediction))
  }
  ELNET.F = function(...){
    cat("Executing Elastic Net...", "\n")
    alphas = seq(0,1,by=0.1)
    foldid = sample(rep_len(1:10, length(phenoTrain)), length(phenoTrain))
    cat("  Rough Tuning Alpha...", "\n")
    genoTrain = genoTrain
    phenoTrain = phenoTrain
    tmp = foreach(i=alphas, .packages="glmnet") %dopar% {
      glmnet::cv.glmnet(genoTrain, phenoTrain, alpha=i, foldid=foldid, 
                        family="gaussian", standardize=FALSE)
    }
    cvm = vector(length=length(alphas), mode="numeric")
    for (i in 1:length(alphas)){
      cvm[i] = min(tmp[[i]]$cvm)
    }
    if(fineTune){
      new.alpha = (which.min(cvm)-1)*0.1
      alpha.min = new.alpha - 0.1
      alpha.max = new.alpha + 0.1
      if(alpha.min<0) {alpha.min = 0}
      if(alpha.max>1) {alpha.max = 1}
      alphas = seq(alpha.min, alpha.max, 0.01)
      cat("  Fine Tuning Alpha...", "\n")
      tmp = foreach(i=alphas, .packages="glmnet") %dopar% {
        glmnet::cv.glmnet(genoTrain, phenoTrain, alpha=i, foldid=foldid, 
                          family="gaussian", standardize=FALSE)
      }
      cvm = vector(length=length(alphas), mode="numeric")
      for (i in 1:length(alphas)){
        cvm[i] = min(tmp[[i]]$cvm)
      }
      model = tmp[[which.min(cvm)]]
      parms = c(((which.min(cvm)-1)*0.01 + alpha.min), model$lambda.min)
    }else{
      model = tmp[[which.min(cvm)]]
      parms = c(((which.min(cvm)-1)*0.1), model$lambda.min)
    }

    cat("     Alpha =  ", parms[1], "\n")
    cat("     Lambda = ", parms[2], "\n")
    model = glmnet::glmnet(genoTrain, phenoTrain, alpha=parms[1], family="gaussian", 
                           standardize=FALSE)
    prediction = predict(model, newx=genoPredict, s=parms[2])
    invisible(as.vector(prediction))
  }
  RF.F = function(...){
    cat("Executing Random Forest...", "\n")
    genoTrain = genoTrain
    phenoTrain = phenoTrain
    model = foreach(ntree=rep(ntree, n.core), .combine=randomForest::combine, 
                    .packages="randomForest")%dopar%{
      randomForest::randomForest(genoTrain, phenoTrain, ntree=ntree, proximity=FALSE)
    }
    prediction = predict(model, genoPredict)
    invisible(as.vector(prediction))
  }
  
  #execute selected methods
  if(any(methods=="RRBLUP")){
    output$RRBLUP = RRBLUP.F()
  }
  if(any(methods=="GAUSS")){
    output$GAUSS = GAUSS.F()
  }
  if(any(methods!="RRBLUP") & any(methods!="GAUSS")){
    cat("Building Workstation Cluster...", "\n")
    cl = parallel::makeCluster(n.core)
    registerDoParallel(cl)
  }
  if(any(methods=="PLSR")){
    output$PLSR = PLSR.F()
  }
  if(any(methods=="ELNET")){
    output$ELNET = ELNET.F()
  }
  if(any(methods=="RF")){
    output$RF = RF.F()
  }
  if(any(methods!="RRBLUP") & any(methods!="GAUSS")){
    parallel::stopCluster(cl)
  }
  if(any(methods=="AVE")){
    output$AVE = apply(as.matrix(output[,-length(names(output))]), 1, mean)
  }
  names(output) = methods
  cat("DONE","\n")
  invisible(output)
}
