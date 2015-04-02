hap2marker <-
  function(files = c("","",""), miss.max=0.7, remove.blank=TRUE, n.core="auto", impute.method=c("none", "median", "mean", "EM", "RF"), ...){
    switch(impute.method[1], "median"=require(randomForest), "mean"=require(rrBLUP), 
           "EM"=require(rrBLUP), "RF"=require(missForest))
    if(n.core=="auto"){
      n.core=parallel::detectCores()
    }
    cat("Reading hapfiles...", "\n")
    
    header = read.table(files[1],colClasses="character",stringsAsFactors=F, 
                        comment.char="", nrows=1)
    line.names = as.character(header[,-(1:11)])
    
    hap1 = scan(files[1],character(), comment.char="", na.strings="N", skip=1, quiet=T)
    dim(hap1) = c(length(header),length(hap1)/length(header))
    hap1 = t(hap1)
    hap1 = hap1[as.numeric(hap1[,11])>=1-miss.max,-c(3:5,7,10:11)]
    hap2 = scan(files[2],character(), comment.char="", na.strings="N", skip=1, quiet=T)
    dim(hap2) = c(length(header),length(hap2)/length(header))
    hap2 = t(hap2)
    hap2 = hap2[as.numeric(hap2[,11])>=1-miss.max,-c(3:5,7,10:11)]
    data = rbind(hap1,hap2)
    rm(hap1,hap2)
    data = data[!duplicated(data[,1:3]),]
    hap3 = scan(files[3],character(), comment.char="", na.strings="N", skip=1, quiet=T)
    dim(hap3) = c(length(header),length(hap3)/length(header))
    hap3 = t(hap3)
    hap3 = hap3[as.numeric(hap3[,11])>=1-miss.max,-c(3:5,7,10:11)]
    data = rbind(data,hap3)
    rm(hap3, header)
    data = data[!duplicated(data[,1:3]),]
    colnames(data) = c("seq","alleles","pos","A","B",line.names)
    if(remove.blank){
      data = data[,!grepl("BLANK", colnames(data), ignore.case=TRUE)]
      line.names = colnames(data)[-(1:5)]
    }
    data[,"pos"] = as.numeric(data[,"pos"]) + 1
    
    cat("Converting to numeric coding...", "\n")
    data[data=="H"]="0"
    recode = function(x){
      alleles=unlist(strsplit(x[2],"/"))
      if(as.numeric(x[4])<as.numeric(x[5])){
        alleles=alleles[c(2,1)] 
      }
      x[2] = paste(alleles[1],alleles[2],sep="/")
      x[x==alleles[1]]="1"
      x[x==alleles[2]]="-1"
      return(x)
    }
    cl = parallel::makeCluster(n.core)
    data = parallel::parApply(cl=cl, data, 1, recode)
    data = t(data)
    
    cat("Removing redundant markers...", "\n")
    data = data[!duplicated(data[,-(1:5)]),]
    tmp.data = data[data[,4]==data[,5],-(4:5)]
    data = data[data[,4]!=data[,5],-(4:5)]
    if(nrow(tmp.data)>0){
      if(nrow(tmp.data)>1){
        tmp.data = t(tmp.data)
        registerDoParallel(cl)
        tmp = foreach(i=1:(ncol(tmp.data)-1), .combine=c) %dopar% {
          i = i
          j = (i+1):ncol(tmp.data)
          test = as.matrix(as.character(as.numeric(tmp.data[-(1:3),i])*(-1)))
          result1 = apply(as.matrix(tmp.data[4:33,j]), 2, identical, test[1:30])
          if(any(result1)){
            j = j[result1]
            result = apply(as.matrix(tmp.data[-(1:3),j]), 2, identical, test)
            if(any(result, na.rm=T)){
              -i
            }else{
              NA
            }
          }else{
            NA
          }
        }
        tmp = na.omit(as.integer(tmp))
        if(length(tmp)>0){
          tmp.data = tmp.data[,tmp]
        }
        rm(tmp)
        tmp.data = t(tmp.data)
      }
      data = rbind(data,tmp.data)
    }
    rm(tmp.data)
    parallel::stopCluster(cl)
    rm(cl)
    rownames(data) = paste("M", 1:(nrow(data)), sep="")
    data = list(markers=data[,-(1:3)], annotations=data[,1:3])
    tmp.dim = dim(data[[1]])
    tmp.name = dimnames(data[[1]])
    data[[1]] = as.numeric(data[[1]])
    dim(data[[1]]) = tmp.dim
    dimnames(data[[1]]) = tmp.name
    data[[1]] = t(data[[1]])
    rm(line.names,tmp.dim,tmp.name)
    
    #impute missing data
    if(impute.method[1]!="none"){
      cat("Performing marker imputation...", "\n")
      if(impute.method[1]=="median"){
        data[[1]] = randomForest::na.roughfix(data[[1]])
      } else if(impute.method[1]=="mean" | impute.method[1]=="EM"){
        if(.Platform$OS.type=="windows") n.core=1
        data[[1]] = rrBLUP::A.mat(data[[1]], impute.method=impute.method[1], n.core=n.core, 
                             return.imputed=TRUE, ...)$imputed
      } else if(impute.method[1]=="RF"){
        cl = parallel::makeCluster(n.core)
        registerDoParallel(cl)
        data[[1]] = missForest::missForest(data[[1]], parallelize="variables", ...)$ximp
        parallel::stopCluster(cl)
        rm(cl)
      } else warning("impute.method not recognized, no imputation performed")
    }
    
    invisible(data)
  }
