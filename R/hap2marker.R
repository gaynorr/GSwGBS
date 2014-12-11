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
  data = data[!duplicated(data[,1:3]),c(-1,-3)]
  colnames(data) = c("alleles","A","B",line.names)
  if(remove.blank){
    data = data[,!grepl("BLANK", colnames(data), ignore.case=TRUE)]
    line.names = colnames(data)[-(1:3)]
  }
  
  cat("Converting to numeric coding...", "\n")
  data[data=="H"]="0"
  recode = function(x){
    alleles=unlist(strsplit(x[1],"/"))
    if(as.numeric(x[2])<as.numeric(x[3])){
      alleles=alleles[c(2,1)] 
    }
    x = x[-1]
    x[x==alleles[1]]="1"
    x[x==alleles[2]]="-1"
    return(as.numeric(x))
  }
  cl = parallel::makeCluster(n.core)
  data = parallel::parApply(cl=cl, data, 1, recode)
  parallel::stopCluster(cl)
  rm(cl)
  data = t(data)
  cat("Removing redundant markers...", "\n")
  tmp.data = data[data[,1]==data[,2],]
  data = data[data[,1]!=data[,2],-(1:2)]
  data = unique(data)
  data = t(data)
  if(nrow(tmp.data)>0){
    tmp.data = tmp.data[,-(1:2)]
    tmp.data = unique(tmp.data)
    if(nrow(tmp.data)>1){
      tmp.data = t(tmp.data)
      cl = parallel::makeCluster(n.core)
      registerDoParallel(cl)
      tmp = foreach(i=1:(ncol(tmp.data)-1), .combine=c) %dopar% {
		i = i
        j = (i+1):ncol(tmp.data)
        test = as.matrix(tmp.data[,i]*(-1))
        result1 = apply(as.matrix(tmp.data[1:30,j]), 2, identical, test[1:30])
        if(any(result1)){
          j = j[result1]
          result = apply(as.matrix(tmp.data[,j]), 2, identical, test)
          if(any(result, na.rm=T)){
            -i
          }else{
            NA
          }
        }else{
          NA
        }
      }
      parallel::stopCluster(cl)
      rm(cl)
      tmp = na.omit(as.integer(tmp))
      if(length(tmp)>0){
        tmp.data = tmp.data[,tmp]
      }
    }
    data = cbind(data,tmp.data)
    rm(tmp)
  }
  rm(tmp.data)
  colnames(data) = paste("M", 1:(ncol(data)), sep="")
  rownames(data) = line.names
  rm(line.names)
  
  #impute missing data
  if(impute.method[1]!="none"){
    cat("Performing marker imputation...", "\n")
    if(impute.method[1]=="median"){
      data = randomForest::na.roughfix(data)
    } else if(impute.method[1]=="mean" | impute.method[1]=="EM"){
      if(.Platform$OS.type=="windows") n.core=1
      data = rrBLUP::A.mat(data, impute.method=impute.method[1], n.core=n.core, 
                           return.imputed=TRUE, ...)$imputed
    } else if(impute.method[1]=="RF"){
      cl = parallel::makeCluster(n.core)
      registerDoParallel(cl)
      data = missForest::missForest(data, parallelize="variables", ...)$ximp
      parallel::stopCluster(cl)
      rm(cl)
    } else warning("impute.method not recognized, no imputation performed")
  }
  
  invisible(data)
}
