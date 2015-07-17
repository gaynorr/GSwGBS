\name{GS.model}
\alias{GS.model}
\title{
Genomic Selection Modeler
}
\description{
Makes genomic selection predictions using a variety of methods implemented in other R packages.
}
\usage{
GS.model(phenoTrain, genoTrain, genoPredict, 
n.core = "auto", ntree = 1000, fineTune = T, 
methods = c("RRBLUP", "GAUSS", "PLSR", "ELNET", "RF"))
}
\arguments{
  \item{phenoTrain}{
Phenotypic data for your training population. Must be a numeric vector listed in the same order as genoTrain.
}
  \item{genoTrain}{
Genotypic data for your training population. Must be a numerically coded (i.e.; 1, 0, -1) in a matrix or data frame with markers as columns and individuals as rows. Missing values are not allowed. Row names must be unique and can't be the same as those used for genoPredict if using method="GAUSS".
}
  \item{genoPredict}{
Genotypic data for your testing population. Formatting is the same as for genoTrain.
}
  \item{n.core}{
Number of cores to use for parallel computing. If n.core="auto", n.core will be set to the number of detected cores. This may not be desirable due to memory requirements and/or overworking of CPU.
}
  \item{ntree}{
For method="RF", the number of trees in the model.
}
  \item{fineTune}{
For method="ELNET", should fine tuning of the alpha parameter be performed.
}
  \item{methods}{
A vector of the methods to use for predictions. If more than one method is given, the average of all given predictions will also be returned as 'AVE'.
}
}
\details{
Uses genotypic and phenotypic data from a training population to make prediction on new lines based on their genotypic data. Designed to take genotypic data in the format outputted by \link{hap2marker}. A variety of methods are used depending on which values are listed as 'methods'.


"RRBLUP": Uses the 'rrBLUP' package to produce predictions based on estimated marker effects (RR-BLUP) determined using a mixed model. See \code{\link[rrBLUP]{mixed.solve}}


"GUASS": Uses the 'rrBLUP' package to produce predictions based on line effects (G-BLUP) determined using a mixed model. A Gaussian kernel is used for the genotypic covariance. See \code{\link[rrBLUP]{kin.blup}}


"PLSR": Uses the 'pls' package to produce predictions using partial least squares regression. The number of components used for prediction is determined with 10-fold cross-validation in the training population to find the lowest value for the bias-corrected cross-validation estimate. See \code{\link[pls]{plsr}}


"ELNET": Uses the 'glmnet' package to produced predictions using a linear model with elastic net regularization. The elastic net mixing parameter (alpha) and tuning parameter (lambda) are set using a grid selection technique that picks values with the lowest mean cross-validation error from a 10-fold cross-validation across a sequence of alpha and lambda values. The sequence of lambda values tested is determined by the glmnet function. The sequence of alpha values ranges from 0 to 1 in increments of 0.1. If fineTune=T, a new sequence of alpha values is used to perform a second stage of grid selection. The new sequence of alpha values ranges from the optimum alpha minus 0.09, with a minimum of 0, to that optimum plus 0.09, with a maximum of 1, in increments of 0.01. See \code{\link[glmnet]{glmnet}}


"RF": Uses the 'randomForest' package to produce predictions using a random forest regression. The number of trees in the forest is equal to 'ntree' and the number of markers sampled at each split is equal to the total number of markers divided by 3. See \code{\link[randomForest]{randomForest}}
}
\value{
Returns a data frame with predictions for lines in 'genoPredict'. Row names come from 'genoPredict'. Variable names are the methods used to make the predictions.
}
\examples{
\dontrun{
geno = matrix(runif(10000,-1,1),100,100)
dimnames(geno) = list(paste("Line",1:100,sep=""), 
                      paste("M",1:100,sep=""))
trait = rnorm(90)
training = as.data.frame(cbind(trait, geno[1:90,]))
testing = as.data.frame(geno[91:100,])
tmp = GS.model(phenoTrain=training[,1], 
               genoTrain=training[,-1], 
               genoPredict=testing, 
               methods="RRBLUP")
}
}