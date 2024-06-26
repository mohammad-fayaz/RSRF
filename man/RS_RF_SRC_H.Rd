% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RSRF_Version_1_0.R
\name{RS_RF_SRC_H}
\alias{RS_RF_SRC_H}
\title{Random Splitting Random Forest SRC for Hybrid Data}
\usage{
RS_RF_SRC_H(
  ResponseVar,
  Covariates,
  RawData,
  Stat,
  min,
  max,
  Params = c(0, 0, 0),
  Distribution,
  Params2 = c(0, 0, 0),
  Distribution2,
  type,
  m = 100,
  k = 500,
  ResponseType = c("Categorical", "Continuous", "Multivariate", "Unsupervised",
    "Survival", "Competing Risk"),
  SplitRule,
  Importance,
  Block.Size,
  Ensemble,
  Proximity,
  RFMethod = c("RandomForest", "Bagging"),
  NUM_Covariates = NUM_Covariates,
  VIP_Algorithm = c("Normal", "Normal_Joint", "RS", "RS_Joint")
)
}
\arguments{
\item{ResponseVar}{The response variable. Add c(rep(0,4),real_response-variable),}

\item{Covariates}{The non-functional covariates. Same name : data.frame(c_1 = c(rep(0,4),real_cov_1),c_2 = c(rep(0,4),real_cov_2)),}

\item{RawData}{The functional data [rows are number of curves][columns are time points]}

\item{Stat}{The summary statistics for the random interval: mean,median,Q25,Q75,sd,max,min,range}

\item{min}{The minimum value of the domain of the functional data.}

\item{max}{The maximum value of the domain of the functional data.}

\item{Params}{The parameters for random number distribution. It includes:
\itemize{
      \item \code{Exponential} : c(rate,0,0)           rate = 1/theta
      \item \code{Normal}      : c(mean , sd,0)
      \item \code{Uniform}     : c(min,max,0)
      }}

\item{Distribution}{The distribution for random number (r*) : Exponential,Normal, Uniform}

\item{Params2}{If type is "Overlap",the parameters for random number for distribution 2:
\itemize{
     \item \code{Exponential} : Params2[1] is rate , Params2[2] is constant c>1. Default = 5
     \item \code{Normal}      : Paramsp2[1] is mu, Paramsp2[2] is sigma,Paramsp2[3] is constant c>1.Default = 5,
     \item \code{Uniform}     : Paramsp2[1] is min, Paramsp2[2] is max,Paramsp2[3] is constant c>1. Default = 5,
     }}

\item{Distribution2}{If type is "Overlap",the distribution  :
Exponential,Normal, Uniform}

\item{type}{Two options: 1-Disjoint, 2-Overlap}

\item{m}{Number of Forests (default = 100)}

\item{k}{Iterations in Each Random Forest (default = 500)}

\item{ResponseType}{The Response Type for model predictions: "Categorical", "Continuous", "Multivariate","Unsupervised","Survival" and "Competing Risk".}

\item{SplitRule}{The Splitting Rule based on the response variable.}

\item{RFMethod}{"Bagging" for (MD-BG) and "RandomForest" for (MD-RF)}

\item{NUM_Covariates}{Number of Covariets. Each functional covariates is one covariates.
For example a model with 3 functional covariates and
2 non-functional covariates have 2+3 = 5 covariates.}

\item{importance}{The variable importance (VIMP) as in rfsrc() from randomForestSRC.}

\item{block.size}{The block.size as in rfsrc() from randomForestSRC.}

\item{ensemble}{The ensemble as in rfsrc() from randomForestSRC.}

\item{bootstrap}{The bootstrap as in rfsrc() from randomForestSRC.}

\item{samptype}{The sampling type as in rfsrc() from randomForestSRC.}

\item{sampsize}{The sampling size as in rfsrc() from randomForestSRC.}

\item{proximity}{The proximity as in rfsrc() from randomForestSRC.}
}
\value{
The output is a list. It contains:
    list[[1]]   ## the summary statistics
    list[[2]]   ## The variable importance
}
\description{
Random splitting random forest for Survival, Regression and Classification (SRC) for hybrid data.
}
\examples{
library(RandomForestMixedData)

 ######## Simulations
#Number of Samples is 100
nSample_1 <- 100
SIM_1_FUNC_DATA <- makeSIMData(nSS = nSample_1, responseType =c("Survival"),Seed= 2400)
par(mfrow=c(2,3))
matplot(t(SIM_1_FUNC_DATA[[1]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[1](t)))
matplot(t(SIM_1_FUNC_DATA[[2]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[2](t)))
matplot(t(SIM_1_FUNC_DATA[[3]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[3](t)))
matplot(t(SIM_1_FUNC_DATA[[4]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[4](t)))
boxplot(SIM_1_FUNC_DATA[[5]][1][-(1:4),][(1:nSample_1/2)], SIM_1_FUNC_DATA[[5]][1][-(1:4),][-(1:nSample_1/2)],col=c("Green","Grey"),main=expression(X[1]))
boxplot(SIM_1_FUNC_DATA[[5]][2][-(1:4),][(1:nSample_1/2)], SIM_1_FUNC_DATA[[5]][2][-(1:4),][-(1:nSample_1/2)],col=c("Green","Grey"),main=expression(X[2]))
######## Model Specification and Running
## Survival
time  <- c(SIM_1_FUNC_DATA[[6]])
event <- c(SIM_1_FUNC_DATA[[7]])
Ys <- data.frame(time,event)

## Functional Covariate Specifications

######### With Functional PCA (FPCA)
FDList <- list()
FDList[[1]] <- UTIL_REC_FPCA(InputFunction = t(SIM_1_FUNC_DATA[[1]]),numberofPCA = 1)
FDList[[2]] <- UTIL_REC_FPCA(InputFunction = t(SIM_1_FUNC_DATA[[2]]),numberofPCA = 1)
FDList[[3]] <- UTIL_REC_FPCA(InputFunction = t(SIM_1_FUNC_DATA[[3]]),numberofPCA = 1)
FDList[[4]] <- UTIL_REC_FPCA(InputFunction = t(SIM_1_FUNC_DATA[[4]]),numberofPCA = 1)

#matplot(t(matrix(unlist(FDList[[1]]),ncol=100)),type="l")
#matplot(t(matrix(unlist(FDList[[2]]),ncol=100)),type="l")
#matplot(t(matrix(unlist(FDList[[3]]),ncol=100)),type="l")
#matplot(t(matrix(unlist(FDList[[4]]),ncol=100)),type="l")

######### Without Functional PCA (FPCA)

FDList <- list()
FDList[[1]] <- SIM_1_FUNC_DATA[[1]]
FDList[[2]] <- SIM_1_FUNC_DATA[[2]]
FDList[[3]] <- SIM_1_FUNC_DATA[[3]]
FDList[[4]] <- SIM_1_FUNC_DATA[[4]]

#matplot(t(matrix(unlist(FDList[[1]]),ncol=100)),type="l")
#matplot(t(matrix(unlist(FDList[[2]]),ncol=100)),type="l")
#matplot(t(matrix(unlist(FDList[[3]]),ncol=100)),type="l")
#matplot(t(matrix(unlist(FDList[[4]]),ncol=100)),type="l")

 ## Parameter Specification for each Functioanl Covariate
Min <- list()
Max <- list()
PARAMS_1 <- list()
DIST_1 <- list()
TYPE <- list()
PARAMS_2 <- list()
STATS <- list()
DIST_2 <- list()

### Functional Covariate 1 Specificaion
Min[[1]] <- 1                    ## minimum time domain
Max[[1]] <- 100                  ## maximum time domain
DIST_1[[1]] <- c("Exponential")  ## The Random splitting distribution
PARAMS_1[[1]] <- c(0.1,0,0)      ## The parameter for the distribution
TYPE[[1]] <- c("Disjoint")       ## Type of splitting
DIST_2[[1]] <- c("")             ## if type is overlap , the distribution of the overlap
PARAMS_2[[1]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
STATS[[1]]<- c("mean")          ## the statistics for each interval

### Functional Covariate 2 Specificaion
Min[[2]] <- 1                    ## minimum time domain
Max[[2]] <- 100                  ## maximum time domain
DIST_1[[2]] <- c("Exponential")  ## The Random splitting distribution
PARAMS_1[[2]] <- c(0.1,0,0)      ## The parameter for the distribution
TYPE[[2]] <- c("Disjoint")       ## Type of splitting
DIST_2[[2]] <- c("")             ## if type is overlap , the distribution of the overlap
PARAMS_2[[2]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
STATS[[2]]<- c("mean")          ## the statistics for each interval

### Functional Covariate 3 Specificaion
Min[[3]] <- 1                    ## minimum time domain
Max[[3]] <- 100                  ## maximum time domain
DIST_1[[3]] <- c("Exponential")  ## The Random splitting distribution
PARAMS_1[[3]] <- c(0.1,0,0)      ## The parameter for the distribution
TYPE[[3]] <- c("Disjoint")       ## Type of splitting
DIST_2[[3]] <- c("")             ## if type is overlap , the distribution of the overlap
PARAMS_2[[3]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
STATS[[3]]<- c("mean")          ## the statistics for each interval

### Functional Covariate 4 Specificaion
Min[[4]] <- 1                    ## minimum time domain
Max[[4]] <- 100                  ## maximum time domain
DIST_1[[4]] <- c("Exponential")  ## The Random splitting distribution
PARAMS_1[[4]] <- c(0.1,0,0)      ## The parameter for the distribution
TYPE[[4]] <- c("Disjoint")       ## Type of splitting
DIST_2[[4]] <- c("")             ## if type is overlap , the distribution of the overlap
PARAMS_2[[4]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
STATS[[4]]<- c("mean")          ## the statistics for each interval

## Non-Functional Covariate Specifications
Covariates_all <- data.frame(c_01= SIM_1_FUNC_DATA[[5]][,1],c_02= SIM_1_FUNC_DATA[[5]][,2])

## General Settings
mtree <- 100
k <- 500

######## Running the model
BG_Res <- RS_RF_SRC_H(
  ResponseVar=Ys,
  Covariates=Covariates_all,
  RawData=FDList,
  Stat=STATS ,
  min=Min,
  max=Max,
  Params=PARAMS_1,
  Distribution=DIST_1,
  Params2=PARAMS_2,
  Distribution2=DIST_2,
  type=TYPE,
  m=500,
  k=k,
  ResponseType="Survival",
  RFMethod = c("Bagging"),
  NUM_Covariates=1,
  Block.Size = 10 ,
  SplitRule = c("logrank"),
  Importance= TRUE,
  Ensemble = c("all"),
  Proximity = TRUE)
}
\references{
\itemize{
 \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
 \item \code{2-} : Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The annals of applied statistics, 2(3), 841-860.
 \item \code{3-} : Ishwaran, H., & Lu, M. (2019). Standard errors and confidence intervals for variable importance in random forest regression, classification, and survival. Statistics in medicine, 38(4), 558-582.
}
}
