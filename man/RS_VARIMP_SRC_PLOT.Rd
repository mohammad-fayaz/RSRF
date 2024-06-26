% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RSRF_Version_1_0.R
\name{RS_VARIMP_SRC_PLOT}
\alias{RS_VARIMP_SRC_PLOT}
\title{The Variable Importance Plot for Survival Random Splitting Random Forest}
\usage{
RS_VARIMP_SRC_PLOT(
  Dataset = RF_OUTPUT,
  VIP_Algorithm = c("RS"),
  COVAR_TYPE = c("Functional", "Non-Functional"),
  COV_NUMBER = 3,
  Type = c("Train", "Test")
)
}
\arguments{
\item{Dataset}{It is `RF_OUTPUT$IMP_VALUE` from `RS_RF_SRC_H`.}

\item{VIP_Algorithm}{Only works for RS and}

\item{COVAR_TYPE}{Type of covariate: `Functional` or `Non-Functional`.}

\item{COV_NUMBER}{Functional covariate index :`1`, `2`, `3` and etc.}

\item{Type}{It has two options `Train` and `Test`.}
}
\value{
The variable importance plots for functional and non-functional covariates.
}
\description{
The Variable Importance Plot for functional and non-functional covariates
}
\details{
The `type` is `Train` or `Test` (Train= in train dataset , 2=in test dataset)
}
\examples{
#RF_Result from previous example.
par(mfrow=c(5,2))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=1,Type=c("Train"))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=1,Type=c("Test"))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=2,Type=c("Train"))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=2,Type=c("Test"))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=3,Type=c("Train"))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=3,Type=c("Test"))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=4,Type=c("Train"))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=4,Type=c("Test"))
RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Non-Functional"), COV_NUMBER=1,Type=c("Train"))
}
\references{
\itemize{
 \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
 \item \code{2-} : Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The annals of applied statistics, 2(3), 841-860.
 \item \code{3-} : Ishwaran, H., & Lu, M. (2019). Standard errors and confidence intervals for variable importance in random forest regression, classification, and survival. Statistics in medicine, 38(4), 558-582.
}
}
