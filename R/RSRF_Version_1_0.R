#' Making Simulation Data
#'
#' Making Simulated Data
#' @param nSS_1           Number of Total Samples in Group 1
#' @param nSS_2           Number of Total Samples in Group 2
#' @param responseType    Type of response ("Categorical","Continuous" and "Survival")
#' @param Seed            Seed number for memory.
#' @details               It makes three functional and two non-functional covariates. The first nSS/2 is case (green color) and second nSS/2 is contorl (grey color).
#' @return                A list with 6 elements.   The 1 to 4 are functional covariets , the 5 is the dataframe for non-functional covariates and  is the response.
#' @references
#' \itemize{
#'  \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
#'  \item \code{2-} : Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The annals of applied statistics, 2(3), 841-860.
#'  \item \code{3-} : Ishwaran, H., & Lu, M. (2019). Standard errors and confidence intervals for variable importance in random forest regression, classification, and survival. Statistics in medicine, 38(4), 558-582.
#' }
#' @examples
#'  #Number of Samples is 100
#'   nSample_1 <- 100
#'   nSample_2 <- 100
#'   SIM_1_FUNC_DATA <- makeSIMData(nSS_1 = nSample_1,nSS_2 = nSample_2, responseType ="Categorical",Seed= 2400)
#'   par(mfrow=c(2,3))
#'   matplot(t(SIM_1_FUNC_DATA[[1]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[1](t)))
#'   matplot(t(SIM_1_FUNC_DATA[[2]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[2](t)))
#'   matplot(t(SIM_1_FUNC_DATA[[3]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[3](t)))
#'   matplot(t(SIM_1_FUNC_DATA[[4]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[4](t)))
#'   boxplot(SIM_1_FUNC_DATA[[5]][1][-(1:4),][(1:nSample_1/2)], SIM_1_FUNC_DATA[[5]][1][-(1:4),][-(1:nSample_1/2)],col=c("Green","Grey"),main=expression(X[1]))
#'   boxplot(SIM_1_FUNC_DATA[[5]][2][-(1:4),][(1:nSample_1/2)], SIM_1_FUNC_DATA[[5]][2][-(1:4),][-(1:nSample_1/2)],col=c("Green","Grey"),main=expression(X[2]))
#'   par()
#' @export
makeSIMData <- function(nSS_1,nSS_2, responseType,Seed){
  
  set.seed(Seed)
  
  #### Define Simuation Functions
  t_1 <- round(seq(0,10*pi,length.out=100),1)
  x_1_1 <- function(x) sin(sqrt(x)) + rnorm(1,mean=0,sd=0.01)
  x_1_2 <- function(x) c(sin(sqrt(x[1:35]))+0.02 + rnorm(1,mean=0,sd=0.01) ,
                         cos((x[36:45]))+ rnorm(1,mean=0,sd=0.01) ,
                         sin(sqrt(x[46:100]))+0.02)+ rnorm(1,mean=0,sd=0.01)
  
  t_2 <- round(seq(0,2*pi,length.out=100),1)
  x_2_1 <- function(x) sin(x) + rnorm(1,mean=0,sd=0.3)
  x_2_2 <- function(x) sin(x)+ rnorm(1,mean=0,sd=0.3)
  
  t_3 <- round(seq(0,20*pi,length.out=100),1)
  x_3_1 <- function(x) sin(x) + rnorm(1,mean=0,sd=0.1)
  x_3_2 <- function(x) c(sin(x[1:30]) +  rnorm(1,mean=0,sd=0.1),
                         sin(x[31:55]) + 2 + rnorm(1,mean=0,sd=0.1) ,
                         sin(x[56:100]) + rnorm(1,mean=0,sd=0.1))
  
  t_4 <- round(seq(0,pi*50,length.out=100),2)
  x_4_1 <- function(x) (sin(x))^3 + rnorm(1,mean=0,sd=0.1)
  x_4_2 <- function(x) c((sin(x[1:10]))^2 + rnorm(1,mean=0,sd=0.1),
                         (sin(x[11:90]))^3 + rnorm(1,mean=0,sd=0.1),
                         (sin(x[91:100]))^2 + rnorm(1,mean=0,sd=0.1))
  
  #### Functional Covariates
  X_1_1 <- x_1_1(x=t_1)
  for( i in 1:(nSS_1-1) ) X_1_1 <- rbind(X_1_1, x_1_1(x=t_1) )
  X_1_2 <- x_1_2(x=t_1)
  for( i in 1:(nSS_2-1) ) X_1_2 <- rbind(X_1_2, x_1_2(x=t_1) )
  X_t_1 <- rbind(X_1_1,X_1_2)
  
  X_2_1 <- x_2_1(x=t_2)
  for( i in 1:(nSS_1-1) ) X_2_1 <- rbind(X_2_1, x_2_1(x=t_2) )
  X_2_2 <- x_2_2(x=t_2)
  for( i in 1:(nSS_2-1) ) X_2_2 <- rbind(X_2_2, x_2_2(x=t_2) )
  X_t_2 <- rbind(X_2_1,X_2_2)
  
  X_3_1 <- x_3_1(x=t_3)
  for( i in 1:(nSS_1-1) ) X_3_1 <- rbind(X_3_1, x_3_1(x=t_3) )
  X_3_2 <- x_3_2(x=t_3)
  for( i in 1:(nSS_2-1)  ) X_3_2 <- rbind(X_3_2, x_3_2(x=t_3) )
  X_t_3 <- rbind(X_3_1,X_3_2)
  
  X_4_1 <- x_4_1(x=t_4)
  for( i in 1:(nSS_1-1) ) X_4_1 <- rbind(X_4_1, x_4_1(x=t_4) )
  X_4_2 <- x_4_2(x=t_4)
  for( i in 1:(nSS_2-1) ) X_4_2 <- rbind(X_4_2, x_4_2(x=t_4) )
  X_t_4 <- rbind(X_4_1,X_4_2)
  
  #### Non-Functional Covariates
  c_01 <- c(0,0,0,0,c(rnorm(n = nSS_1 , mean = 10.5 , sd = 0.2 ), rnorm(n = nSS_2 , mean = 10.7 , sd = 0.2 )))
  c_02 <- c(0,0,0,0,c(rnorm(n = nSS_1 , mean = 60 , sd = 20), rnorm(n = nSS_2 , mean = 1  , sd = 15)))
  Covariates_all  <- data.frame(c_01,c_02)
  
  if(responseType =="Categorical"){
    Yres <- c(rep(0,4),rep(0,nSS_1),rep(1,nSS_2))
    output <- list (X_t_1,X_t_2,X_t_3,X_t_4,Covariates_all,Yres)
  }
  if(responseType =="Continuous"){
    Yres <- c(rep(0,4),rnorm(nSS_1, mean  = 2, sd = 4),rnorm(nSS_2, mean  = 8, sd = 12))
    output <- list (X_t_1,X_t_2,X_t_3,X_t_4,Covariates_all,Yres)
  }
  if(responseType =="Survival"){
    Time <- c(rep(0,4),rweibull(n = nSS_1, shape = 2,scale=4),rweibull(n = nSS_2, shape = 18,scale=6))
    Event <- c(rep(0,4),rep(0,nSS_1),rep(1,nSS_2))
    output <- list (X_t_1,X_t_2,X_t_3,X_t_4,Covariates_all,Time,Event)
  }
  
  return(output)
}

#' Split the Domain 
#'
#' It splits the domain (interval) from  minimum to maximum by the random numbers generated from Normal, Exponential or Uniform distributions.
#' @param min            The minimum value of the domain.
#' @param max            The maximum value of the domain.
#' @param Params         The parameters of the distribution. 
#'                           In the "Exponential" is (rate, 0,0)
#'                           In the "Uniform" is  (min,max,0)
#'                           In the "Normal"  is  (Mu, Sigma,0)
#' @param Distribution   The Distribution name. ("Exponential", "Uniform" and "Normal").
#' @param Params2        The parameters of the second distribution, if type is "Overlap".
#'                           In the "Exponential" is (rate, Constant,0), Constant default is 5. Final rate is rate*Constant.
#'                           In the "Uniform" is  (min,max,Constant), Constant default is 5. Final min and max are min/Constant and max/Constant.
#'                           In the "Normal"  is  (Mu, Sigma,Constant) , Constant default is 5. Final sd is sd/Constant.
#' @param Distribution2  The Distribution name. ("Exponential", "Uniform" and "Normal"), if type is "Overlap".
#' @param type           The type of splitting, "Disjoint" or "Overlap". 
#' @return It returns a list:
#' @return                num: the number of splits.
#' @return                Rs : the random number.
#' @return                lr : the minimum value of the split.
#' @return                ur : the maximum value of the split.
#' @references
#' \itemize{
#'  \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
#' }
#' @examples
#' ## Split a randomly a domain from 0 to 10 by random numbers from exponential distribution with rate = 1 and disjoint.
#' set.seed(1234)
#' Split_Domain(min=0,max=10, Params=c(1,0,0), Distribution = "Exponential",type="Disjoint")
#' ## Split a randomly a domain from 0 to 10 by random numbers from exponential distribution with rate = 1 and Overlap. The the overlap is also obtained from random numbers from the Exponential distribution with the rate parameter equal to 1 and divided by the constant equal to 5  .
#' set.seed(1234)
#' Split_Domain(min=0,max=10, Params=c(1,0,0), Distribution = "Exponential",Params2=c(1,5,0),Distribution2="Exponential" , type="Overlap")
#' @export
Split_Domain <- function(min,max,Params=c(0,0,0),Distribution,Params2=c(0,0,0),Distribution2,type)
{
  
  num <- Rs <- lr <- ur <- 0
  output <- data.frame(num,Rs,lr, ur)
  if(Distribution == "Exponential")  { Fiter <- 2*(max*Params[1]) }
  if(Distribution == "Uniform")      { Fiter <- 2*(max*(Params[2]-Params[1])) }
  if(Distribution == "Normal")       { Fiter <- 2*(max*Params[1]+max*Params[2])}
  
  for(i in 1:Fiter){
    if(Distribution == "Exponential"){ Rs <-  rexp(n=1 ,rate = Params[1]) }   ## Paramsp[1] is rate = 1/theta}
    if(Distribution == "Uniform")    { Rs <-  runif(n=1,min  = Params[1] , max = Params[2])} ## Params[1] is minimum and Params[2] is maximum. 
    if(Distribution == "Normal")     {
      RRnorm <- -1
      while(RRnorm < 0) {RRnorm <-  rnorm(n=1, mean = Params[1],sd = Params[2]) } ## Params[1] is mean and Params[2] is sd. 
      Rs <- RRnorm
    }
    
    if(type=="Disjoint"){
      
      if(i == 1){
        minV = min;
        output[i,] = c(i,round(Rs,3),round(minV,3),round(minV+Rs,3))
      }
      if (i >1)
      {
        minV = as.numeric(output[i-1,4])
        if(max > round(minV+Rs,3)){
          output[i,] = c(i,round(Rs,3),round(minV,3),round(minV+Rs,3))
        } else {
          output[i,] = c(i,round(max-round(minV,3),3),round(minV,3),max)
          
          break
        }
      }
    }
    
    if(type=="Overlap"){
      Rs2 <-0
      if(Distribution2 == "Exponential"){ Rs2 <-  rexp(n=1,rate = Params2[1] * Params2[2]) }   ## Paramsp2[1] is rate , Params2[2] is constant c>1 Tuning Default = 5
      if(Distribution2 == "Uniform")    { Rs2 <-  runif(n=1,min=Params2[1]/Params2[3],max=Params2[2]/Params2[3])}  ## Params2[1] is min, Params2[2] is max, Params2[3] is constant c>1 Default = 5,
      if(Distribution2 == "Normal")     {
        RRnorm2 <- -1
        while(RRnorm2 < 0) {RRnorm2 <-  rnorm(n=1,mean = Params2[1],sd = Params2[2]/Params2[3])  }  ## Paramsp2[1] is mu, Params2[2] is sigma, Params2[3] is constant c>1 Default = 5,
        Rs2 <- RRnorm2
      }
      
      
      if(i == 1){
        minV = min;
        output[i,] = c(i,round(Rs,3),round(minV,3),round(minV+Rs,3))
      }
      if (i >1)
      {
        minV = as.numeric(output[i-1,4])
        if(max > round(minV+Rs,3)){
          output[i,] = c(i,round(Rs,3),round(minV-Rs2,3),round(minV+Rs,3))
        } else {
          output[i,] = c(i,round(max-round(minV,3),3),round(minV,3),max)
          
          break
        }
      }
      
    }
    
    
  }
  output
}


#' Mapping the functional data into summary statistics by the random splits.  
#'
#' It maps the functional data by the random splits with the summary statistics. 
#' @param RawData        The functional data that defines in the domain. 
#' @param Stat           The summary statistics for each split including "mean","median","Q25","Q75","sd","max","min","range".
#' @param min            The minimum value of the domain.
#' @param max            The maximum value of the domain.
#' @param Params         The parameters of the distribution. 
#'                           In the "Exponential" is (rate, 0,0)
#'                           In the "Uniform" is  (min,max,0)
#'                           In the "Normal"  is  (Mu, Sigma,0)
#' @param Distribution   The Distribution name. ("Exponential", "Uniform" and "Normal").
#' @param Params2        The parameters of the second distribution, if type is "Overlap".
#'                           In the "Exponential" is (rate, Constant,0), Constant default is 5. Final rate is rate*Constant.
#'                           In the "Uniform" is  (min,max,Constant), Constant default is 5. Final min and max are min/Constant and max/Constant.
#'                           In the "Normal"  is  (Mu, Sigma,Constant) , Constant default is 5. Final sd is sd/Constant.
#' @param Distribution2  The Distribution name. ("Exponential", "Uniform" and "Normal"), if type is "Overlap".
#' @param type           The type of splitting, "Disjoint" or "Overlap". 
#' @return It returns a list:
#' @return                num: the number of splits.
#' @return                Rs : the random number.
#' @return                lr : the minimum value of the split.
#' @return                ur : the maximum value of the split.
#' @references
#' \itemize{
#'  \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
#' }
#' @examples
#' nSample_1 <- 10  ##Number of samples Group 1
#' nSample_2 <- 5   ##Number of samples Group 2
#' SIM_1_FUNC_DATA <- makeSIMData(nSS_1 = nSample_1,nSS_2 =nSample_2, responseType ="Categorical",Seed= 2400)
#' matplot(t(SIM_1_FUNC_DATA[[3]]),type="l")  ### Plot functional data number 3.
#' ### Mapping the functioanl data into the mean with the random splits from exponential distribtion and rate = 0.05.
#' set.seed(123)
#' MAP_Mean <- Mapping(RawData=SIM_1_FUNC_DATA[[3]],Stat=c("mean"), min=0,max=100,Params=c(0.05,0,0),Distribution= "Exponential",Params2=c(0,0,0),Distribution2,type="Disjoint")
#' set.seed(123)
#' MAP_Median <- Mapping(RawData=SIM_1_FUNC_DATA[[3]],Stat=c("median"), min=0,max=100,Params=c(0.05,0,0),Distribution= "Exponential",Params2=c(0,0,0),Distribution2,type="Disjoint")
#' set.seed(123)
#' MAP_Q25 <- Mapping(RawData=SIM_1_FUNC_DATA[[3]],Stat=c("Q25"), min=0,max=100,Params=c(0.05,0,0),Distribution= "Exponential",Params2=c(0,0,0),Distribution2,type="Disjoint")
#' set.seed(123)
#' MAP_Q75 <- Mapping(RawData=SIM_1_FUNC_DATA[[3]],Stat=c("Q75"), min=0,max=100,Params=c(0.05,0,0),Distribution= "Exponential",Params2=c(0,0,0),Distribution2,type="Disjoint")
#' set.seed(123)
#' MAP_sd <- Mapping(RawData=SIM_1_FUNC_DATA[[3]],Stat=c("sd"), min=0,max=100,Params=c(0.05,0,0),Distribution= "Exponential",Params2=c(0,0,0),Distribution2,type="Disjoint")
#' set.seed(123)
#' MAP_max <- Mapping(RawData=SIM_1_FUNC_DATA[[3]],Stat=c("max"), min=0,max=100,Params=c(0.05,0,0),Distribution= "Exponential",Params2=c(0,0,0),Distribution2,type="Disjoint")
#' set.seed(123)
#' MAP_min <- Mapping(RawData=SIM_1_FUNC_DATA[[3]],Stat=c("min"), min=0,max=100,Params=c(0.05,0,0),Distribution= "Exponential",Params2=c(0,0,0),Distribution2,type="Disjoint")
#' set.seed(123)
#' MAP_range <- Mapping(RawData=SIM_1_FUNC_DATA[[3]],Stat=c("range"), min=0,max=100,Params=c(0.05,0,0),Distribution= "Exponential",Params2=c(0,0,0),Distribution2,type="Disjoint")
#' par(mfrow=c(3,3))
#' matplot(t(SIM_1_FUNC_DATA[[3]])[,15],type="l",xlab = "Time", main="Real Curves")
#' MINS = min(t(SIM_1_FUNC_DATA[[3]])[,15])
#' MAXS = max(t(SIM_1_FUNC_DATA[[3]])[,15])
#' plot(rnorm(1), ylim=c(MINS,MAXS), xlim=c(0,100),type="n",ylab="",main="Mean")
#' for(i in 1:9){
#'   lines(y = c(MAP_Mean[[1]]$Curve_15[i],MAP_Mean[[1]]$Curve_15[i]) , x = c(MAP_Mean[[1]]$lr[i], MAP_Mean[[1]]$ur[i]))
#' }
#' plot(rnorm(1), ylim=c(MINS,MAXS), xlim=c(0,100),type="n",ylab="",main="Median")
#' for(i in 1:9){
#'   lines(y = c(MAP_Median[[1]]$Curve_15[i],MAP_Median[[1]]$Curve_15[i]) , x = c(MAP_Median[[1]]$lr[i], MAP_Median[[1]]$ur[i]))
#' }
#' plot(rnorm(1), ylim=c(MINS,MAXS), xlim=c(0,100),type="n",ylab="",main="Q25")
#' for(i in 1:9){
#'   lines(y = c(MAP_Q25[[1]]$Curve_15[i],MAP_Q25[[1]]$Curve_15[i]) , x = c(MAP_Q25[[1]]$lr[i], MAP_Q25[[1]]$ur[i]))
#' }
#' plot(rnorm(1), ylim=c(MINS,MAXS), xlim=c(0,100),type="n",ylab="",main="Q75")
#' for(i in 1:9){
#'   lines(y = c(MAP_Q75[[1]]$Curve_15[i],MAP_Q75[[1]]$Curve_15[i]) , x = c(MAP_Q75[[1]]$lr[i], MAP_Q75[[1]]$ur[i]))
#' }
#' plot(rnorm(1), ylim=c(MINS,MAXS), xlim=c(0,100),type="n",ylab="",main="SD")
#' for(i in 1:9){
#'   lines(y = c(MAP_sd[[1]]$Curve_15[i],MAP_sd[[1]]$Curve_15[i]) , x = c(MAP_sd[[1]]$lr[i], MAP_sd[[1]]$ur[i]))
#' }
#' plot(rnorm(1), ylim=c(MINS,MAXS), xlim=c(0,100),type="n",ylab="",main="Max")
#' for(i in 1:9){
#'   lines(y = c(MAP_max[[1]]$Curve_15[i],MAP_max[[1]]$Curve_15[i]) , x = c(MAP_max[[1]]$lr[i], MAP_max[[1]]$ur[i]))
#' }
#' plot(rnorm(1), ylim=c(MINS,MAXS), xlim=c(0,100),type="n",ylab="",main="Min")
#' for(i in 1:9){
#'   lines(y = c(MAP_min[[1]]$Curve_15[i],MAP_min[[1]]$Curve_15[i]) , x = c(MAP_min[[1]]$lr[i], MAP_min[[1]]$ur[i]))
#' }
#' plot(rnorm(1), ylim=c(MINS,MAXS), xlim=c(0,100),type="n",ylab="",main="Range")
#' for(i in 1:9){
#'   lines(y = c(MAP_range[[1]]$Curve_15[i],MAP_range[[1]]$Curve_15[i]) , x = c(MAP_range[[1]]$lr[i], MAP_range[[1]]$ur[i]))
#' }
#' @export
Mapping <- function(RawData,Stat, min,max,Params=c(0,0,0),Distribution,Params2=c(0,0,0),Distribution2,type){
  NPoints <- dim(RawData)[2]
  NCurves <- dim(RawData)[1]
  index   <- seq(from=min,to=max,by=1)[1:NPoints]
  Split_Plan <- Split_Domain(min=min,max=max,Params=Params,Distribution=Distribution ,Params2=Params2,Distribution2=Distribution2,type)
  result <- t(Split_Plan)
  
  for(i in 1:NCurves){
    
    TheCurve <- RawData[i,]
    res <- data.frame(TheCurve,index)
    stat<-0
    
    for(j in 1:dim(Split_Plan)[1]){
      
      SelectedInterval <- res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1]
      
      
      
      if(Stat== "mean" )   { stat[j]   <- mean(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],na.rm =TRUE)}
      
      if(Stat=="median")   { stat[j]   <- median(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],na.rm =TRUE)}
      
      if(Stat=="Q25")      { stat[j]   <- quantile(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],0.25,na.rm =TRUE)}
      
      if(Stat=="Q75")      { stat[j]   <- quantile(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],0.75,na.rm =TRUE)}
      
      if(Stat=="sd")      { stat[j]   <- sd(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],na.rm =TRUE)}
      
      if(Stat=="max")      { stat[j]   <- max(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],na.rm =TRUE)}
      
      if(Stat=="min")      { stat[j]   <- min(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],na.rm =TRUE)}
      
      if(Stat=="range")    { stat[j]   <- max(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],na.rm =TRUE) - min(res[which(res$index >= Split_Plan$lr[j] & res$index < Split_Plan$ur[j]),1],na.rm =TRUE)}
      
    }
    result <- rbind(result,stat)
    
  }
  dim(result)
  result <- data.frame(t(result))
  names(result)[seq(from = 5 ,to = 4 + (NCurves),by=1)] <- paste("Curve_",seq(from=1,to=NCurves,by=1),sep="")
  dim(result)
  
  ### Reordering Data
  result <- result[,c(1:4,seq(from = 5 ,to = 4 + (NCurves*1),by=1))]
  
  result <- result[complete.cases(result),]
  dim(result)
  
  fres <- list(result)
  fres
}


#' Making Final Data set with non-functional covariate, random split functional covariates and response values.   
#'
#' It maps the functional data by the random splits with the summary statistics. 
#' @param ResponseVar    The response variable (Continuous, Categorical or Survival ) 
#' @param Covariates     The non-functional Covariates.
#' @param RawData        The functional data that defines in the domain. 
#' @param Stat           The summary statistics for each split including "mean","median","Q25","Q75","sd","max","min","range".
#' @param min            The minimum value of the domain.
#' @param max            The maximum value of the domain.
#' @param Params         The parameters of the distribution. 
#'                           In the "Exponential" is (rate, 0,0)
#'                           In the "Uniform" is  (min,max,0)
#'                           In the "Normal"  is  (Mu, Sigma,0)
#' @param Distribution   The Distribution name. ("Exponential", "Uniform" and "Normal").
#' @param Params2        The parameters of the second distribution, if type is "Overlap".
#'                           In the "Exponential" is (rate, Constant,0), Constant default is 5. Final rate is rate*Constant.
#'                           In the "Uniform" is  (min,max,Constant), Constant default is 5. Final min and max are min/Constant and max/Constant.
#'                           In the "Normal"  is  (Mu, Sigma,Constant) , Constant default is 5. Final sd is sd/Constant.
#' @param Distribution2  The Distribution name. ("Exponential", "Uniform" and "Normal"), if type is "Overlap".
#' @param type           The type of splitting, "Disjoint" or "Overlap". 
#' @return It returns a list:
#' @return                num: the number of splits.
#' @return                Rs : the random number.
#' @return                lr : the minimum value of the split.
#' @return                ur : the maximum value of the split.
#' @references
#' \itemize{
#'  \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
#' }
#' @examples
#'
#' @export
Complete_DS <- function(ResponseVar,Covariates,RawData,Stat, min,max,Params=c(0,0,0),Distribution,Params2=c(0,0,0),Distribution2,type)
{
  NCurves <- dim(RawData)[1]
  Map_Res0 <- Mapping(RawData=RawData,Stat=Stat, min=min,max=max,Params  = Params,Distribution  =Distribution,
                      Params2 = Params2,Distribution2 =Distribution2,
                      type = type)
  Map_Res <- data.frame(t(Map_Res0[[1]]))
  
  # No covariates
  RNames <- row.names(Map_Res)
  ResponseVar <- data.frame(ResponseVar)
  result <- data.frame(ResponseVar,Covariates,Map_Res )
  
  row.names(result) <- RNames
  if(names(result)[2] == "Covariates") {
    names(result)[2] = "C1"
  }
  
  # result <- data.frame(ResponseVar,Map_Res,TotDev )
  return(result)
  
}


#' Random Splitting Random Forest SRC for Hybrid Data
#'
#' Random splitting random forest for Survival, Regression and Classification (SRC) for hybrid data.
#' @param ResponseVar  The response variable. Add c(rep(0,4),real_response-variable),
#' @param Covariates   The non-functional covariates. Same name : data.frame(c_1 = c(rep(0,4),real_cov_1),c_2 = c(rep(0,4),real_cov_2)),
#' @param RawData      The functional data [rows are number of curves][columns are time points]
#' @param Stat         The summary statistics for the random interval: mean,median,Q25,Q75,sd,max,min,range
#' @param min          The minimum value of the domain of the functional data.
#' @param max          The maximum value of the domain of the functional data.
#' @param Params       The parameters for random number distribution. It includes:
#'                     \itemize{
#'                           \item \code{Exponential} : c(rate,0,0)           rate = 1/theta
#'                           \item \code{Normal}      : c(mean , sd,0)
#'                           \item \code{Uniform}     : c(min,max,0)
#'                           }
#' @param Distribution  The distribution for random number (r*) : Exponential,Normal, Uniform
#' @param Params2       If type is "Overlap",the parameters for random number for distribution 2:
#'                      \itemize{
#'                           \item \code{Exponential} : Params2[1] is rate , Params2[2] is constant c>1. Default = 5
#'                           \item \code{Normal}      : Paramsp2[1] is mu, Paramsp2[2] is sigma,Paramsp2[3] is constant c>1.Default = 5,
#'                           \item \code{Uniform}     : Paramsp2[1] is min, Paramsp2[2] is max,Paramsp2[3] is constant c>1. Default = 5,
#'                           }
#' @param Distribution2 If type is "Overlap",the distribution  :
#'                      Exponential,Normal, Uniform
#' @param type          Two options: 1-Disjoint, 2-Overlap
#' @param m             Number of Forests (default = 100)
#' @param k             Iterations in Each Random Forest (default = 500)
#' @param ResponseType  The Response Type for model predictions: "Categorical", "Continuous", "Multivariate","Unsupervised","Survival" and "Competing Risk".
#' @param SplitRule     The Splitting Rule based on the response variable.
#' @param importance    The variable importance (VIMP) as in rfsrc() from randomForestSRC.
#' @param block.size    The block.size as in rfsrc() from randomForestSRC.
#' @param ensemble      The ensemble as in rfsrc() from randomForestSRC.
#' @param bootstrap     The bootstrap as in rfsrc() from randomForestSRC.
#' @param samptype      The sampling type as in rfsrc() from randomForestSRC.
#' @param sampsize      The sampling size as in rfsrc() from randomForestSRC.
#' @param proximity     The proximity as in rfsrc() from randomForestSRC.
#' @param RFMethod      "Bagging" for (MD-BG) and "RandomForest" for (MD-RF)
#' @param NUM_Covariates Number of Covariets. Each functional covariates is one covariates.
#'                   For example a model with 3 functional covariates and
#'                   2 non-functional covariates have 2+3 = 5 covariates.
#' @return  The output is a list. It contains:
#'     list[[1]]   ## the summary statistics
#'     list[[2]]   ## The variable importance
#' @references
#' \itemize{
#'  \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
#'  \item \code{2-} : Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The annals of applied statistics, 2(3), 841-860.
#'  \item \code{3-} : Ishwaran, H., & Lu, M. (2019). Standard errors and confidence intervals for variable importance in random forest regression, classification, and survival. Statistics in medicine, 38(4), 558-582.
#' }
#' @examples
#' library(RandomForestMixedData)
#'
#'  ######## Simulations
#' #Number of Samples is 100
#' nSample_1 <- 100
#' SIM_1_FUNC_DATA <- makeSIMData(nSS = nSample_1, responseType =c("Survival"),Seed= 2400)
#' par(mfrow=c(2,3))
#' matplot(t(SIM_1_FUNC_DATA[[1]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[1](t)))
#' matplot(t(SIM_1_FUNC_DATA[[2]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[2](t)))
#' matplot(t(SIM_1_FUNC_DATA[[3]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[3](t)))
#' matplot(t(SIM_1_FUNC_DATA[[4]]),type="l",col=c(rep("Green",nSample_1/2),rep("Grey",nSample_1/2)),xlab = expression(t),ylab = expression(X(t)) ,main= expression(X[4](t)))
#' boxplot(SIM_1_FUNC_DATA[[5]][1][-(1:4),][(1:nSample_1/2)], SIM_1_FUNC_DATA[[5]][1][-(1:4),][-(1:nSample_1/2)],col=c("Green","Grey"),main=expression(X[1]))
#' boxplot(SIM_1_FUNC_DATA[[5]][2][-(1:4),][(1:nSample_1/2)], SIM_1_FUNC_DATA[[5]][2][-(1:4),][-(1:nSample_1/2)],col=c("Green","Grey"),main=expression(X[2]))

#' ######## Model Specification and Running
#' ## Survival
#' time  <- c(SIM_1_FUNC_DATA[[6]])
#' event <- c(SIM_1_FUNC_DATA[[7]])
#' Ys <- data.frame(time,event)
#'
#' ## Functional Covariate Specifications
#'
#' ######### With Functional PCA (FPCA)
#' FDList <- list()
#' FDList[[1]] <- UTIL_REC_FPCA(InputFunction = t(SIM_1_FUNC_DATA[[1]]),numberofPCA = 1)
#' FDList[[2]] <- UTIL_REC_FPCA(InputFunction = t(SIM_1_FUNC_DATA[[2]]),numberofPCA = 1)
#' FDList[[3]] <- UTIL_REC_FPCA(InputFunction = t(SIM_1_FUNC_DATA[[3]]),numberofPCA = 1)
#' FDList[[4]] <- UTIL_REC_FPCA(InputFunction = t(SIM_1_FUNC_DATA[[4]]),numberofPCA = 1)
#'
#' #matplot(t(matrix(unlist(FDList[[1]]),ncol=100)),type="l")
#' #matplot(t(matrix(unlist(FDList[[2]]),ncol=100)),type="l")
#' #matplot(t(matrix(unlist(FDList[[3]]),ncol=100)),type="l")
#' #matplot(t(matrix(unlist(FDList[[4]]),ncol=100)),type="l")
#'
#' ######### Without Functional PCA (FPCA)
#'
#' FDList <- list()
#' FDList[[1]] <- SIM_1_FUNC_DATA[[1]]
#' FDList[[2]] <- SIM_1_FUNC_DATA[[2]]
#' FDList[[3]] <- SIM_1_FUNC_DATA[[3]]
#' FDList[[4]] <- SIM_1_FUNC_DATA[[4]]
#'
#' #matplot(t(matrix(unlist(FDList[[1]]),ncol=100)),type="l")
#' #matplot(t(matrix(unlist(FDList[[2]]),ncol=100)),type="l")
#' #matplot(t(matrix(unlist(FDList[[3]]),ncol=100)),type="l")
#' #matplot(t(matrix(unlist(FDList[[4]]),ncol=100)),type="l")
#'
#'  ## Parameter Specification for each Functioanl Covariate
#' Min <- list()
#' Max <- list()
#' PARAMS_1 <- list()
#' DIST_1 <- list()
#' TYPE <- list()
#' PARAMS_2 <- list()
#' STATS <- list()
#' DIST_2 <- list()
#'
#' ### Functional Covariate 1 Specificaion
#' Min[[1]] <- 1                    ## minimum time domain
#' Max[[1]] <- 100                  ## maximum time domain
#' DIST_1[[1]] <- c("Exponential")  ## The Random splitting distribution
#' PARAMS_1[[1]] <- c(0.1,0,0)      ## The parameter for the distribution
#' TYPE[[1]] <- c("Disjoint")       ## Type of splitting
#' DIST_2[[1]] <- c("")             ## if type is overlap , the distribution of the overlap
#' PARAMS_2[[1]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
#' STATS[[1]]<- c("mean")          ## the statistics for each interval
#'
#' ### Functional Covariate 2 Specificaion
#' Min[[2]] <- 1                    ## minimum time domain
#' Max[[2]] <- 100                  ## maximum time domain
#' DIST_1[[2]] <- c("Exponential")  ## The Random splitting distribution
#' PARAMS_1[[2]] <- c(0.1,0,0)      ## The parameter for the distribution
#' TYPE[[2]] <- c("Disjoint")       ## Type of splitting
#' DIST_2[[2]] <- c("")             ## if type is overlap , the distribution of the overlap
#' PARAMS_2[[2]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
#' STATS[[2]]<- c("mean")          ## the statistics for each interval
#'
#' ### Functional Covariate 3 Specificaion
#' Min[[3]] <- 1                    ## minimum time domain
#' Max[[3]] <- 100                  ## maximum time domain
#' DIST_1[[3]] <- c("Exponential")  ## The Random splitting distribution
#' PARAMS_1[[3]] <- c(0.1,0,0)      ## The parameter for the distribution
#' TYPE[[3]] <- c("Disjoint")       ## Type of splitting
#' DIST_2[[3]] <- c("")             ## if type is overlap , the distribution of the overlap
#' PARAMS_2[[3]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
#' STATS[[3]]<- c("mean")          ## the statistics for each interval
#'
#' ### Functional Covariate 4 Specificaion
#' Min[[4]] <- 1                    ## minimum time domain
#' Max[[4]] <- 100                  ## maximum time domain
#' DIST_1[[4]] <- c("Exponential")  ## The Random splitting distribution
#' PARAMS_1[[4]] <- c(0.1,0,0)      ## The parameter for the distribution
#' TYPE[[4]] <- c("Disjoint")       ## Type of splitting
#' DIST_2[[4]] <- c("")             ## if type is overlap , the distribution of the overlap
#' PARAMS_2[[4]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
#' STATS[[4]]<- c("mean")          ## the statistics for each interval
#'
#' ## Non-Functional Covariate Specifications
#' Covariates_all <- data.frame(c_01= SIM_1_FUNC_DATA[[5]][,1],c_02= SIM_1_FUNC_DATA[[5]][,2])
#'
#' ## General Settings
#' mtree <- 100
#' k <- 500
#'
#' ######## Running the model
#' BG_Res <- RS_RF_SRC_H(
#'   ResponseVar=Ys,
#'   Covariates=Covariates_all,
#'   RawData=FDList,
#'   Stat=STATS ,
#'   min=Min,
#'   max=Max,
#'   Params=PARAMS_1,
#'   Distribution=DIST_1,
#'   Params2=PARAMS_2,
#'   Distribution2=DIST_2,
#'   type=TYPE,
#'   m=500,
#'   k=k,
#'   ResponseType="Survival",
#'   RFMethod = c("Bagging"),
#'   NUM_Covariates=1,
#'   Block.Size = 10 ,
#'   SplitRule = c("logrank"),
#'   Importance= TRUE,
#'   Ensemble = c("all"),
#'   Proximity = TRUE)
#' @export
RS_RF_SRC_H <- function(ResponseVar,Covariates,RawData,Stat, min,max,Params=c(0,0,0),Distribution,
                        Params2=c(0,0,0),Distribution2,type,m=100,k=500,
                        ResponseType=c("Categorical", "Continuous", "Multivariate","Unsupervised","Survival","Competing Risk"),
                        SplitRule, Importance,Block.Size, Ensemble,Proximity,
                        RFMethod = c("RandomForest","Bagging"),NUM_Covariates=NUM_Covariates,
                        VIP_Algorithm = c("Normal","Normal_Joint","RS","RS_Joint")){
  #### Calling Libraries
  library(rpart)
  library(data.table)
  library(caret)
  library(e1071)
  library(randomForest)
  library(progress)
  library(randomForestSRC)
  library(survival)
  
  #   m                      ## Number of Forests (default = 100)
  #   NUM_Covariates         ## Number of Covariets
  #   k                     ## Iterations in Each Random Forest (default = 500)
  #   ResponseType           ## The Response Type for model predictions
  
  ####### Making Output Structure
  if(ResponseType == c("Categorical")) result_VIM_F_all <- data.frame(iteration =0,VarF=0,type=0,VIP.all=0,VIP.0=0,VIP.1=0)
  if(ResponseType == c("Continuous"))  result_VIM_F_all <- data.frame(iteration =0,VarF=0,type=0,VIP=0)
  if(ResponseType == c("Survival"))    result_VIM_F_all <- data.frame(iteration =0,VarF=0,type=0,VIP=0)
  
  IMP_Final <-data.frame(iteration = 0, VARS = 0 , VIMP_Train = 0, VIMP_Test = 0, num = 0 ,Rs = 0 , lr = 0, ur =0 )
  
  ## Elapsed Time
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = m,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)      # Width of the progress bar
  ## Learning Steps
  for(h in 1:m){
    #h=1
    Tree_Final_OutPut <- list()
    preds_list <- list()
    
    ####### Step 1 : Curve Splitting Method (RS, RS_Joint, Normal, Normal_Joint)
    
    if (VIP_Algorithm == c("RS")  | VIP_Algorithm == c("RS_Joint")){
      #### Curve Splitting
      NumCurves <- length(RawData)
      if(NumCurves >= 1){
        Comple_DS<- list()
        if(NumCurves < 10){
          CurveNames <- paste("x_0",seq(from=1,to=NumCurves,by=1),".",sep="")
        }else if(NumCurves  >= 10){
          CurveNames_1 <- paste("x_0",seq(from=1,to=9,by=1),".",sep="")
          CurveNames_2 <- paste("x_",seq(from=10,to=NumCurves,by=1),".",sep="")
          CurveNames <- c(CurveNames_1,CurveNames_2)
        }
        
        for(f in 1:NumCurves){
          CurveReady <- Complete_DS(ResponseVar=ResponseVar,Covariates=Covariates,RawData=RawData[[f]],Stat=Stat[[f]],min=Min[[f]],max=Max[[f]],Params=Params[[f]],Distribution=Distribution[[f]],Params2=Params2[[f]],Distribution2=Distribution2[[f]],type=type)
          nameCurve <- names(CurveReady)
          names(CurveReady) <- gsub("X",paste(CurveNames[f]),nameCurve)
          if(f ==1){
            #Comple_DS[[f]]  <- CurveReady[,-dim(CurveReady)[2]]
            Comple_DS[[f]]  <- CurveReady
          }else if(f>1){
            Comple_DS[[f]] <- CurveReady[,grep("x_",names(CurveReady))]  }
        }
      } else{
        CurveReady <- Complete_DS(ResponseVar=ResponseVar,Covariates=Covariates, RawData=RawData[[f]],Stat=Stat[[f]], min=Min[[f]],max=Max[[f]],Params=Params[[f]],Distribution=Distribution[[f]],Params2=Params2[[f]],Distribution2=Distribution2[[f]],type=type)
        Comple_DS  <- CurveReady[,-dim(CurveReady)[2]]
      }
      #print("1_Part 1")
      
      if(ResponseType != "Survival"){
        ### Response
        ResponseVar_INDX  <- which(names(Comple_DS[[1]]) == "ResponseVar")
        ResponseVar <- Comple_DS[[1]][ResponseVar_INDX]
      }else{
        ResponseVar_INDX  <- which(names(Comple_DS[[1]]) == "time" | names(Comple_DS[[1]]) == "event" )
        ResponseVar <- Comple_DS[[1]][ResponseVar_INDX]
      }
      ### Covariates
      Covariates_INDX   <- grep("c_",names(Comple_DS[[1]]))
      Covariates_Value  <- Comple_DS[[1]][grep("c_",names(Comple_DS[[1]]))]
      
      ### Functional
      number_non_f_covariates <- length(Covariates_INDX)
      number_f_covariates     <- length(Comple_DS)
      number_total_covariates <- number_non_f_covariates + number_f_covariates
      
      
      if(RFMethod == "RandomForest"){
        number_covariates_used <- NUM_Covariates # Input
        number_use_covariates <- sort(sample(number_total_covariates,NUM_Covariates,replace = FALSE))
      }
      if(RFMethod == "Bagging"){
        number_covariates_used <- number_total_covariates
        number_use_covariates <- sort(sample(number_total_covariates,number_total_covariates,replace = FALSE))
      }
      
      #print("2_Part 2")
      
      used_non_f_covariates <- used_f_covariates <- used_non_f_covariates_all <- used_f_covariates_all <- 0
      for(i in 1:number_covariates_used)
      {
        if(number_use_covariates[i] <= number_non_f_covariates)
        {
          used_non_f_covariates <- number_use_covariates[i]
          used_non_f_covariates_all <- c(used_non_f_covariates_all,used_non_f_covariates)
        }else{
          used_f_covariates <- number_use_covariates[i]
          used_f_covariates_all <- c(used_f_covariates_all,used_f_covariates)
        }
      }
      
      if(length(used_non_f_covariates_all) >1)  used_non_f_covariates_all <- used_non_f_covariates_all[-1]
      if(length(used_f_covariates_all) >1)      used_f_covariates_all <- used_f_covariates_all[-1]
      
      ####### No Functional Covariate
      if(length(used_f_covariates_all)  == 1 &  used_f_covariates_all == 0)
      {
        used_f_covariates_all <- 0
      }else{
        used_f_covariates_all <- used_f_covariates_all - number_non_f_covariates
      }
      
      Covariates_Value_Used <- Covariates_Value[used_non_f_covariates_all]
      
      if(used_f_covariates_all > 0 )
      {
        ufcl <- length(used_f_covariates_all)
        fdat_all <- 0
        for(i in 1:ufcl){
          if(used_f_covariates_all[i] == 1)
          {
            fdat_all <- data.frame(fdat_all,Comple_DS[[1]][grep("x_01.",names(Comple_DS[[1]]))])
          }else{
            fdat_all <- data.frame(fdat_all,Comple_DS[[used_f_covariates_all[i]]])
          }
        }
        fdat_all <- fdat_all[,-1]
        Comple_DS <- data.frame(ResponseVar,Covariates_Value_Used,fdat_all)
      } else {
        Comple_DS <- data.frame(ResponseVar,Covariates_Value_Used)
      }
      
      AddData   <- Comple_DS[1:4,]
      InputData <- Comple_DS[-c(1:4),]
    }
    
    
    if (VIP_Algorithm == c("Normal")  | VIP_Algorithm == c("Normal_Joint")){
      Comple_DS <- RawData
      if(h == 1){
        if(ResponseType == "Survival"){
          ResponseVar <- ResponseVar[-c(1:4),]
        }else{
          ResponseVar <- ResponseVar[-c(1:4)]
        }
        # print(dim(ResponseVar))
        #Covariates_Value <- Covariates_all[-c(1:4),]
      }else{
        ResponseVar <- ResponseVar
        # Covariates_Value <- Covariates_all
      }
      Covariates_Value <- Covariates_all[-c(1:4),]
      #ResponseVar <- ResponseVar
      
      number_non_f_covariates <- dim(Covariates_all)[2]
      number_f_covariates     <- length(Comple_DS)
      number_total_covariates <- number_non_f_covariates + number_f_covariates
      
      
      if(RFMethod == "RandomForest"){
        number_covariates_used <- NUM_Covariates # Input
        number_use_covariates <- sort(sample(number_total_covariates,NUM_Covariates,replace = FALSE))
      }
      if(RFMethod == "Bagging"){
        number_covariates_used <- number_total_covariates
        number_use_covariates <- sort(sample(number_total_covariates,number_total_covariates,replace = FALSE))
      }
      
      used_non_f_covariates <- used_f_covariates <- used_non_f_covariates_all <- used_f_covariates_all <- 0
      for(i in 1:number_covariates_used){
        if(number_use_covariates[i] <= number_non_f_covariates)
        {
          used_non_f_covariates <- number_use_covariates[i]
          used_non_f_covariates_all <- c(used_non_f_covariates_all,used_non_f_covariates)
        }else{
          used_f_covariates <- number_use_covariates[i]
          used_f_covariates_all <- c(used_f_covariates_all,used_f_covariates)
        }
      }
      
      if(length(used_non_f_covariates_all) >1)  used_non_f_covariates_all <- used_non_f_covariates_all[-1]
      if(length(used_f_covariates_all) >1)      used_f_covariates_all <- used_f_covariates_all[-1]
      
      if(length(used_f_covariates_all)  == 1 &  sum(used_f_covariates_all) == 0)
      {
        used_f_covariates_all <- 0
      }else{
        used_f_covariates_all <- used_f_covariates_all - number_non_f_covariates
      }
      
      Covariates_Value_Used <- Covariates_Value[used_non_f_covariates_all]
      
      if(length(used_f_covariates_all) > 0 & used_f_covariates_all>0 )
      {
        ufcl <- length(used_f_covariates_all)
        fdat_all <- 0
        for(i in 1:ufcl){
          if(used_f_covariates_all[i] == 1)
          {
            fdat_all <- data.frame(fdat_all,Comple_DS[[1]])
            mxnum <- dim(fdat_all)[2]-1
            fdat_all <- fdat_all[,-1]
            names(fdat_all) <- paste("x_01.",1:mxnum,sep="")
          }else{
            if(i < 10 ){
              res_1 <- data.frame(Comple_DS[[used_f_covariates_all[i]]])
              names(res_1) <- paste("x_0",used_f_covariates_all[i],".",1:dim(res_1)[2],sep="")
              fdat_all <- data.frame(fdat_all,res_1)
            }else{
              res_1 <- data.frame(Comple_DS[[used_f_covariates_all[i]]])
              names(res_1) <- paste("x_",used_f_covariates_all[i],".",1:dim(res_1)[2],sep="")
              fdat_all <- data.frame(fdat_all,res_1)
            }
          }
        }
        Comple_DS <- data.frame(ResponseVar,Covariates_Value_Used,fdat_all)
      } else {
        Comple_DS <- data.frame(ResponseVar,Covariates_Value_Used)
      }
      #AddData   <- Comple_DS[1:4,]
      InputData <- Comple_DS
    }
    pb$tick()
    
    #### Train and Test Dataset
    if(ResponseType != "Survival"){
      Train_Rows <- sample(nrow(InputData),nrow(InputData)*0.7,replace = FALSE)
      Train_DS   <- InputData[Train_Rows,]
      Test_DS    <- InputData[-Train_Rows,]
      x_train = Train_DS[,-1]
      y_train = Train_DS[,1]
      x_test  = Test_DS[,-1]
      y_test  = Test_DS[,1]
      Train_DS   <- data.frame(Train_DS)
      Test_DS    <- data.frame(Test_DS)
    }else{
      Train_Rows <- sample(nrow(InputData),nrow(InputData)*0.7,replace = FALSE)
      Train_DS <- as.matrix(InputData[Train_Rows,])
      Test_DS  <- as.matrix(InputData[-Train_Rows,])
      x_train = as.matrix(Train_DS[,-c(1,2)])
      y_train = as.matrix(Train_DS[,c(1,2)])
      x_test  = as.matrix(Test_DS[,-c(1,2)])
      y_test  = as.matrix(Test_DS[,c(1,2)])
      Train_DS <- data.frame(Train_DS)
      Test_DS <- data.frame(Test_DS)
    }
    
    
    if(ResponseType == "Categorical")
    {
      while(length(unique(as.factor(y_train))) != length(unique(as.factor(y_test))))
      {
        Train_Rows <- sample(nrow(InputData),nrow(InputData)*0.7,replace = FALSE)
        Train_DS   <- InputData[Train_Rows,]
        Test_DS    <- InputData[-Train_Rows,]
        x_train = Train_DS[,-1]
        y_train = Train_DS[,1]
        x_test  = Test_DS[,-1]
        y_test  = Test_DS[,1]
        Train_DS   <- data.frame(Train_DS)
        Test_DS    <- data.frame(Test_DS)
      }
      y_train <- as.factor(y_train)
      y_test <- as.factor(y_test)
    }
    
    
    
    if(RFMethod == c("RandomForest"))
    {
      mtry_selected <- NUM_Covariates
    }
    if(RFMethod == c("Bagging"))
    {
      mtry_selected <- dim(x_train)[2]
    }
    
    #### Random Forest (RF-SRC)
    if(ResponseType == "Categorical"){
      rf_model <- rfsrc(ResponseVar ~ . ,
                        data = Train_DS ,
                        ntree = k,
                        mtry=mtry_selected,
                        block.size = Block.Size ,
                        splitrule  = SplitRule,
                        importance = Importance,
                        ensemble = Ensemble,
                        proximity = Proximity,
                        perf.type = "brier")
    }
    
    if(ResponseType == "Continuous"){
      rf_model <- rfsrc(ResponseVar ~ . ,
                        data = Train_DS ,
                        ntree = k,
                        mtry=mtry_selected,
                        block.size = Block.Size ,
                        splitrule  = SplitRule,
                        importance = Importance,
                        ensemble = Ensemble,
                        proximity = Proximity )
    }
    
    if(ResponseType == "Survival")
    {
      rf_model <- rfsrc(Surv(time,event) ~ . ,
                        data = Train_DS ,
                        ntree = k,
                        mtry=mtry_selected,
                        block.size = Block.Size,
                        splitrule  = SplitRule,
                        importance = Importance,
                        ensemble = Ensemble,
                        proximity = Proximity )
    }
    
    
    ### Prediction
    rf_prd <- predict(rf_model,
                      newdata  = Test_DS ,
                      importance =Importance,
                      block.size = Block.Size,
                      ensemble = Ensemble,
                      outcome = c("test"))
    
    Time_Model <- rf_model$ctime.external
    Time_Pred  <- rf_prd$ctime.external
    
    ####### Prediction Error
    if(ResponseType == "Survival")
    {
      ### Harrell’s C-index (1 minus concordance)
      Harrell_C_Index <- get.cindex(time = Train_DS$time, censoring = Train_DS$event, predicted = rf_model$predicted.oob)
      
      MSE_TRAIN <- (sum(as.numeric(rf_model$predicted) - y_train[,1])^2)/length(y_train[,1])
      RSQ_TRAIN <-  cor(as.numeric(rf_model$predicted) , y = as.numeric(y_train[,1]))
      
      MSE_TEST  <- (sum(as.numeric(rf_prd$predicted) - y_test[,1])^2)/length(y_test[,1])
      RSQ_TEST  <- cor(as.numeric(rf_prd$predicted) , y_test[,1])
      
      Values <- c(Harrell_C_Index,MSE_TRAIN,RSQ_TRAIN,MSE_TEST,RSQ_TEST)
      lables <- c(rep("Estimate",5))
      DT_VAL <- c("General",rep(c("Train"),2),rep(c("Test"),2))
      STAT_VAL <- c(rep(c("Harrell_C_Index"),1),rep(c("MSE"),1),rep(c("RSQ"),1),rep(c("MSE"),1),rep(c("RSQ"),1))
      iters <- rep(h,5)
      
      ##### Brier Scores
      Brier_TRAIN_KM    <- get.brier.survival(rf_model, cens.mode="km")
      Brier_TRAIN_KM_Score <- Brier_TRAIN_KM$brier.score
      Brier_TRAIN_KM_CRPS <- Brier_TRAIN_KM$crps
      Brier_TRAIN_KM_CRPS_STD <- Brier_TRAIN_KM$crps.std
      
      Brier_TRAIN_RFSRC    <- get.brier.survival(rf_model, cens.mode="rfsrc")
      Brier_TRAIN_RFSRC_Score <- Brier_TRAIN_KM$brier.score
      Brier_TRAIN_RFSRC_CRPS <- Brier_TRAIN_KM$crps
      Brier_TRAIN_RFSRC_CRPS_STD <- Brier_TRAIN_KM$crps.std
      
      Brier_TEST_KM    <- get.brier.survival(rf_prd, cens.mode="km")
      Brier_TEST_KM_Score <- Brier_TRAIN_KM$brier.score
      Brier_TEST_KM_CRPS <- Brier_TRAIN_KM$crps
      Brier_TEST_KM_CRPS_STD <- Brier_TRAIN_KM$crps.std
      
      Brier_TEST_RFSRC    <- get.brier.survival(rf_prd, cens.mode="rfsrc")
      Brier_TEST_RFSRC_Score <- Brier_TRAIN_KM$brier.score
      Brier_TEST_RFSRC_CRPS <- Brier_TRAIN_KM$crps
      Brier_TEST_RFSRC_CRPS_STD <- Brier_TRAIN_KM$crps.std
      
      
      
      
      if(h ==1){
        result_All <- data.frame(iters=0,DT_VAL=0,STAT_VAL=0,lables=0,Values=0 )
        IMP_Final <- data.frame(iteration=0,VARS=0,VIMP_Train=0,VIMP_Test=0,num=0,Rs=0,lr=0,ur=0)
        Brier_DS_ALL <-data.frame(iters =0 , time = 0,
                                  Brier_TRAIN_KM_Score = 0, Brier_TRAIN_RFSRC_Score = 0,
                                  Brier_TEST_KM_Score = 0,Brier_TEST_RFSRC_Score=0,
                                  Brier_TRAIN_KM_CRPS = 0 , Brier_TRAIN_KM_CRPS_STD = 0 ,
                                  Brier_TRAIN_RFSRC_CRPS = 0 , Brier_TRAIN_RFSRC_CRPS_STD = 0 ,
                                  Brier_TEST_KM_CRPS = 0 , Brier_TEST_KM_CRPS_STD = 0 ,
                                  Brier_TEST_RFSRC_CRPS= 0 , Brier_TEST_RFSRC_CRPS_STD = 0 )
      }
      
      result_1 <- data.frame(iters,DT_VAL,STAT_VAL,lables,Values )
      result_All <- rbind(result_All,result_1)
      ####
      Brier_DS <- data.frame(iters=h,
                             time = Brier_TRAIN_KM_Score[,1],
                             Brier_TRAIN_KM_Score = Brier_TRAIN_KM_Score[,2],
                             Brier_TRAIN_RFSRC_Score = Brier_TRAIN_RFSRC_Score[,2],
                             Brier_TEST_KM_Score = Brier_TEST_KM_Score[,2],
                             Brier_TEST_RFSRC_Score=Brier_TEST_RFSRC_Score[,2],
                             Brier_TRAIN_KM_CRPS = Brier_TRAIN_KM_CRPS,Brier_TRAIN_KM_CRPS_STD = Brier_TRAIN_KM_CRPS_STD,
                             Brier_TRAIN_RFSRC_CRPS = Brier_TRAIN_RFSRC_CRPS, Brier_TRAIN_RFSRC_CRPS_STD=Brier_TRAIN_RFSRC_CRPS_STD,
                             Brier_TEST_KM_CRPS = Brier_TEST_KM_CRPS, Brier_TEST_KM_CRPS_STD=Brier_TEST_KM_CRPS_STD,
                             Brier_TEST_RFSRC_CRPS=Brier_TEST_RFSRC_CRPS, Brier_TEST_RFSRC_CRPS_STD=Brier_TEST_RFSRC_CRPS_STD)
      Brier_DS_ALL <- rbind(Brier_DS_ALL,Brier_DS)
    }
    
    
    #print("4_Part 4")
    
    if(ResponseType == "Continuous")
    {
      MSE_TRAIN <- (sum(as.numeric(rf_model$predicted) - y_train)^2)/length(y_train)
      RSQ_TRAIN <- cor(as.numeric(rf_model$predicted) , y_train)
      
      MSE_TEST  <- (sum(as.numeric(rf_prd$predicted) - y_test)^2)/length(y_test)
      RSQ_TEST  <- cor(as.numeric(rf_prd$predicted) , y_test)
      
      Values <- c(MSE_TRAIN,RSQ_TRAIN,MSE_TEST,RSQ_TEST)
      lables <- c(rep("Estimate",4))
      DT_VAL <- c(rep(c("Train"),2),rep(c("Test"),2))
      STAT_VAL <- c(rep(c("MSE"),1),rep(c("RSQ"),1),rep(c("MSE"),1),rep(c("RSQ"),1))
      iters <- rep(h,4)
      
      if(h ==1){
        result_All <- data.frame(iters=0,DT_VAL=0,STAT_VAL=0,lables=0,Values=0 )
        IMP_Final <- data.frame(iteration=0,VARS=0,VIMP_Train=0,VIMP_Test=0,num=0,Rs=0,lr=0,ur=0)
      }
      
      result_1 <- data.frame(iters,DT_VAL,STAT_VAL,lables,Values )
      result_All <- rbind(result_All,result_1)
    }
    
    if(ResponseType == "Categorical")
    {
      ### AUC and Brier Error
      AUC_Train <- get.auc(Train_DS$ResponseVar, rf_model$predicted.oob)
      Brier_Error_Train <- get.brier.error(Train_DS$ResponseVar, rf_model$predicted.oob)
      Conf_Train <- get.confusion(Train_DS$ResponseVar, rf_model$predicted.oob)
      MisClass_error_Train <- get.misclass.error(Train_DS$ResponseVar, rf_model$predicted.oob)
      
      AUC_Test <- get.auc(Test_DS$ResponseVar, rf_prd$predicted)
      Brier_Error_Test <- get.brier.error(Test_DS$ResponseVar, rf_prd$predicted)
      Conf_Test <- get.confusion(Test_DS$ResponseVar, rf_prd$predicted)
      MisClass_error_Test <- get.misclass.error(Test_DS$ResponseVar, rf_prd$predicted)
      
      Train_PRED <- rep(1,length(rf_model$predicted[,1]))
      Train_PRED[which(rf_model$predicted[,1] >rf_model$predicted[,2] )] <- 0
      Test_PRED <- rep(1,length(rf_prd$predicted[,1]))
      Test_PRED[which(rf_prd$predicted[,1] >rf_prd$predicted[,2] )] <- 0
      Train_CM <- confusionMatrix(factor(Train_PRED), factor(Train_DS$ResponseVar))
      Test_CM  <- confusionMatrix(factor(Test_PRED), factor(Test_DS$ResponseVar))
      CONF_TRAIN <- data.frame(iter= h ,label="Train", ROWS=row.names(Train_CM$table),Train_CM$table)
      CONF_TEST <-  data.frame(iter= h ,label="Test", ROWS=row.names(Test_CM$table),Test_CM$table)
      result_1_CATE <- rbind(CONF_TRAIN, CONF_TEST)
      
      #### CM Statistics
      Train_DF <- data.frame(ref=y_train,pred=Train_PRED)
      Train_CM   <-confusionMatrix(data=as.factor(Train_DF$pred),reference=as.factor(Train_DF$ref),positive	=c("1"))
      Train_Result <- rbind(t(t(Train_CM$overall)),t(t(Train_CM$byClass)))
      Train_Result <- data.frame(iter=h ,name="Train", Statistics =rownames(Train_Result),Values =as.numeric(Train_Result[,1]) )
      Train_Result <-rbind(Train_Result,
                           data.frame(iter=h ,name="Train", Statistics ="Brier Error",Values =as.numeric(Brier_Error_Train)),
                           data.frame(iter=h ,name="Train", Statistics ="AUC",Values =as.numeric(AUC_Train) ),
                           data.frame(iter=h ,name="Train", Statistics ="MisClass_error_1",Values =as.numeric(MisClass_error_Train[1]) ),
                           data.frame(iter=h ,name="Train", Statistics ="MisClass_error_2",Values =as.numeric(MisClass_error_Train[2]) ))
      
      Test_DF <- data.frame(ref=y_test,pred=Test_PRED)
      Test_CM <- confusionMatrix(data=as.factor(Test_DF$pred),reference=as.factor(Test_DF$ref),positive	=c("1"))
      Test_Result <- rbind(t(t(Test_CM$overall)),t(t(Test_CM$byClass)))
      Test_Result <- data.frame(iter=h ,name="Test", Statistics =rownames(Test_Result),Values =as.numeric(Test_Result[,1]) )
      Test_Result <-rbind(Test_Result,
                          data.frame(iter=h ,name="Test", Statistics ="Brier Error",Values =as.numeric(Brier_Error_Test)),
                          data.frame(iter=h ,name="Test", Statistics ="AUC",Values =as.numeric(AUC_Test) ),
                          data.frame(iter=h ,name="Test", Statistics ="MisClass_error_1",Values =as.numeric(MisClass_error_Test[1]) ),
                          data.frame(iter=h ,name="Test", Statistics ="MisClass_error_2",Values =as.numeric(MisClass_error_Test[2]) ))
      
      cStats_result <- rbind(Train_Result,Test_Result)
      
      if(h ==1){
        cTable_All<- result_1_CATE[1,]
        cTable_All$iter<-0
        cStats_result_ALL <-cStats_result[1,]
        cStats_result_ALL$iter <- 0
        IMP_Final <- data.frame(iteration=0,VARS=0,VIMP_Train=0,VIMP_Test=0,num=0,Rs=0,lr=0,ur=0)
      }
      
      cTable_All <- rbind(cTable_All,result_1_CATE)
      cStats_result_ALL <- rbind(cStats_result_ALL,cStats_result)
    }
    
    
    ############# Variable Importance Plot
    if(VIP_Algorithm == "RS"){
      if(ResponseType == "Categorical"){
        ### Variable Importance
        result_VIM_F_all <- data.frame(i=0,
                                       VarF=0,
                                       type=0,
                                       VIP.all=0,
                                       VIP.0=0,
                                       VIP.1=0)
        IMP_Type_Train <- rf_model$importance[,1]
        IMP_Type_Test  <- rf_prd$importance[,1]
      }else{
        IMP_Type_Train <- rf_model$importance
        IMP_Type_Test  <- rf_prd$importance
      }
      
      if(is.null(IMP_Type_Test) == FALSE){
        
        IMP_Var_DS <- data.frame(VARS       = names(IMP_Type_Train),
                                 VIMP_Train = IMP_Type_Train,
                                 VIMP_Test  = IMP_Type_Test)
      }else{
        IMP_Var_DS <- data.frame(VARS       = names(IMP_Type_Train),
                                 VIMP_Train = IMP_Type_Train,
                                 VIMP_Test  = 0)
      }
      
      if(ResponseType != "Survival"){
        AddData_DS <- data.frame(VARS=row.names(t(AddData)),t(AddData))[-1,]
      }else{
        AddData_DS <- data.frame(VARS=row.names(t(AddData)),t(AddData))[-c(1:2),]
      }
      
      IMP_Var_DS$VARS  <- AddData_DS$VARS
      IMP_Merge <- merge(x=IMP_Var_DS, y = AddData_DS , by = "VARS")
      names(IMP_Merge) <- c("VARS","VIMP_Train","VIMP_Test","num","Rs","lr","ur")
      IMP_Merge <- data.frame(iteration=h,IMP_Merge)
      IMP_Final <- rbind(IMP_Final,IMP_Merge)
      #print(paste("Number of Forest: ",h))
    }
    if(VIP_Algorithm == "RS_Joint"){
      #if(ResponseType == "Categorical"){
      
      Covariates_name <- rf_model$xvar.names
      result_VIM_F  <- 0
      result_VIM_NF <- 0
      
      Num_Non_FUNC <- Covariates_name[grep(paste("c_",sep=""),Covariates_name)]
      Num_FUNC <- Covariates_name[grep(paste("x_",sep=""),Covariates_name)]
      
      if(length(Num_FUNC) > 0){
        
        for(i in 1:length(used_f_covariates_all))
        {
          if(i < 10){
            VIM_F <- data.frame(iteration = h , VarF = paste("x_0",used_f_covariates_all,sep="")[i] ,type="Functional", VIP = vimp(rf_model, Covariates_name[grep(paste("x_0",used_f_covariates_all[i],sep=""),Covariates_name)],joint=TRUE)$importance)
            result_VIM_F <- rbind(result_VIM_F,VIM_F)
          }
          if(i >= 10){
            VIM_F <- data.frame(iteration = h , VarF = paste("x_0",used_f_covariates_all,sep="")[i] ,type="Functional", VIP = vimp(rf_model, Covariates_name[grep(paste("x_",used_f_covariates_all[i],sep=""),Covariates_name)],joint=TRUE)$importance)
            result_VIM_F <- rbind(result_VIM_F,VIM_F)
          }
        }
      }else{
        result_VIM_F <- 0
      }
      
      if(length(Num_Non_FUNC) > 0){
        
        for(i in 1:length(used_non_f_covariates_all))
        {
          if(i < 10){
            VIM_NF <- data.frame(iteration = h , VarF = paste("c_0",used_non_f_covariates_all,sep="")[i] ,type="Non-Functional", VIP = vimp(rf_model, Covariates_name[grep(paste("c_0",used_non_f_covariates_all[i],sep=""),Covariates_name)],joint=TRUE)$importance)
            result_VIM_NF <- rbind(result_VIM_NF,VIM_NF)
          }
          if(i >= 10){
            VIM_NF <- data.frame(iteration = h , VarF = paste("c_0",used_non_f_covariates_all,sep="")[i] ,type="Non-Functional", VIP = vimp(rf_model, Covariates_name[grep(paste("c_",used_non_f_covariates_all[i],sep=""),Covariates_name)],joint=TRUE)$importance)
            result_VIM_NF <- rbind(result_VIM_NF,VIM_NF)
          }
        }
      }else{
        result_VIM_NF <- 0
      }
      
      result_VIM_FF <- rbind(result_VIM_F,result_VIM_NF)
      result_VIM_F_all <- rbind(result_VIM_F_all,result_VIM_FF)
      #  }
    }
    
    
    if(VIP_Algorithm == "Normal"){
      if(ResponseType == "Categorical"){
        ### Variable Importance
        IMP_Type_Train <- rf_model$importance[,1]
        IMP_Type_Test  <- rf_prd$importance[,1]
      }else{
        IMP_Type_Train <- rf_model$importance
        IMP_Type_Test  <- rf_prd$importance
      }
      
      if(is.null(IMP_Type_Test)== FALSE){
        
        IMP_Final_c <- data.frame(iteration = h  ,
                                  VARS = names(IMP_Type_Train) ,
                                  VIMP_Train = IMP_Type_Train,
                                  VIMP_Test = IMP_Type_Test,
                                  num = 0 ,
                                  Rs = 0 ,
                                  lr = 0,
                                  ur =0 )
      }else{
        IMP_Final_c <- data.frame(iteration = h  ,
                                  VARS = names(IMP_Type_Train) ,
                                  VIMP_Train = IMP_Type_Train,
                                  VIMP_Test = 0,
                                  num = 0 ,
                                  Rs = 0 ,
                                  lr = 0,
                                  ur =0 )
      }
      IMP_Final <- rbind(IMP_Final_c,IMP_Final)
      
    }
    if(VIP_Algorithm == "Normal_Joint"){
      #  if(ResponseType == "Categorical"){
      
      Covariates_name <- rf_model$xvar.names
      result_VIM_F  <- 0
      result_VIM_NF <- 0
      
      Num_Non_FUNC <- Covariates_name[grep(paste("c_",sep=""),Covariates_name)]
      Num_FUNC <- Covariates_name[grep(paste("x_",sep=""),Covariates_name)]
      
      if(length(Num_FUNC) > 0){
        
        for(i in 1:length(used_f_covariates_all))
        {
          if(i < 10){
            VIM_F <- data.frame(iteration = h , VarF = paste("x_0",used_f_covariates_all,sep="")[i] ,type="Functional", VIP = vimp(rf_model, Covariates_name[grep(paste("x_0",used_f_covariates_all[i],sep=""),Covariates_name)],joint=TRUE)$importance)
            result_VIM_F <- rbind(result_VIM_F,VIM_F)
          }
          if(i >= 10){
            VIM_F <- data.frame(iteration = h , VarF = paste("x_0",used_f_covariates_all,sep="")[i] ,type="Functional", VIP = vimp(rf_model, Covariates_name[grep(paste("x_",used_f_covariates_all[i],sep=""),Covariates_name)],joint=TRUE)$importance)
            result_VIM_F <- rbind(result_VIM_F,VIM_F)
          }
        }
      }else{
        result_VIM_F <- 0
      }
      
      if(length(Num_Non_FUNC) > 0){
        
        for(i in 1:length(used_non_f_covariates_all))
        {
          if(i < 10){
            VIM_NF <- data.frame(iteration = h , VarF = paste("c_0",used_non_f_covariates_all,sep="")[i] ,type="Non-Functional", VIP = vimp(rf_model, Covariates_name[grep(paste("c_0",used_non_f_covariates_all[i],sep=""),Covariates_name)],joint=TRUE)$importance)
            result_VIM_NF <- rbind(result_VIM_NF,VIM_NF)
          }
          if(i >= 10){
            VIM_NF <- data.frame(iteration = h , VarF = paste("c_0",used_non_f_covariates_all,sep="")[i] ,type="Non-Functional", VIP = vimp(rf_model, Covariates_name[grep(paste("c_",used_non_f_covariates_all[i],sep=""),Covariates_name)],joint=TRUE)$importance)
            result_VIM_NF <- rbind(result_VIM_NF,VIM_NF)
          }
        }
      }else{
        result_VIM_NF <- 0
      }
      
      result_VIM_FF <- rbind(result_VIM_F,result_VIM_NF)
      result_VIM_F_all <- rbind(result_VIM_F_all,result_VIM_FF)
      #  }
      
    }
    #print(h)
    ######## Finishing the Loop
  }
  
  if(VIP_Algorithm == "RS" | VIP_Algorithm == "Normal")
  {
    IMP_Final$COV_Type <- as.character(substr(IMP_Final$VARS,start=1,stop=1))
    IMP_Final$COV_NUM <- as.numeric(substr(IMP_Final$VARS,start=3,stop=4))
    IMP_Final$SplitNumberF <- as.numeric(substr(IMP_Final$VARS,start=6,stop=10))
    IMP_Final <- IMP_Final[order(IMP_Final$iteration,IMP_Final$COV_Type,IMP_Final$COV_NUM,IMP_Final$SplitNumberF),]
    
    
    if(ResponseType != "Categorical")
    {
      outputs <- list(Result = result_All[-1,] , IMP_VALUE = IMP_Final[-1,])
    }else{
      cStats_result_ALL <- cStats_result_ALL[-1,]
      cTable_All <- cTable_All[-1,]
      result_All <- list(CTable = cTable_All,CStats = cStats_result_ALL)
      outputs <- list(Result = result_All , IMP_VALUE = IMP_Final[-1,])
    }
  }
  
  if(VIP_Algorithm == "RS_Joint" | VIP_Algorithm == "Normal_Joint")
  {
    result_VIM_F_0 <- result_VIM_F_all
    DimSize <- dim(result_VIM_F_0)[1]
    result_VIM_F_0$COV_Type <- 0
    
    FuncDS <- subset(result_VIM_F_0,type==c("Functional"))
    if(dim(FuncDS)[1] > 0 ){
      FuncDS$COV_Type <-  FuncDS$VarF
      #FuncDS$COV_NUM <- used_f_covariates_all
    }
    
    NonFuncDS <- subset(result_VIM_F_0,type==c("Non-Functional"))
    if(dim(NonFuncDS)[1] > 0 ){
      NonFuncDS$COV_Type <- NonFuncDS$VarF
      #NonFuncDS$COV_NUM <-  used_non_f_covariates_all
    }
    IMP_Final <- rbind(NonFuncDS ,FuncDS )
    IMP_Final <- IMP_Final[order(IMP_Final$iteration),]
    
    if(ResponseType != "Categorical")
    {
      outputs <- list(Result = result_All[-1,] , IMP_VALUE = IMP_Final)
    }else{
      cStats_result_ALL <- cStats_result_ALL[-1,]
      cTable_All <- cTable_All[-1,]
      result_All <- list(CTable = cTable_All,CStats = cStats_result_ALL)
      outputs <- list(Result = result_All , IMP_VALUE = IMP_Final)
    }
  }
  
  if(ResponseType == "Survival")
  {
    outputs <- list(Result = outputs$Result,
                    IMP_VALUE = outputs$IMP_VALUE,
                    Brier_DS = Brier_DS_ALL)
  }
  return(outputs)
}







#' Summary of Random Splitting Random Forest SRC for Hybrid Data
#'
#' The summary of Bagging (MD-BG)  and Random Forest  (MD-RF)
#' @param Dataset         The 'RF_OUTPUT$Result' from `RS_RF_SRC_H`.
#' @param Outcome_Type    The output type : `Continuous`, `Categorical` and `Survival`.
#' @details               It summarizes the random forest output.
#' @return                It produces the model performance based on the type of outcome.
#'                        \item Categorical
#'                        \item \item Accuracy (Mean, SD, Q1, Q3) for Train and Test.
#'                        \item \item Sensitivity (Mean, SD, Q1, Q3) for Train and Test.
#'                        \item \item Specificity (Mean, SD, Q1, Q3) for Train and Test.
#'                        \item Continuous
#'                        \item \item MSE (Mean, SD) for Train and Test.
#'                        \item \item RQR (Mean, SD) for Train and Test.
#'                        \item Survival
#'                        \item \item MSE (Mean, SD) for Train and Test.
#'                        \item \item RQR (Mean, SD) for Train and Test.
#'                        \item \item Harrell_C_Index (Mean, SD) for Train and Test.
#' @references
#' \itemize{
#'  \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
#'  \item \code{2-} : Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The annals of applied statistics, 2(3), 841-860.
#'  \item \code{3-} : Ishwaran, H., & Lu, M. (2019). Standard errors and confidence intervals for variable importance in random forest regression, classification, and survival. Statistics in medicine, 38(4), 558-582.
#' }
#' @examples
#'  #Simulation
#'   nSample_1 <- nSample_2 <- 100
#'   SIM_1_FUNC_DATA <- makeSIMData(nSS_1 = nSample_1,nSS_2 = nSample_2, responseType =c("Survival"),Seed= 2400)
#'   ## Survival
#'   time  <- c(SIM_1_FUNC_DATA[[6]])
#'   event <- c(SIM_1_FUNC_DATA[[7]])
#'   Ys <- data.frame(time,event)
#'    FDList <- list()
#'    FDList[[1]] <- SIM_1_FUNC_DATA[[1]]
#'    FDList[[2]] <- SIM_1_FUNC_DATA[[2]]
#'    FDList[[3]] <- SIM_1_FUNC_DATA[[3]]
#'    FDList[[4]] <- SIM_1_FUNC_DATA[[4]]
#'    ## Parameter Specification for each Functioanl Covariate
#'    Min <- list()
#'    Max <- list()
#'    PARAMS_1 <- list()
#'    DIST_1 <- list()
#'    TYPE <- list()
#'    PARAMS_2 <- list()
#'    STATS <- list()
#'    DIST_2 <- list()
#'    ### Functional Covariate 1 Specificaion
#'    Min[[1]] <- 1                    ## minimum time domain
#'    Max[[1]] <- 100                  ## maximum time domain
#'    DIST_1[[1]] <- c("Exponential")  ## The Random splitting distribution
#'    PARAMS_1[[1]] <- c(0.1,0,0)      ## The parameter for the distribution
#'    TYPE[[1]] <- c("Disjoint")       ## Type of splitting
#'    DIST_2[[1]] <- c("")             ## if type is overlap , the distribution of the overlap
#'    PARAMS_2[[1]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
#'    STATS[[1]]<- c("mean")          ## the statistics for each interval
#'    ### Functional Covariate 2 Specificaion
#'    Min[[2]] <- 1                    ## minimum time domain
#'    Max[[2]] <- 100                  ## maximum time domain
#'    DIST_1[[2]] <- c("Exponential")  ## The Random splitting distribution
#'    PARAMS_1[[2]] <- c(0.1,0,0)      ## The parameter for the distribution
#'    TYPE[[2]] <- c("Disjoint")       ## Type of splitting
#'    DIST_2[[2]] <- c("")             ## if type is overlap , the distribution of the overlap
#'    PARAMS_2[[2]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
#'    STATS[[2]]<- c("mean")          ## the statistics for each interval
#'    ### Functional Covariate 3 Specificaion
#'    Min[[3]] <- 1                    ## minimum time domain
#'    Max[[3]] <- 100                  ## maximum time domain
#'    DIST_1[[3]] <- c("Exponential")  ## The Random splitting distribution
#'    PARAMS_1[[3]] <- c(0.1,0,0)      ## The parameter for the distribution
#'    TYPE[[3]] <- c("Disjoint")       ## Type of splitting
#'    DIST_2[[3]] <- c("")             ## if type is overlap , the distribution of the overlap
#'    PARAMS_2[[3]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
#'    STATS[[3]]<- c("mean")          ## the statistics for each interval
#'    ### Functional Covariate 4 Specificaion
#'    Min[[4]] <- 1                    ## minimum time domain
#'    Max[[4]] <- 100                  ## maximum time domain
#'    DIST_1[[4]] <- c("Exponential")  ## The Random splitting distribution
#'    PARAMS_1[[4]] <- c(0.1,0,0)      ## The parameter for the distribution
#'    TYPE[[4]] <- c("Disjoint")       ## Type of splitting
#'    DIST_2[[4]] <- c("")             ## if type is overlap , the distribution of the overlap
#'    PARAMS_2[[4]] <- c(0,0,0)        ## if type is overlap , the parameter of the distribution
#'    STATS[[4]]<- c("mean")          ## the statistics for each interval
#'    ## Non-Functional Covariate Specifications
#'    Covariates_all <- data.frame(c_01= SIM_1_FUNC_DATA[[5]][,1],c_02= SIM_1_FUNC_DATA[[5]][,2])
#'    ## Categorical Response 
#'    BG_Res <- RS_RF_SRC_H(
#'    ResponseVar=as.factor(Ys[,2]),
#'    Covariates=Covariates_all,
#'    RawData=FDList,
#'    Stat=STATS ,
#'    min=Min,
#'    max=Max,
#'    Params=PARAMS_1,
#'    Distribution=DIST_1,
#'    Params2=PARAMS_2,
#'    Distribution2=DIST_2,
#'    type=TYPE,
#'    m=100,
#'    k=50,
#'    ResponseType="Categorical",
#'    RFMethod = c("Bagging"),
#'    NUM_Covariates=1,
#'    Block.Size = 10 ,
#'    SplitRule = c("gini"),
#'    Importance= TRUE,
#'    Ensemble = c("all"),
#'    Proximity = TRUE,
#'    VIP_Algorithm = c("RS"))
#'    ### Summary 
#'    RS_Summary(Dataset =BG_Res$Result , Outcome_Type = c("Categorical"))
#'    
#'    ## Continuous Response 
#'    BG_Res <- RS_RF_SRC_H(
#'    ResponseVar=Ys[,1],
#'    Covariates=Covariates_all,
#'    RawData=FDList,
#'    Stat=STATS ,
#'    min=Min,
#'    max=Max,
#'    Params=PARAMS_1,
#'    Distribution=DIST_1,
#'    Params2=PARAMS_2,
#'    Distribution2=DIST_2,
#'    type=TYPE,
#'    m=100,
#'    k=50,
#'    ResponseType="Continuous",
#'    RFMethod = c("Bagging"),
#'    NUM_Covariates=1,
#'    Block.Size = 10 ,
#'    SplitRule = c("mse"),
#'    Importance= TRUE,
#'    Ensemble = c("all"),
#'    Proximity = TRUE,
#'    VIP_Algorithm = c("RS"))
#'    ### Summary 
#'    RS_Summary(Dataset =BG_Res$Result , Outcome_Type = c("Continuous"))
#'   
#'    ## Survival Response 
#'    BG_Res <- RS_RF_SRC_H(
#'    ResponseVar=Ys,
#'    Covariates=Covariates_all,
#'    RawData=FDList,
#'    Stat=STATS ,
#'    min=Min,
#'    max=Max,
#'    Params=PARAMS_1,
#'    Distribution=DIST_1,
#'    Params2=PARAMS_2,
#'    Distribution2=DIST_2,
#'    type=TYPE,
#'    m=100,
#'    k=50,
#'    ResponseType="Survival",
#'    RFMethod = c("Bagging"),
#'    NUM_Covariates=1,
#'    Block.Size = 10 ,
#'    SplitRule = c("logrank"),
#'    Importance= TRUE,
#'    Ensemble = c("all"),
#'    Proximity = TRUE,
#'    VIP_Algorithm = c("RS"))
#'    ### Summary 
#'    RS_Summary(Dataset =BG_Res$Result , Outcome_Type = c("Survival"))
#' @export
RS_Summary <- function(Dataset =RF_OUTPUT, Outcome_Type = c("Continuous","Categorical","Survival"))
{
  if(Outcome_Type == c("Continuous"))
  {
    Train_Data_MSE <- subset(Dataset,  DT_VAL == "Train" & STAT_VAL == "MSE" & lables == "Estimate" ,select=c("Values") )
    Train_Data_RSQ <- subset(Dataset,  DT_VAL == "Train" & STAT_VAL == "RSQ" & lables == "Estimate" ,select=c("Values") )
    Test_Data_MSE  <- subset(Dataset,  DT_VAL == "Test"  & STAT_VAL == "MSE" & lables == "Estimate" ,select=c("Values") )
    Test_Data_RSQ  <- subset(Dataset,  DT_VAL == "Test"  & STAT_VAL == "RSQ" & lables == "Estimate" ,select=c("Values") )
    
    Train_MSE_Mean <- mean(unlist(Train_Data_MSE),na.rm = TRUE)
    Train_MSE_SD <- sd(unlist(Train_Data_MSE),na.rm = TRUE)
    Train_RSQ_Mean <- mean(unlist(Train_Data_RSQ),na.rm = TRUE)
    Train_RSQ_SD <- sd(unlist(Train_Data_RSQ),na.rm = TRUE)
    
    Test_MSE_Mean <- mean(unlist(Test_Data_MSE),na.rm = TRUE)
    Test_MSE_SD <- sd(unlist(Test_Data_MSE),na.rm = TRUE)
    Test_RSQ_Mean <- mean(unlist(Test_Data_RSQ),na.rm = TRUE)
    Test_RSQ_SD <- sd(unlist(Test_Data_RSQ),na.rm = TRUE)
    
    Final_result <- data.frame(Train_MSE_Mean,Train_MSE_SD,
                               Train_RSQ_Mean,Train_RSQ_SD,
                               Test_MSE_Mean,Test_MSE_SD,Test_RSQ_Mean,Test_RSQ_SD)
    
    return(Final_result)
    
  }
  
  if(Outcome_Type == c("Categorical"))
  {
    
    Train_ACC  <- as.numeric(unlist(subset(Dataset$CStats,   name == "Train"  & Statistics == "Accuracy"  ,select=c("Values") )))
    Test_ACC   <- as.numeric(unlist(subset(Dataset$CStats,   name == "Test"  & Statistics == "Accuracy"  ,select=c("Values") )))
    Train_SEN  <- as.numeric(unlist(subset(Dataset$CStats,   name == "Train"  & Statistics == "Sensitivity"  ,select=c("Values") )))
    Test_SEN   <- as.numeric(unlist(subset(Dataset$CStats,   name == "Test"  & Statistics == "Sensitivity"  ,select=c("Values") )))
    Train_SPE  <- as.numeric(unlist(subset(Dataset$CStats,   name == "Train"  & Statistics == "Specificity"  ,select=c("Values") )))
    Test_SPE   <- as.numeric(unlist(subset(Dataset$CStats,   name == "Test"  & Statistics == "Specificity"  ,select=c("Values") )))
    
    Train_ACC_Mean <- mean(Train_ACC,na.rm = TRUE)
    Train_ACC_Q1   <- quantile(Train_ACC,c(0.25),na.rm = TRUE)
    Train_ACC_Q3   <- quantile(Train_ACC,c(0.75),na.rm = TRUE)
    Train_ACC_sd   <- sd(Train_ACC,na.rm = TRUE)
    Test_ACC_Mean  <- mean(Test_ACC,na.rm = TRUE)
    Test_ACC_Q1    <-  quantile(Test_ACC,c(0.25),na.rm = TRUE)
    Test_ACC_Q3    <-  quantile(Test_ACC,c(0.75),na.rm = TRUE)
    Test_ACC_sd    <- sd(Test_ACC,na.rm = TRUE)
    
    Train_SEN_Mean <- mean(Train_SEN,na.rm = TRUE)
    Train_SEN_Q1   <- quantile(Train_SEN,c(0.25),na.rm = TRUE)
    Train_SEN_Q3   <- quantile(Train_SEN,c(0.75),na.rm = TRUE)
    Train_SEN_sd   <- sd(Train_SEN,na.rm = TRUE)
    Test_SEN_Mean  <- mean(Test_SEN,na.rm = TRUE)
    Test_SEN_Q1    <- quantile(Test_SEN,c(0.25),na.rm = TRUE)
    Test_SEN_Q3    <- quantile(Test_SEN,c(0.75),na.rm = TRUE)
    Test_SEN_sd    <- sd(Test_SEN,na.rm = TRUE)
    
    Train_SPE_Mean <- mean(Train_SPE,na.rm = TRUE)
    Train_SPE_Q1   <- quantile(Train_SPE,c(0.25),na.rm = TRUE)
    Train_SPE_Q3   <- quantile(Train_SPE,c(0.75),na.rm = TRUE)
    Train_SPE_sd   <- sd(Train_SPE,na.rm = TRUE)
    Test_SPE_Mean  <- mean(Test_SPE,na.rm = TRUE)
    Test_SPE_Q1    <- quantile(Test_SPE,c(0.25),na.rm = TRUE)
    Test_SPE_Q3    <- quantile(Test_SPE,c(0.75),na.rm = TRUE)
    Test_SPE_sd    <- sd(Test_SPE,na.rm = TRUE)
    
    Final_result <- data.frame(Train_ACC_Mean,Train_ACC_sd,Train_ACC_Q1,Train_ACC_Q3,
                               Test_ACC_Mean,Test_ACC_sd,Test_ACC_Q1,Test_ACC_Q3,
                               Train_SEN_Mean,Train_SEN_sd,Train_SEN_Q1,Train_SEN_Q3,
                               Test_SEN_Mean,Test_SEN_sd,Test_SEN_Q1,Test_SEN_Q3,
                               Train_SPE_Mean,Train_SPE_sd,Train_SPE_Q1,Train_SPE_Q3,
                               Test_SPE_Mean,Test_SPE_sd,Test_SPE_Q1,Test_SPE_Q3)
    
    return(Final_result)
    
  }
  
  if(Outcome_Type == c("Survival"))
  {
    Train_Data_MSE  <- subset(Dataset,  DT_VAL == "Train" & STAT_VAL == "MSE" & lables == "Estimate" ,select=c("Values") )
    Train_Data_RSQ  <- subset(Dataset,  DT_VAL == "Train" & STAT_VAL == "RSQ" & lables == "Estimate" ,select=c("Values") )
    Test_Data_MSE   <- subset(Dataset,  DT_VAL == "Test"  & STAT_VAL == "MSE" & lables == "Estimate" ,select=c("Values") )
    Test_Data_RSQ   <- subset(Dataset,  DT_VAL == "Test"  & STAT_VAL == "RSQ" & lables == "Estimate" ,select=c("Values") )
    Harrell_C_Index <- subset(Dataset,  DT_VAL == "General"  & STAT_VAL == "Harrell_C_Index" & lables == "Estimate" ,select=c("Values") )
    
    
    Train_MSE_Mean <- mean(unlist(Train_Data_MSE),na.rm = TRUE)
    Train_MSE_SD <- sd(unlist(Train_Data_MSE),na.rm = TRUE)
    Train_RSQ_Mean <- mean(unlist(Train_Data_RSQ),na.rm = TRUE)
    Train_RSQ_SD <- sd(unlist(Train_Data_RSQ),na.rm = TRUE)
    
    Test_MSE_Mean <- mean(unlist(Test_Data_MSE),na.rm = TRUE)
    Test_MSE_SD <- sd(unlist(Test_Data_MSE),na.rm = TRUE)
    Test_RSQ_Mean <- mean(unlist(Test_Data_RSQ),na.rm = TRUE)
    Test_RSQ_SD <- sd(unlist(Test_Data_RSQ),na.rm = TRUE)
    
    Harrell_C_Index_Mean <- mean(Harrell_C_Index$Values)
    Harrell_C_Index_SD  <- sd(Harrell_C_Index$Values)
    
    Final_result <- data.frame(Train_MSE_Mean,Train_MSE_SD,
                               Train_RSQ_Mean,Train_RSQ_SD,
                               Test_MSE_Mean,Test_MSE_SD,Test_RSQ_Mean,Test_RSQ_SD,
                               Harrell_C_Index_Mean,Harrell_C_Index_SD)
    
    return(Final_result)
    
  }
}



#' The Variable Importance Plot for Survival Random Splitting Random Forest
#'
#' The Variable Importance Plot for functional and non-functional covariates
#' @param  Dataset           It is `RF_OUTPUT$IMP_VALUE` from `RS_RF_SRC_H`.
#' @param  VIP_Algorithm     Only works for RS and 
#' @param  COVAR_TYPE        Type of covariate: `Functional` or `Non-Functional`.
#' @param  COV_NUMBER        Functional covariate index :`1`, `2`, `3` and etc.
#' @param  Type              It has two options `Train` and `Test`.
#' @details                  The `type` is `Train` or `Test` (Train= in train dataset , 2=in test dataset)
#' @return                   The variable importance plots for functional and non-functional covariates.
#' @references
#' \itemize{
#'  \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
#'  \item \code{2-} : Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The annals of applied statistics, 2(3), 841-860.
#'  \item \code{3-} : Ishwaran, H., & Lu, M. (2019). Standard errors and confidence intervals for variable importance in random forest regression, classification, and survival. Statistics in medicine, 38(4), 558-582.
#' }
#' @examples
#' #RF_Result from previous example.
#' par(mfrow=c(5,2))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=1,Type=c("Train"))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=1,Type=c("Test"))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=2,Type=c("Train"))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=2,Type=c("Test"))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=3,Type=c("Train"))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=3,Type=c("Test"))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=4,Type=c("Train"))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Functional"), COV_NUMBER=4,Type=c("Test"))
#' RS_VARIMP_SRC_PLOT(Dataset = BG_Res$IMP_VALUE,VIP_Algorithm  =  c("RS"), COVAR_TYPE = c("Non-Functional"), COV_NUMBER=1,Type=c("Train"))
#' @export
RS_VARIMP_SRC_PLOT <- function(Dataset =RF_OUTPUT , VIP_Algorithm  =  c("RS"),COVAR_TYPE = c("Functional","Non-Functional"),COV_NUMBER=3, Type=c("Train","Test"))
{
  if(VIP_Algorithm  ==  c("RS") )
  {
    if(COVAR_TYPE == c("Functional"))
    {
      Analysis_Data <- subset(Dataset, COV_Type == "x" & COV_NUM == COV_NUMBER )
      Domain_T <- seq(from=min(Analysis_Data$lr),to=max(Analysis_Data$lr),by=0.01)
      Iter_name <- unique(Analysis_Data$iteration)
      Iter_number <- length(Iter_name)
      IMP_Matrix <- matrix(0,ncol=Iter_number+1,nrow = length(Domain_T))
      IMP_Matrix[,1] <- Domain_T
      
      for(i in 1:Iter_number)
      {
        Iter_Data <- subset(Analysis_Data, iteration == Iter_name[i])
        NSplit <- dim(Iter_Data)[1]
        
        for(j in 1:NSplit)
        {
          lr_filter <- Iter_Data[j,]$lr
          ur_filter <- Iter_Data[j,]$ur
          rn_filter <- which(IMP_Matrix[,1]>= lr_filter & IMP_Matrix[,1]<ur_filter )
          if(Type==c("Train")) IMP_Matrix[rn_filter,(i+1)] <- Iter_Data[j,]$VIMP_Train
          if(Type==c("Test")) IMP_Matrix[rn_filter,(i+1)] <- Iter_Data[j,]$VIMP_Test
        }
      }
      
      VAR_IMP_Curve_VALUE  <- rowMeans(IMP_Matrix[,-1])
      VAR_IMP_Curve_50  <- apply(IMP_Matrix[,-1],1,quantile, probs=c(0.50),na.rm=TRUE)
      VAR_IMP_Curve_25     <- apply(IMP_Matrix[,-1],1,quantile, probs=c(0.25),na.rm=TRUE)
      VAR_IMP_Curve_75     <- apply(IMP_Matrix[,-1],1,quantile, probs=c(0.75),na.rm=TRUE)
      VAR_IMP_Curve_Domain <- IMP_Matrix[,1]
      result<- data.frame(VAR_IMP_Curve_Domain,VAR_IMP_Curve_VALUE,VAR_IMP_Curve_50,VAR_IMP_Curve_25,VAR_IMP_Curve_75)
      ymax <- max(VAR_IMP_Curve_75) + 1*sd(VAR_IMP_Curve_75)
      ymin <- min(VAR_IMP_Curve_25) - 1*sd(VAR_IMP_Curve_25)
      plot(x=result$VAR_IMP_Curve_Domain, y=result$VAR_IMP_Curve_VALUE,type="l",ylab="",xlab="Domain",main=paste("The VI Plot: Curve",COV_NUMBER),sub=paste(Type),lty=1,lwd=2,col="Black",ylim=c(ymin,ymax))
      lines(x=result$VAR_IMP_Curve_Domain, y=result$VAR_IMP_Curve_50,type="l",col="Red",lwd=2,lty=3)
      lines(x=result$VAR_IMP_Curve_Domain, y=result$VAR_IMP_Curve_25,type="l",col="Blue",lwd=1,lty=2)
      lines(x=result$VAR_IMP_Curve_Domain, y=result$VAR_IMP_Curve_75,type="l",col="Blue",lwd=1,lty=2)
      abline(h=0)
      
    }
    
    if(COVAR_TYPE == c("Non-Functional"))
    {
      Analysis_Data <- subset(Dataset, COV_Type == "c" )
      
      
      if(Type==c("Train"))
      {
        Mean_Value     <- aggregate(Analysis_Data$VIMP_Train, list(Analysis_Data$COV_NUM), mean)
        Quantile_value <- data.frame(aggregate(Analysis_Data$VIMP_Train, list(Analysis_Data$COV_NUM),quantile, probs=c(0.25,0.50,0.75)))
        Result <- data.frame(Quantile_value$Group.1,Quantile_value$x,AVG=Mean_Value$x)
        names(Result) <- c("Covariates","Q_25","Q_50","Q_75","Mean")
      }
      if(Type==c("Test"))
      {
        Mean_Value     <- aggregate(Analysis_Data$VIMP_Test, list(Analysis_Data$COV_NUM), mean)
        Quantile_value <- data.frame(aggregate(Analysis_Data$VIMP_Test, list(Analysis_Data$COV_NUM),quantile, probs=c(0.25,0.50,0.75)))
        Result <- data.frame(Quantile_value$Group.1,Quantile_value$x,AVG=Mean_Value$x)
        names(Result) <- c("Covariates","Q_25","Q_50","Q_75","Mean")
      }
      
      ymin <- min(Result$Q_25) - 1*sd(Result$Q_25)
      ymax <- max(Result$Q_75) + 1*sd(Result$Q_75)
      
      plot(x=factor(Result$Covariates), y = Result$Mean, ylim=c(ymin, ymax),ylab="",xlab="Covariates",main="Varibale Importance Plot",sub=Type)
      plot(x=factor(Result$Covariates), y = Result$Q_50,col="Blue",lty=2,add=TRUE)
      plot(x=factor(Result$Covariates), y = Result$Q_25,col="Blue",lty=3,add=TRUE)
      plot(x=factor(Result$Covariates), y = Result$Q_75,col="Blue",lty=3,add=TRUE)
      abline(h=0)
    }
  }else{
    print("Under Development!")
  }
}



#' Plot the Brier Score
#'
#' Plot the Brier Score for the Survival Response based on the two methods Kaplan-Meier censoring distribution (KM) and RSF censoring distribution (RFSRC).
#' @param  Dataset           It is `RF_OUTPUT$Brier_DS` from `RS_RF_SRC_H`.
#' @param  method            It has two options: "KM","RFSRC".
#' @param  dataDivid         It has two options `Train` and `Test`.
#' @details                  The `type` is `Train` or `Test` (Train= in train dataset , 2=in test dataset)
#' @return                   The Brier score plot and the summary statistics. 
#'                           the mean and SD of The continuous rank probability score (CRPS) (min, Q1, median, mean, Q3 and max).
#' @references
#' \itemize{
#'  \item \code{1-} : Möller, A., Tutz, G., & Gertheiss, J. (2016). Random forests for functional covariates. Journal of Chemometrics, 30(12), 715-725.
#'  \item \code{2-} : Ishwaran, H., Kogalur, U. B., Blackstone, E. H., & Lauer, M. S. (2008). Random survival forests. The annals of applied statistics, 2(3), 841-860.
#'  \item \code{3-} : Ishwaran, H., & Lu, M. (2019). Standard errors and confidence intervals for variable importance in random forest regression, classification, and survival. Statistics in medicine, 38(4), 558-582.
#' }
#' @examples
#' #RF_Result from previous example.
#' par(mfrow=c(2,2))
#' RS_Brier_Surv(dataset = BG_Res$Brier_DS,  method = c("KM"), dataDivid = c("TRAIN"))
#' RS_Brier_Surv(dataset = BG_Res$Brier_DS,  method = c("KM"), dataDivid = c("TEST"))  
#' RS_Brier_Surv(dataset = BG_Res$Brier_DS,  method = c("RFSRC"), dataDivid = c("TRAIN"))
#' RS_Brier_Surv(dataset = BG_Res$Brier_DS,  method = c("RFSRC"), dataDivid = c("TEST")) 
#' @export
RS_Brier_Surv <- function(dataset = DS,  method = c("KM","RFSRC"), dataDivid = c("TRAIN","TEST"))
{
  dataset <- dataset[-1,]
  TimeCol <- dataset$time
  names_ds <- names(dataset)
  Stati_names <- names_ds[grepl(method,names_ds)]
  select_name <- Stati_names[grepl(dataDivid,Stati_names)]
  Dataset_Select <- dataset[select_name]
  DS_Fin <- data.frame(TimeCol,Dataset_Select)
  DS_Fin <- data.table(DS_Fin)
  names(DS_Fin) <- c("Time","Score","CRPS","CRSP_STD")
  plot(DS_Fin[,mean(Score),by=Time], xlab = c("Time"), ylab = paste(dataDivid, "-" , method),lwd=3,main="Average")
  CRPS_SUM     <- summary(DS_Fin$CRPS)
  CRPS_SUM_STD <- summary(DS_Fin$CRSP_STD)
  result <- rbind(CRPS_SUM, CRPS_SUM_STD)
  return(result)
}





