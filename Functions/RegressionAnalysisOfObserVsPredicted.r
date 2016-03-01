
require('smatr') # these test if the slope = siginificantly different to 1 and the intercept is sig diff to 0 


# this function does regression analysis of the performance of predictive model
# see  Pineiro et al. 2008 ,Ecological modelling 216 (2008) 316-322 and refs therein & Nash and Sutcliffe (1970)  &   Moriasi et al 2007
# the plot is presented with Observed (Y) in the y-axis vs. predicted (YHat)in the x-axis (OP) regressions
do.lm <- function(x=NULL,  y=x$observations, yhat=x$predictions, doLegend=F, method = "OLS", alpha = 0.05, Ylab=NULL, YhatLab=NULL, ...) {
 
  if(!is.null(x)) {
    y <- x[,1]             # observations
    yhat <- x[,2]          # predictions
  }
  
  if(any(is.null(c(Ylab, YhatLab)))) {
  YhatLab <-  "Predicted"
  Ylab <- "Observed"
  } 
 
  # remove NA values from either obs or cals
  BothRealNumbers <- !is.na(y) & !is.na(yhat)
  y <- y[BothRealNumbers]
  yhat <- yhat[BothRealNumbers]
  n <- length(yhat)
  
  if(n <= 1) {
    r2.reg <- NA
    NSE  <- NA
    RMSD <- NA
    MAD <- NA
    SSPE <- NA
    BETA <- NA
    Interc<-NA
    bias <- NA
    pbias  <- NA
    Ubias <- NA
    Ubeta <- NA
    Ue  <- NA
    RSR <- NA
    SigSlope  <- NA
    SigInter <- NA
    slope  <- NA
    inter  <- NA
    Pslope  <- NA
    Pinter  <- NA
  } else {
    SSPE <- sum((y-yhat)^2)  # sum of squared predictive errors a measure of goodness of fit.
    OBS <- mean(y)
    PRE <- mean(yhat)
    bias <- mean(y-yhat)     # bias measures the average tendency of the predicted values to be larger or smaller than the observed
                             # Positive values indicate the model underestimates and negative values indicate overestimation bias (see Moriasi et al 2007)
    # Percent bias (PBIAS) measures the average tendency of the simulated data to be larger or smaller than their observed counterparts (Gupta et al., 1999 in Moriasi et al 2007).
    pbias <- (sum(y - yhat))*100 / sum(y)
    regress <- lm(y~yhat)
    BETA <- as.numeric(regress$coefficients[2])    # the slope of the regression
    Interc <- as.numeric(regress$coefficients[1])    # the intercept of the regression
    est <- as.numeric(predict(regress))            # the fitted values
    # these are the constituent components of the bias known as Theil's partial inequality coefficients (Ubias, Uslope and Uerror),
    # see refs in Pineiro et al. 2008 ,Ecological modelling 216 (2008) 316-322
    Ubias <- (n*((OBS-PRE)^2))/SSPE                    # bias proportion
    Ubeta <- (((BETA-1)^2) *sum((yhat-PRE)^2))/SSPE    # slope proportion CONSISTENCY
    Ue <- (sum((est-y)^2))/SSPE
    MySlopeTest <- slope.test(y=y, x=yhat, test.value = 1, data=NULL, method = method, alpha = alpha, intercept = TRUE)
    MyElevTest  <- elev.test (y=y, x=yhat, test.value = 0,            method = method, alpha = alpha)
    CIslope <- as.numeric(MySlopeTest$ci)
    CIinter <- as.numeric(MyElevTest$a.ci)
    Pslope <- as.numeric(MySlopeTest$p)
    Pinter <- as.numeric(MyElevTest$p)
    slope = ifelse(CIslope[1]<1 & CIslope[2]>1, "NotSigDiffTo1",  "SigDifTto1")
    inter = ifelse(CIinter[1]<0 & CIinter[2]>0, "NotSigDiffTo0",  "SigDifTto0")
    r2.reg <- summary(regress)$r.squared
    NSE <- 100*(1 -  sum((y-yhat)^2)/sum((y-OBS)^2))  # Nash-Sutcliffe efficiency (NSE): The Nash-Sutcliffe efficiency (Nash and Sutcliffe, 1970).
                                                # NSE indicates how well the plot of observed versus simulated/predicted  data fits the 1:1 line.
                                                # this is what Piñeiro et al., 2008 refer to (confusingly) as r2.
    RMSD <- sqrt(1/(n-1) * (sum((yhat-y)^2)))
    MAD <- median(abs(yhat-y))
 
    RSR <- sqrt(sum((y - yhat)^2) ) / sqrt(sum((y - OBS)^2) )
 
     minT <- min(c(y,yhat, (max(c(y,yhat)))))  # to make pretty plots
     maxT<- max(c(y,yhat,(min(c(y,yhat)))))
 
     plot(y=y, x=yhat, ylab=Ylab, xlab=YhatLab, cex.lab=1.3, #bty="n",  Observed (Y) in the y-axis vs. predicted (YHat)in the x-axis (OP) regressions
         xlim = c(minT, maxT), ylim= c(minT, maxT), ...)
     lines(x = c(minT, maxT), y= c(minT, maxT), type = "l", lwd=2, lty=2, col= "red")  # the one to one line (1:1 line)
     points(x=yhat, y=est, type="l", lwd=2, lty=1, col="blue" )  # plot the regression
     if(doLegend == T) legend("bottomright",  legend=c(paste("NSE = ", round(NSE,2), "%"), paste("RMSD =",round(RMSD,2)), paste("Bias =", round(bias,2))), cex = 1.1, bty="n", inset =0.0001)
 
  }
  out <-   data.frame(n = n,
                r2.reg = round(r2.reg,2),
                NSE = round(NSE,2),
                RMSD = round(RMSD,3),
                      Interc = round(Interc,3),
                MAD =  round(MAD,3),
                SSPE = round(SSPE,2),
                BETA = round(BETA,2),
                bias = round(bias,3),
                pbias = round(pbias, 3),
                Ubias = round(Ubias,3),
                Ubeta =  round(Ubeta,3),
                Ue = round(Ue,2),
                RSR = round(RSR, 3),
                SigSlope = slope,
                SigInter = inter,
                Pslope = round(Pslope,3),
                Pinter = round(Pinter,3))
 
  return(out)
}
#########

