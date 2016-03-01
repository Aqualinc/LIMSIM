# Tons Functions for RF models (June 2012)
# Several functions to assist with using Random Forest models, 
# the functions are particularly aimed at analysis of multivariate and univariate spatial data (where objects have spatial coordinates)

# MapVar                =	 Map a variable
# map.class             =  Map a classification
# classify.domain       =  define a classification from an input dataset at levels of detail specified
# resp.mod.fit          =  fits individual RF models to each colom of a response dataframe (representing taxa or any other response)
# jackRF                =  jack Knife proceedure for individual Random Forest models 
# var.order.rf          =  puts the predictors in the decreasing order mportance
# BVpartialAnal         =  evaluates a bivariate partial dependence for variables in a random forest model 
# BVPartPlot            =  plots  a bivariate partial dependence for variables in a random forest model     
# ProbPartialPlot       =  partial dependence for variables in a classification model where the y scale is probability (not logistic)       
# CalcPartialDependence =  computes data for multiple PartialDependence plots (with multiple response)
# PlotPartialDependence =  PLOTs the multiple PartialDependence plots  from the returned object from function CalcPartialDependence
# MultipleImpPlot       =  Produces multiple importance plots from a list of RF models  and a table of the same
# varSelRF.reg          =  modified version of varSelRF from package::varSelRF, variable selection for continuous response variables (regression)
# varSelRF.RegPlot      =  plots a varSelRF.reg object
# varSelRFBoot.reg      =  modifed  varSelRFBoot to boot the varSelRF.reg function


require('randomForest')

var.order.rf <- function(x, MyType=1) { # x is a RF model, MyTYpe=1 -> MeanDecreaseAccuracy, MyType = 2 -> mean decrease in node impurity
# puts the predictors in the decreasing order mportance based on Mean Decrease Gini index (classification) or IncNodePurity (regression)
if(x$type == "classification" & MyType == 1) measure <- "MeanDecreaseAccuracy"
if(x$type == "regression" & MyType == 1) measure <- "%IncMSE"
if(x$type == "classification" & MyType == 2) measure <- "MeanDecreaseGini"
if(x$type == "regression" & MyType == 2) measure <- "IncNodePurity"
import <- importance(x)[,measure]
imp <- sort(import, decreasing = T)
imp.df <- as.data.frame(matrix(nrow=length(imp), ncol=2))
names(imp.df) <- c("Variable", paste("Importance -", measure))
imp.df[,1] <- names(imp)
imp.df[,2] <- as.numeric(imp)
imp.df
}

##################################################################################################### 
#  jackRF jack Knife proceedure
##################################################################################################### 
 
jackRF <- function(MyYVar, MyXPreds, ...) {
    n.sites <- nrow(MyXPreds) 
    results <- as.data.frame(matrix(NA, nrow=n.sites, ncol=2))
    names(results) <- c("observed", "predicted")
    for (i in 1:n.sites) {
        hold <- 1:n.sites == i
        this.RF <- randomForest(x = MyXPreds[!hold, ], y=MyYVar[!hold])
        this.pred <- predict(this.RF, newdata=MyXPreds[hold, ], type="response")
        results[i, ]   <- c(MyYVar[hold] ,this.pred) 
        cat("row= ", i, "observed = ", MyYVar[hold] ,"predicted=", this.pred,  "\n" )           
    }
    return(results)
}


# ------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------

  #getS3method("partialPlot","randomForest") <- STARTING CODE FROM HERE
BVpartialAnal <- function (x, pred.data, x.var, y.var, which.class=F, plot = TRUE, add = FALSE, w,
    n.pt = min(length(unique(pred.data[, xname])), 51), rug = TRUE,
    xlab = deparse(substitute(x.var)), ylab = "", main = paste("Partial Dependence on",
        deparse(substitute(x.var))), ...)
{
 #this function evaluates a bivariate partial dependence for variables in a random forest model
# implemented for continuous predictors and categorical (CLASSES) and regression (continuous) response


  sol <- c(call = match.call())
  class(sol) <- "RFpartplotBV"
  classRF <- x$type != "regression"  # is the RF model Classification or Regression?
  if (is.null(x$forest))
        stop("The randomForest object must contain the forest", "\n")

    x.var <- substitute(x.var)
    y.var <- substitute(y.var)

    xname <- if (is.character(x.var))  { x.var
   } else {
        if (is.name(x.var))
            deparse(x.var)
        else {
            eval(x.var)
        }
    }

    yname <- if (is.character(y.var))  { y.var
    } else {
        if (is.name(y.var))
            deparse(y.var)
        else {
            eval(y.var)
        }
    }


    xv <- pred.data[, xname]    # variable 1 being analysed
    yv <- pred.data[, yname]    # variable 2 being analysed
    n <- nrow(pred.data)
    
if(classRF==TRUE) {

           if (which.class==F) {
            classes <- colnames(x$votes)
            N.classes <- length(classes)
           }   else {
           N.classes<-1  # only one class
              if (missing(w))              # determine the focus class or focus on 1st class
                  w <- rep(1, n)
              if (classRF) {
                  if (missing(which.class)) {
                      focus <- 1
                  }
                  else {
                      focus <- charmatch(which.class, colnames(x$votes))
                      if (is.na(focus))
                          stop(which.class, "is not one of the class labels.")     # determine the focus class or focus on 1st class
                  }
              }
             }



       x.pt <- seq(min(xv), max(xv), length = n.pt)   # values to test??
       y.pt <- seq(min(yv), max(yv), length = n.pt)

      vector("list", length=N.classes)

    z.values <-   array(data = NA, dim = c(n.pt, n.pt, N.classes), dimnames = NULL)    # response array


         # y.pt <- numeric(n.pt*n.pt)

        for (j in 1:n.pt) {
        for (i in 1:n.pt) {      # for each value in x.pt

            x.data <- pred.data
            x.data[, xname]<-rep(x.pt[i], nrow(pred.data))    # replace the
            x.data[, yname]<- rep(y.pt[j], nrow(pred.data))   # replace the values of the variable to be tesed with the values in x.pt

           # if (classRF) {

                pr <- predict(x, x.data, type = "prob")
##########################################################################################
# NOTE i SIMPLIFIED AND MODIFIED HERE!!  This does not allow prior weighting of the probabilties
# the original univariate function does(USING argument "w" )
##########################################################################################

               for (k in 1:N.classes) {                      # do for all classes
                z.values[i,j,k] <- mean(pr[,k])
                }  # end loop for each class
                
            } # end the i loop
            } # end the j loop
            

if (which.class==F) {            
sol$class.names <- classes
}
sol$focus <- which.class
 }  # end the classification RF model analysis

 if(classRF==FALSE) {  # for a regression RF model

       x.pt <- seq(min(xv), max(xv), length = n.pt)   # values to test??
       y.pt <- seq(min(yv), max(yv), length = n.pt)

    z.values <-   array(data = NA, dim = c(n.pt, n.pt, 1), dimnames = NULL)    # response array

        for (j in 1:n.pt) {
        for (i in 1:n.pt) {      # for each value in x.pt

            x.data <- pred.data
            x.data[, xname]<-rep(x.pt[i], nrow(pred.data))    # replace the
            x.data[, yname]<- rep(y.pt[j], nrow(pred.data))   # replace the values of the variable to be tesed with the values in x.pt

                pr <- predict(x, x.data, type = "response")
                z.values[i,j,1] <- mean(pr)

         } # end the i loop
         } # end the j loop
            
} # end the regression RF model analysis

    sol$x.pt <-  x.pt
    sol$y.pt <-  y.pt
    sol$z.values <-  z.values
    sol$xname <- xname
    sol$yname <- yname
    sol$classRF <- x$type # keep the type of model for plotting
    sol
}
# ------------------------------------------------------------------------------------------

  BVPartPlot <-  function(x,              # the partial dependance analysis from BVpartialAnal
                          class.plot=1,          #  the class to be plotted (these are stored in array z.values
                          theta=140,          #  rotation
                          phi=15,          #  elevation
                          ticktype="simple",         # type of axies ticke can be  "simple" or "detailed"
                          MyMain= "BV plot")        #   
                          {

 # plot function for the bivariate partial dependance analysis
 require("graphics")
 classRF <- x$classRF 
  
     if(classRF != "regression") {
     MyZ <- x$z.values[,,class.plot]
     MyZlab <- paste("Partial probability of class ", class.plot) 
     } else {
     MyZ <-  x$z.values[,,1]
     MyZlab <- "Partial dependence of response"  
     }                                     

persp(x = x$x.pt, y = x$y.pt, z=MyZ , theta=theta, phi=phi, ticktype=ticktype, main=MyMain, xlab=x$xname, ylab=x$yname, zlab=MyZlab)
     
    } # end plot

####################################################################################################################################################
# UNIvariate Partial Plots where response is in PROBABILTIE NOT LOGIT scale
####################################################################################################################################################
#
ProbPartialPlot <- function (x, pred.data, x.var, which.class, w, plot = TRUE, add = FALSE,   # response is in PROBABILTIE NOT LOGIT scale
    n.pt = min(length(unique(pred.data[, xname])), 51), rug = TRUE,
    xlab = deparse(substitute(x.var)), ylab = "", main = paste("Partial Dependence on",
        deparse(substitute(x.var))), ...)
        
        # for classification this returns the probabilty not the logit.
{
    classRF <- x$type != "regression"
    if (is.null(x$forest))
        stop("The randomForest object must contain the forest.\n")

        
        if(is.character(x.var)) {
         xname <- x.var 
         } else {        
    x.var <- substitute(x.var)
    }
    xname <- if (is.character(x.var))
        x.var
    else {
        if (is.name(x.var))
            deparse(x.var)
        else {
            eval(x.var)
        }
    }
    xv <- pred.data[, xname]
    n <- nrow(pred.data)
    if (missing(w))
        w <- rep(1, n)
    if (classRF) {
        if (missing(which.class)) {
            focus <- 1
        }
        else {
            focus <- charmatch(which.class, colnames(x$votes))
            if (is.na(focus))
                stop(which.class, "is not one of the class labels.")
        }
    }
    if (is.factor(xv) && !is.ordered(xv)) {
        x.pt <- levels(xv)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
            if (classRF) {            
                pr <- predict(x, x.data, type = "prob")
                y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] >
                  0, pr[, focus], 1)) - rowMeans(log(ifelse(pr >
                  0, pr, 1))), w, na.rm = TRUE)
            }
            else y.pt[i] <- weighted.mean(predict(x, x.data),
                w, na.rm = TRUE)
        }
        if (add) {
            points(1:length(x.pt), y.pt, type = "h", lwd = 2,
                ...)
        }
        else {
            if (plot)
                barplot(y.pt, width = rep(1, length(y.pt)), col = "blue",
                  xlab = xlab, ylab = ylab, main = main, names.arg = x.pt,
                  ...)
        }
    }
    else {
        if (is.ordered(xv))
            xv <- as.numeric(xv)
        x.pt <- seq(min(xv), max(xv), length = n.pt)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- rep(x.pt[i], n)
            

            if (classRF) {
                pr <- predict(x, x.data, type = "prob")
              # doctored from here
                  LOGIT <- weighted.mean(log(ifelse(pr[, focus] ==
                  0, 1, pr[, focus])) - rowMeans(log(ifelse(pr ==
                  0, 1, pr))), w, na.rm = TRUE)
                        
                y.pt[i]  <- exp(LOGIT)/(1+exp(LOGIT))   # return the probabilty by inverse Logit 
            }
            else {
                y.pt[i] <- weighted.mean(predict(x, x.data),
                  w, na.rm = TRUE)
            }
        }
        if (add) {
            lines(x.pt, y.pt, ...)
        }
        else {
            if (plot)
                plot(x.pt, y.pt, type = "l", xlab = xlab, ylab = ylab,
                  main = main, ...)
        }
        if (rug && plot) {
            if (n.pt > 10) {
                rug(quantile(xv, seq(0.1, 0.9, by = 0.1)), side = 1)
            }
            else {
                rug(unique(xv, side = 1))
            }
        }
    }
    
    if (n.pt > 10) {   # to keep the RUG data
                ForRug <-  quantile(xv, seq(0.1, 0.9, by = 0.1))
            }
            else {
                ForRug <-  unique(xv, side = 1)
                }
           
    invisible(list(x = x.pt, y = y.pt, ForRug = ForRug, Xlab=xlab, Ylab=ylab, Main = main))
}



######################################################################################################
# this function simply modified version of partialPlt to ensure a name (character) can be given to specify a variable  
MyPartialPlot <- function (x, pred.data, x.var, which.class, w, plot = TRUE, add = FALSE, 
    n.pt = min(length(unique(pred.data[, xname])), 51), rug = TRUE, 
    xlab = deparse(substitute(x.var)), ylab = "", main = paste("Partial Dependence on", 
        deparse(substitute(x.var))), ...) 
{
    classRF <- x$type != "regression"
    if (is.null(x$forest)) 
        stop("The randomForest object must contain the forest.\n")
        
if(is.character(x.var)) {                           # the function was modified here 
         xname <- x.var                             # the function was modified here 
         } else {                                   # the function was modified here 
    x.var <- substitute(x.var)                      # the function was modified here 
    }
    xname <- if (is.character(x.var))
        x.var
    else {
        if (is.name(x.var)) 
            deparse(x.var)
        else {
            eval(x.var)
        }
    }
    xv <- pred.data[, xname]
    n <- nrow(pred.data)
    if (missing(w)) 
        w <- rep(1, n)
    if (classRF) {
        if (missing(which.class)) {
            focus <- 1
        }
        else {
            focus <- charmatch(which.class, colnames(x$votes))
            if (is.na(focus)) 
                stop(which.class, "is not one of the class labels.")
        }
    }
    if (is.factor(xv) && !is.ordered(xv)) {
        x.pt <- levels(xv)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- factor(rep(x.pt[i], n), levels = x.pt)
            if (classRF) {
                pr <- predict(x, x.data, type = "prob")
                y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] > 
                  0, pr[, focus], .Machine$double.eps)) - rowMeans(log(ifelse(pr > 
                  0, pr, .Machine$double.eps))), w, na.rm = TRUE)
            }
            else y.pt[i] <- weighted.mean(predict(x, x.data), 
                w, na.rm = TRUE)
        }
        if (add) {
            points(1:length(x.pt), y.pt, type = "h", lwd = 2, 
                ...)
        }
        else {
            if (plot) 
                barplot(y.pt, width = rep(1, length(y.pt)), col = "blue", 
                  xlab = xlab, ylab = ylab, main = main, names.arg = x.pt, 
                  ...)
        }
    }
    else {
        if (is.ordered(xv)) 
            xv <- as.numeric(xv)
        x.pt <- seq(min(xv), max(xv), length = n.pt)
        y.pt <- numeric(length(x.pt))
        for (i in seq(along = x.pt)) {
            x.data <- pred.data
            x.data[, xname] <- rep(x.pt[i], n)
            if (classRF) {
                pr <- predict(x, x.data, type = "prob")
                y.pt[i] <- weighted.mean(log(ifelse(pr[, focus] == 
                  0, .Machine$double.eps, pr[, focus])) - rowMeans(log(ifelse(pr == 
                  0, .Machine$double.eps, pr))), w, na.rm = TRUE)
            }
            else {
                y.pt[i] <- weighted.mean(predict(x, x.data), 
                  w, na.rm = TRUE)
            }
        }
        if (add) {
            lines(x.pt, y.pt, ...)
        }
        else {
            if (plot) 
                plot(x.pt, y.pt, type = "l", xlab = xlab, ylab = ylab, 
                  main = main, ...)
        }
        if (rug && plot) {
            if (n.pt > 10) {
                rug(quantile(xv, seq(0.1, 0.9, by = 0.1)), side = 1)
            }
            else {
                rug(unique(xv, side = 1))
            }
        }
    }
    
      if (n.pt > 10) {   # to keep the RUG data
                ForRug <-  quantile(xv, seq(0.1, 0.9, by = 0.1))
            }
            else {
                ForRug <-  unique(xv, side = 1)
                } 
                
    invisible(list(x = x.pt, y = y.pt, ForRug = ForRug, Xlab=xlab, Ylab=ylab, Main = main))
}


#################################################################################################
#################################################################################################

###########################################################################################################

MultipleImpPlot <- function (Model, sort = TRUE,  NoVars= 10,   # X is either a list of models or a single model object
                              ImpType = 2, class = NULL, scale = TRUE, main = deparse(substitute(x)), ...)
# NoVars = max to include on the plot    
# MultipleImpPlot =   Produces multiple importance plots from a list of RF models
{

    if(class(Model)=="list") {
                MyNames <- names(Model)
                names(MyNames) <- MyNames
        } else {
                Model <- list(Model)         # put a single model into a list of 1
                MyNames <- "Response"
                names(MyNames) <- MyNames
                names(Model) <-  MyNames
        }
        
Nmods <- length(MyNames)
n.var = min(NoVars, nrow(Model[[1]]$importance))
VarNames <- row.names(Model[[1]]$importance)

    if (!inherits(Model[[1]], "randomForest"))
        stop("This function only works for objects of class `randomForest'")
        
         if (Nmods > 1) {
        op <- par(mfrow = c(1, Nmods), mar = c(1, 1, 1, 1)) # mgp = c(2, 0.8, 0), oma = c(0, 0, 2, 0), no.readonly = TRUE

        on.exit(par(op))
    }
    
  imp <- lapply(Model, function(x) importance(x, class = class, scale = scale, type = ImpType)) # list of model importances
  OrderByModel <- list("vector")
        
for (d in 1:Nmods) {
        ord <- if (sort)
            rev(order(imp[[d]], decreasing = TRUE)[1:n.var])
        else 1:n.var
       if (colnames(imp[[d]]) %in% c("IncNodePurity",  "MeanDecreaseGini")) {
           xmin <-   0
           }  else {
            xmin <- min(imp[[d]][ord])
           }           
      OrderByModel[[d]] <-  match(VarNames, var.order.rf(Model[[d]], MyType=ImpType)[,1])
      
      dotchart(imp[[d]][ord], xlab = colnames(imp[[d]]), labels = row.names(imp[[d]])[ord],
                  main = MyNames[d], xlim = c(xmin, max(imp[[d]])), ...)
                  
}
    
    
TabOrderByModel <-  data.frame(VarNames, OrderByModel)
names(TabOrderByModel) <- c("Variable",MyNames)   
   
invisible(TabOrderByModel)  # return the tabulated orders
   
}

#######################################################################################
# this function calculates data for multiple PartialDependence plots  
#######################################################################################

CalcPartialDependence <- function(Model,              # this needs to be either a list of models or a single model
                                  MyPreds = Preds,   # the predictors
                                  varOrder = NULL,    # specify the variable plotting order (if NULL will be based on importance
                                  ImpType = 2,        #  the importance type
                                  Nvars = 8,          # how many variables to compute Partial Dependence for
                                  Myn.pt = 11,        # number of points to compute on varible gradients
                                  MyClass = 1,        # the class if the model is a classification model, if class= "ALL" plot responses for all classses
                                  Range=T)            # range standardise the response (range 0:1 for variable with greatest importance in each modeled response
{                                                     # if Range=F and classification model, the response scale is probability

# function returns a list with each element being PP data for each variable specified by either varOrder OR Nvars
# each list element contains the data needed to draw a PP plot for each model contained in the list of models (Model)
# the PP data consists x, y, and ForRug
# Y is the marginal response for regression or marginal probabilities for classification
  
  sol <- c(call = match.call())
  class(sol) <- "CalcPartialDependence"

        if(class(Model)=="list") {
                MyNames <- names(Model)
                names(MyNames) <- MyNames
        } else {
                Model <- list(Model)         # put a single model into a list of 1
                MyNames <- "Response"
                names(MyNames) <- MyNames
                names(Model) <-  MyNames
        }

        if(ImpType == 1 & ncol(importance(Model[[1]]))==1) stop("cannot have accuracy importance if forest does not contain importance measures")
        # if the variable are not specified use the Nvar most important for the first model
        if(is.null(varOrder)) varOrder <- var.order.rf(Model[[1]], MyType=ImpType)[1:Nvars, 1]
        
        names(varOrder) <- varOrder
        OUT <- vector("list", Nvars)
        names(OUT) <- varOrder
                                                             
        VarImps <- var.order.rf(Model[[1]], MyType=ImpType)
        VarImportances <-  VarImps[match(varOrder, VarImps$Variable),2] # get the importances for the variables
        names(VarImportances) <-  varOrder
               
if(Model[[1]]$type == "classification" & MyClass == "ALL") {
        
        TheClasses <- names(table(Model[[1]]$predicted))
        names(TheClasses) <- TheClasses 
        
        for (i in 1:Nvars) {
        cat("This Var ", varOrder[i], "\n")          
          OUT[[i]] <- lapply(TheClasses, function(x) ProbPartialPlot(Model[[1]], pred.data=MyPreds,  x.var=varOrder[i],
          which.class=x, plot = F, n.pt = Myn.pt))  # note this is a modified version of the partialPlot function
          }
          
    if(Range == T) {
        require(vegan)    
                 
          ClassYRanges <- vector("list", length(TheClasses))   # get the range of the PP response for each class
          names(ClassYRanges) <- TheClasses
          
          for (e in 1:length(TheClasses)) {          # for each CLASS
                  ThisClass <- TheClasses[[e]]                                             
          for (k in 1:Nvars) {                    # for each variable                     
                 ThisVar <- varOrder[k]          
              ClassYRanges[[ThisClass]] <- c(ClassYRanges[[ThisClass]], OUT[[ThisVar]][[ThisClass]]$y)
              cat("Class ", ThisClass, "Var ", ThisVar,   "\n") 
              }}
         
          MyRange <- lapply(ClassYRanges, range)  
         
         for (e in 1:length(TheClasses)) {          # for each CLASS
                  ThisClass <- TheClasses[[e]]
         for (k in 1:Nvars) {                    # for each variable
                 ThisVar <- varOrder[k]
                                                                                                                                  
                  cat("Range Adj response of class", ThisClass, "and ", ThisVar, "with range ", round(range(OUT[[ThisVar]][[ThisClass]]$y),2), 
                  "to have max range ", round(unlist(MyRange[ThisClass]),2),   "\n")
                                                                                         
                  OUT[[ThisVar]][[ThisClass]]$y <- as.numeric(decostand(matrix(OUT[[ThisVar]][[ThisClass]]$y), "range", 
                  MARGIN=2,  range.global=as.matrix(unlist(MyRange[ThisClass]))))
                  
                 }
          } 
  } # end range standardise         
} else { # end if class= "ALL" model all classses       
                           
    if(Model[[1]]$type == "classification") {
          for (i in 1:Nvars) {
          cat("This Variable ", varOrder[i], "\n")
          OUT[[i]] <- lapply(MyNames, function(x) ProbPartialPlot(Model[[x]], pred.data=MyPreds,  x.var=varOrder[i],
          which.class=MyClass, plot = F, n.pt = Myn.pt))  # note this is a modified version of the partialPlot function
          }
    } else {
          for (i in 1:Nvars) {
          cat("This Variable ", varOrder[i], "\n")
          OUT[[i]] <- lapply(MyNames, function(x) MyPartialPlot(Model[[x]], pred.data=MyPreds,  x.var=varOrder[i],
          plot = F, n.pt = Myn.pt))    # note this is a modified version of the partialPlot function TOO!!
          }
    }
          
    if(Range == T) {
        require(vegan)        
        # first find the most important predictor (and its partial plot range in each model)
         ModelImps <- lapply(Model, var.order.rf, MyType=ImpType)
         MostImpVar <- lapply(ModelImps, function(x) x$Variable)
         
          VarRanges <- vector("list", length(MyNames))
          names(VarRanges) <- MyNames
              
         if(Model[[1]]$type == "classification") {   # calculate the range of the PP for ALL variables in each model              
                  for (i in 1:length(MyNames)) {           
                  VarRanges[[i]] <- lapply(MostImpVar[[i]], function(x) range(ProbPartialPlot(Model[[i]], pred.data=MyPreds,  x.var=x,
                  which.class=MyClass, plot = F, n.pt = Myn.pt)$y))  # note this is a modified version of the partialPlot function
              }
         } else {              
                  for (i in 1:length(MyNames)) {           
                  VarRanges[[i]] <- lapply(MostImpVar[[i]], function(x) range(MyPartialPlot(Model[[i]], pred.data=MyPreds,  x.var=x,
                  which.class=MyClass, plot = F, n.pt = Myn.pt)$y)) 
              }
              }
              
          MyRange <- lapply(VarRanges, range)  
         
          for (e in 1:length(MyNames)) {          # for each response
                  ThisModel <- MyNames[[e]]
          for (k in 1:Nvars) {                    # for each variable
                 ThisVar <- varOrder[k]
                 
                  cat("Range Adj ", ThisModel, "and ", ThisVar, "with range ", round(range(OUT[[ThisVar]][[ThisModel]]$y),2), 
                  "to have max range ", round(unlist(MyRange[ThisModel]),2),   "\n")
                                                                            
                  OUT[[ThisVar]][[ThisModel]]$y <- as.numeric(decostand(matrix(OUT[[ThisVar]][[ThisModel]]$y), "range", 
                  MARGIN=2,  range.global=as.matrix(unlist(MyRange[ThisModel]))))
                  
                 }
          } 
    }  # end range adjustment
}          

    sol$classRF <- Model[[1]]$type # keep the type of model for plotting
    sol$VarImportances <- VarImportances
    sol$OUT <- OUT
    sol
}
#######################################################################################
# this function PLOTs the multiple PartialDependence plots  from the returned object from function CalcPartialDependence
#######################################################################################
#######################################################################################
# this function calculates data for multiple PartialDependence plots  
# INCLUDING Confidence intervals
#######################################################################################

CalcPartialDependenceCI <- function(AllPreds = Preds,   # the predictors
                                    MyResp = YFrame,    # the responses
                                    MyType = "ab",      # "ab" for a classification anything else for regression.
                                  MyVarOrder = NULL,    # specify the variable plotting order (if NULL will be based on importance
                                  ImpType = 2,        #  the importance type
                                  Nvars = 8,          # how many variables to compute Partial Dependence for
                                  MyN.pt = 10,        # number of points to compute on varible gradients
                                  ThisClass = 1,        # the class if the model is a classification model
                                  nboot=10,
                                  MyCIs = c(0.025, 0.975),
                                  MyRange=F)            # range standardise the response (range 0:1 for variable with greatest importnce in model
{
# function returns a list with each element being PP data for each variable specified by either varOrder OR Nvars
# each list element contains the data needed to draw a PP plot for each model contained in the list of models (Model)
# the PP data consists x, y, and ForRug
# Y is the marginal response for regression or marginal probabilities for classification

  sol <- c(call = match.call())
  class(sol) <- "CalcPartialDependenceCI"
  
 if(!is.data.frame(MyResp)) {     # make sure response is a DF
  MyResp <- as.data.frame(MyResp) 
  names(MyResp) <- "Response"
  }
    
  # set up some DataFrames to store the Y outputs over nboot times.  
  if(MyType == "ab" & ThisClass == "ALL") {        #either modelling several responses OR several levels of a single (factor) response  
  MyNames <- levels(factor(MyResp[,1]))
  names(MyNames) <- MyNames
  Nmodels <- length(MyNames)
  } else { 
  Nmodels <- ncol(MyResp)
  MyNames <- names(MyResp)
  names(MyNames) <- MyNames
  } 
   
 # if the variable are not specified use the Nvar most important for the first model
if(is.null(MyVarOrder))  {
BaseModel <- resp.mod.fit(resp=data.frame(MyResp[,1]), pred=AllPreds, Occ=0, Type = MyType, importance=T) 
MyVarOrder <- var.order.rf(BaseModel[[1]], MyType=ImpType)[1:Nvars, 1]
}

names(MyVarOrder) <- MyVarOrder

         # make list of lists to store the output
         MakeList <- function(x, myN = length(MyNames)) {
          mylist <- vector(mode = "list", length = myN)
          names(mylist) <- MyNames
          return(mylist)
        }

BootOUT <-  lapply(1:Nvars, function(x) lapply(MakeList(), function(x) as.data.frame(matrix(nrow=nboot, ncol=MyN.pt))))
names(BootOUT) <- MyVarOrder 

for (b in 1:nboot) {
    cat(".........................................................................", "\n")
    cat("This Boot", b,  "\n")
  
    PICK  <- sample(1:nrow(AllPreds), size=nrow(AllPreds),  replace = TRUE)              # bootstrap sample to define a model
    ThisResp <- data.frame(MyResp[PICK, ])    
    Model <- resp.mod.fit(resp=ThisResp, pred=AllPreds[PICK, ], Occ=20, Type = MyType, importance=T)
    #Model <- lapply(Model$TaxaModels, function(x) return(x))

       ComputeBoot <- CalcPartialDependence(Model,              # this needs to be either a list of models or a single model
                                        MyPreds = AllPreds,   # the predictors
                                        varOrder = MyVarOrder,    # specify the variable plotting order (if NULL will be based on importance
                                        ImpType = 2,        #  the importance type
                                        Nvars = length(MyVarOrder),          # how many variables to compute Partial Dependence for
                                        Myn.pt = MyN.pt,        # number of points to compute on varible gradients
                                        MyClass = ThisClass,        # the class if the model is a classification model
                                        Range=MyRange) 
      
      # store the Y outputs here each n itteration.
      for (e in 1:length(MyNames)) {          # for each response
                        ThisModel <- MyNames[e]
                for (k in 1:Nvars) {                    # for each variable
                       ThisVar <- MyVarOrder[k]                                                                          
                        BootOUT[[ThisVar]][[ThisModel]][b,] <-  ComputeBoot$OUT [[ThisVar]][[ThisModel]]$y
                        }             
      
      }
} # end loop for each boot

# process the dataframe above here to extract the median and the CIs
        for (e in 1:length(MyNames)) {          # for each response
                          ThisModel <- MyNames[e]
                  for (k in 1:Nvars) {                    # for each variable
                         ThisVar <- MyVarOrder[k]                                                                                           
                          ComputeBoot$OUT[[ThisVar]][[ThisModel]]$y <-  apply(BootOUT[[ThisVar]][[ThisModel]], 2, median)
                          ComputeBoot$OUT[[ThisVar]][[ThisModel]]$uCI <- apply(BootOUT[[ThisVar]][[ThisModel]], 2,  quantile, probs=MyCIs[1])
                          ComputeBoot$OUT[[ThisVar]][[ThisModel]]$lCI <-   apply(BootOUT[[ThisVar]][[ThisModel]], 2,  quantile, probs=MyCIs[2])                  
                          }  
                 }


    sol$classRF <- Model[[1]]$type # keep the type of model for plotting
    sol$OUT <- ComputeBoot$OUT
    sol
}
#######################################################################################


#######################################################################################
# this function PLOTs the multiple PartialDependence plots  from the returned object from 
# function CalcPartialDependence   AND  CalcPartialDependenceCI
#######################################################################################

PlotPartialDependence <- function(PP,                # this is the returned list from function CalcPartialDependence
                                  Ylim = NULL,       # set limits of Y axis
                                  LogVars = NULL,    # which vars to plot on a log scale
                                  MyLayout = c(2,4), # layout rows cols
                                  LabAxes =c(1,5),   # LABEL axis on which plots ?
                                  Ylab = "Marginal Response",
                                  Pos = "topleft",   # position of legend bow in first plot
                                  Lwd = 1,           # line weight
                                  MyCols = NULL,
                                  MySpan = 0.2,       # smoothing parameter (see function loess)
                                  ShowImportance = T, # include the actual importance meaures on the X  axis lable
                                  RoundImp = 2,       # significant digits to round the importance measure to
                                  ...)               # passed t0 plot
                                     
{

 
VarImportances <- PP$VarImportances
varOrder <- names(PP$OUT)
Nvars <- length(varOrder)
Nplots <- min(MyLayout[1]*MyLayout[2], Nvars)
ResponseNames <- names(PP$OUT[varOrder[1]][[1]])

if (is.null(MyCols))  MyCols <- 1:length(names(PP$OUT)) # automatically assign colours

if(!is.null(PP$OUT[[varOrder[1]]][[ResponseNames[1]]]$uCI)) {   # if there are CIs then plot them
PlotCI <- TRUE
} else {
PlotCI <- FALSE}
  

if(is.null(Ylim)) {
AllY <- vector("numeric")   # get limits for plot Y 

for (k in 1:Nvars) {
       ThisVar <- varOrder[k]
       PlotData <- PP$OUT[[ThisVar]]  # for each variable
    for (e in 1:length(ResponseNames)) {          # for each response    
    if(PlotCI) {
    AllY <- c(AllY, PlotData[[ResponseNames[e]]]$uCI)
    } else {  
 AllY <- c(AllY, PlotData[[ResponseNames[e]]]$y)
 }  
       }
}       
       Ylim <- range(AllY)       
} 

x11(width=17, height=10)                      # ‘c(bottom, left, top, right)     
par(mfrow=MyLayout, cex.lab=1.3, cex.axis=1, oma=c(0,3,0,0), mar=c(4, 1.5, 1, 0.5) + 0.1)
  
  for (i in 1:Nplots) {
       ThisVar <- varOrder[i]
       PlotData <- PP$OUT[[ThisVar]]
       ResponseNames <- names(PlotData)
    
         # use a smoother to simplify the shape of the partial plots 
         MyLoess <- loess(PlotData[[1]]$y ~ PlotData[[1]]$x, span=MySpan)
         MyLoessY <-  MyLoess$fitted
         
         if(ShowImportance == T) {       # include importance or not?
         ThisXlab <- paste(ThisVar, " (", round(VarImportances[i], RoundImp), ")", sep="")
         } else {
         ThisXlab <-  ThisVar
         }
                  
           # Log axis for some variables
       if(!is.na(match(ThisVar, LogVars))) {

         plot(PlotData[[1]]$x, MyLoessY, type="l", ylim=Ylim, lty=1, yaxt="n",
         xlab=ThisXlab, log="x", ylab="", lwd=Lwd, col=MyCols[1], ...)

        } else {
        plot(PlotData[[1]]$x, MyLoessY, type="l", ylim=Ylim, lty=1, yaxt="n",
        xlab=ThisXlab,  ylab="", lwd=Lwd, col=MyCols[1], ...)
        }       
              if(PlotCI) {
              points(PlotData[[1]]$x, PlotData[[1]]$uCI,type="l", col="grey75", lty=3) 
              points(PlotData[[1]]$x, PlotData[[1]]$lCI,type="l", col="grey75", lty=3) 
              }
        
        if(length(ResponseNames) > 1) {
           for(j in 2:length(ResponseNames)) {
           # use a smoother to simplify the shape of the partial plots 
            MyLoess <- loess(PlotData[[j]]$y ~ PlotData[[j]]$x, span=MySpan)
            MyLoessY <-  MyLoess$fitted
            
            points(PlotData[[j]]$x, MyLoessY, type="l", lty=j, lwd=Lwd, col=MyCols[j], ...)
              
              if(PlotCI) {
              points(PlotData[[j]]$x, PlotData[[j]]$uCI,type="l", col="grey75", lty=3) 
              points(PlotData[[j]]$x, PlotData[[j]]$lCI,type="l", col="grey75", lty=3) 
              }
        }     
        }
            
      rug(as.numeric(quantile(PlotData[[1]]$ForRug, probs = seq(0, 1, 0.1), na.rm = TRUE)), side=1)

      if(!is.na(match(i, LabAxes))) {  # if this plot is one to label the axes
      #title(ylab=Ylab)
      # 1=below, 2=left, 3=above and 4=right.
      axis(side=2)
      }

     # Add legend  to first plot
     #VarTaxaImp <- sapply(ModelImport, function(x)  return(match(Impvar.names[i], x$Variable)))
     if(i==1) legend(Pos,  ResponseNames, col = MyCols, lty = 1:length(ResponseNames), lwd=Lwd, inset=0.03,
     pch = c(-1, -1), merge = TRUE,  bty="n" )
  }
 mtext(text=Ylab, side = 2, line= 1, outer=T) 
}
###########################################################################################################

######################################################################################################################
 
resp.mod.fit  <- function(resp, pred, Occ=20, Type = "ab", ...) {
# Occ = percent occupancy of individual taxa to be modelled.
# Type = pa will produce P/A (classification) models

if(nrow(resp)!=nrow(pred)) stop("Wrong data dimensions")
           
if(Type == "pa") {    # convert to PA data
require(vegan)
    resp <- decostand(resp, "pa")
    PerOccupy <- apply(resp, 2, function(x) sum(x>0)/nrow(resp)*100)
    TaxaToModel <- PerOccupy >= Occ
    NMods <- sum(TaxaToModel)
    TaxaMods <- colnames(resp)[TaxaToModel]
    RespToUse <- resp[,TaxaToModel]
} else {
    NMods <- ncol(resp)
    TaxaMods <- colnames(resp)
}
forest.mods <- vector("list", NMods) #  keep each forest model
names(forest.mods) <- TaxaMods
cat("........................................................", "\n")
cat("Fitting Random Forest models for ", NMods, "responses", "\n")

if(Type == "pa") {

      for (i in 1:NMods) {
          cat("........................................................", "\n")
          cat("Fitting presence/absence model to response ",i, ":", TaxaMods[i], "\n")
              forest.mods[[i]] <- randomForest(x=pred, y=as.factor(as.numeric(RespToUse[,i]>0)), ntree=500, ...) # store this model
          }
} else {
        for (i in 1:NMods) {
        cat("........................................................", "\n")
        cat("Fitting abundance model to response: ",i, ":", TaxaMods[i], "\n")
            forest.mods[[i]] <- randomForest(x=pred, y= resp[,i], ntree=500, ...) # store this model
}
}
return(forest.mods)
} # end resp.mod.fit

######################################################################################################################


MapVar <- function(x=domain, divs=10, Var=1, coords, res = 30000, reverse = F, pos="topleft", ...) {
#This function maps variables in domain (variable is specified by Var). The coordinates of the 
#points (classification units/entities) must be specified by cords. The gradient in the variable 
#is shown by colors from red to blue(the number of which are defined by divs). 
          
 ob.name <- names(x)[Var]  #name of Var

   if(nrow(coords)>res) {
              pick <- sample(1:nrow(coords), size=res, replace = F)
              } else {
              pick<-c(1:nrow(coords))
              }
              coords.s <- coords[pick,]
              domain.var <- x[pick,Var]
  
  HeatCol <- rainbow(divs, alpha = 1, start=0, end=4/6)
 if(reverse) HeatCol <- rev(rainbow(divs, alpha = 1, start=0, end=4/6))

cut.var <- cut(domain.var, breaks=divs, labels=F)
NamesCutVar <- names(table(cut.var)) 
div.mean <-  round(tapply(X=domain.var, INDEX=as.factor(cut.var), FUN = mean),2)
OrdVar <- names(sort(table(cut.var), decreasing=T))  # put in decreasing order by size of class 

   plot(coords.s, col="white")
   for (i in 1:divs) points(coords.s[cut.var==OrdVar[i],],col=HeatCol[match(OrdVar[i], NamesCutVar)], pch=15, ...) 
   
legend(pos, legend=div.mean, title=paste("Divisions of ", ob.name, sep=""),
 text.col = "black",  pch = 15 , bg = 'gray90', col=HeatCol,
            inset = .01) 
            
} # end MapVar

varSelRF.reg  <- function (
                           xdata,          # A data frame or matrix, with subjects/cases in rows and variables in columns. NAs not allowed.
                           Resp,           # The dependent variable; must be a continuous numeric
                           c.sd = 1,       # The factor that multiplies the sd. to decide on stopping the tierations or choosing the final solution. See reference for details.
                           nboot=10,       # number of booot straps to use to compute the mean and standard error of the prediction errors
                           mtryFactor = 1, # The multiplication factor of ?{number.of.variables} for the number of variables to use for the ntry argument of randomForest.
                           ntree = 5000,   # The number of trees to use for the first forest; same as ntree for randomForest.
                           ntreeIterat = 2000, # The number of trees to use (ntree of randomForest) for all additional forests.
                           vars.drop.num = NULL,  # The number of variables to exclude at each iteration.
                           vars.drop.frac = 0.2,  # The fraction of variables, from those in the previous forest, to exclude at each iteration.
                           whole.range = TRUE,    # 	If TRUE continue dropping variables until a forest with only two variables is built, and choose the best model from the complete series of models.
                           recompute.var.imp = FALSE,  # If TRUE recompute variable importances at each new iteration. (not recomended)
                           verbose = TRUE,             # Give more information about what is being done.
                           returnFirstForest = TRUE,   # An (optional) object of class randomForest previously fitted
                           keep.forest = FALSE) {       # Same argument as in randomForest function. If the forest is kept, it will be returned as part of the "rf.model" component of the output.
    
# modified version of varSelRF from package::varSelRF
# this will do variable selection for continuous response variables (regression)
# there is also modified plotting and boot strap (to smooth variation over runs (bootstraps)
# in this function the mean error and its standard error are computed by running bootstrap samples of the models at each itteration of the variable
# reduction process. The number of bootstraps are specified by nboot. The mean and se of the OOB error is calculated over the bootstrapped OOB error 
# for each bootsrap.
    
    if (is.factor(Resp))
        stop("Response should be a continuous varaiable")
    if ((is.null(vars.drop.num) & is.null(vars.drop.frac)) |
        (!is.null(vars.drop.num) & !is.null(vars.drop.frac)))
        stop("One (and only one) of vars.drop.frac and vars.drop.num must be NULL and the other set")
    max.num.steps <- dim(xdata)[2]
    num.subjects <- dim(xdata)[1]
    if (is.null(colnames(xdata)))
        colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep = "")
    n.vars <- vars <- OOB.rf <- OOB.sd <- rep(NA, max.num.steps)
       mtry <- floor(sqrt(ncol(xdata)) * mtryFactor)
          
oobError <- function(rf) {   # modified to return the MSE
        ooo <- mean((rf$predicted - rf$y)^2)
        return(ooo)
    }
    
# this function boots the RF model and returns the OOB error from each boot strapped model        
bootrf <- function(nboot=nboot, ...) {  # args passed to randomForest 
                         oobboot <- 1:nboot
                         names(oobboot) <- 1:nboot
                         cat("Computing OOB error", "\n")
                         for (i in 1:nboot) {
                                        boot.Sample <- sample(1:num.subjects, size=num.subjects, replace = T)
                                        thisRFboot <- randomForest(x = xdata[boot.Sample, ] , y = Resp[boot.Sample ], ntree = ntree,
                                                                    mtry = mtry, importance = F, keep.forest = F)
                                       oobboot[i] <- oobError(thisRFboot)
                                       cat(i, " ")
                                       }
                                       return(oobboot)
                                       }                     
      
 rf <- randomForest(x = xdata, y = Resp, ntree = ntree,   # the initial model with all variables
            mtry = mtry, importance = TRUE, keep.forest = keep.forest)
        
    if (returnFirstForest)
        FirstForest <- rf
    else FirstForest <- NULL
        
    initial.RFboot  <-  bootrf(nboot=nboot)   # determine mean and se of OOB error
    m.iterated.ob.error <- m.initial.ob.error <- mean(initial.RFboot)
    sd.iterated.ob.error <- sd.initial.ob.error <-  sd(initial.RFboot)/sqrt(nboot)

    if (verbose) {
        print(paste("Initial OOB error: mean = ", round(m.initial.ob.error,
            4), "; se = ", round(sd.initial.ob.error, 4), sep = "")) }
   
    importances <- importance(rf, type = 1, scale = FALSE)
    selected.vars <- order(importances, decreasing = TRUE)
    ordered.importances <- importances[selected.vars]
    initialImportances <- importances
    initialOrderedImportances <- ordered.importances
    j <- 1
    n.vars[j] <- dim(xdata)[2]
    vars[j] <- paste(colnames(xdata), collapse = " + ")
    OOB.rf[j] <- m.iterated.ob.error
    OOB.sd[j] <- sd.iterated.ob.error
    var.simplify <- TRUE
    
          while (var.simplify) {
              
              gc()

              last.rf <- rf
              last.vars <- selected.vars
              previous.m.error <- m.iterated.ob.error
              previous.sd.error <- sd.iterated.ob.error
              if (length(selected.vars) <= 2) {
                  var.simplify <- FALSE
                  break
              }
              if (recompute.var.imp & (j > 1)) {
                  importances <- importance(rf, type = 1, scale = FALSE)
                  tmp.order <- order(importances, decreasing = TRUE)
                  selected.vars <- selected.vars[tmp.order]
                  ordered.importances <- importances[tmp.order]
              }
              num.vars <- length(selected.vars)
              if (is.null(vars.drop.num))
                  vars.drop <- round(num.vars * vars.drop.frac)
              else vars.drop <- vars.drop.num
              if (num.vars >= (vars.drop + 2)) {
                  selected.vars <- selected.vars[1:(num.vars - vars.drop)]
                  ordered.importances <- ordered.importances[1:(num.vars -
                      vars.drop)]
              }
              else {
                  selected.vars <- selected.vars[1:2]
                  ordered.importances <- ordered.importances[1:2]
              }
              if ((length(selected.vars) < 2) | (any(selected.vars <
                  1))) {
                  var.simplify <- FALSE
                  break
              }
              mtry <- floor(sqrt(length(selected.vars)) * mtryFactor)
              if (mtry > length(selected.vars))
                  mtry <- length(selected.vars)
              if (recompute.var.imp)
                  rf <- randomForest(x = xdata[, selected.vars], y = Resp,
                      importance = TRUE, ntree = ntree, mtry = mtry,
                      keep.forest = keep.forest)
              else rf <- randomForest(x = xdata[, selected.vars], y = Resp,
                  importance = FALSE, ntree = ntreeIterat, mtry = mtry,
                  keep.forest = keep.forest)
                  
                  iterated.RFboot  <-  bootrf(nboot=nboot, x = xdata[, selected.vars], y = Resp,
                                       ntree = ntree, mtry = mtry, importance = F, keep.forest = F)
                   m.iterated.ob.error <-  mean(iterated.RFboot)                     
                   sd.iterated.ob.error <- sd(iterated.RFboot)/sqrt(nboot)

              if (verbose) {
                  print(paste("..... iteration ", j, "; OOB error: mean = ",
                      round(m.iterated.ob.error, 4), "; se = ", round(sd.iterated.ob.error,
                        4), "; num. vars = ", length(selected.vars),
                      sep = ""))
              }
              j <- j + 1
              n.vars[j] <- length(selected.vars)
              vars[j] <- paste(colnames(xdata)[selected.vars], collapse = " + ")
              OOB.rf[j] <- m.iterated.ob.error
              OOB.sd[j] <- sd.iterated.ob.error
              
              if (!whole.range & ((m.iterated.ob.error > (m.initial.ob.error +
                  c.sd * sd.initial.ob.error)) | (m.iterated.ob.error >
                  (previous.m.error + c.sd * previous.sd.error))))
                  var.simplify <- FALSE
          }    
    
    if (!whole.range) {
        if (!is.null(colnames(xdata)))
            selected.vars <- sort(colnames(xdata)[last.vars])
        else selected.vars <- last.vars
        out <- list(selec.history = data.frame(Number.Variables = n.vars,
            Vars.in.Forest = vars, OOB = OOB.rf, sd.OOB = OOB.sd)[1:j,
            ], rf.model = last.rf, selected.vars = selected.vars,
            selected.model = paste(selected.vars, collapse = " + "),
            best.model.nvars = length(selected.vars), initialImportances = initialImportances,
            initialOrderedImportances = initialOrderedImportances,
            ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor,
            firstForest = FirstForest)
        class(out) <- "varSelRF.reg"
        return(out)
    }
    else {
         
        n.vars <- n.vars[1:j]
        vars <- vars[1:j]
        OOB.rf <- OOB.rf[1:j]
        OOB.sd <- OOB.sd[1:j]
        
        # determine best model 
        min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]
        best.pos <- which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <=
            min.oob.ci)])]
            
        selected.vars <- sort(unlist(strsplit(vars[best.pos],
            " + ", fixed = TRUE)))
        out <- list(selec.history = data.frame(Number.Variables = n.vars, 
            Vars.in.Forest = vars, OOB = OOB.rf, se.OOB = OOB.sd), c.sd = c.sd, 
            rf.model = NA, selected.vars = selected.vars, selected.model = paste(selected.vars,
                collapse = " + "), best.model.nvars = n.vars[best.pos],
            initialImportances = initialImportances, initialOrderedImportances = initialOrderedImportances,
            ntree = ntree, ntreeIterat = ntreeIterat, mtryFactor = mtryFactor,
            firstForest = FirstForest)
        class(out) <- "varSelRF.reg"
        return(out)
    }
}

# ploting function  for  varSelRF.Reg
# getS3method("plot","varSelRF")
varSelRF.RegPlot <- function (x, nvar = NULL, ...)
{

        op <- par(las = 1)
          on.exit(par(op))

    if (is.null(nvar))
        nvar <- min(30, length(x$initialOrderedImportances))
        
        c.sd <- x$c.sd

  x11()   # plot the initial importances

myImp <- as.vector(x$initialImportances)
names(myImp) <- row.names(x$initialImportances)
myImp <- sort(myImp, decreasing=F)
dotchart(myImp, main = "Initial importances", xlab = "Importances (unscaled)")

# plot the best model 
n.vars <- x$selec.history$Number.Variables
OOB.rf  <-  x$selec.history$OOB
OOB.sd <-  x$selec.history$se.OOB 
 
   # determine best model 
   

min.oob.ci <- min(OOB.rf) + c.sd * OOB.sd[which.min(OOB.rf)]  # min oob PLUS its standard error
best.pos <- which(OOB.rf <= min.oob.ci)[which.min(n.vars[which(OOB.rf <= min.oob.ci)])]
best.model.nvars <- n.vars[best.pos]  # model whose mean OOB error is within (c.sd * OOB.sd) of minimum error model
           
seU <-  x$selec.history$OOB + x$selec.history$se.OOB
seL <-  x$selec.history$OOB - x$selec.history$se.OOB

x11()
        ylim <- range(c(x$selec.history$OOB + 2 * x$selec.history$se.OOB, x$selec.history$OOB - 2 * x$selec.history$se.OOB))
        plot(x$selec.history$Number.Variables, x$selec.history$OOB,
            type = "b", xlab = "Number of variables used", ylab = "OOB error", lwd = 2, 
            sub = "Most parsimonious model whose mean OOB error is within (c.sd * SE of OOB) of model with minimum error", 
            log = "x", ylim = ylim, cex.sub = 0.7, ...)
        abline(v=best.model.nvars,  col="green")
        abline(h= min.oob.ci, col="red")
        arrows(x$selec.history$Number.Variables, seU, x$selec.history$Number.Variables, seL, length=.05,angle=90, code=3, col="grey")  # add SE errors bars to plot
                
  legend("topright",  c("Mean OOB error", "SE of OOB error", paste("Minimum OOB error +", c.sd, "* se"), 
  paste("Minimum predictors = ", best.model.nvars)),
   col = c("black","gray70", "red", "green"), text.col = "black", lwd =c(2, 1,1,1), lty = c(1, 1, 1), pch = c(1, -1, -1,-1),
           merge = TRUE, bg = 'white', inset = .05)

}

 # modifed  varSelRFBoot to boot the varSelRF.reg function     # how many boot straps runs of the variable selection proceedure (i.e. varSelRFBoot)
 # greatly modified NO cluster functionality, and simplified overall solution 
varSelRFBoot.reg <- function (xdata,   # A data frame or matrix, with subjects/cases in rows and variables in columns. NAs not allowed.
                              Resp ,    # response (continuous variable)
                              c.sd = 1, 
                              mtryFactor = 1, 
                              ntree = 500, 
                              Nboot=10,     # number of booot straps to use to compute the mean and standard error of the prediction errors
                              ntreeIterat = 500, 
                              vars.drop.frac = 0.2, 
                              bootnumber = 200,  # how many boot straps runs of the variable selection proceedure (i.e. varSelRFBoot)
                              whole.range = TRUE, 
                              recompute.var.imp = FALSE, 
                              srf = NULL, 
                              verbose = TRUE,
                              ...)
{

    if (is.null(colnames(xdata)))
        colnames(xdata) <- paste("v", 1:dim(xdata)[2], sep = "")
              if (!is.null(srf)) {
                  if (class(srf) != "varSelRF")
                      stop("srf must be the results of a run of varSelRF")
                  n.ntree <- srf$ntree
                  n.ntreeIterat <- srf$ntreeIterat
                  n.mtryFactor <- srf$mtryFactor
                  if ((n.ntree != ntree) | (n.mtryFactor != mtryFactor) |
                      (n.ntreeIterat != ntreeIterat))
                      warning("Using as ntree and mtryFactor the parameters obtained from srf",
                          immediate. = TRUE)
                  ntree <- n.ntree
                  mtryFactor <- n.mtryFactor
                  rm(n.ntree, n.mtryFactor)
                  all.data.run <- srf
              }  else {       
        all.data.run <- varSelRF.reg(
                                     Resp = Resp,       # response data
                                     xdata = xdata,     # predictors
                                     nboot=Nboot,       # number of booot straps to use to compute the mean and standard error of the prediction errors
                                     c.sd = c.sd, 
                                     mtryFactor = mtryFactor, 
                                     ntree = ntree,
                                     ntreeIterat = ntreeIterat, 
                                     vars.drop.frac = vars.drop.frac, 
                                     whole.range = whole.range, 
                                     recompute.var.imp = recompute.var.imp, 
                                     verbose = T)
    }
    columns.data <- which(colnames(xdata) %in% all.data.run$selected.vars)
    all.data.rf.mtry <- floor(mtryFactor * sqrt(length(columns.data)))
    if (all.data.rf.mtry > length(columns.data))
        all.data.rf.mtry <- length(columns.data)
    all.data.rf.predict <- randomForest(y = Resp, x = xdata[,
        columns.data], ntree = all.data.run$ntree, mtry = all.data.rf.mtry,
        xtest = xdata[, columns.data], ytest = Resp, keep.forest = FALSE)
    full.pred <- all.data.rf.predict$test$predicted
    all.data.selected.vars <- all.data.run$selected.vars
    all.data.best.model.nvars <- all.data.run$best.model.nvars
    N <- length(Resp)
    solution.sizes <- rep(NA, bootnumber)
    overlap.with.full <- rep(NA, bootnumber)
    vars.in.solutions <- vector()

          bootTrainTest <- function(dummy, c.sd, mtryFactor, ntree,
              ntreeIterat, whole.range, recompute.var.imp, ...) {
              N <- length(ClassTheCluster)
              sample.again <- TRUE
              while (sample.again) {
                  bootsample <- unlist(tapply(1:N, ClassTheCluster,
                      function(x) sample(x, size = length(x), replace = TRUE)))
                  nobootsample <- setdiff(1:N, bootsample)
                  if (!length(nobootsample))
                      sample.again <- TRUE
                  else sample.again <- FALSE
              }
              train.data <- xdataTheCluster[bootsample, , drop = FALSE]
              test.data <- xdataTheCluster[nobootsample, , drop = FALSE]
              train.class <- ClassTheCluster[bootsample]
              test.class <- ClassTheCluster[nobootsample]
              boot.run <- varSelRF.reg(Resp = train.class, xdata = train.data,
                  c.sd = c.sd, mtryFactor = mtryFactor, ntree = ntree,
                  ntreeIterat = ntreeIterat, whole.range = whole.range,
                  recompute.var.imp = recompute.var.imp, vars.drop.frac = vars.drop.frac)
              output.cl <- list()
              output.cl$best.model.nvars <- boot.run$best.model.nvars
              output.cl$selected.model <- boot.run$selected.model
              output.cl$selected.vars <- boot.run$selected.vars
              output.cl$nobootsample <- nobootsample
              output.cl$bootsample <- bootsample
              output.cl$initialImportances <- boot.run$initialImportances
              output.cl$initialOrderedImportances <- boot.run$initialOrderedImportances
              output.cl$selec.history <- boot.run$selec.history
              boot.col.data <- which(colnames(xdataTheCluster) %in%
                  boot.run$selected.vars)
              run.test.mtry <- floor(mtryFactor * sqrt(length(boot.col.data)))
              if (run.test.mtry > length(boot.col.data))
                  run.test.mtry <- length(boot.col.data)
              boot.run.test <- randomForest(y = train.class, x = train.data[,
                  boot.col.data, drop = FALSE], ntree = boot.run$ntree,
                  mtry = run.test.mtry, keep.forest = FALSE, xtest = test.data[,
                      boot.col.data, drop = FALSE])$test
              output.cl$class.pred.array <- boot.run.test$predicted
              output.cl$prob.pred.array <- boot.run.test$votes
              rm(boot.run.test)
              gc()
              return(output.cl)
          }

        boot.runs <- list()
        xdataTheCluster <- xdata
        ClassTheCluster <- Resp
        cat("\n",      "Running bootstrap iterations", "\n")
        for (nboot in 1:bootnumber) {
            cat("bootstrap number", nboot, "\n")
            boot.runs[[nboot]] <- bootTrainTest(nboot=Nboot, c.sd = c.sd,
                mtryFactor = mtryFactor, ntree = ntree, ntreeIterat = ntreeIterat,
                whole.range = whole.range, recompute.var.imp = recompute.var.imp)
        }
        cat("\n")

    solutions <- unlist(lapply(boot.runs, function(z) {
        paste(sort(z$selected.vars), collapse = " + ") }))

    vars.in.solutions <- lapply(boot.runs, function(z) z$selected.vars)
    solution.sizes <- unlist(lapply(boot.runs, function(z) z$best.model.nvars))

    # deleted some original code for classification solution from here -- simplified results

    Median.Solution <-  median(solution.sizes) # median of No. variable per solution

      Unique.vars.in.solution <- unique(unlist(vars.in.solutions))
      names(Unique.vars.in.solution) <- Unique.vars.in.solution
      # the number of times each varaible was selected
      No.times.var.selected <- unlist(lapply(Unique.vars.in.solution, function(x) sum(x == unlist(vars.in.solutions))))
      Final.Vars <- names(sort(No.times.var.selected, decreasing=T)[1:Median.Solution])
      
       FinalReducedModel <- randomForest(y = Resp, x = xdata[,Final.Vars], importance = T, keep.forest = T)      # response data
                                     
    out <- list(FinalReducedModel=FinalReducedModel, Final.Vars=Final.Vars,  No.times.var.selected= No.times.var.selected,
        number.of.bootsamples = bootnumber,
        all.data.randomForest = all.data.rf.predict, all.data.vars = all.data.selected.vars,
        all.data.run = all.data.run,
         number.of.vars = solution.sizes,
         all.vars.in.solutions = vars.in.solutions,
        all.solutions = solutions,  allBootRuns = boot.runs)
    class(out) <- "varSelRFBoot.rep"
    return(out)
}


##############################################################
# function: map.class 
##############################################################
 map.class <- function(classes, level, coords, res=20000, pos=NULL,...)  {

 # function to map a classification
 ob.name <- deparse(substitute(classes))
 tclass <-  as.factor(classes[,level])
 n.class <- length(which(as.vector((table(tclass) >0))))
 Tclass <-  table(tclass) 
 OrdClass <- names(sort(Tclass, decreasing=T))  # put in decreasing order by size of class
 PropClass <- round(100*(Tclass/length(tclass)),1) # proportion of domain asigned to each class
 class.labs <- as.factor(names(Tclass))
 class.names <-  names(Tclass)
 col.c <- rep(c(1:8), length.out=n.class)
 pch.c <- rep(c(rep(15,8),rep(16,8), rep(17,8)), length.out=n.class)

 if(nrow(coords)>res) {
              pick <- sample(1:nrow(coords), size=res, replace = F)
              } else {
              pick<-c(1:nrow(coords))
              }
              coords.s <- coords[pick,]
              tclass.s <- tclass[pick] 
              
 x11()
 plot(coords.s, col="white",xlab = "east", ylab = "north", ...)
     for(i in 1:n.class) {
     points(coords.s[tclass.s==OrdClass[i], ], col=col.c[match(OrdClass[i], class.names)], pch=pch.c[i], cex=0.6)
     }

  if(!is.null(pos)) {
     legend(pos, legend=paste(class.names, " (", PropClass, ")", sep=""), col = col.c, title="Class (% domain)",
     text.col = "black",  pch = pch.c, bg = 'gray90',  inset = .01) 
     }

 } # end map classification function
 
MapVar <- function(x=domain, divs=10, Var=1, coords, res = 30000, reverse = F, pos="topleft", ...) {
#This function maps variables in domain (variable is specified by Var). The coordinates of the 
#points (classification units/entities) must be specified by cords. The gradient in the variable 
#is shown by colors from red to blue(the number of which are defined by divs). 
          
 ob.name <- names(x)[Var]  #name of Var

   if(nrow(coords)>res) {
              pick <- sample(1:nrow(coords), size=res, replace = F)
              } else {
              pick<-c(1:nrow(coords))
              }
              coords.s <- coords[pick,]
              domain.var <- x[pick,Var]
  
  HeatCol <- rainbow(divs, alpha = 1, start=0, end=4/6)
 if(reverse) HeatCol <- rev(rainbow(divs, alpha = 1, start=0, end=4/6))

cut.var <- cut(domain.var, breaks=divs, labels=F)
NamesCutVar <- names(table(cut.var)) 
div.mean <-  round(tapply(X=domain.var, INDEX=as.factor(cut.var), FUN = mean),2)
OrdVar <- names(sort(table(cut.var), decreasing=T))  # put in decreasing order by size of class 

par(mar= c(1, 1, 2, 1) +  0.1)
   plot(coords.s, col="white", xaxt="n", yaxt="n", ...)
   for (i in 1:divs) points(coords.s[cut.var==OrdVar[i],],col=HeatCol[match(OrdVar[i], NamesCutVar)], pch=15) 
   
legend(pos, legend=div.mean, title=paste("Divisions of ", ob.name, sep=""),
 text.col = "black",  pch = 15 , bg = 'gray90', col=HeatCol,
            inset = .01) 
            
} # end MapVar

######################################################################################################################
######################################################################################################################

classify.domain <- function(x, MyMethod="manhattan", clusters="pam", level=seq(5,50,5), large = 1500 ) {
  # this function will define a classification from an input dataset at level of detail specified
  # x needs to be a matrix of transformed variables representing the domain
  # the function will choose the best clustering method based on the size of domain.
  require(cluster)
  require(vegan)

  ClassDomain <- as.data.frame(matrix(nrow=nrow(x), ncol=length(level)))  # class assignmments
  row.names(ClassDomain) <- row.names(x)
  names(ClassDomain) <- paste("Lev_", as.character(level), sep="")
   
  if(nrow(x) > large) {
    cat("........................................................", "\n")
    cat("Defining classification level", max(level), "by first Partioning using Clara","\n")
    domain.part <- clara(x, k=max(level), metric = MyMethod) # fast partitioning using clara
    part.grp <- as.data.frame(matrix(nrow = max(level), ncol = length(level)))  # make a blank matrix to store cluster membership at each level of classification detail
    names(part.grp) <- paste("Lev_", level, sep="")
    part.grp[,length(level)] <- as.factor(1:max(level))
    ClassDomain[,length(level)] <- as.factor(domain.part$clustering)

    if(clusters!="pam") {
    part.hclust <- hclust(vegdist(domain.part$medoids, method=MyMethod), method=clusters)   # now hierarchical clustering       
   
      for(i in 1:(length(level)-1)) { 
      cat("Defining classification level", level[i],"using", clusters,"\n") 
      part.grp[,i] <- as.factor(cutree(part.hclust, k = level[i]))  # define clusters at level i
      ClassDomain[,i] <- as.factor(part.grp[match(as.vector(domain.part$clustering), part.grp[,length(level)]), i]) # store the site to class assignments
      }
     
     } else {
     for(i in 1:(length(level)-1)) {
     cat("Defining classification level", level[i],"using Clara", "\n")       
     #part.grp[,i] <- pam(x=vegdist(domain.part$medoids, MyMethod), k = level[i])$clustering  # use PAM to define the lower level classes
     #ClassDomain[,i] <- as.factor(part.grp[match(as.vector(domain.part$clustering), part.grp[,length(level)]), i]) # store the site to class assignments
     ClassDomain[,i] <- as.factor(clara(x, k=level[i], metric = MyMethod)$clustering)
     }
     }

  } else {  # if the data set is not large
    cat("........................................................", "\n")
    cat("Defining classification small","\n")

    domainDis <-  vegdist(x, method=MyMethod) # nonextended RESPONSE distances

    if(clusters!="pam") {
      domain.hclus <- hclust(domainDis, method = clusters) # use hierachical clustering
    }

    for(i in 1:length(level)) {
      cat("........................................................", "\n")
      cat("Defining classification level ", level[i],"\n")
      if(clusters=="pam") {
        domain.grp <- as.factor(pam(domainDis, k=level[i])$clustering)
      } else {
        domain.grp <- as.factor(cutree(domain.hclus, k = level[i]))
      }
      ClassDomain[,i] <- domain.grp # store the site to class assignments
    }
  }

return(ClassDomain)
}
################################################################################### 