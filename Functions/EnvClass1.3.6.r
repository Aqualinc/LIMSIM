# EnvClass Version 1.3.6, 27 March 2012, improved MapVar and map.class
# EnvClass Version 1.3.5, 12 Sept 2011, improved mapping functions map.class and  MapVar
# EnvClass Version 1.3.4, 12 Sept 2011, improved taxa.class.fit function
# EnvClass Version 1.3.3, 16 March 2011, improved the classify.domain function 
# EnvClass Version 1.3.2, 8 March 2011, modify forest class to include Kmeans clustering! improved the classify.domain function 
# EnvClass Version 1.3.1, 19 October 2010, improved the dismod.transform function, minor changes to classify.domain function
# EnvClass Version 1.3, 21 AUGUST 2010, incuded predict first classify leter functions - taxa.class.fit/predict
# EnvClass Version 1.2.6, 7 AUGUST 2010, mods to dismod.fit including new dismod.plot and changes to making predictors
# EnvClass Version 1.2.4, 28 June 2010, mods to dismod.plot to plot the functions fitted to each predictor 
# fitted in a dissimilarity model returned by dismod.fit. The names of the fitted varaibles and their transformations are now seprated
# by an undersore in the fit.preds item of the object returned by dismod.fit 
# EnvClass Version 1.2.3, 15 June 2010, fix error using tree in function evaluate
# EnvClass Version 1.2.2, 2 July 2009, changes to previous version 1.2.1
# added the ... operator to meandist, anosift, homog, mansift and partial.mansift in order 
# to pass additional arguments to vegdist 


# List of all Classification Functions (The script contains other functions called by these functions) 	
#Classification Testing and Evaluating Functions	
      #Map a classification	                      map.class
      #Mean Dissimilarity Analysis              	meandist
      #Mean Dissimilarity Dendrogram            	mddendro
      #Homogeneity                                Analysis	homog
      #Anosift – ANOSIM	                          anosift
      #Mansift – Mantel test	                    mansift
      #Partial.Mansift – Partial Mantel test	    Partial.Mansift 
      #Enchain Classification Analyses	          enchain
      #Plot an Enchain Analysis	                  plot.enchain
      #Evaluate Classifications	                  evaluate
      #Map an Evaluate Analysis	                  map.evaluate
      #Plot an Evaluate Analysis	                plot.evaluate
#Classification defining Functions	
      #Map a variable	                            MapVar 
      #Describe a classification	                DescribeClass
      #Define a classification	                  ClassifyDomain
      #Fit a Matrix Correlation Model	            matcor.fit 
      #Predict from a Matrix Correlation Model	  matcor.predict 
      #Fit a Dissimilarity Model	                dismod.fit 
      #Plot a Dissimilarity Model	                dismod.plot (dismod)
      #Predict from a Dissimilarity Model	        dismod.predict (dismod)
      #Cassify from a dissimilarity model	        dismod.class (dismod)
      #Transform predictor space using a dissimilarity model	                           dismod.transform(dismod)
      #Fit and predict a classification based on multiple linear discriminant analysis	 lda.class
      #Fit and predict a tree based classification	                                     tree.class
      #Fit and predict a Random Forest based classification	                             forest.class
      #Fit a predict first and classify later based classification	                     taxa.class.fit  
      #Fit a predict first and classify later based classification	                     taxa.class.predict  


require("vegan")


# EnvClass Version 1.2.1, 22 May 2009, fixes bugs in previous version 1.2
# added code to function 'evaluate' and 'enchain' to catch illogical specification of NumClass
# added code to function 'evaluate' to warn that the CART model that defines 
# the georegions cannot produce as many regions as there are classes being tested. 
# and set the 'minsize' argument default to 2 so that the CART model fits the data perfectly i.e. produces
# as many regions as there are classes being tested.
# made the 'minsize' argument in tree.class = 2 (default) this allows the tree to perfectly fit the data 
# set the row.names of the argument 'domain' to be the row.names of the 'domain.classes' item of objects returned by:
# tree.class, forest.class, dismod.class, lda.class.

meandist <- function(testdata, groups, method="euclidean", mins=2, permutations=F, nboot=F, boot.perc=0.9, cboot=0.95, ...)
{
# Meandist calculates the mean distance (or dissimilarity) statistic for a multivariate dataset
# where cases are grouped into several groups, such as a classification.
# Meandist has been adapted from the "meansim" statistic (Van Sickle 1997) using functions from
# the VEGAN package and this package is loaded in order to use the function.
# The meandist statistic is the mean between group distance minus the mean within distance distance.
# The mean distances for each class are weighted by the number of cases in the group (n) divided by
# the total number of cases (N). A permutation test is made on the data by permuting the cases
# between the groups while holding n for each group constant. Meandist works with the vegdist
# function and therefore uses distance (dissimilarity) measures, which is the reverse of VanSickle's
# original implementation which uses similarities, however the statistics and p-values are equivalent.

# The function allows the minimum number of cases in each group to be set with the argument "mins".
# This reduces the possibilty that extreme values in groups with small numbers of cases will affect the statistic.
# The function returns several components including the meandist statistic, the mean within class distance
# in each group, the mean between group distance.

# The function also calculates confidence intervals for the meandist statistic using a bootstrap procedure.
# Bootstrapping with replacement is not effective for distance matrices because the distance between a sample and
# its replicate is zero. The method of Goslee and Urban (2007) is therefore used where bootstrap samples are
# constructed by randomly selecting some proportion of the data (set by argument  boot.perc) without replacement. The default value for boot.perc is 90% as recommended by  Goslee and Urban (2007). Meandist is calculated for some number of replicates of such samples
# (set with the argument nboot) and the 5th and 95th confidence intervals are calculated from the distribution of values of meandist.

# NOTE that "enchain" is an additional utility that calculates mean distance statistics for several levels of classification detail (or classifications)

# Van Sickle, J. 1997. Using mean similarity dendrograms to evaluate classifications.
# Journal of Agricultural, Biological and Environmental Statistics 2:370-388.

# Goslee, S. C., and D. L. Urban. 2007. The ecodist Package for Dissimilarity-based Analysis of Ecological Data.
# Journal of Statistical Software 22

require("vegan")

  #function tests
  if(mins<2) stop("mins must be >1")
  if(!is.factor(groups)) stop("groups must be a factor vector")
  if(boot.perc>1 | boot.perc<0) stop("p.perc must by a value between 0 and 1") #NB boot.perc is considered as 0 if FALSE

  if(class(testdata)!= "dist") {
  if(nrow(testdata)!=length(groups)) stop("Wrong data dimensions")
  } else {
  if(length(testdata) != (length(groups)*(length(groups)-1))/2)   stop("Wrong data dimensions")
}   # end tests
  

  #inner-functions
  matched <- function(irow, icol, grouping){#  extracting within and between
    grouping[irow] == grouping[icol]
  }

  CompSeq<-function(Thisirow,Thisicol,Thisgrouping,Thisx){
    within <- matched(Thisirow, Thisicol, Thisgrouping)  # function above - assigns  distances to either within or between
    cl.vec <- rep("Between", length(Thisx))  # makes a vector of length N (= no. distances)all entries are "between"
    take <- as.numeric(Thisirow[within])   # take is a vector of L equal to number of withins
    cl.vec[within] <- levels(Thisgrouping)[Thisgrouping[take]]
    cl.vec <- factor(cl.vec, levels = c("Between", levels(Thisgrouping)))   # factor vector of between or within class categories for each distance measure!!
    # add  code here to compute MEANDIST
    group.num <- as.vector(table(Thisgrouping))  # numeric vector of No. sites by group(class)
    ThisN<-sum(group.num)
    group.means <- tapply(Thisx, cl.vec, mean)   # calculate mean distances by group
    between <- group.means[1]   # extract the mean between groups distance from group.means
    within <- group.means[2:length(group.means)]  # extract all the within group distances from the group.means vector
    weighted.within.mean <- sum(within*(group.num/ThisN),na.rm=T) # calculate the weighted average within group distance among groups with more than 1 member
    meandist <- between - weighted.within.mean       # statistic for Similarities
    stats<-list(within=within,meandist=as.numeric(meandist),cl.vec=cl.vec,weighted.within.mean=weighted.within.mean,between=between)
    return(stats)
  }

  as.randtest<-function (sim, obs, alter = c("greater", "less", "two-sided"),call = match.call()){
    res <- list(sim = sim, obs = obs)
    res$alter <- match.arg(alter)
    res$rep <- length(sim)
    res$expvar <- c(Std.Obs = (res$obs - mean(sim))/sd(sim),
        Expectation = mean(sim), Variance = var(sim))
    if (res$alter == "greater") {
        res$pvalue <- (sum(sim >= obs) + 1)/(length(sim) + 1)
    }
    else if (res$alter == "less") {
        res$pvalue <- (sum(sim <= obs) + 1)/(length(sim) + 1)
    }
    else if (res$alter == "two-sided") {
        sim0 <- abs(sim - mean(sim))
        obs0 <- abs(obs - mean(sim))
        res$pvalue <- (sum(sim0 >= obs0) + 1)/(length(sim) +
            1)
    }
    res$call <- call
    class(res) <- "randtest"
    return(res)
  }
#end of inner-functions

  # remove cases in classes with n fewer than "mins"
  grps<-groups
  levels(grps)<-table(groups)>=mins
  grouping<-factor(groups[grps==T])
  classes.tested <- length(table(grouping))
     
  
    if(class(testdata)!= "dist") {
        testset<-testdata[grps==T,] #eliminate non selected cases from  testdata input data matrix
        dis <- vegdist(testset, method = method, ...) # similarity/distance matrix for retained sites
        x <- as.dist(dis)  #  the distance matrix
        } else {
        dis <- testdata
        x<-as.dist(as.matrix(testdata)[grps==T,grps==T]) #eliminate non selected cases from distance testdata input distance matrix
        } 
    
  N <- attributes(x)$Size  # number of cases in the dataset
  irow <- as.vector(as.dist(row(matrix(nrow = N, ncol = N))))  # from here to break original ANOSIM code
  icol <- as.vector(as.dist(col(matrix(nrow = N, ncol = N))))

  temp<-CompSeq(irow,icol,grouping,x)
  sol <- c(call = match.call())
  class(sol) <- "meandist"

  sol$dissimilarity <- attr(dis, "method")
  sol$mins <- mins
  sol$nperm <- permutations
  sol$N  <- N
  sol$classes.tested <- classes.tested
  meandist<-temp$meandist
  sol$statistic <- meandist
  sol$class.vec <- temp$cl.vec   # factor vector of between or within class categories for each distance measure!!
  sol$within <- temp$within # vector of withins by class
  sol$weighted.within.mean <- temp$weighted.within.mean
  sol$between <- temp$between
 

  if(nboot) {
    boot.results<-sapply(numeric(nboot),function(y){
    boot.N<-ceiling(boot.perc*N)#return the 'top' integer of p.perc*N
    boot.test <- sample(N, boot.N)#sample only p.perc of the values
    sboot<-sort(boot.test)#there is no need to test the group effect here ! that's 'permutations' job
    boot.grouping<-factor(grouping[sboot])
    #elimination of groups with only one member
    boot.grps<-boot.grouping
    levels(boot.grps)<-table(boot.grps)>1
    boot.grouping<-factor(boot.grouping[boot.grps==T])
    sboot<-sboot[boot.grps==T]
    boot.N<-length(boot.grouping)
    #end of the elimination of groups with only one member
    boot.x<-as.dist(as.matrix(x)[sboot,sboot]) #eliminate non selected cases from x
    boot.irow <- as.vector(as.dist(row(matrix(nrow = boot.N, ncol = boot.N))))
    boot.icol <- as.vector(as.dist(col(matrix(nrow = boot.N, ncol = boot.N))))
    boot.temp<-CompSeq(boot.irow,boot.icol,boot.grouping,boot.x)
    xx<-boot.temp$meandist   #store iteration  statistic in vector.
    })

    sol$boot<-boot.results
    pval <- (1 - cboot)/2
    sol$CI<-quantile(boot.results, probs=c(pval,1-pval), na.rm = TRUE, names=F)
    }
  if (permutations) {
    perm.results<-sapply(numeric(permutations),function(y){
      p <- sample(N, N)
      p.grouping <- grouping[p] #  this permutes the groups but keeps the number of cases in each group contstant
      p.temp<-CompSeq(irow,icol,p.grouping,x)
      y<-p.temp$meandist   #store itteration  statistic in vector
      })
    sol$perm<-as.randtest(perm.results,meandist,alter="greater")
    sol$p <- sol$perm$pvalue
    sol$perms <- sol$perm$sim
    }

  sol
}  # end meandist

mddendro <- function(meandist.object, main=NULL, xlim=NULL, xlab=NULL)
{
 # The function plots mean dissimilarity dendrograms (Van Sickle 1997) from a meandist object. 
 #Van Sickle's dendrograms were based on similarity and thus, plots using measures of association 
 #such as Bray Curtis will be the mirror image of these. However, mean distance dendrograms also enable 
 #distance based association measures to be used such as Euclidean distance (see Snelder et al. 2005).
# Van Sickle, J. 1997. Using mean similarity dendrograms to evaluate classifications.
# Journal of Agricultural, Biological and Environmental Statistics 2:370-388.

#Snelder, T.H., Woods, R. and Biggs, B.J.F., 2005. Improved eco-hydrological classification of rivers. 
# River Research and Applications, 21: 609-628.      

  between <-  as.numeric(meandist.object$between)
  within <- as.numeric(meandist.object$within)
  
  plotpos <- seq(1,length(within))
  class.labs <- names(meandist.object$within)
  
  ob.name <- deparse(substitute(meandist.object))
  
  if (is.null(main)) main <- paste(ob.name, "mean distance dendrogram")
  if (is.null(xlab)) xlab <- paste(meandist.object$dissimilarity, "Distance", sep=" ")
  
  if (is.null(xlim)) {  
  x0 <- min(within, between)
  x1 <-max(within, between)
  xd<- max(between-x0, x1-between)
   xlim<-c(between-xd -(diff(range(within,between))*0.25), between+xd+(diff(range(within,between))*0.25))
  }
  
  plot(range(within), c(0,(length(within)+1)), col="white", bty = 'n',cex.lab = 1.2, yaxt="n",
  ylab=" ", xlab=xlab, xlim=xlim,
  main=main)
  points( c(between, between),c(0,(length(within)+1)) , type="l", col = "black", lwd = 1)

  for (i in 1:length(within)) arrows(between,plotpos[i], within[i],plotpos[i], length=.05,angle=90,code=3,  lwd = 1)
     
  lable.posH <- as.numeric(within>between)
  lable.posL <- as.numeric(within<between)
      
  lab.pos <- apply(cbind( within*lable.posH + lable.posH*(diff(range(within,between))*0.2),   # compute lable positions
  within*lable.posL - lable.posL*(diff(range(within,between))*0.2) ) , 1, max)

  text(x=lab.pos, y= plotpos, labels= class.labs, cex=1.2)  
  } # end function


homog <- function(testdata, groups, method="euclidean", mins=2, permutations=F, nboot=F, boot.perc=0.9, cboot=0.95, ...)
{

# Function homog calculates homogeneity statistics (Bedward et al 1992) for a single level of classification detail

# Homogeneity is  1-W/T, where W is the mean distance between sites within the same group and
# T is the mean dissimilarity across the entire dataset (without
# any grouping). For any given set of evaluation sites, T remains the same
# regardless of the classification being tested. Homogeneity ranges
# between 0 (no discrimination between classes) and 1 (perfect discrimination).

# The homogeneity statistic is calculated for multivariate data
# where cases are grouped into several groups, such as a classification.
# homog uses functions from the VEGAN package and this package is loaded in order to use the function.

# The homog statistic is 1 - the mean within group distance divided the mean distance across all cases in the dataset
# The homog statistic  is calculated with and without weighting. The function returns "statistic" for the weighted version and Homo for the unweghted version.
# If weighted, the mean within class distances are multiplied by the number of cases in the class (n) divided by the total number of cases (N).

# A permutation test is made on the data by permuting the cases
# between the groups while holding n for each group constant.

# The function allows the minimum number of cases in each group to be set with the argument "mins".
# This reduces the possibilty that extreme values in groups with small numbers of cases will affect the statistic.
# The function returns several components

# The function also calculates confidence intervals for the homo statistic (weigthed version) using a bootstrap procedure.
# Bootstrapping with replacement is not effective for distance matrices because the distance between a sample and
# its replicate is zero. The method of Goslee and Urban (2007) is therefore used where bootstrap samples are
# constructed by randomly selecting some proportion of the data (set by argument  boot.perc) without replacement. The default value for boot.perc is 90% as recommended by  Goslee and Urban (2007). Meandist is calculated for some number of replicates of such samples
# (set with the argument nboot) and the 5th and 95th confidence intervals are calculated from the distribution of values of the statistic.

# NOTE that "enchain" is an additional utility that calculates homogeneity statistics for several levels of classification detail (or classifications)

# Bedward, M. et al (1992) Homogeneity analysis: assessing the utility of classifications and maps
# of natural resources. Australian Journal of Ecology 17:133-139.

# Goslee, S. C., and D. L. Urban. 2007. The ecodist Package for Dissimilarity-based Analysis of Ecological Data.
# Journal of Statistical Software 22

require("vegan")

  sol <- c(call = match.call())
  class(sol) <- "homog"

  #function tests
  if(mins<2) stop("mins must be >1")
  if(!is.factor(groups)) stop("groups must be a factor vector")
  if(boot.perc>1 | boot.perc<0) stop("p.perc must by a value between 0 and 1") #NB boot.perc is considered as 0 if FALSE

  if(class(testdata)!= "dist") {
  if(nrow(testdata)!=length(groups)) stop("Wrong data dimensions")
  } else {
  if(length(testdata) != (length(groups)*(length(groups)-1))/2)   stop("Wrong data dimensions")
  }   # end tests


            #inner-functions
          matched <- function(irow, icol, grouping){   #  extracting within and between
            grouping[irow] == grouping[icol]
          }

          CompSeq<-function(Thisirow,Thisicol,Thisgrouping,Thisx){
            within <- matched(Thisirow, Thisicol, Thisgrouping)  # function above - assigns  distances to either within or between
            cl.vec <- rep("Between", length(Thisx))  # makes a vector of length N (= no. distances)all entries are "between"
            take <- as.numeric(Thisirow[within])   # take is a vector of L equal to number of withins
            cl.vec[within] <- levels(Thisgrouping)[Thisgrouping[take]]
            cl.vec <- factor(cl.vec, levels = c("Between", levels(Thisgrouping)))   # factor vector of between or within class categories for each distance measure!!
            # add  code here to compute MEANDIST
            group.num <- as.vector(table(Thisgrouping))  # numeric vector of No. sites by group(class)
            ThisN<-sum(group.num)
            group.means <- tapply(Thisx, cl.vec, mean)   # calculate mean distances by group
            between <- group.means[1]   # extract the mean between groups distance from group.means
            within <- group.means[2:length(group.means)]  # extract all the within group distances from the group.means vector
            weighted.within.mean <- sum(within*(group.num/ThisN),na.rm=T) # calculate the weighted average within group distance among groups with more than 1 member
            meandist <- between - weighted.within.mean       # statistic for Similarities
            stats<-list(within=within,meandist=as.numeric(meandist),cl.vec=cl.vec,weighted.within.mean=weighted.within.mean,between=between)
            return(stats)
          }

          as.randtest<-function (sim, obs, alter = c("greater", "less", "two-sided"),call = match.call()){
            res <- list(sim = sim, obs = obs)
            res$alter <- match.arg(alter)
            res$rep <- length(sim)
            res$expvar <- c(Std.Obs = (res$obs - mean(sim))/sd(sim),
                Expectation = mean(sim), Variance = var(sim))
            if (res$alter == "greater") {
                res$pvalue <- (sum(sim >= obs) + 1)/(length(sim) + 1)
            }
            else if (res$alter == "less") {
                res$pvalue <- (sum(sim <= obs) + 1)/(length(sim) + 1)
            }
            else if (res$alter == "two-sided") {
                sim0 <- abs(sim - mean(sim))
                obs0 <- abs(obs - mean(sim))
                res$pvalue <- (sum(sim0 >= obs0) + 1)/(length(sim) +
                    1)
            }
            res$call <- call
            class(res) <- "randtest"
            return(res)
          }
        #end of inner-functions


  if(class(testdata)!= "dist") {
  dis <- vegdist(testdata, method = method, ...) # distance matrix for ALL sites
  } else {
  dis <- testdata
  }
  Te <- mean(as.vector(dis))  # mean dissimilarity all sites
  
  
     # remove cases in classes with n fewer than "mins"
  grps<-groups
  levels(grps)<-table(groups)>=mins
  grouping<-factor(groups[grps==T])
  classes.tested <- length(table(grouping))
  x<-as.dist(as.matrix(dis)[grps==T,grps==T]) #eliminate non selected cases from distance testdata input distance matrix


  N <- attributes(x)$Size  # number of cases in the dataset
  irow <- as.vector(as.dist(row(matrix(nrow = N, ncol = N))))  #
  icol <- as.vector(as.dist(col(matrix(nrow = N, ncol = N))))

  temp<-CompSeq(irow,icol,grouping,x)

  sol$dissimilarity <- attr(dis, "method")
  sol$mins <- mins
  sol$nperm <- permutations
  sol$N  <- N
  sol$classes.tested <- classes.tested
  sol$within <- temp$within # vector of withins by class
  sol$weighted.within <- temp$weighted.within.mean
  sol$between <- temp$between
  sol$statistic <- 1- (temp$weighted.within.mean/Te)
  sol$homo <- 1- ((mean(temp$within))/Te)
  sol$Te <- Te

  if(nboot) {
    boot.results<-sapply(numeric(nboot),function(y){
    boot.N<-ceiling(boot.perc*N)#return the 'top' integer of p.perc*N
    boot.test <- sample(N, boot.N)#sample only p.perc of the values
    sboot<-sort(boot.test)#there is no need to test the group effect here ! that's 'permutations' job
    boot.grouping<-factor(grouping[sboot])
    #elimination of groups with only one member
    boot.grps<-boot.grouping
    levels(boot.grps)<-table(boot.grps)>1
    boot.grouping<-factor(boot.grouping[boot.grps==T])
    sboot<-sboot[boot.grps==T]
    boot.N<-length(boot.grouping)
    #end of the elimination of groups with only one member
    boot.x<-as.dist(as.matrix(x)[sboot,sboot]) #eliminate non selected cases from x
    boot.irow <- as.vector(as.dist(row(matrix(nrow = boot.N, ncol = boot.N))))
    boot.icol <- as.vector(as.dist(col(matrix(nrow = boot.N, ncol = boot.N))))
    boot.temp<-CompSeq(boot.irow,boot.icol,boot.grouping,boot.x)
    xx<-1- ((mean(boot.temp$weighted.within.mean))/Te)   #store iteration  statistic in vector.
    })

    sol$boot<-boot.results
    pval <- (1 - cboot)/2
    sol$CI<-quantile(boot.results, probs=c(pval,1-pval), na.rm = TRUE, names=F)
    }

  if (permutations) {
   y.wt<-vector("numeric", length=permutations)
   y <-vector("numeric", length=permutations)
  for (j in 1:permutations)  {
      p <- sample(N, N)
      p.grouping <- grouping[p] #  this permutes the groups but keeps the number of cases in each group contstant
      p.temp<-CompSeq(irow,icol,p.grouping,x)
      y.wt[j] <- 1- (p.temp$weighted.within.mean/Te)   #store itteration  statistic in vector
      y[j] <- 1- ((mean(p.temp$within))/Te)                     #store itteration  statistic in vector
      }
    sol$perms <- y.wt
    sol$rand.homo <- y
    sol$p <- (sum(as.numeric(y.wt>sol$statistic))+1)/(permutations+1)  # compute p!
    }

  sol
  }    # end function homog


anosift <- function(testdata, groups, method="euclidean", mins=2, permutations=F, nboot=F, boot.perc=0.9, cboot=0.95, ...)
{

# This function performs an ANOSIM test after removing cases in classes with fewer replicates than "mins".
# The function also calculates confidence intervals for the homo statistic (weigthed version) using a bootstrap procedure.
# Bootstrapping with replacement is not effective for distance matrices because the distance between a sample and
# its replicate is zero. The method of Goslee and Urban (2007) is therefore used where bootstrap samples are
# constructed by randomly selecting some proportion of the data (set by argument  boot.perc) without replacement. The default value for boot.perc is 90% as recommended by  Goslee and Urban (2007). Meandist is calculated for some number of replicates of such samples
# (set with the argument nboot) and the 5th and 95th confidence intervals are calculated from the distribution of values of the statistic.

# testdata can be a distance object or raw dataframe

# NOTE that "enchain" is an additional utility that calculates ANOSIM R statistics for several levels of classification detail (or classifications)




sol <- c(call = match.call())
class(sol) <- "anosift"

#function tests
  if(mins<2) stop("mins must be >1")
  if(!is.factor(groups)) stop("groups must be a factor vector")
  if(boot.perc>1 | boot.perc<0) stop("p.perc must by a value between 0 and 1") #NB boot.perc is considered as 0 if FALSE

  if(class(testdata)!= "dist") {
  if(nrow(testdata)!=length(groups)) stop("Wrong data dimensions")
  } else {
  if(length(testdata) != (length(groups)*(length(groups)-1))/2)   stop("Wrong data dimensions")
  }  # end tests

  # remove cases in groups with replicates  fewer than "mins"
  grps<-groups
  levels(grps)<-table(groups)>=mins
  grouping<-factor(groups[grps==T])
  classes.tested <- length(table(grouping))

  if(class(testdata)!= "dist") {
  testset<-testdata[grps==T,] #eliminate non selected cases from  testdata input data matrix
  dis <- vegdist(testset, method = method, ...) # distance matrix for retained sites
  x <- as.dist(dis)  #  the distance matrix
  } else {
  dis <- testdata
  x<-as.dist(as.matrix(testdata)[grps==T,grps==T]) #eliminate non selected cases from distance testdata input distance matrix
  }

  N <- attributes(x)$Size  # number of cases in the dataset

  test <- anosim(x, grouping, permutations=permutations)


    if(nboot) {
    boot.results<-sapply(numeric(nboot),function(y){
    boot.N<-ceiling(boot.perc*N)#return the 'top' integer of p.perc*N
    boot.test <- sample(N, boot.N)#sample only p.perc of the values
    sboot<-sort(boot.test)#there is no need to test the group effect here ! that's 'permutations' job
    boot.grouping<-factor(grouping[sboot])
    #elimination of groups with only one member
    boot.grps<-boot.grouping
    levels(boot.grps)<-table(boot.grps)>1
    boot.grouping<-factor(boot.grouping[boot.grps==T])
    sboot<-sboot[boot.grps==T]
    boot.N<-length(boot.grouping)
    #end of the elimination of groups with only one member
    boot.x<-as.dist(as.matrix(x)[sboot,sboot]) #eliminate non selected cases from x

    boot.temp<-anosim(boot.x,boot.grouping,permutations=0)
    xx<- boot.temp$statistic   #store iteration  statistic in vector.
    })

    sol$boot<-boot.results
    pval <- (1 - cboot)/2
    sol$CI<-quantile(boot.results, probs=c(pval,1-pval), na.rm = TRUE, names=F)
    }

 sol$test
 sol$statistic <- test$statistic
 sol$classes.tested <- classes.tested
 sol$N  <- N
 sol$permutations <- test$permutations
 sol$perms <- test$perm
 sol$p <- test$signif
 sol
 }
 
 

mansift<- function(testdata, groups, method="euclidean", mins=2, permutations=F, nboot=F, boot.perc=0.9, cboot=0.95, ...)
{

 # This function performs a  mantel test on a classification by converting the class categories to a distance matrix (Zar. 1984)

   sol <- c(call = match.call())
  class(sol) <- "mansift"

 #function tests
  if(mins<2) stop("mins must be >1")
  if(!is.factor(groups)) stop("groups must be a factor vector")
  if(boot.perc>1 | boot.perc<0) stop("p.perc must by a value between 0 and 1") #NB boot.perc is considered as 0 if FALSE

  if(class(testdata)!= "dist") {
  if(nrow(testdata)!=length(groups)) stop("Wrong data dimensions")
  } else {
  if(length(testdata) != (length(groups)*(length(groups)-1))/2)   stop("Wrong data dimensions")
  }     # end tests

  grps<-groups   # remove cases in classes with n fewer than "mins"
  levels(grps)<-table(groups)>=mins
  grouping<-factor(groups[grps==T])
  classes.tested <- length(table(grouping))
  names.classes.tested <- names(table(grouping))

  if(class(testdata)!= "dist") {
  testset<-testdata[grps==T,]     #eliminate non selected cases from  testdata input data matrix
  x.dis <- vegdist(testset, method = method, ...) # distance matrix for retained sites
  } else {
  dis <- testdata
  x.dis <- as.dist(as.matrix(testdata)[grps==T,grps==T]) #eliminate non selected cases from distance testdata input distance matrix
  }

    N<- attributes(x.dis)$Size  # number of cases in the dataset

  c.dummy <- matrix(nrow=attributes(x.dis)$Size, ncol=classes.tested)   # convert classes for selected sites to  dummy variables
for (i in 1:classes.tested) {
  c.dummy[,i] <- ifelse(grouping == names.classes.tested[i], 1, 0)
  }
  c.dis<- vegdist(c.dummy, method="bray")

    # use mantel from vegan!
  vegan.mantel <- function(...)
   get("mantel", grep("package:vegan$", search()))(...)

    result <- vegan.mantel(x.dis, c.dis, method="pearson", permutations=permutations)

   if(nboot) {
    boot.results<-sapply(numeric(nboot),function(y){
    boot.N<-ceiling(boot.perc*N)#return the 'top' integer of p.perc*N
    boot.test <- sample(N, boot.N)#sample only p.perc of the values
    sboot<-sort(boot.test)#there is no need to test the group effect here ! that's 'permutations' job
    boot.grouping<-factor(grouping[sboot])
    #elimination of groups with only one member
    boot.grps<-boot.grouping
    levels(boot.grps)<-table(boot.grps)>1
    boot.grouping<-factor(boot.grouping[boot.grps==T])
    sboot<-sboot[boot.grps==T]
    boot.N<-length(boot.grouping)
    #end of the elimination of groups with only one member
    boot.x<-as.dist(as.matrix(x.dis)[sboot,sboot]) #eliminate non selected cases from x

    boot.temp <-  vegan.mantel(boot.x, as.dist(as.matrix(c.dis)[sboot,sboot]), method="pearson", permutations=0)

    xx<- boot.temp$statistic   #store iteration  statistic in vector.
    })

   sol$boot<-boot.results
    pval <- (1 - cboot)/2
    sol$CI<-quantile(boot.results, probs=c(pval,1-pval), na.rm = TRUE, names=F)
    }


  sol$dissimilarity <- attr(x.dis, "method")
  sol$mins <- mins
  sol$nperm <- permutations
  sol$N  <- N
  sol$classes.tested <- classes.tested
  sol$statistic <- as.numeric(result$statistic)
  sol$perms <- result$perm  #A vector of permuted values.
  sol$p <- as.numeric(result$signif)
  sol$result # return the object from the vegdist mantel function

  sol
  } # end mansift 

partial.mansift<- function(testdata, groups, coords, method="euclidean", mins=2, permutations=F, nboot=F, boot.perc=0.9, cboot=0.95, ...)
{
# function "p.mansift" performs a PARTIAL mantel test on a classification by converting the class categories to a distance matrix (Zar. 1984)
# coordinates are provided by argument "coords" and the effect of intersite distance is partialled out.
require("ecodist")

   sol <- c(call = match.call())
  class(sol) <- "p.mansift"

 #function tests
  if(mins<2) stop("mins must be >1")
  if(!is.factor(groups)) stop("groups must be a factor vector")
  if(boot.perc>1 | boot.perc<0) stop("p.perc must by a value between 0 and 1") #NB boot.perc is considered as 0 if FALSE

  if(class(testdata)!= "dist") {
  if(nrow(testdata)!=length(groups)) stop("Wrong data dimensions")
  } else {
  if(length(testdata) != (length(groups)*(length(groups)-1))/2)   stop("Wrong data dimensions")
  }

 if (!is.null(coords)) {       # perform checks on coords
    if(class(coords)!= "dist") {
  if(nrow(coords)!=length(groups)) stop("Wrong data dimensions for coords")
  } else {
  if(length(coords) != (length(groups)*(length(groups)-1))/2)   stop("Wrong data dimensions for coords")
  }
}   # end tests


  grps<-groups   # remove cases in classes with n fewer than "mins"
  levels(grps)<-table(groups)>=mins
  grouping<-factor(groups[grps==T])
  classes.tested <- length(table(grouping))
  names.classes.tested <- names(table(grouping))

  if(class(testdata)!= "dist") {
  testset<-testdata[grps==T,]     #eliminate non selected cases from  testdata input data matrix
  x.dis <- vegdist(testset, method = method, ...) # similarity/distance matrix for retained sites
  } else {
  dis <- testdata
  x.dis <- as.dist(as.matrix(testdata)[grps==T,grps==T]) #eliminate non selected cases from distance testdata input distance matrix
  }

    N<- attributes(x.dis)$Size  # number of cases in the dataset


  c.dummy <- matrix(nrow=attributes(x.dis)$Size, ncol=classes.tested)   # convert classes for selected sites to  dummy variables
for (i in 1:classes.tested) {
  c.dummy[,i] <- ifelse(grouping == names.classes.tested[i], 1, 0)
  }
  c.dis<- vegdist(c.dummy, method="bray")

    if(class(coords)!= "dist") {   #  coordinates - perform regression and retain  residuals for testing.
    coords.testset<-coords[grps==T,] #eliminate non selected cases from  testdata input data matrix
    coords.dis <- vegdist(coords.testset, method = "euclidean") # geographic distance matrix for retained sites
    g.dis <- as.dist(coords.dis)  #  the geographic distance matrix
    } else {
    coords.dis <- coords
    g.dis <- as.dist(as.matrix(coords)[grps==T,grps==T]) #eliminate non selected cases from coords input distance matrix
    }

    pmgram.response <- pmgram(x.dis, g.dis, nperm=0, resids=T)
   x.dis<- rewrap(pmgram.response$resids)  # residuals of response on distance

    pmgram.class <-  pmgram(c.dis, g.dis, nperm=0, resids=T)
   c.dis <-  rewrap(pmgram.class$resids)      # residuals of class on distance

   sol$pmgram.response <- pmgram.response
   sol$pmgram.class <- pmgram.class

    # use mantel from vegan!
  vegan.mantel <- function(...)
   get("mantel", grep("package:vegan$", search()))(...)

    result <- vegan.mantel(x.dis, c.dis, method="pearson", permutations=permutations)

   if(nboot) {
    boot.results<-sapply(numeric(nboot),function(y){
    boot.N<-ceiling(boot.perc*N)#return the 'top' integer of p.perc*N
    boot.test <- sample(N, boot.N)#sample only p.perc of the values
    sboot<-sort(boot.test)#there is no need to test the group effect here ! that's 'permutations' job
    boot.grouping<-factor(grouping[sboot])
    #elimination of groups with only one member
    boot.grps<-boot.grouping
    levels(boot.grps)<-table(boot.grps)>1
    boot.grouping<-factor(boot.grouping[boot.grps==T])
    sboot<-sboot[boot.grps==T]
    boot.N<-length(boot.grouping)
    #end of the elimination of groups with only one member
    boot.x<-as.dist(as.matrix(x.dis)[sboot,sboot]) #eliminate non selected cases from x

    boot.temp <-  vegan.mantel(boot.x, as.dist(as.matrix(c.dis)[sboot,sboot]), method="pearson", permutations=0)

    xx<- boot.temp$statistic   #store iteration  statistic in vector.
    })

   sol$boot<-boot.results
    pval <- (1 - cboot)/2
    sol$CI<-quantile(boot.results, probs=c(pval,1-pval), na.rm = TRUE, names=F)
    }


  sol$dissimilarity <- attr(x.dis, "method")
  sol$mins <- mins
  sol$nperm <- permutations
  sol$N  <- N
  sol$classes.tested <- classes.tested
  sol$statistic <- as.numeric(result$statistic)
  sol$perms <- result$perm  #A vector of permuted values.
  sol$p <- as.numeric(result$signif)
  sol$result # return the object from the vegdist mantel function

  sol
  } # end partial.mansift
  
enchain <- function(testdata, classes, method="euclidean", NumClass=NULL, stat=meandist, mins=2, permutations=F, nboot=F, boot.perc=0.9, cboot=0.95, coords=NULL)  {
# "enchain" is a utility that calculates either ANOSIM, mean dissimilarity analysis or homogeneity analysis or Mantel test statistics for
# several levels of classification detail (or classifications)
# the argument "classes" is a dataframe or matrix whose columns are the class memebr assignments for the test sites
# NumClass is a vector whose values specify the number of classes for each classification specified by coloms in "classes
# the order of values in NumClass must match the corresponding order of classifications repreneted by "classes"
# If NumClass is not specified the number of classes is derived from the colomn itself
# if statistic is the partial mantel test (p.mansift) then coords are expected.
 
 require("vegan")
 
  sol <- c(call = match.call())
  class(sol) <- "enchain"
  
  # check for illogical NumClass
if (!is.null(NumClass)) {
rep.classes <- apply(classes, 2, function(x) length(table(x)))
  if(length(NumClass)!= length(rep.classes)) {
  cat("The number of classifications (or classification levels) specified by the 'NumClass' argument does not", "\n",
  "correspond to the number of classifications specified by the 'classes' argument.", "\n")
  stop() }
  if(sum(as.numeric(NumClass<rep.classes))>0) {
   bads <- as.character(which(NumClass<rep.classes))
   cat("The number of classes specified by the 'NumClass' argument is exceeded by the number of classes specified in the 'classes'", "\n",
   "argument for columns specified below.", "\n",
   "Bad columns are ", print(bads, sep=" "), "\n")
  stop()
  }
}

classes.m <- as.matrix(classes)
No.levels <- ncol(classes.m)
N.cases <- nrow(classes.m)

results <- as.data.frame(matrix(nrow=No.levels, ncol=7))
row.names(results) <- names(classes)
names(results) <- c("N", "classes.tested", "Statistic", "CI_L", "CI_U", "NumClasses", "p")
             
if(permutations!=F) {
perm.statistic <- matrix(nrow=permutations, ncol=No.levels)
}

  if(class(testdata)!= "dist") {
  testdata <- vegdist(testdata, method = method)  # similarity/distance matrix for retained sites
  }
  
for (i in 1:No.levels) {    # do for each classification being tested

if (!is.null(NumClass)) {
K <- NumClass[i]
} else {
K <- length(table(classes.m[,i]))  # number of classes
keep.NumClass <- length(table(classes.m[,i]))
}
groups <- as.factor(classes.m[,i])


classes.testable <- sum(as.numeric(table(groups) >= mins))

 if(classes.testable>1) {

    if (!is.null(coords)) {
    this.lev <- stat(testdata=testdata, groups, coords=coords, method=p.mansift, mins=mins, permutations=permutations, nboot=nboot, boot.perc=boot.perc, cboot=cboot)
    } else {
    this.lev <- stat(testdata=testdata, groups, method=method, mins=mins, permutations=permutations, nboot=nboot, boot.perc=boot.perc, cboot=cboot)
    }


  results[i,1] <- this.lev$N
  results[i,2]<- this.lev$classes.tested
  results[i,3] <- this.lev$statistic

  if(nboot) {
  results[i,4:5] <- this.lev$CI
  }

  if (!is.null(NumClass)) {
 results[i,6] <-  NumClass[i]
 } else {
 results[i,6] <- keep.NumClass
 }
    
  if(permutations!=F) {
 perm.statistic[,i] <- this.lev$perms
 results[i,7]  <- this.lev$p
 }

cat("completed test", i, "of", names(classes[i]), "-class level", "statistic =", this.lev$statistic, "\n")
} else {
 results[i,1] <- 0
 results[i,2]<- "NT"
 results[i,3] <- "NT"
 cat("Column", i, "being", names(classes[i]), "has insufficient replicates to test - try changing mins", "\n")
 }


 }  # end class loop

 if(permutations!=F) {
 sol$perm.statistic <- perm.statistic
 }
 
 if (!is.null(coords)) {
 sol$pmgram.response<-this.lev$pmgram.response
 }

 sol$stat <- class(this.lev)
 sol$mins <- mins
 sol$results <- results
 sol$No.levels <- No.levels
 sol$permutations <- permutations
 sol$N.cases <- N.cases
 sol$nboot <- nboot

  sol
} # end enchain function 


plot.enchain <- function(enchain.object1, enchain.object2=NULL, enchain.object3=NULL, pos="bottomright", xmax=NULL, ymin=NULL) {
 # plot function for  "enchain"  objects
 marge <- 0.05    # adds margin to maximum values on Y axis

if (!is.null(ymin))  boty<-ymin else boty<-0 
    
if (!is.null(enchain.object2))   {

 if((enchain.object1$stat!= enchain.object2$stat) & 
     (enchain.object1$stat != "mansift" & enchain.object2$stat != "p.mansift")  &
     (enchain.object2$stat!= "mansift" & enchain.object1$stat != "p.mansift") 
     )
     stop ("enchain objects must have the same statistics")

 ob.name1 <- deparse(substitute(enchain.object1))
 ob.name2 <- deparse(substitute(enchain.object2))


 if(enchain.object1$stat=="meandist") {tit<-paste(ob.name1, "and", ob.name2); subt <-paste("Tested with Mean Dissimilarity Analysis"); ylable="CS"}
 if(enchain.object1$stat=="anosift") {tit<-paste(ob.name1, "and", ob.name2); subt <-paste("Tested with ANOSIM Global R")  ; ylable="ANOSIM R"}
 if(enchain.object1$stat=="homog") {tit<-paste(ob.name1, "and", ob.name2); subt <-paste("Tested with Homogeneity Analysis") ; ylable="Weighted Homogeneity 1-(W/T)"}
 if(enchain.object1$stat=="mansift") {tit<-paste(ob.name1, "and", ob.name2); subt <-paste("Tested with Mantel Test") ; ylable="Mantel R"}
 if(enchain.object1$stat=="p.mansift") {tit<-paste(ob.name1, "and", ob.name2); subt <-paste("Tested with Mantel Test") ; ylable="Mantel R"}


 good1 <-  enchain.object1$results[1] > 0              # cut off rows for which tests could not be made due to insufficient data
 results.good1 <-  enchain.object1$results[good1, ]
 perm.stat.good1 <- enchain.object1$perm.statistic[, good1]

 good2 <-  enchain.object2$results[1] > 0              # cut off rows for which tests could not be made due to insufficient data
 results.good2 <-  enchain.object2$results[good2, ]
 perm.stat.good2 <- enchain.object2$perm.statistic[, good2]
    
    
  if(enchain.object1$nboot & enchain.object2$nboot) {     
  yrange <- as.numeric(c(min(boty,results.good1$Statistic, results.good2$Statistic,results.good1$CI_L, results.good2$CI_L),
  max(results.good1$CI_U, results.good2$CI_U)))
  ylimt <- yrange +  yrange*c(-marge, marge)
  }  else {
  yrange = as.numeric(c(min(boty,results.good1$Statistic, results.good2$Statistic),
  max(0,results.good1$Statistic, results.good2$Statistic)))
  ylimt <- yrange +  yrange*c(-marge, marge)
  }
  if (!is.null(xmax)) {
  xlimt<-c(0, xmax) 
  } else {
  xlimt <- c(0, max(results.good1$NumClass, results.good2$NumClass,  na.rm = TRUE))
  }
   
if (!is.null(enchain.object3))   {  # set up for enchain.object3

 if ((enchain.object2$stat!= enchain.object3$stat) & 
     (enchain.object2$stat != "mansift" & enchain.object3$stat != "p.mansift")  &
     (enchain.object3$stat!= "mansift" & enchain.object2$stat != "p.mansift") 
     )
     stop ("enchain objects must have the same statistics")
     
 ob.name3 <- deparse(substitute(enchain.object3))   
  tit <- paste(ob.name1, "and", ob.name2, "and", ob.name3)
  
 good3 <-  enchain.object3$results[1] > 0              # cut off rows for which tests could not be made due to insufficient data
 results.good3 <-  enchain.object3$results[good3, ]
 perm.stat.good3 <- enchain.object3$perm.statistic[, good3]
  
   if(enchain.object1$nboot & enchain.object2$nboot & enchain.object3$nboot) {
  yrange <- as.numeric(c(min(boty,results.good1$Statistic, results.good2$Statistic, results.good3$Statistic, results.good1$CI_L, results.good2$CI_L, results.good3$CI_L),
  max(results.good1$CI_U, results.good2$CI_U, results.good3$CI_U)))
  ylimt <- yrange +  yrange*c(-marge, marge)
  }  else {
  yrange = as.numeric(c(min(boty,results.good1$Statistic, results.good2$Statistic, results.good3$Statistic),
  max(0,results.good1$Statistic, results.good2$Statistic, results.good2$Statistic)))
  ylimt <- yrange +  yrange*c(-marge, marge)
  }
 
   if (!is.null(xmax)) {
  xlimt<-c(0, xmax) 
  } else {
  xlimt <- c(0, max(results.good1$NumClass, results.good2$NumClass, results.good3$NumClass,  na.rm = TRUE))
  }
 }   # end of enchain.object3 set up
 
 x11()
  par(mfrow=c(1,1))
  
 plot(results.good1$NumClass, results.good1$Statistic, type="b", col="red", lwd=2,
 ylim=ylimt, xlim=xlimt,
   xlab ="Classification Level (Number of Classes)", ylab=ylable)
   points(results.good2$NumClass, results.good2$Statistic, type="b", col="blue", lwd=2)
   
   title(main=tit,  sub=subt,
   cex.main = 1,   font.main= 4, col.main= "blue",
      cex.sub = 0.75, font.sub = 3, col.sub = "red")
      
    if (enchain.object1$nboot & enchain.object2$nboot) {
 arrows(results.good1$NumClass, results.good1$CI_L, results.good1$NumClass, results.good1$CI_U,
 length=.05,angle=90,code=3, col="gray50")

  arrows(results.good2$NumClass, results.good2$CI_L, results.good2$NumClass, results.good2$CI_U,
 length=.05,angle=90,code=3, col="gray50")
 }    
      
      
          if (!is.null(enchain.object3))   {
             points(results.good3$NumClass, results.good3$Statistic, type="b", col="green", lwd=2)
          
              if (enchain.object3$nboot & enchain.object2$nboot & enchain.object1$nboot) {
                 arrows(results.good3$NumClass, results.good3$CI_L, results.good3$NumClass, results.good3$CI_U,
                 length=.05,angle=90,code=3, col="gray50")
                 }
                 
               legend(x=pos,  c(ob.name1, ob.name2, ob.name3), col = c("red","blue", "green"),
               text.col = "black", lty = c(1, 1,1), pch = c(-1, -1,-1), lwd=2,
               merge = TRUE, bg = 'gray90', inset = .05)                    
             }  else {
      
     legend(x=pos,  c(ob.name1, ob.name2), col = c("red","blue"),
     text.col = "black", lty = c(1, 1), pch = c(-1, -1), lwd=2,
     merge = TRUE, bg = 'gray90', inset = .05)
     }

 
 }   else  {  # end of multiple enchain object plotting 
 
ob.name1 <- deparse(substitute(enchain.object1))

if(enchain.object1$stat=="meandist") {tit<-paste(ob.name1); subt <-paste("Tested with Mean Dissimilarity Analysis"); ylable="CS"}
 if(enchain.object1$stat=="anosift") {tit<-paste(ob.name1); subt <-paste("Tested with ANOSIM Global R")  ; ylable="ANOSIM R"}
 if(enchain.object1$stat=="homog") {tit<-paste(ob.name1); subt <-paste("Tested with Homogeneity Analysis") ; ylable="Weighted Homogeneity 1-(W/T)"}
 if(enchain.object1$stat=="mansift") {tit<-paste(ob.name1); subt <-paste("Tested with Mantel Test") ; ylable="Mantel R"}
 if(enchain.object1$stat=="p.mansift") {tit<-paste(ob.name1); subt <-paste("Tested with Partial Mantel Test") ; ylable="Mantel R"}

 good <-  enchain.object1$results[1] > 0              # cut off rows for which tests could not be made due to insufficient data
 results.good <-  enchain.object1$results[good, ]
 perm.stat.good <- enchain.object1$perm.statistic[, good]

    if (!is.null(xmax)) {  
  xlimt <- c(0, xmax) 
  } else {
  xlimt <- c(0, max(results.good$NumClass,na.rm = TRUE))
  } 

  if(enchain.object1$nboot) {
  yrange <- as.numeric(c(min(boty,results.good$Statistic, perm.stat.good,results.good$CI_L),
  max(results.good$CI_U)))
  ylimt <- yrange +  yrange*c(-marge, marge)
 }  else {
 yrange <- as.numeric(c(min(boty,results.good$Statistic, perm.stat.good),
  max(results.good$Statistic)))
  ylimt <- yrange +  yrange*c(-marge, marge)
  }
      
if(enchain.object1$permutations==F) {
 x11()
     split.screen(c(2,1))        # split display into two screens
     split.screen(c(1,2), screen = 2)
     screen(1)
       
plot(results.good$NumClass, results.good$Statistic, type="b", col="red", lwd=2,
ylim=ylimt, xlim=xlimt,
main=tit,
xlab ="Classification Level (Number of Classes)", ylab=ylable)


} else {
x11()
     split.screen(c(2,1))        # split display into two screens
     split.screen(c(1,2), screen = 2)
     screen(1)

plot(results.good$NumClass, results.good$Statistic, type="b", col="red", lwd=2,
ylim=ylimt, xlim=xlimt, 
main=tit,
xlab ="Classification Level (Number of Classes)", ylab=ylable)

meds<-vector("numeric", length=nrow(results.good))

for (k in 1:nrow(results.good)) {   # for each classification tested
points(rep(results.good$NumClass[k], times=nrow(perm.stat.good)), perm.stat.good[,k], type="p", col="grey40")
meds[k] <- median(perm.stat.good[,k])
}

points(results.good$NumClass, meds, type="l", lty=2)

}

if(enchain.object1$nboot) {
 arrows(results.good$NumClass, results.good$CI_L, results.good$NumClass, results.good$CI_U,
 length=.05,angle=90,code=3, col="gray50")
 }

  screen(3)
plot(results.good$NumClass, (results.good$N/enchain.object1$N.cases)*100, type="b", col="red", lwd=2,
ylim=c(0,100),
main="Sites Participating",
xlab ="Classification Level", ylab="Proportion (%)" )

 screen(4)
plot(results.good$NumClass, (as.numeric(results.good$classes.tested)/results.good$NumClasses)*100, type="b", col="red", lwd=2,
ylim=c(0,100),
main="Classes Tested",
xlab ="Classification Level", ylab="Proportion (%)" )

} # end else
} # end plot.Enchain function

evaluate <- function(
                    testdata, # the test datset can be either raw sites by variables dataframe or distance onbject
                    classes,  # matrix of sites by classes with columns representing diffrent levels of detail or classifications to be tested
                    coords=NULL,  # coordinates for sites if distance clusters and geographic regions are to be evaluated
                    NumClass=NULL,  # vector with actual number of c lasses for each grouping reprsented in "classes"
                    stat=meandist,  # the test statistic to be used, one of "anosift", "meandist" or "homog" or "mansift"
                    theta=c(1,4,5,6), # rotations of the coordinates to be tested
                    method="euclidean", # method to be used to compute test dat distances
                    clusters="pam",    # clustering method to be used
                    mins=2,   # minimum number of sites to include a class in the analysis
                    nboot=F,  # number of bootstrap samples to be taken to compute CI
                    minsize = 10,   #  control argument for tree, default of   default of 2 produces a tree that fits the data perfectly
                    georegions=FALSE  # are georegions to be fitted and tested
){

# Function evaluate performs test of classification performance compared with a close to optimal contiguous geographic partition
# testdata can be raw data or a distance matrix!
# permutations are not applied but bootsrapping can be used to estimate CIs for Statistics
# coords is optional, if included the test will include simple distance based clusters and will fit geographic regions using classification trees fitted using X and Y cooords
# rotation of cordinate space will be tested by angles specified in Theta  <- 2=90°, 4=45°, 6=30°, 15=15°
# the a post clusters are based on hierarchical classification of testdaya distances using linkage method = "clusters"
# NumClass is a vector whose values specify the number of classes for each classification specified by coloms in "classes
# the order of values in NumClass must match the corresponding order of classifications repreneted by "classes"
# If NumClass is not specified the number of classes is derived from the colomn itself
# The comparison statistic is chosen by the user and is one of "anosift", "meandist" or "homog" or "mansift"

require(vegan)
require(stats)
require(cluster)

if (georegions) {
require(tree)
}

sol <- c(call = match.call())
  class(sol) <- "eval"

classes.m <- as.matrix(classes)
No.classifns <- ncol(classes.m)
N.cases <- nrow(classes.m)

# check for illogical NumClass
if (!is.null(NumClass)) {
rep.classes <- apply(classes, 2, function(x) length(table(x)))
  if(length(NumClass)!= length(rep.classes)) {
  cat("The number of classifications (or classification levels) specified by the 'NumClass' argument does not", "\n",
  "correspond to the number of classifications specified by the 'classes' argument.", "\n")
  stop() }
  if(sum(as.numeric(NumClass<rep.classes))>0) {
   bads <- as.character(which(NumClass<rep.classes))
   cat("The number of classes specified by the 'NumClass' argument is exceeded by the number of classes specified in the 'classes'", "\n",
   "argument for columns specified below.", "\n",
   "Bad columns are ", print(bads, sep=" "), "\n")
  stop()
  }
}

results <- as.data.frame(matrix(nrow=No.classifns, ncol=12))
row.names(results) <- names(classes)
names(results) <- c("No.Classes", "Stat_Class", "SE_Class", "% tested", "Stat_Apost", "SE_Apost", "Stat_Dist", "SE_Dist", "Stat_Geo","SE_Geo", "No.Geo.Regions", "thetamax")

result.rot <- as.data.frame(matrix(NA, nrow=No.classifns, ncol = length(theta))) # keep results of coordinate rotation tests
row.names(result.rot) <- names(classes)
names(result.rot) <- as.character(theta)

  if(class(testdata)!= "dist") {
  testdata.d <- vegdist(testdata, method = method) # similarity/distance matrix for retained sites
  } else {
  testdata.d  <- testdata
  }

  if(clusters!="pam") {
 apost <- hclust(testdata.d, method = clusters)
 }

 if (!is.null(coords)) {
  coords.d <-  dist(coords, method="euclidean")

         if(clusters!="pam") {
 dist.clust <- hclust(coords.d, method = clusters)
 }
   }

    

apost.vecs <- as.data.frame(matrix(nrow=N.cases, ncol=No.classifns))
geo.vecs <- as.data.frame(matrix(nrow=N.cases, ncol=No.classifns))
dist.vecs <- as.data.frame(matrix(nrow=N.cases, ncol=No.classifns))

for (i in 1:No.classifns) {    # do for each classification being tested

if (!is.null(NumClass)) {
K <- NumClass[i]
} else {
K <- length(table(classes.m[,i]))  # number of classes
}
class.groups <- as.factor(classes.m[,i])

  if (sum(as.numeric(table(class.groups)>=mins))<2) {
  results[i,1:3] <- 0
  cat("insufficent replicates to test this classification (level)", i, "classification", "\n")
  } else {
  results[i,1] <- K
test.class <- stat(testdata.d, class.groups, method=method, mins=mins, permutations=F, nboot=nboot, boot.perc=0.9, cboot=0.95)
results[i,2] <- test.class$statistic
if(nboot) {
class.boot <- test.class$boot
results[i,3] <- sd(class.boot)/sqrt(nboot)   # calculate standard error
}
results[i,4] <- round((test.class$classes.tested/K)*100, digits = 0)
}

if(clusters=="pam") {
apost <- pam(testdata.d, k = K)
apost.groups <- as.factor(apost$clustering)

  if (!is.null(coords)) {
  dist.clust <-  pam(dist(coords, method="euclidean"), k = K)
  dist.groups <-  as.factor(dist.clust$clustering)
  dist.vecs[i] <- dist.groups
  }

 } else {
apost.groups <- as.factor(cutree(apost, k = K))
   if (!is.null(coords)) {
      dist.groups <-  as.factor(cutree(dist.clust, k = K))
      dist.vecs[i] <- dist.groups
      }
     }
apost.vecs[i] <- as.vector(apost.groups)

  if (sum(as.numeric(table(apost.groups)>=mins))<2) {
  results[i,5] <- 0
  cat("insufficent replicates in the apost classification to test the", i, "classification (level)")
  } else {

    # keep the apost groupings
test.apost <- stat(testdata.d, apost.groups, method=method, mins=mins, permutations=F, nboot=nboot, boot.perc=0.9, cboot=0.95)
results[i,5] <- test.apost$statistic
if(nboot) {
apost.boot <- test.apost$boot
results[i,6] <- sd(apost.boot)/sqrt(nboot)   # calculate standard error
}
}

    if (!is.null(coords)) {
# test  distance clusters
        test.dist <- stat(testdata.d, dist.groups, method=method, mins=mins, permutations=F, nboot=nboot, boot.perc=0.9, cboot=0.95)
        results[i,7] <- test.dist$statistic
        if(nboot) {
dist.boot <- test.dist$boot
results[i,8] <- sd(dist.boot)/sqrt(nboot)   # calculate standard error
}
}

if (georegions) {

        for (j in 1:length(theta)) {

          rot.angle <- pi/theta[j]
          rotate<-matrix(c(cos(rot.angle), -sin(rot.angle), sin(rot.angle), cos(rot.angle)),2,2)
          coords.rot <-as.matrix(coords)%*%rotate         
          geo.tree <- tree(apost.groups ~ coords.rot[,1] + coords.rot[,2],control = tree.control(nobs=N.cases, minsize = minsize, mindev = 0))
          geo.pruned <- prune.tree(geo.tree, k = NULL, best = K)    #   NOTE THIS WILL ONLY WORK FOR SUFFICIENT NUMBERS OF CASES
          temp.geo <- as.factor(as.vector(geo.pruned$where))
          
              # safety move
          if(length(table(temp.geo))<K) {
          cat("WARNING : Cannot produce the same naumber of optimal rectangular geographic regions with as there are test classes")
          cat("Try reducing the argument; minsize", "\n") 
          }

          if (sum(as.numeric(table(temp.geo)>=mins))<2) {
            result.rot[i,j] <- -999
            cat("insufficent replicates to test this classification (level)", i, "and rotation", j,  "geo", "\n")
            } else {
                  test.geo <- stat(testdata.d, temp.geo, method=method, mins=mins, permutations=F, nboot=F)
                  result.rot[i,j] <- results[i,2]/test.geo$statistic
                  }
                  }  # end rotation trials

       if(max(result.rot[i, ])==-999) {
       cat("insufficent replicates in the geo classification to test the", i, "classification (level)")
       } else {

    best.rot<- match(max(result.rot[i, ]), result.rot[i, ])
          rot.angle <- pi/theta[best.rot]
          rotate<-matrix(c(cos(rot.angle), -sin(rot.angle), sin(rot.angle), cos(rot.angle)),2,2)
          coords.rot <-as.matrix(coords)%*%rotate        
          geo.tree <- tree(apost.groups ~ coords.rot[,1] + coords.rot[,2], control = tree.control(nobs=N.cases, minsize = minsize, mindev = 0)) # problems may arise here if tree depth exceeded
   if (length(table(geo.tree$where)) < K) {
   cat("WARNING; the CART model that defines the georegions cannot produce as many regions as there are classes being tested", "\n")
    }
          geo.pruned <- prune.tree(geo.tree, k = NULL, best = K)    # NOTE  THIS WILL ONLY WORK FOR SUUFICIENT NUMBERS OF CASES
          geo.groups <- as.factor(as.vector(geo.pruned$where))
          geo.vecs[i] <- geo.groups
                  cat("using", length(table(geo.groups)), "tree based clusters with coords rotated by" ,pi/theta[best.rot], "radians as geographic regions to test", K, "classes", "\n")
                  test.geo <- stat(testdata.d, geo.groups, method=method, mins=mins, permutations=F, nboot=nboot, boot.perc=0.9, cboot=0.95)

        results[i,9] <- test.geo$statistic
          if(nboot) {
          results[i,10] <- sd(test.geo$boot)/sqrt(nboot)   # calculate standard error
          }
        results[i,11] <-   length(table(geo.groups))
        results[i,12] <-  theta[best.rot]
        } # end else
        } # end geographic regions

}   # end for loop for each classification
sol$results <- results
sol$geo.vecs <- geo.vecs
sol$apost.vecs <- apost.vecs
sol$dist.vecs <- dist.vecs
sol$coords <- coords
sol$nboot <- nboot
sol$results.rotate <- result.rot
sol$stat <- class(test.class)
sol$clusters <- clusters
sol$georegions <- georegions
sol
}  # end function


plot.evaluate <- function(evaluate.object, pos="right") {
# This function plots the statistics from an "evaluate" object
 good <-  evaluate.object$results[1] > 0              # cut off rows for which tests could not be made due to insufficient data
 results.good <-  evaluate.object$results[good, ]

      ob.name <- deparse(substitute(evaluate.object))
      x11()
      
      if  (is.null(evaluate.object$coords)) {
       plot(results.good[,1], results.good[,5]/results.good[,2], type="l", col="blue",
       lwd = 2,
       ylim=c(min(0,na.omit(results.good[,5]/results.good[,2])), max(na.omit(results.good[,5]/results.good[,2]))),
      xlab="Classification Level (No. of classes)",
      ylab=paste("Statistic for alternative classification/test classification"),
      main= paste(ob.name, "evaluated using", evaluate.object$stat,"and" ,evaluate.object$clusters))
        abline(h=1,  col = "grey70")
      legend(pos,  c("a posteriori classification"), col = c("blue"),
                  text.col = "black", lwd = 2, pch = -1,
                  merge = TRUE,  inset = .05)
      } else {
      
      if (evaluate.object$georegions==TRUE)    {
      ylimt=c(min(c(na.omit(results.good[,7]/results.good[,2]),na.omit(results.good[,9]/results.good[,2]))),
      max(c(na.omit(results.good[,5]/results.good[,2]),na.omit(results.good[,7]/results.good[,2]), na.omit(results.good[,9]/results.good[,2]))))
      } else {
      ylimt=c(min(c(na.omit(results.good[,7]/results.good[,2]))),
      max(c(na.omit(results.good[,5]/results.good[,2]),na.omit(results.good[,7]/results.good[,2]))))
      }

      plot(results.good[,1], results.good[,5]/results.good[,2], type="l", col="blue",
      ylim=ylimt , lwd = 2,
      xlab="Classification Level (No. of classes)",
      ylab=paste("Statistic for alternative classification / test classification"),
      main= paste(ob.name, ": Evaluated using", evaluate.object$stat,"and" ,evaluate.object$clusters))
      points(results.good[,1], results.good[,7]/results.good[,2], type="l", col="red", lwd = 2)
      points(results.good[,1], results.good[,9]/results.good[,2], type="l", col="green", lwd = 2)
      abline(h=1,  col = "grey70")
      
       if (evaluate.object$georegions==TRUE)    {
      legend(title= "Alternative classifications", x=pos,  c("a posteriori classification", "distance clusters", "optimal rectangular geographic regions"), col = c("blue", "red", "green"),
                  text.col = "black", lwd = c(2, 2, 2), pch = c(-1, -1, -1),
                  merge = TRUE,  inset = .05)
       } else {
       legend(pos,  c("a posteriori classification", "distance clusters", "optimal rectangular geographic regions"), col = c("blue", "red", "green"),
                  text.col = "black", lwd = c(2, 2, 2), pch = c(-1, -1, -1),
                  merge = TRUE,  inset = .05)
                  }
                  }
       } # end function




map.evaluate <- function(evaluate.object, coords, class.level, newcoords=NULL, res=30000) {     #
 # This function produces maps from an "evaluate" object
if(is.null(evaluate.object$coords)) stop ("Cannot plot map if coordinates were not used in function evaluate")
coords <- evaluate.object$coords
ob.name <- deparse(substitute(evaluate.object))

      i<-class.level

        geo.groups <- evaluate.object$geo.vecs[,i]
        apost.groups <- evaluate.object$apost.vecs[,i]
        dist.groups <- evaluate.object$dist.vecs[,i]

              geo.classes <- as.numeric(names(table(geo.groups)))
              colour = seq(from = 1, to = length(geo.classes), by = 1)
              sym<- rep(c(19:25), times=length(geo.classes))


              dist.classes <- as.numeric(names(table(dist.groups)))
              apost.classes <- as.numeric(names(table(apost.groups)))

              x11()

              plot(coords[ ,1], coords[ ,2], col="white", xlab = "east", ylab = "north",
              main = paste(length(geo.classes), "Optimal rectangular geographic regions"),
              sub=ob.name)
              for (c in 1:length(geo.classes)) {
              points(coords[geo.groups==geo.classes[c],1], coords[geo.groups==geo.classes[c],2], col=colour[c], pch=sym[c])
              }

              x11()

                      plot(coords[ ,1], coords[ ,2], col="white", xlab = "east", ylab = "north",
                      main = paste(length(apost.classes), "A posterori groups"),
                      sub=ob.name)
                      for (c in 1:length(dist.classes)) {
                      points(coords[apost.groups==c,1], coords[apost.groups==c,2], col=colour[c], pch=sym[c])
                      }
                      
              x11()

                      plot(coords[ ,1], coords[ ,2], col="white", xlab = "east", ylab = "north",
                      main = paste(length(dist.classes), "Distance clusters"),
                      sub=ob.name)
                      for (c in 1:length(geo.classes)) {
                      points(coords[dist.groups==c,1], coords[dist.groups==c,2], col=colour[c], pch=sym[c])
                      }
                      

    if(!is.null(newcoords)) {
require(tree)

minsize = 10   # Control parameter for the classification tree used in constructing the contiguous geographic regions.

sol <- c(call = match.call())
  class(sol) <- "eval"

 if(is.null(evaluate.object$coords)) stop ("Cannot plot map if coordinates were not used in function evaluate")
ob.name <- deparse(substitute(evaluate.object))

   if(nrow(newcoords)>res) {
              pick <- sample(1:nrow(newcoords), size=res, replace = F)
              } else {
              pick<-c(1:nrow(newcoords))
              }
              newcoords.s <- newcoords[pick,]

      i<-class.level

  apost.groups <-as.factor(evaluate.object$apost.vecs[,i])
  K <- length(table(apost.groups))
  N.cases <- length(apost.groups)

          rot.angle <- pi/evaluate.object$results$thetamax[i]
          rotate<-matrix(c(cos(rot.angle), -sin(rot.angle), sin(rot.angle), cos(rot.angle)),2,2)
          coords.rot <-as.data.frame(as.matrix(evaluate.object$coords)%*%rotate)
          geo.tree <- tree(apost.groups ~ ., coords.rot, control = tree.control(nobs=N.cases, minsize = minsize, mindev = 0))
          geo.pruned <- prune.tree(geo.tree, k = NULL, best = K)    #   NOTE THIS WILL ONLY WORK FOR SUUFICIENT NUMBERS OF CASES

          temp.geo <- as.factor(as.vector(geo.pruned$where))  # note the where groups need to now be modelled!!!!
          geo.tree.where <- tree(temp.geo ~ ., coords.rot, control = tree.control(nobs=N.cases, minsize = minsize, mindev = 0))

          rotnewcoords <- as.data.frame(as.matrix(newcoords.s)%*%rotate)
          names(rotnewcoords) <- c("V1", "V2")

          predictions <- predict.tree(geo.tree.where, newdata=rotnewcoords, type="class")

    sol$predictions <- predictions   # keep these for elsewhere

              where.classes <- as.numeric(names(table(predictions)))
              colour = rep(c(1:8), length.out=length(where.classes)) # from=1, to=8, )

              # open graphics device for CV output

              x11()

              plot(c(min(newcoords[ ,1]), max(newcoords[ ,1])),c(min(newcoords[ ,2]), max(newcoords[ ,2])), col="white", xlab = "east", ylab = "north",
              main = paste(evaluate.object$results$No.Classes[i], "Optimal rectangular geographic regions"),
              sub=ob.name)
              for (c in 1:length(where.classes)) {
              points(newcoords.s[predictions==where.classes[c],1], newcoords.s[predictions==where.classes[c],2], col=colour[c], pch=20, cex=1)
              }
              sol
              }
             }  # end map.evaluate function

#################################################################################################################################################################
# dismod functions 
###############################################################################################################################################################

# function to cbind all the elements in a list of data frames
# works even when all NA's
Doug.cbind.list <- function (dfs)
{
    data.frame(do.call("cbind", lapply(dfs, function(df) df)))
}
#########

  slide.x <-  function(x)  {   # function slides vector x with min values <= ZERO to allow transformation
     mn <- min(x)
   if(mn <= 0) x.slid <- x + (abs(mn) + 1)  
   else x.slid <- x
   return(x.slid)
    } # end slide.x function

   slide <-  function(x)  {   # function slides ALL predictors with min values <= ZERO to allow transformation  
    x.slid <- as.data.frame(apply(x, 2, slide.x))
    return(x.slid)
     } # end slide function

#################################################################################################################################################################

     var.transformer <- function(v, tranf = tran) {
              v.tran <-  as.matrix(as.vector(decostand(v, method="range")))
              # seven transformations                     
              if(any(tranf%in%"CU")) v.tran <- cbind(v.tran, as.vector(decostand(v^3, method="range")))         
              if(any(tranf%in%"SQ")) v.tran <- cbind(v.tran, as.vector(decostand(v^2, method="range")))
              if(any(tranf%in%"R2")) v.tran <- cbind(v.tran, as.vector(decostand(v^0.5, method="range")))
              if(any(tranf%in%"R4")) v.tran <- cbind(v.tran, as.vector(decostand(v^0.25, method="range")))
              if(any(tranf%in%"L10")) v.tran <- cbind(v.tran, as.vector(decostand(log10(v), method="range")))
              if(any(tranf%in%"LN")) v.tran <- cbind(v.tran, as.vector(decostand(log(v), method="range")))             
              if(any(tranf%in%"2LOG")) v.tran <- cbind(v.tran, as.vector(decostand(log10(slide.x(log10(v))), method="range")))  # note that log10 of numbers smaller than 1 are -ve and need to be slid (to be >0) to be logged the 2nd time.                             
              labs <- vector("character") # lables  paste names+trans
              labs[1] <- names(v)
              if(any(tranf%in%"CU")) labs <- c(labs, paste(names(v) ,"CU", sep="_"))
              if(any(tranf%in%"SQ")) labs <- c(labs, paste(names(v) ,"SQ", sep="_"))
              if(any(tranf%in%"R2")) labs <- c(labs, paste(names(v) ,"R2", sep="_"))
              if(any(tranf%in%"R4")) labs <- c(labs, paste(names(v) ,"R4", sep="_"))
              if(any(tranf%in%"L10")) labs <- c(labs, paste(names(v) ,"L10", sep="_"))
              if(any(tranf%in%"LN")) labs <- c(labs, paste(names(v) ,"LN", sep="_"))
              if(any(tranf%in%"2LOG")) labs <- c(labs, paste(names(v) ,"2LOG", sep="_"))
              v.tran <- as.data.frame(v.tran)
              names(v.tran) <- labs
              return(v.tran)
              }
  
 trans <-  function(x, tranB=tran, PLOT=T) {   # function FIRST makes a chosen set of transformations of the predictors AND THEN range standardises
          x.tran <- as.data.frame(matrix(nrow=nrow(x), ncol=0))  # blank DF
          M <- ncol(x)  # number of variables
          if(min(x) <= 0) stop("cannot transform variables less than or equal to 0")
for (i in 1:M) x.tran <- cbind.data.frame(x.tran,  var.transformer(v=x[i], tranf = tranB))
return(x.tran)              
              } # end trans function
                           
#################################################################################################################################################################
# function returns a new distance matrix comprising a rewrapped distance vector
###############################################################################################################################################################              
              
               rewrap <- function(vec) {   # function returns a new distance matrix comprising a rewrapped distance vector
              ndis <- length(vec) # number of distances in vector
              np <- (1 + sqrt(1 + 8*ndis))/2  # np = No. of sites - number of distances = n(n-1)/2; solve for positive root of n
              # Convert vec, the vector of distances, to D.tri, an upper triangular matrix.
              D.tri <- matrix(nrow=np, ncol=np)
              i <- 1
              n1 <- 1
              n2 <- i*(np) - 1
              while(i <= np - 1)
              {
              	D.tri[i,(i + 1):np] <- vec[n1:n2]
              	n1 <- n2 + 1
              	n2 <- n2 + (np - i - 1)
              i <- (i + 1)
              }
              # Convert D.tri, the upper triangular matrix into D.sym, a symmetric matrix.
              D.sym <- matrix(nrow=np, ncol=np)
              i <- 1
              while(i <= np)
              {
              	D.sym[i,i] <- 0
              	j <- (i + 1)
              	while(j <= np)
              	{
              		D.sym[j,i] <- D.tri[i,j]
              		D.sym[i,j] <- D.tri[i,j]
              		j <- (j + 1)
              	}
              i <- (i + 1)
              }
              as.dist(D.sym)
              }   # end of rewrap  function
                
#################################################################################################################################################################
# make preds combines the functions  slide.x , slide, var.transformer, trans
###############################################################################################################################################################
                                                        
 makepreds <- function(pred, geo=F, tranA=tran) {   # this function makes the complete set of  distance pairs to be used as predictors          
          if (geo==TRUE) {   # if geographic variables are included, euclidean distances between sites are included including 5 transformed versions
          cat("........................................................", "\n")
          cat("Including geographic distances as predictors", "\n")
          coords <- pred[ ,1:2] #  coordinates
          pred <- pred[ ,3:ncol(pred)]
          geo.dis <- geodist(coords, tranB=tranA)
          }                                                   
       pred.slid <- slide(pred)   # note predictors are "slid" to allow transformations of negative numbers
       pred.trans <- trans(pred.slid, tranB = tranA)   # for each predictor compute all transformation options and range standardise (min=O, max=1)              
       # predictors are converted to distance pairs (manhattan distance of transformed and range standardised predictors)
          pred.dis <-  as.data.frame(apply(pred.trans, 2, function(x) as.vector(dist(x, method="manhattan")))) # convert to disssimilarities      
          if (geo==TRUE) pred.dis <- cbind.data.frame(geo.dis, pred.dis)          
           return(pred.dis)
           } # end makepreds function

#################################################################################################################################################################
# geodist is a function to make geographic distance pairs  
###############################################################################################################################################################
      
  geodist <- function(coords, tranB=tran)  {  # compute geographical distances including transformed versions
              geo.dist <- as.data.frame(as.vector(dist(coords, method="euclidean")))
              names(geo.dist) <- "geodist" 
              geo.d <- var.transformer(v=geo.dist+1, tranf = tranB) # note sites may have the same cordinates (e.g., due to lack of spatil resolution) thus log(x+1) to avoid NaN                    
              return(geo.d)
              }    #end geodist function   
                 
#################################################################################################################################################################
# dismod.fit function  
###############################################################################################################################################################

dismod.fit <- function(resp, pred, method="bray", totP=10, toolong = NULL, geo=F, crit=0, tran=c("CU", "SQ", "R2","R4", "L10", "LN", "2LOG"), folds=10, repCV=1)   ##fits a dissimilarity model
{
# This function fits a dissimilarity model
if(nrow(resp)!=nrow(pred)) stop("Wrong data dimensions")

  sol <- c(call = match.call())
  class(sol) <- "dismod"
  sol$pred <- pred
  sol$resp <- resp

 # initiate timing call
 z1 <- unclass(Sys.time())

          resp.nonex <-  vegdist(resp, method=method) #nonextended RESPONSE distances

          if(!is.null(toolong)) {
          cat("........................................................", "\n")
          cat("Extending distances.....", "\n")
          resp.ini.d <-  stepacross(resp.nonex, path = "extended", toolong) # compute extended RESPONSE distances with cutoff default = 0.9
          resp.ini.vec <- as.vector(resp.ini.d)  # response distances as vector
          } else {
          resp.ini.d <- resp.nonex
          resp.ini.vec <- as.vector(resp.nonex)  # NON-EXTENDED response distances as vector
          }

          pred.ini.d <- dist(decostand(pred, method="range"),method="manhattan")   # this is the range std manhattan using all predictors
          R.ini<-cor(as.vector(pred.ini.d), as.vector(resp.ini.d), method="pearson") # do initial test with input predictors
          R.ini.stat <- R.ini
          cat("........................................................", "\n")         
          cat("Untuned R using all predictors", R.ini.stat, "\n")
          sol$untunedR  <- R.ini.stat      #save result


          N <-nrow(pred)  # number of cases in the dataset
          Ndists <- (N*(N-1))/2     # no distances in the complete distance matrices
          check.var <- apply(pred, 2, function(x) sd(x))   # check variables have some variation and remove any variables with no variation
          pred.use <-  pred[, check.var > 0]
          M <- ncol(pred.use) # number of variables
                                          
pred.dis <- makepreds(pred.use, geo=geo, tran=tran) # make distance pairs of all variables plus transformations
     
# START CV
               
RepCV <- vector("list", repCV) # store output from repeats of CV 
repX <- 1                    
while (repX <= repCV)  { # control flow for repeat CV loops 

allR <- data.frame(matrix(nrow = totP, ncol = 0))   #make blank df for results of each fold
testR <- vector(mode = "numeric", length = folds) # this stores the testing R values from each fold
coeffs.cv <- vector("list", folds)
 
fold.vec <- sample(rep(c(1:folds), length.out=N))  # ramdom slpit of N cases into 10 folds

# open graphics device for CV output
         x11()
         par(mfrow=c(3,4))

    for (f in 1:folds)  {      # OUTER LOOP TO MAKE CV FOLDS

      # Select fold from pred.dis and resp.ini.d
      preds.f <- as.numeric(fold.vec != f) # new vector 1 if included in fold, 0 if in holdout set
      N.f <- sum(preds.f)  # number of cases in the fold dataset
      Ndists.f <- (N.f*(N.f-1))/2     # no distances in the fold distance matrix
      N.h <- N - N.f
      Ndists.ho <- (N.h*(N.h-1))/2     # no distances in the holdout distance matrix

      pick.vec <- as.vector(dist(preds.f, method="canberra"))  # any distance NOT equal to ZERO is to be EXCLUDED
      pick.vec <- replace(pick.vec, is.na(pick.vec), 10) # replace NAs to be replaced by 10 - these are the holdout distances

      resp.vec <- resp.ini.vec[pick.vec == 0]  # initial response
      resp.y <- resp.vec # intialize response vector

      if (sd(resp.vec) == 0) {         # safety move here could have response subset with no variabilty, hence distances will be all 0
      cat("this fold has no variabilty, skip to next subset", "\n")
      bestvarcoeff <- NULL
      testR <- rep(NA, times=folds)
      } else
      {   # ENTER MIDDLE loop to itteratively add best variable to model up to totP = total number of variables to add and TEST

          bestvar <-vector(mode = "numeric", length = totP)       #store list of vars giving best model
          bestvarcoeff <- vector(mode = "numeric", length = totP)  # store best var coeffs
          testR <- vector(mode = "numeric", length = totP)       # store test results
          fitR <- vector(mode = "numeric", length = totP)       # store R for fitted model


                 for (j in 1:totP)  {      # MIDDLE LOOP TO ADD VARIABLES TO MODEL  j=NO.predictors to add

                        varcor <- apply(pred.dis[pick.vec == 0, ], 2, function(x) { if (sd(x)==0) 0 else cor(x, resp.y)})  # safety move in case of no variabilty in predictor fold
                        varcor.pos <- as.numeric(varcor>0)  # vector to remove all NEGATIVE r
                        modvar <- as.vector(varcor*varcor.pos)  # remove variables with -ve corellation
                        bestvar[j] <- as.numeric(match(max(modvar), modvar))  #  best var for j-th itteration

                        # fit a linear model here

                        if(j==1) {  # if more than one variable is fitted need to recalculate coeffs
                            var.m <- pred.dis[pick.vec == 0, bestvar[j]] # best variable itteration 1
                            mod <- lm(resp.y ~ var.m)
                            bestvarcoeff[j] <-  as.numeric(mod$coefficients[2])
                            resp.y <- mod$residuals   # residuals of best model to date
                            }
                        else {
                            bestvarcoeff <- vector(mode = "numeric", length = totP)  # reset the best var coeffs to ZERO this step calcs a new set
                            var.m <- pred.dis[pick.vec == 0, bestvar[1:j]]  # put selected variables to date in dataframe
                            mod <- lm(resp.vec ~ .,var.m) # refit model
                            mod.coef <- mod$coefficients[2:(j+1)]  # these are the fitted coefs
                            coef.pos <- mod.coef * as.numeric(mod.coef>0)  # if the coeffs are +ve then they are retained else = 0
                            coef.pos <- replace(coef.pos, is.na(coef.pos), 0) # IF THERE ARE pred.dis with COR = 1 function lm will produce coeffs = NA, replace NAs by 0
                            bestvarcoeff[1:j] <- coef.pos[1:j]  # term may be eliminated if coeffs are less than 0
                            resp.y <- mod$residuals  # residuals of best model to date
                            } # end of coefficient refitting

                                fitR[j] <- cor(mod$fitted.values, resp.vec)  #  best R (fitted variable) for j-th itteration


                                      # test fitted model with holdout fold
                                      test.preds <- as.matrix(pred.dis[pick.vec == 10, bestvar[1:j]])   # test fold - columns to select are stored in bestvar
                                      test.preds.d.vec <- test.preds %*% bestvarcoeff[1:j]  #  distances = addition of selected pred.dis*coeff
                                      test.resp.vec <- resp.ini.vec[pick.vec == 10] # take rows and colomns from extended distance matric to test this fold

                                      if (sd(test.resp.vec) == 0) {       # safety move here, small set of sites can have all zero distances
                                      R <- 0
                                      } else {
                                      R <- as.numeric(cor(test.preds.d.vec, test.resp.vec, method="pearson"))
                                      }
                                      testR[j] <- R #save results
          }  #end middle loop
          }  # end if

          coeffs.cv[[f]]<- bestvarcoeff
          allR<-cbind.data.frame(allR,testR)
          p<-c(1:totP)

          # plot results of CV fold
          fit.max <- max(fitR)
          test.max <- max(testR)
          y.max <- max(fit.max, test.max)
          y.min <- min(fit.max, test.max,0)

          plot(p, fitR, type="b", xlab ="Number of Predictors", ylab="Correlation", ylim = c(y.min, y.max), col="blue",
          main=(paste("CV fold", f)))
          lines(p, testR, type="b", col="red")
          legend("bottomright",  c("fitted R", "test R"), col = c("blue", "red"),
            text.col = "black", lty = c(2, 2), pch = c(-1, -1),
            merge = TRUE, bg = 'gray90', inset = .05)

     } # end outer loop
          cat("........................................................", "\n")     
          cat("Completed CV number ", repX, "\n")
          cat("........................................................", "\n")              
          RepCV[[repX]] <- allR
          repX <- repX + 1   
 } # control flow repeat CV   
                  
  aveR <- apply(Doug.cbind.list(RepCV), 1,  mean, na.rm = TRUE) # find mean and SD for each combine results of all CV repeats
  se <- apply(Doug.cbind.list(RepCV), 1, function(x) sd(x, na.rm = TRUE)/sqrt(folds*repCV))  # calculate standard error

    seU <- aveR+(se) #
    seL <- aveR-(se)
    maxJ <- max(aveR) - crit*(max(aveR) - seL[match(max(aveR), aveR)]) # calculates the minimum value of the justifiable mean R
    predsJ <- match(1, as.numeric(aveR >= maxJ))  # finds index of first item in 'mean' greater than maxJ (i.e. the justifiable no. of preds)

   x11()           # plot final CV results
   upper <- max(seU)
   lower <-min(seL)
   plot(p, aveR, type="l",  xlab ="Number of Predictors", ylab="Correlation", main="Final CV results", ylim=c(lower,upper), col="red", lwd = 2) #overplot the se lines?
   arrows(p,seU,p,seL,length=.05,angle=90,code=3, col="gray70")
   abline(h=maxJ, v=predsJ, col = "black")
   legend("bottomright",  c("Mean R", "SE","Justified No. of predictors"), col = c("red", "gray70", "black"),
           text.col = "black", lty = c(2, 2), pch = c(-1, -1),
           merge = TRUE, bg = 'gray90', inset = .05)

   # fit final model

    cat("fitting final model with", predsJ, "predictors", "\n")

      bestvar.fin <-vector(mode = "numeric", length = predsJ)       #store list of vars for FINAL model
      bestvarcoeff.fin <- vector(mode = "numeric", length = predsJ)  # store  var coeffs for FINAL model
      finalR <- vector(mode = "numeric", length = predsJ)       # store test results for FINAL MODEL
      bestvar.names.fin <-vector(mode ="character", length=predsJ)

      resp.fin.y <- resp.ini.vec  # initialize response distances

     # final loop replicates the CV step above and could be made into a function!

            for (h in 1:predsJ)  {      # FINAL model OUTER LOOP TO ADD VARIABLES TO MODEL  predsJ = NO.predictors to add


                        varcor.fin <- apply(pred.dis,2, function(x) cor(x, resp.fin.y))
                        varcor.pos.fin <- as.numeric(varcor.fin>0)  # vector to remove all NEGATIVE r
                        modvar.fin <- as.vector(varcor.fin*varcor.pos.fin)  # remove variables with -ve corellation
                        bestvar.fin[h] <- as.numeric(match(max(modvar.fin), modvar.fin))  #  best var for h-th itteration
                        bestvar.names.fin[h] <- as.character(names(pred.dis[match(max(modvar.fin), modvar.fin)]))

                        # fit a linear model here

                        if(h == 1) {  # if more than one variable is fitted need to recalculate coeffs
                            var.fin.m <- pred.dis[, bestvar.fin[h]] # best variable itteration 1
                            mod.fin <- lm(resp.ini.vec ~ var.fin.m)
                            bestvarcoeff.fin[h] <-  as.numeric(mod.fin$coefficients[2])
                            resp.fin.y <- mod.fin$residuals  # residuals of best model to date
                            }
                            else {
                            bestvarcoeff.fin <- vector(mode = "numeric", length = predsJ)  # reset the best var coeffs to ZERO this step calcs a new set
                            var.fin.m <- pred.dis[, bestvar.fin[1:h]]  # put h-th selected variables into matrix
                            mod.fin <- lm(resp.ini.vec ~ .,var.fin.m) # refit model
                            mod.fin.coef <- mod.fin$coefficients[2:(h+1)]  # these are the fitted coeffs
                            coef.fin.pos <- mod.fin.coef * as.numeric(mod.fin.coef>0)  # if the coeffs are +ve then they are retained else = 0
                            coef.fin.pos <- replace(coef.fin.pos, is.na(coef.fin.pos), 0) # IF THERE ARE pred.dis with COR = 1 function lm will produce coeffs = NA, replace NAs by 0
                            bestvarcoeff.fin[1:h] <- coef.fin.pos[1:h]  # term may be eliminated if coeffs are less than 0
                            resp.fin.y <- mod.fin$residuals  # residuals of best model to date
                            } # end of coefficient refitting

                             }  #end FINAL model OUTER loop

      cat("........................................................", "\n")         
      cat("final model - selected variables", "\n")
      print(bestvar.names.fin)
      cat("........................................................", "\n")
      cat("coefficients for final variables", "\n")
      print(bestvarcoeff.fin)
      cat("........................................................", "\n")         
      cat("Final Model R (CV)", aveR[predsJ], "\n")
      cat("Final Model R (Fitted)", cor(mod.fin$fitted.values, resp.ini.vec), "\n")

z2 <- unclass(Sys.time())
elapsed.time.minutes <- round((z2 - z1)/ 60,2)  #calculate the total elapsed time

      sol$R.cv <- aveR[predsJ]
      sol$R.fit <- cor(mod.fin$fitted.values, resp.ini.vec)
      sol$fit.preds <- bestvar.names.fin
      sol$fit.coeffs <- bestvarcoeff.fin
      sol$intercept <- mod.fin$coefficients[1]
      sol$geo <- geo
      sol$method <- method
      sol$toolong <- toolong
      sol$tran <- tran
      sol$RepCV <- RepCV
      sol$elapsed.time.minutes = elapsed.time.minutes
      sol
}     #end of dismod.fit function

#################################################################################################################################################################
# dismod.plot function  
###############################################################################################################################################################

 dismod.plot  <- function(x, plot.layout = c(3, ceiling(length(x$fit.preds)/3)), keep=FALSE)   {
  # this function plots a disimilarity model
 require(vegan)
  sol <- c(call = match.call())
  class(sol) <- "dismod"

  ob.name <- deparse(substitute(x))
          resp.nonex <-  vegdist(x$resp, method=x$method) #calculated RESPONSE distances

    if(!is.null(x$toolong)) {
          cat("Extending distances.....", "\n")
          resp.d <-  stepacross(resp.nonex, path = "extended", toolong=x$toolong) # compute extended RESPONSE distances
          resp.vec <- as.vector(resp.d)  # response distances as vector
          } else {
          resp.vec <- as.vector(resp.nonex) }

           pred.dis <- dismod.transform(x, x$pred, distance=TRUE)   # predictor distances fitted by model
           mod.resp.vec <- x$intercept +  as.vector(pred.dis)   # modeled response as a vector

           minval <- min(mod.resp.vec)
           maxval <-max(mod.resp.vec)
           minvalx <- pred.dis[match(minval,mod.resp.vec)]
           maxvalx <- pred.dis[match(maxval,mod.resp.vec)]
           resplo <-c(minval, maxval)
           xplo <-c(minvalx, maxvalx)
           ylimT <- c(0, max(resplo, resp.vec)*1.1)

           x11()
             plot(as.vector(pred.dis), resp.vec, xlab ="Predictor dissimilarity", ylab="Response dissimilarity", col="gray50",
             main=paste("Plot of dissimilarity model ", ob.name), pch=".", cex=1, ylim=ylimT)
             lines(xplo, resplo , col="red", type="l", lwd = 2 )
              legend("bottomright",  c("Observed response", "Fitted response"), col = c("gray50", "red"),
               text.col = "black",  pch = c(20, 20), bg = 'gray90', inset = .05)

  # This second part of this function illustrates
  # the transformations and weights that are aplied to the predictors
     
      allFitPrds.list <- lapply(x$fit.preds, DougSplit, split = "_")  # cannot use sapply because not all vars have transformation

allFitPrds <- data.frame(nrow=length(x$fit.preds), ncol=2)
names(allFitPrds) <- c("Variable", "Transformation")
      for (i in 1:length(x$fit.preds)) {
       allFitPrds[i,1] <- allFitPrds.list[[i]][1]
       allFitPrds[i,2] <- allFitPrds.list[[i]][2]
      }
allFitPrds[,2] <- replace(allFitPrds[,2], is.na(allFitPrds[,2]), "Not transformed")

pred.faked <- as.data.frame(matrix(rep(seq(0,1, length.out = 20), nrow(allFitPrds)), ncol=nrow(allFitPrds), nrow=20))
names(pred.faked) <- allFitPrds[,1]        

pred.faked.slid <- slide(pred.faked)    # slide then transform the variables
pred.faked.trans <-  trans(pred.faked.slid, tran=x$tran)[x$fit.preds]

 pred.faked.tranWeight <- as.data.frame(matrix(NA, ncol=length(x$fit.preds), nrow=20))   # apply the model coefficients to the transformed variables
 names(pred.faked.tranWeight) <- x$fit.preds
 for (a in 1:length(x$fit.preds))  pred.faked.tranWeight[,a] <- pred.faked.trans[,a] * unlist(x$fit.coeffs)[a]

 y.lim <- c(0, max(pred.faked.tranWeight))  # max value for y-axis      

x11(width = 12, height = 15); par(mfrow=plot.layout)

for (j in 1:length(x$fit.preds)) {    # plot       
    plot(pred.faked[,j], pred.faked.tranWeight[,j] ,
    type="l", col="blue", lwd=2, xlim=c(0,1),  bty="n",
    xlab=paste(names(pred.faked)[j]),  ylim = y.lim,
    ylab=paste( "Transformed ", names(pred.faked)[j]))
    legend("topleft",paste(allFitPrds[j,2]), bty="n", text.col="red")
    }

          if(keep) {
          sol$pred.dis <- rewrap(pred.dis)      # predictor distance matrix
          sol$resp.dis <- rewrap(resp.vec)           # response as a distance matri
          sol$mod.resp.dis <- rewrap(mod.resp.vec)  # modeled response as a distance matrix
          sol$illus.predraw <- pred.faked
          sol$illus.predtranwt <- pred.faked.tranWeight
          sol
          }
             } # end of dismod.plot function
             
#################################################################################################################################################################
# dismod.predict and dismod.transform functions 
###############################################################################################################################################################

dismod.predict  <- function(dismod.object, domain){
   # This function predicts response distances using a disimilarity model and a set of predictor values for new sites
  sol <- c(call = match.call())
  class(sol) <- "dismod"
newpreds <- domain
 mod.env.dist <- dismod.transform(dismod.object, newpreds, distance=TRUE)
 resp.dist.pred <- dismod.object$intercept +   mod.env.dist
 sol$mod.response <- rewrap(resp.dist.pred)  # modeled response as distance matrix
 sol$pred.dis <- rewrap(mod.env.dist)         # predictors as distance matrix
 sol
 } # end function

#########################################################################################################

DougSplit <- function(myStr, ...) {  # string split function 
        strsplit(myStr, ...)[[1]] }
        
        
dismod.transform  <- function(dismod.object, domain, distance=FALSE)   {
# this function takes a disimilarity model and a set of predictor values for new sites and transforms
# this predictor space to that fitted by the model
# outputs either a dataframe of the converted variables or a distance matrix object

    newpreds <- domain
    allFitPrds.list <- lapply(dismod.object$fit.preds, DougSplit, split = "_")  # cannot use sapply because not all vars have transformation

allFitPrds <- data.frame(nrow=length(dismod.object$fit.preds), ncol=2) # Df with variable + transformation 
names(allFitPrds) <- c("Variable", "Transformation")
      for (i in 1:length(dismod.object$fit.preds)) {
       allFitPrds[i,1] <- allFitPrds.list[[i]][1]
       allFitPrds[i,2] <- allFitPrds.list[[i]][2]
      }


SelPreds <- subset(domain, select = allFitPrds[,1])  # select only the fitted variables
SlidePreds <- slide(SelPreds) # slide ALL predictors with min values <= ZERO to allow transformation
TranfPreds <- SlidePreds
 
 for (j in 1:length(dismod.object$fit.preds)) {
     if(!is.na(allFitPrds[j,2])) {
     TranfPreds[,j] <- var.transformer(v=SlidePreds[,j], tranf = allFitPrds[j,2])[,2]  # if the fitted variable was transformed, range std and transform
     } else {
     TranfPreds[,j] <-  as.vector(decostand(SlidePreds[,j], method="range"))  # if the fitted variable was not transformed, simply range std
     }
 }

 if (distance==TRUE) {
        pred.dis <-  as.data.frame(apply(TranfPreds, 2, function(x) as.vector(dist(x, method="manhattan")))) # convert to disssimilarities 
        resp.dis <- as.matrix(pred.dis) %*% unlist(dismod.object$fit.coeffs)  # apply model variables weigthed by the coefficeints
        dist.mat <- rewrap(resp.dis)
        dist.mat  # returns the distance matrix

      }  else  { # convert variables to range standardised and transformed versions, weighted by the model coefficients

     if (dismod.object$geo==T)  stop ("cannot produce transformed version of variables that include geographic distance - distance must be TRUE")
     
     else  {
      sel.preds.w <- t(t(TranfPreds)* dismod.object$fit.coeffs)  # transformed variables weigthed by the coefficeints
      as.data.frame(sel.preds.w) # returns the transformed variables weigthed by the coefficeints
           } #end else
           }
      } #end function

#################################################################################################################################################################
# dismod.class function  
###############################################################################################################################################################
     
dismod.class <- function(dismod.object, levs=seq(2,20,1), clusters="ward", domain)  {
# This function fits and predicts a classification from a disimilarity model
# using the training set to define classes
require(vegan)
require(cluster)

  sol <- c(call = match.call())
  class(sol) <- "dismod.class"

  train.cases <- nrow(dismod.object$pred)    # number of cases in the training set

  if(dismod.object$geo==TRUE) {
 if(clusters=="pam") {
 stop("cannot use PAM if dissimilarity model includes geographic distance", "\n")
   }
   }

pred.tran <- dismod.transform(dismod.object, dismod.object$pred, distance=TRUE)   # distance matrix transformed according to model
if(clusters!="pam") type.hclust <- hclust(pred.tran, method = clusters)  # hierachical clustering to define classes

train.types <- as.data.frame(matrix(nrow=train.cases, ncol=length(levs))) # DF to store classes at each level for training cases
type.medoids <-vector("list")  # to store medoids (or pred means) of each type at each level
for (k in 1:length(levs)) {
    if(clusters=="pam") {
    types <- pam(pred.tran, k = levs[k], metric="manhattan")     # make types at each level
    train.types[,k]<- types$clustering                           # assign cases to types
    type.medoids[[k]]<- types$medoids                             # keep medoids of types
} else {
train.types[,k] <- as.vector(cutree(type.hclust, k=levs[k]))  # assignment of training cases to types at each level

  type.means <- as.data.frame(matrix(NA, nrow=levs[k], ncol(dismod.object$pred)))
  for(m in 1:levs[k])  {
     type.means[m,] <- apply(dismod.object$pred[train.types[,k]==m,], 2,  mean)# mean values of preds for classes
     }
  type.medoids[[k]]<- type.means #  store mean values of preds for classes at each level
  names(type.medoids[[k]]) <- names(domain)
   }
 } # end make types and medoids
domain.types <- as.data.frame(matrix(nrow=nrow(domain), ncol=length(levs)))   # now assign new cases to types
row.names(domain.types) <- row.names(domain) 

 for (i in 1:nrow(domain)) { # for each case in domain
  for (k in 1:length(levs)) {  # for each level
  if(clusters=="pam") {
  centroids <- dismod.object$pred[type.medoids[[k]], ]  # the centroids of the types (defined from training set)
  pred.centroids <- rbind(domain[i,], centroids)   #  bind case i to centroids
  dist.mat <- as.matrix(dismod.transform(dismod.object, domain=pred.centroids , distance=TRUE)) # transform as for dis model
  dist.type <- dist.mat[2:nrow(dist.mat),1] # the first column are distances of case to centroids
  domain.types[i,k] <- match(min(dist.type), dist.type)  # closest type to case i in domain
 } else  {
  centroids <- type.medoids[[k]]  # the type.mean of the types (defined from training set)
  pred.centroids <- rbind(domain[i,], centroids)   #  bind case i to centroids
  dist.mat<- as.matrix(dismod.transform(dismod.object, domain=pred.centroids, distance=TRUE))
  dist.type <- dist.mat[2:nrow(dist.mat),1] # the first column are distances of case to centroids
  domain.types[i,k] <- names(dist.type)[match(min(dist.type), dist.type)]  # closest type to case i in domain
  }  # end else
  }  # end levels
  cat("processed case" , i, "\n")
    } # end of cases cycle

  names(domain.types) <- paste("Type_", as.character(levs), sep="")    # name the columns
  names(train.types) <- paste("Type_", as.character(levs), sep="")    # name the columns

sol$domain.class <- domain.types
sol$medoids <- type.medoids
sol$train.class <- train.types
sol
}   # end dismod.class function




dismod.fit.X <- function(resp, pred, method="bray", totP=10, toolong = 0.90, min.and.max=NULL, crit=0, tran=c("SQ", "R2", "L10", "LN"))   ##fits a dissimilarity model 
   {
  #this version is slower but uses less memory and can fit larger models 
  
  # the following functions are all used in the dismod.fit.X functions
     
             manhat <- function(x) {    # function computes RANGE standardised manhattan distances
         mn <-    unlist(apply(t(x),1,min)) #computes ALL min values
         mx <-    unlist(apply(t(x),1,max))#computes ALL max values
         diffmaxmin<- mx-mn # compute ALL differences
         diffmaxmin[diffmaxmin==0]<-1 #when diffmaxmin became 1 when max=min. avoid NaN.
         x.temp <- t(abs(mn - t(x))/diffmaxmin)
         dist(x.temp, method="manhattan")   #end manhat function
         }

         man <- function(x) {           #  function computes manhattan distances  AND RETURNS AS VECTOR
          as.vector(dist(x, method="manhattan"))   
         }                              #end manhat function


   rstand <-  function(x, minmax.var = NULL) { # function to range standardise individual predictors
      if (is.null(minmax.var)) {
                r.std <- abs(min(x) - x)/((max(x)-min(x)))
       }  else {
                mn <- min(minmax.var)
                mx <- max(minmax.var)
                r.std <- abs(mn - x)/(mx-mn)
                }
       r.std
       } # end rstand function


 slide <-  function(x, min.and.max=NULL)  {   # function slides predictors with min values <= ZERO to allow transformation
  x.slide <- as.data.frame(matrix(nrow=nrow(x), ncol=ncol(x)))
  M <- ncol(x)  # number of variables
  if (is.null(min.and.max)) {
    for (i in 1:M) {
      var.temp <- as.vector(x[i])
      x.slide[i] <- var.temp
      mn <- min(var.temp)
        if(mn <= 0) {   #
        var.mod <- var.temp + (abs(mn) + 1)
        x.slide[i] <- var.mod
        }
    }
  }
  else {
  for (i in 1:M) {
      var.temp <- as.vector(x[,i])
      x.slide[i] <- var.temp
      mn <- min(min.and.max[,i])
        if(mn <= 0) {
        var.mod <- var.temp + (abs(mn) + 1)
        x.slide[i] <- var.mod
        }
    }
  }
  names(x.slide) <- names(x)
  x.slide
  } # end slide function

  trans <-  function(x, min.and.max=NULL, tran=tran)  {   # function FIRST makes a chosen set of transformations of the predictors AND THEN range standardises
          x.tran <- as.data.frame(matrix(nrow=nrow(x), ncol=0))  # blank DF
          M <- ncol(x)  # number of variables
          if(min(x) <= 0) cat("cannot transform variables less than or equal to 0")

          if (is.null(min.and.max)) {
            for (i in 1:M) {
              var.tran <- matrix(nrow=nrow(x), ncol=0)
              y <- as.vector(x[i])
            
              var.tran <-  cbind(var.tran, rstand(y)) 
              if(any(tran%in%"CU")) var.tran <- cbind(var.tran, rstand(y^3))
              if(any(tran%in%"SQ")) var.tran <- cbind(var.tran, rstand(y^2))
              if(any(tran%in%"R2")) var.tran <- cbind(var.tran, rstand(y^0.5))
              if(any(tran%in%"R4")) var.tran <- cbind(var.tran, rstand(y^0.25)) 
              if(any(tran%in%"L10")) var.tran <- cbind(var.tran, rstand(log10(y)))
              if(any(tran%in%"LN")) var.tran <- cbind(var.tran, rstand(log(y)))   # six transformations
                
              labs <- vector("character") # lables  paste names+trans
              labs[1] <- names(x[i])
              if(any(tran%in%"CU")) labs <- c(labs, paste((names(x[i])) ,"CU", sep=""))
              if(any(tran%in%"SQ")) labs <- c(labs, paste((names(x[i])) ,"SQ", sep=""))
              if(any(tran%in%"R2")) labs <- c(labs, paste((names(x[i])) ,"RT2", sep=""))
              if(any(tran%in%"R4")) labs <- c(labs, paste((names(x[i])) ,"RT4", sep=""))
              if(any(tran%in%"L10")) labs <- c(labs, paste((names(x[i])) ,"L10", sep=""))
              if(any(tran%in%"LN")) labs <- c(labs, paste((names(x[i])) ,"LN", sep=""))
              var.tran <- as.data.frame(var.tran)
              names(var.tran) <- labs
              x.tran <- cbind.data.frame(x.tran, var.tran)
              }
              } else {
              min.and.max <- as.data.frame(min.and.max)
              for (i in 1:M) {
              var.tran <- matrix(nrow=nrow(x), ncol=0) # blank matrix
              y <- as.vector(x[i])
              var.tran <-  cbind(var.tran, rstand(y, minmax.var = slide(min.and.max[i])))
              if(any(tran%in%"CU")) var.tran <- cbind(var.tran, rstand((y^3), minmax.var= slide(min.and.max[i])^3))
              if(any(tran%in%"SQ")) var.tran <- cbind(var.tran, rstand((y^2), minmax.var= slide(min.and.max[i])^2))
              if(any(tran%in%"R2")) var.tran <- cbind(var.tran, rstand((y^0.5), minmax.var= slide(min.and.max[i])^0.5))
              if(any(tran%in%"R4")) var.tran <- cbind(var.tran, rstand((y^0.25), minmax.var= slide(min.and.max[i])^0.25))
              if(any(tran%in%"L10")) var.tran <-cbind(var.tran, rstand((log10(y)), minmax.var= log10(slide(min.and.max[i]))))
              if(any(tran%in%"LN")) var.tran <- cbind(var.tran, rstand((log(y)), minmax.var= log(slide(min.and.max[i]))))
              
              labs <- vector("character") # lables  paste names+trans
              labs[1] <- names(x[i]) 
              if(any(tran%in%"CU")) labs <- c(labs, paste((names(x[i])) ,"CU", sep=""))
              if(any(tran%in%"SQ")) labs <- c(labs, paste((names(x[i])) ,"SQ", sep=""))
              if(any(tran%in%"R2")) labs <- c(labs,paste((names(x[i])) ,"RT2", sep=""))
              if(any(tran%in%"R4")) labs <- c(labs, paste((names(x[i])) ,"RT4", sep=""))
              if(any(tran%in%"L10")) labs <- c(labs,paste((names(x[i])) ,"L10", sep=""))
              if(any(tran%in%"LN")) labs <- c(labs,paste((names(x[i])) ,"LN", sep=""))
              var.tran <- as.data.frame(var.tran)
              names(var.tran) <- labs
              x.tran <- cbind.data.frame(x.tran, var.tran)
              }
              }
              x.tran
              } # end trans function


 rewrap <- function(vec) {   # function returns a new distance matrix comprising a rewrapped distance vector
              ndis <- length(vec) # number of distances in vector
              np <- (1 + sqrt(1 + 8*ndis))/2  # np = No. of sites - number of distances = n(n-1)/2; solve for positive root of n
              # Convert vec, the vector of distances, to D.tri, an upper triangular matrix.
              D.tri <- matrix(nrow=np, ncol=np)
              i <- 1
              n1 <- 1
              n2 <- i*(np) - 1
              while(i <= np - 1)
              {
              	D.tri[i,(i + 1):np] <- vec[n1:n2]
              	n1 <- n2 + 1
              	n2 <- n2 + (np - i - 1)
              i <- (i + 1)
              }
              # Convert D.tri, the upper triangular matrix into D.sym, a symmetric matrix.
              D.sym <- matrix(nrow=np, ncol=np)
              i <- 1
              while(i <= np)
              {
              	D.sym[i,i] <- 0
              	j <- (i + 1)
              	while(j <= np)
              	{
              		D.sym[j,i] <- D.tri[i,j]
              		D.sym[i,j] <- D.tri[i,j]
              		j <- (j + 1)
              	}
              i <- (i + 1)
              }
              as.dist(D.sym)
              }   # end of rewrap  function 
      
         
            geodist <- function(coords, min.and.max=NULL, tran=tran)  {   #  compute geographical distances including transformed versions
              geo.d <- as.vector(dist(coords, method="euclidean"))
              #geo.tran <- matrix(nrow=Ndists, ncol=0)
                    if (is.null(min.and.max)) {
                    geo.tran <- rstand(geo.d)
                    if(any(tran%in%"CU"))  geo.tran <- cbind(geo.tran, rstand(geo.d^3))   # six transformations
                    if(any(tran%in%"SQ"))  geo.tran <- cbind(geo.tran, rstand(geo.d^2))
                    if(any(tran%in%"R2"))  geo.tran <- cbind(geo.tran, rstand(geo.d^0.5))
                    if(any(tran%in%"R4"))  geo.tran <- cbind(geo.tran, rstand(geo.d^0.25))
                    if(any(tran%in%"L10")) geo.tran <- cbind(geo.tran, rstand(log10(geo.d+1))) # note sites may have the same cordinates (e.g., due to lack of spatil resolution) thus log(x+1) to avoid NaN
                    if(any(tran%in%"LN"))  geo.tran <- cbind(geo.tran, rstand(log(geo.d+1)))
                    } else {
              max.dis <- as.vector(dist(min.and.max[, 1:2], method="euclidean"))  # note function needs "coord.minmax"
              min.dis <- min(geo.d)
              minmax.var=c(min.dis, max.dis)
              geo.tran <- rstand(geo.d, minmax.var)
              if(any(tran%in%"CU"))  geo.tran <- cbind(geo.tran, rstand(geo.d^3, (minmax.var)^3 ))    # six transformations
              if(any(tran%in%"SQ"))  geo.tran <- cbind(geo.tran, rstand(geo.d^2, (minmax.var)^2 ))
              if(any(tran%in%"R2"))  geo.tran <- cbind(geo.tran,rstand(geo.d^0.5, (minmax.var)^0.5))
              if(any(tran%in%"R4"))  geo.tran <- cbind(geo.tran, rstand(geo.d^0.25, (minmax.var)^0.25))
              if(any(tran%in%"L10")) geo.tran <- cbind(geo.tran, rstand(log10(geo.d+1), log10(minmax.var)))  # note sites may have the same cordinates (e.g., due to lack of spatil resolution) thus log(x+1) to avoid NaN
              if(any(tran%in%"LN"))  geo.tran <- cbind(geo.tran, rstand(log(geo.d+1), log(minmax.var)))
              }
              labs <- vector("character") # lables  paste names+trans
              labs[1] <- paste("geodist")
              if(any(tran%in%"CU"))  labs <- cbind(labs, paste("geodist" ,"CU", sep=""))
              if(any(tran%in%"SQ"))  labs <- cbind(labs, paste("geodist" ,"SQ", sep=""))
              if(any(tran%in%"R2"))  labs <- cbind(labs, paste("geodist" ,"R2", sep=""))
              if(any(tran%in%"R4"))  labs <- cbind(labs, paste("geodist" ,"R4", sep=""))
              if(any(tran%in%"L10")) labs <- cbind(labs, paste("geodist" ,"L10", sep=""))
              if(any(tran%in%"LN"))  labs <- cbind(labs, paste("geodist" ,"LN", sep=""))
              geo.dis <- as.data.frame(geo.tran)
              names(geo.dis) <- labs
              geo.dis
              }    #end geodist function   
                                      
makepreds.k <- function(pred, min.and.max=NULL, geo=F, tran=tran) {  # this function makes the complete set of  distance pairs to be used as predictors
          N <-nrow(pred)  # number of cases in the dataset
          minmaxmod <- data.frame(matrix(nrow=2, ncol=0))
                   
          pred.slid <- slide(pred, min.and.max)   # note predictors are "slid" to allow transformations of negative numbers
          pred.k <- trans(pred.slid, min.and.max=min.and.max, tran=tran)   # for each predictor compute all transformation options and range standardise (min=O, max=1)

          pred.k
            } # end makepreds function  



if(nrow(resp)!=nrow(pred)) stop("Wrong data dimensions")
 geo=F 
require(vegan)

  sol <- c(call = match.call())
  class(sol) <- "dismod"
  sol$pred <- pred
  sol$resp <- resp
              
       
          resp.nonex <-  vegdist(resp, method=method) #nonextended RESPONSE distances 
          cat("Extending distances.....", "\n")
          resp.ini.d <-  stepacross(resp.nonex, path = "extended", toolong) # compute extended RESPONSE distances with cutoff default = 0.9
          resp.ini.vec <- as.vector(resp.ini.d)  # response distances as vector
          
           
          pred.ini.d <- manhat(pred)  # this is the range std manhattan using all predictors 
          R.ini<-mantel(pred.ini.d, resp.ini.d, method="pearson", permutations=0) # do initial test with input predictors
          R.ini.stat <- R.ini$statistic        
          cat("Untuned R using all predictors", R.ini.stat, "\n")
          sol$untunedR  <- R.ini.stat      #save result 

       
          N <-nrow(pred)  # number of cases in the dataset
          Ndists <- (N*(N-1))/2     # no distances in the complete distance matrices
          check.var <- apply(pred, 2, function(x) sd(x))   # check variables have some variation and remove any variables with no variation
          pred.use <-  pred[, check.var > 0]
          M <- ncol(pred.use) # number of variables
 
pred.k <- makepreds.k(pred.use, min.and.max=min.and.max, geo=geo, tran=tran)    # make distance pairs of all variables plus transformations
                 
# START CV
   
allR <- data.frame(matrix(nrow = totP, ncol = 0))   #make blank df for results of each fold
testR <- vector(mode = "numeric", length = 10) # this stores the testing R values from each fold 
coeffs.cv <- vector("list", 10)

fold.vec <- sample(rep(c(1:10), length.out=N))  # ramdom slpit of N cases into 10 folds 

# open graphics device for CV output
         x11()
         par(mfrow=c(3,4))
    
    fold.fit.preds <- vector("list", length=10)  # save details of each fold model
    fold.fit.coeffs <- vector("list", length=10) # save details of each fold model
         
    for (f in 1:10)  {      # OUTER LOOP TO MAKE CV FOLDS
            
      # Select fold from pred.dis and resp.ini.d
      preds.f <- as.numeric(fold.vec != f) # new vector 1 if included in fold, 0 if in holdout set
      N.f <- sum(preds.f)  # number of cases in the fold dataset
      Ndists.f <- (N.f*(N.f-1))/2     # no distances in the fold distance matrix
      N.h <- N - N.f
      Ndists.ho <- (N.h*(N.h-1))/2     # no distances in the holdout distance matrix
                
      pick.vec <- as.vector(dist(preds.f, method="canberra"))  # any distance NOT equal to ZERO is to be EXCLUDED                 
      pick.vec <- replace(pick.vec, is.na(pick.vec), 10) # replace NAs to be replaced by 10 - these are the holdout distances
      
                 
      resp.vec <- resp.ini.vec[pick.vec == 0]  # initial response 
      resp.y <- resp.vec # intialize response vector
     
      
      if (sd(resp.vec) == 0) {         # safety move here could have response subset with no variabilty, hence distances will be all 0
      cat("this fold has no variabilty, skip to next subset", "\n")
      bestvarcoeff <- NULL
      testR <- rep(NA, times=10)
      } else     
      {   # ENTER MIDDLE loop to itteratively add best variable to model up to totP = total number of variables to add and TEST
                                                            
          bestvar <-vector(mode = "numeric", length = totP)       #store list of vars giving best model
          bestvarcoeff <- vector(mode = "numeric", length = totP)  # store best var coeffs
          testR <- vector(mode = "numeric", length = totP)       # store test results
          fitR <- vector(mode = "numeric", length = totP)       # store R for fitted model

                
                 for (j in 1:totP)  {      # MIDDLE LOOP TO ADD VARIABLES TO MODEL  j=NO.predictors to add         
                              
                        varcor <- vector("numeric", ncol(pred.k))
                        for (z in 1:ncol(pred.k)) {                                        
                        p.q <- man(pred.k[fold.vec != f, z ])
                        varcor[z] <- if (sd(p.q)==0) 0 else cor(p.q, resp.y)  # safety move in case of no variabilty in predictor fold                                                                                            
                       rm(p.q)
                        }
                        
                        varcor.pos <- as.numeric(varcor>0)  # vector to remove all NEGATIVE r 
                        modvar <- as.vector(varcor*varcor.pos)  # remove variables with -ve corellation                                                                      
                        bestvar[j] <- as.numeric(match(max(modvar), modvar))  #  best var for j-th itteration
                       
                        # fit a linear model here
                                                
                        if(j==1) {  # if more than one variable is fitted need to recalculate coeffs
                            var.m <- man(pred.k[fold.vec != f, bestvar[j]]) # best variable itteration 1                                              
                            mod <- lm(resp.y ~ var.m)                            
                            bestvarcoeff[j] <-  as.numeric(mod$coefficients[2])
                            resp.y <- mod$residuals   # residuals of best model to date                                                      
                            rm(var.m)
                            }
                        else {                      
                            bestvarcoeff <- vector(mode = "numeric", length = totP)  # reset the best var coeffs to ZERO this step calcs a new set                                                       
                            
                            var.m <-  as.data.frame(matrix(NA,  nrow = Ndists.f, ncol = 0)) 
                            for(b in 1:j) {
                            var.m <- cbind.data.frame(var.m, man(pred.k[fold.vec != f, bestvar[b]]))  # put selected variables AS DISTANCES into date in dataframe
                            names(var.m) <- as.character(1:b)
                            }
                            
                                                                                                                                                                                          
                            mod <- lm(resp.vec ~ .,var.m) # refit model                           
                            mod.coef <- mod$coefficients[2:(j+1)]  # these are the fitted coefs 
                            coef.pos <- mod.coef * as.numeric(mod.coef>0)  # if the coeffs are +ve then they are retained else = 0                                                    
                            coef.pos <- replace(coef.pos, is.na(coef.pos), 0) # IF THERE ARE pred.dis with COR = 1 function lm will produce coeffs = NA, replace NAs by 0 
                            bestvarcoeff[1:j] <- coef.pos[1:j]  # term may be eliminated if coeffs are less than 0                   
                            resp.y <- mod$residuals  # residuals of best model to date
                            rm(var.m)
                            } # end of coefficient refitting
                            
                                fitR[j] <- cor(mod$fitted.values, resp.vec)  #  best R (fitted variable) for j-th itteration
                                      
                                               
                                                             
                                      # test fitted model with holdout fold                                      
                                      test.preds <- as.matrix(pred.k[fold.vec == f, bestvar[1:j]])   # test fold - columns to select are stored in bestvar
                                      test.preds.w <- t(bestvarcoeff[1:j]* t(test.preds))  #  distances = addition of selected pred.dis*coeff 
                                      test.preds.d.vec <- man(test.preds.w)                                                                                                         
                                      test.resp.vec <- resp.ini.vec[pick.vec == 10] # take rows and colomns from extended distance matric to test this fold
                                      
                                      if (sd(test.resp.vec) == 0) {       # safety move here, small set of sites can have all zero distances
                                      R <- 0
                                      } else {                                        
                                      R <- as.numeric(cor(test.preds.d.vec, test.resp.vec, method="pearson"))
                                      }
                                      testR[j] <- R #save results 
                                      rm(test.preds.d.vec)
                                      
                                                                                                    
          }  #end middle loop   
          }  # end if               
                  
          coeffs.cv[[f]]<- bestvarcoeff                       
          allR<-cbind.data.frame(allR,testR)
          p<-c(1:totP)
          
          # plot results of CV fold
          fit.max <- max(fitR)
          test.max <- max(testR)
          y.max <- max(fit.max, test.max)
          y.min <- min(fit.max, test.max,0)
               
          plot(p, fitR, type="b", xlab ="Number of Predictors", ylab="Correlation", ylim = c(y.min, y.max), col="blue",
          main=(paste("CV fold", f)))
          lines(p, testR, type="b", col="red")          
          legend("bottomright",  c("fitted R", "test R"), col = c("blue", "red"),
            text.col = "black", lty = c(2, 2), pch = c(-1, -1),
            merge = TRUE, bg = 'gray90', inset = .05)        
      
       
    fold.fit.preds[[f]]  <- bestvar  # save model details of each fold model
    fold.fit.coeffs[[f]] <- bestvarcoeff 
                        
     } # end outer loop
     

   aveR <- vector(mode = "numeric", length = totP) # evaluate mean R over 10 folds for each of totP varaiables 
   se <-   vector(mode = "numeric", length = totP)  
   for (c in 1:totP)  {   # find mean and SD for each fold
              
              aveR[c] <- mean(as.numeric(allR[c,], na.rm = TRUE))     # calculate mean of the 10 CV trials for No. vars 1 to, strip NAs 
              se[c]  <- sd(as.numeric(allR[c,]), na.rm = TRUE)/sqrt(10)   # calculate standard error,  strip NAs
    }

    seU <- aveR+(se) # 
    seL <- aveR-(se)
    maxJ <- max(aveR) - crit*(max(aveR) - seL[match(max(aveR), aveR)]) # calculates the minimum value of the justifiable mean R 
    predsJ <- match(1, as.numeric(aveR >= maxJ))  # finds index of first item in 'mean' greater than maxJ (i.e. the justifiable no. of preds)        
    
   x11()           # plot final CV results
   upper <- max(seU)
   lower <-min(seL)
   plot(p, aveR, type="l",  xlab ="Number of Predictors", ylab="Correlation", main="Final CV results", ylim=c(lower,upper), col="red", lwd = 2) #overplot the se lines?
   arrows(p,seU,p,seL,length=.05,angle=90,code=3, col="gray70")
    abline(h=maxJ, v=predsJ, col = "black")
   legend("bottomright",  c("Mean R", "SE","Justified No. of predictors"), col = c("red", "gray70", "black"),
           text.col = "black", lty = c(2, 2), pch = c(-1, -1),
           merge = TRUE, bg = 'gray90', inset = .05)
           
    sol$R.mod <- aveR[predsJ]            
   
   # fit final model                                                                                                                                                                          
                                                                                                                                                                                                                                                                                                                                                                                            
    cat("fitting final model with", predsJ, "predictors", "\n")                                                                                                                               
                                                                                                                                                                                              
      bestvar.fin <-vector(mode = "numeric", length = predsJ)       #store list of vars for FINAL model                                                                                       
      bestvarcoeff.fin <- vector(mode = "numeric", length = predsJ)  # store  var coeffs for FINAL model                                                                                      
      finalR <- vector(mode = "numeric", length = predsJ)       # store test results for FINAL MODEL                                                                                          
      bestvar.names.fin <-vector(mode ="character", length=predsJ)                                                                                                                            
                                                                                                                                                                                                                                                                                                              
      resp.y <- resp.ini.vec  # initialize response distances                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
     # final loop replicates the CV step above and could be made into a function!                                                                               
            
            for (h in 1:predsJ)  {  # FINAL model OUTER LOOP TO ADD VARIABLES TO MODEL  predsJ = NO.predictors to add         

                        
                        varcor <- vector("numeric", ncol(pred.k))
                        for (z in 1:ncol(pred.k)) {                                        
                        p.q <- man(pred.k[, z])
                        varcor[z] <- if (sd(p.q)==0) 0 else cor(p.q, resp.y)  # safety move in case of no variabilty in predictor fold                                                                                            
                        }
                                                                                                                  
                        varcor.pos <- as.numeric(varcor>0)  # vector to remove all NEGATIVE r 
                        modvar.fin <- as.vector(varcor*varcor.pos)  # remove variables with -ve corellation                                                                      
                        bestvar.fin[h] <- as.numeric(match(max(modvar.fin), modvar.fin))  #  best var for h-th itteration
                        bestvar.names.fin[h] <- as.character(names(pred.k[match(max(modvar.fin), modvar.fin)]))                        
                          
                          # fit a linear model here
                                                
                        if(h==1) {  # if more than one variable is fitted need to recalculate coeffs
                            var.m <- man(pred.k[, bestvar.fin[h]]) # best variable itteration 1                                              
                            mod <- lm(resp.y ~ var.m)                            
                            bestvarcoeff.fin[h] <-  as.numeric(mod$coefficients[2])
                            resp.y <- mod$residuals   # residuals of best model to date                                                      
                            }
                        else {                      
                            bestvarcoeff.fin <- vector(mode = "numeric", length = predsJ)  # reset the best var coeffs to ZERO this step calcs a new set                                                       
                            
                            var.m <-  as.data.frame(matrix(NA,  nrow = Ndists, ncol = 0)) 
                            for(b in 1:h) {
                            var.m <- cbind.data.frame(var.m, man(pred.k[ , bestvar.fin[b]]))  # put selected variables AS DISTANCES into date in dataframe
                            names(var.m) <- as.character(1:b)
                            }
                            
                                                                                                                                                                                         
                            mod <- lm(resp.ini.vec ~ .,var.m) # refit model                           
                            mod.coef <- mod$coefficients[2:(h+1)]  # these are the fitted coefs 
                            coef.pos <- mod.coef * as.numeric(mod.coef>0)  # if the coeffs are +ve then they are retained else = 0                                                    
                            coef.pos <- replace(coef.pos, is.na(coef.pos), 0) # IF THERE ARE pred.dis with COR = 1 function lm will produce coeffs = NA, replace NAs by 0 
                            bestvarcoeff.fin[1:h] <- coef.pos[1:h]  # term may be eliminated if coeffs are less than 0                   
                            resp.y <- mod$residuals  # residuals of best model to date
                            } # end of coefficient refitting                                                                                                                                                                                                                                                                                                                                                                                                                                 
                            }  #end FINAL model OUTER loop
                             
                              
       cat("final model - selected variables", "\n")
      print(bestvar.names.fin)
      cat("coefficients for final variables", "\n")
      print(bestvarcoeff.fin)
      cat("Final Model R (CV)", aveR[predsJ], "\n")
      cat("Final Model R (Fitted)", cor(mod$fitted.values, resp.ini.vec), "\n") 
      
      sol$R.cv <- aveR[predsJ]
      sol$R.fit <- cor(mod$fitted.values, resp.ini.vec)
      sol$fit.preds <- bestvar.names.fin
      sol$fit.coeffs <- bestvarcoeff.fin
      sol$intercept <- mod$coefficients[1]
      sol$min.and.max <- min.and.max
      sol$geo <- geo
      sol$method <- method
      sol$toolong <- toolong
      sol$tran <- tran
      sol$allR <-allR
       sol$fold.fit.preds  <- fold.fit.preds # save model details of each fold model
       sol$fold.fit.coeffs <- fold.fit.coeffs       
      sol
}     #end of dismod.fit.X function

tree.class <- function(resp, pred, level=seq(5,100,10), method="canberra", toolong=0.9, clusters="pam", domain, divs=10, minsize = 2)  {
# fit a tree based classification  using a sample provided by pred and resp to a domain provided by newdata
if(nrow(resp)!=nrow(pred)) stop("Wrong data dimensions")

require(vegan)
require(stats)
require(tree)
require(cluster)

N.cases <- nrow(resp)

class.frame <- as.data.frame(matrix(nrow = nrow(domain), ncol = 0))

sol <- c(call = match.call())
class(sol) <- "treeclass"

      resp.d <-  vegdist(resp, method=method) #nonextended RESPONSE distances
      cat("Extending distances.....", "\n")
      resp.dx <-  stepacross(resp.d, path = "extended", toolong) # compute extended RESPONSE distances with cutoff default = 0.9

x=level[length(level)] #  number of classes at the highest level

if(clusters=="pam") {
g.x <- pam(resp.dx, k = x)
x.groups <- as.factor(g.x$clustering)
 } else {
g.x <- hclust(resp.dx, method = clusters)  # NB if groups have less than minsize they will not be modelled by tree
x.groups <- as.factor(cutree(g.x, k = x))
}
rm(resp.d)
rm(resp.dx)

   x.tree <- tree(x.groups ~., pred, control = tree.control(nobs=N.cases, minsize = minsize, mindev = 0))

pick <- rep(1:divs, length.out=nrow(domain))   #  subdivide data for predictions
class.frame <-  as.data.frame(matrix(NA, nrow(domain), length(level)))  # df for results (class prediction)
row.names(class.frame) <- row.names(domain)

   for (j in 1:length(level)) {
          x.pruned <- prune.tree(x.tree, k = NULL, best = level[j])    #
          cat("fitted class number", j, "\n")

        for (d in 1:divs) {
              class.frame[pick==d, j] <- predict.tree(x.pruned, newdata=domain[pick==d, ], type="where")
             cat("preds made for subset", d, "level", j, "\n")
              }
          }
# evaluate exact number of classes at each level of detail - note this differes from "levels" as tree is pruned
number.classes <- vector("numeric", length= ncol(class.frame))
 for(i in 1:ncol(class.frame)) {
number.classes[i] <- length(table(as.vector(class.frame[,i])))
}
          names(class.frame) <- paste("Lev_", as.character(number.classes), sep="")    # name the columns

    sol$domain.class <- class.frame   # keep these for else where
    sol$tree <- x.tree # keep the tree
    sol$number.classes <- number.classes # vector with number of classes at each level of classification detail
    sol

 } # end tree.class
##############################################################
# function: forest.class 
##############################################################
 
forest.class <- function(resp, pred, method="bray", level=seq(5,100,10),  toolong=NULL, clusters="pam", domain, divs=10, keep=F, ntree=500)  {
# fit random forest based partitions using a sample provided by pred and resp to a domain

if(nrow(resp)!=nrow(pred)) stop("Wrong data dimensions")
require(vegan)
require(stats)
require(randomForest)
require(cluster)
N.cases <- nrow(resp)
sol <- c(call = match.call())
class(sol) <- "forestclass"

      resp.d <-  vegdist(resp, method=method) #nonextended RESPONSE distances
      if(!is.null(toolong)) {
      cat("........................................................", "\n")
      cat("Extending distances.....", "\n")
      resp.dx <-  stepacross(resp.d, path = "extended", toolong) # compute extended RESPONSE distances with cutoff default = 0.9
      } else {
      resp.dx <- resp.d
      }
      rm(resp.d) # tidy up

if(keep) {
forest.mods <- vector("list", length=length(level)) #  keep each forest model
}
class.assign <-matrix(nrow=nrow(resp), ncol=length(level))  # keep class assignmments

if(clusters!="pam" & clusters != "Kmean") {
resp.clus <- hclust(resp.dx, method = clusters)   #
} #


pick <- rep(1:divs, length.out=nrow(domain))   #
class.frame <-  as.data.frame(matrix(NA, nrow=nrow(domain), ncol=length(level)))
row.names(class.frame) <- row.names(domain)

for(j in 1:length(level)) {

 if(clusters=="pam") {
resp.grp.l <- as.factor(pam(resp.dx, k=level[j])$clustering)
 }  else {
    if (clusters=="Kmean") {
    resp.grp.l <- as.factor(kmeans(resp, centers=level[j])$cluster)
    } else {
    resp.grp.l <- as.factor(cutree(resp.clus, k = level[j]))
    }
    }
class.assign[,j] <- resp.grp.l # store the site to class assignments

if(keep) {
forest <- randomForest(x=pred, y=resp.grp.l, importance=TRUE, keep.forest=TRUE, ntree=ntree)
forest.mods[[j]] <- assign(paste("forest.", j, sep=""), forest)   # store this model

} else {
forest <- randomForest(x=pred, y=resp.grp.l, ntree=ntree)
}

         for (d in 1:divs) {
                  class.frame[pick==d, j] <- predict(forest, newdata=domain[pick==d, ])  # add prediction to frame
                  cat("predictions made for subset", d, "level", j, "\n")
                  }
  }

 names(class.frame) <- paste("Lev_", as.character(level), sep="")

    sol$domain.class <- class.frame
    sol$train.class <- class.assign
    if (keep) {
    sol$forest.mods <- forest.mods # keep the trees
    sol$train.class <- class.assign
    }
    sol

 } # end forest.class

##############################################################
# function: map.class 
##############################################################
 map.class <- function(classes, level, coords, res=20000, pos=NULL, RapidPlot = F, ...)  {

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
              
par(mar= c(1, 1, 2, 1) +  0.1)
plot(coords.s, col="white", xaxt="n", yaxt="n", ...)  # this is faster but may obscure classes by overplotting 
if(RapidPlot == T) {
    for(i in 1:n.class) {
     points(coords.s[tclass.s==OrdClass[i], ], col=col.c[match(OrdClass[i], class.names)], pch=pch.c[i], cex=0.6)
    }
} else {     
     for(i in 1:nrow(coords.s)) {
       points(coords.s[i, ], col=col.c[match(tclass.s[i], class.names)], pch=pch.c[match(tclass.s[i], class.names)], cex=0.6)
     }
}
     
  if(!is.null(pos)) {
     legend(pos, legend=paste(class.names, " (", PropClass, ")", sep=""), col = col.c, title="Class (% domain)",
     text.col = "black",  pch = pch.c, bg = 'gray90',  inset = .01) 
     }

 } # end map classification function
 

matcor.fit <- function(resp, pred, method="bray", toolong = 0.90, depth=5)  { ## fits a matrix correlation model (methodof Dodkins et al. 2005)
# This function fits a matrix correlation model
if(nrow(resp)!=nrow(pred)) stop("Wrong data dimensions")

require(vegan)

  sol <- c(call = match.call())
  class(sol) <- "matcor.mod"

if(depth >= ncol(pred)) depth <- ncol(pred)

resp.nonex <-  vegdist(resp, method=method) #nonextended RESPONSE distances
cat("Extending distances.....", "\n")
resp.ex <-  stepacross(resp.nonex, path = "extended", toolong) # compute extended RESPONSE distances with cutoff default = 0.9
resp.dis <- as.vector(resp.ex)  # response distances as vector

npred<-ncol(pred)   # number of predictors

# make vars
 allpred <- make.matcor.preds(pred)

  pred.pos <- array(data = c(1:(npred*4)), dim = c(ncol(pred), 4), dimnames = list(names(pred), c(2:5))) # index matrix for vars and subs

 name.best <- vector("character",depth)
 subs.best <- vector("character",depth)
 vars.used <- vector("numeric", depth) # vector of  vars included in order of columns of pred names
 rho <- vector("numeric", depth) # vector of  vars included in order of columns of pred names
 keep.vars<-as.data.frame(matrix(NA, nrow=nrow(pred), ncol=depth))

 for (v in 1:depth) {

 if(v==1) {
   res.cor<- vector("numeric") # vector of correlation results
     for (f in 1:4) {
     for (i in 1:npred) {    # this finds best single variable PLUS number of divs
     pred.dis<- as.vector(dist(allpred[,i,f], method="manhattan"))
     if (sd(pred.dis)==0) res.cor<-c(res.cor,0) else         # safety move in case no variation in predictor
     res.cor <- c(res.cor, cor(x=pred.dis, y = resp.dis, method = c("spearman")) )
     }
     }
     best <- match(max(res.cor), res.cor)

what.var <-   which(pred.pos==best,arr.ind=T)  # locate variable in the position matrix
vars.used[v] <- what.var[1]
name.best[v] <- names(pred[what.var[1]])   # name of best variable
subs.best[v] <- (what.var[2])+1            # number of subdivisons of best variable
keep.vars[,v] <- as.vector(allpred[,what.var[1], what.var[2]]) # keep the used variables
rho[v] <-max(res.cor)

allpred.m <- allpred[,-(vars.used), ]     # modified variables array without selected variables
no.pred.kept <-dim(allpred.m)[2]
pred.pos.m <- array(data = c(1:(no.pred.kept*4)), dim = c(no.pred.kept, 4), dimnames = list(names(pred[-(vars.used)]), c(2:5))) # index matrix for vars and subs

     } else {
     res.cor<- vector("numeric") # reset vector of correlation results

     for (f in 1:4) {
     for (i in 1:ncol(allpred.m)) {    # this finds best single variable PLUS number of divs
     vars.inc <- keep.vars[,1:(v-1)]
     
       if(no.pred.kept==1) {  # this is a safety move in case only 1 variable remains
       pred.dis<- as.vector(dist(cbind.data.frame(vars.inc, allpred.m[,i]), method="manhattan"))
        } else {
     pred.dis<- as.vector(dist(cbind.data.frame(vars.inc, allpred.m[,i,f]), method="manhattan"))
     }
     res.cor <- c(res.cor, cor(x=pred.dis, y = resp.dis, method = c("spearman")) )
     }
     }
     best <- match(max(res.cor), res.cor)

what.var <-   which(pred.pos.m==best,arr.ind=T)  # locate variable in the position matrix
name.best[v] <- row.names(what.var)   # name of best variable
subs.best[v] <- (what.var[2])+1            # number of subdivisions of best variable
vars.used[v] <-  match(name.best[v], names(pred))
rho[v] <-max(res.cor)

allpred.m <- allpred[,-(vars.used), ]     # modified variables array without selected variables
no.pred.kept<- npred-sum(as.numeric(vars.used>0))    # kept=remaining!
pred.pos.m <- array(data = c(1:(no.pred.kept*4)), dim = c(no.pred.kept, 4), dimnames = list(names(pred[-(vars.used)]), c(2:5))) # index matrix for vars and subs
}

cat("Tested variable", v, "selected", name.best[v], "with Rho =", rho[v], "\n")   # out put progress information

     } # end fitting v loop

max.depth.just <- match(max(rho), rho)   # max value of Rho reached
# save the subdisivion thresholds
variable.sub<-vector("list", length=max.depth.just)

for (x in 1:max.depth.just)  {
variable<-pred[ ,vars.used[x]]
variable.sub[[x]] <- as.vector(quantile(variable, probs=c(seq(from=1/as.numeric(subs.best[x]), to=0.999, by = 1/as.numeric(subs.best[x])))))
}

sol$rho <- rho                          # results raw, Rho
sol$vars <-  name.best             # the names of all the predictors fitted - some may not be justified by increasing Rho
sol$var.subs  <- subs.best  # the number of subdivisons of predictors fitted - some may not be justified by increa

sol$fit.preds <-  name.best[1:max.depth.just]        # the fitted predictors justified by increasing rho value
sol$fit.pred.Nsubs <- subs.best[1:max.depth.just]   #  number of subdivisons of fitted predictors
sol$fit.pred.pos <- vars.used[1:max.depth.just]     # colomn position in pred matrix of the fitted predictors
sol$fit.pred.subs <- variable.sub                   # list containing the boundaries of the subdivisions of each fitted predictor
sol
   } # end fit function
   
  # inner function matcor.mod
 make.matcor.preds <- function(pred)  { # makes the dummy variables # matrix correlation model
    # prepare a set of dummy variables for the variables in 2, 3  4 and 5 categories
    divset.2 <- as.data.frame(matrix(NA, nrow=nrow(pred), ncol=ncol(pred)))
    names(divset.2)<-paste(names(pred))#,2,sep="")
   for (i in 1:ncol(pred)) {
   med <- median(pred[,i])
   divset.2[,i] <- as.numeric(pred[,i]>med)
   }
   divset.3 <- as.data.frame(matrix(NA, nrow=nrow(pred), ncol=ncol(pred)))
   names(divset.3)<-paste(names(pred))#, "3", sep="")
   for (i in 1:ncol(pred)) {
   thirds <- as.vector(quantile(pred[,i], probs=c(0.333, 0.666)))
   divset.3[pred[,i]<=thirds[1], i] <- 1
   divset.3[(pred[,i] > thirds[1]) & (pred[,i]< thirds[2]), i] <- 2
   divset.3[(pred[,i] >= thirds[2]), i] <- 3
   }
   divset.4 <- as.data.frame(matrix(NA, nrow=nrow(pred), ncol=ncol(pred)))
   names(divset.4)<-paste(names(pred))#, "4", sep="")
   for (i in 1:ncol(pred)) {
   forths <- as.vector(quantile(pred[,i], probs=c(0.25, 0.5, 0.75)))
   divset.4[pred[,i]<=forths[1], i] <- 1
   divset.4[(pred[,i] > forths[1]) & (pred[,i]< forths[2]), i] <- 2
   divset.4[(pred[,i] >= forths[2]) & (pred[,i]< forths[3]), i] <- 3
   divset.4[(pred[,i] >= forths[3]), i] <- 4
   }
   divset.5 <- as.data.frame(matrix(NA, nrow=nrow(pred), ncol=ncol(pred)))
   names(divset.5)<-paste(names(pred))#, "5", sep="")
   for (i in 1:ncol(pred)) {
   fifths <- as.vector(quantile(pred[,i], probs=c(0.2, 0.4, 0.6, 0.8)))
   divset.5[pred[,i]<=fifths[1], i] <- 1
   divset.5[(pred[,i] > fifths[1]) & (pred[,i]< fifths[2]), i] <- 2
   divset.5[(pred[,i] >= fifths[2]) & (pred[,i]< fifths[3]), i] <- 3
   divset.5[(pred[,i] >= fifths[3]) & (pred[,i]< fifths[4]), i] <- 4
   divset.5[(pred[,i] >= fifths[4]), i] <- 5
   }
   allpred <- array(data = NA, dim = c(nrow(pred), ncol(pred), 4), dimnames = list(row.names(pred), names(pred), c(2:5)))
   for (f in 1:4)   allpred[,,f]<- as.matrix(get(paste("divset",(f+1),sep=".")))    
   return(allpred)
   }  # end make.matcor.preds function
   
   
matcor.predict <- function(matcor.object, domain)  { ## predicts to new cases from a matrix correaltion model  
 
  sol <- c(call = match.call())
  class(sol) <- "matcor.pred"
 
 # domain must contain the same variables (with variable names) as the model pred data but order not important.
 justified <- length(matcor.object$fit.preds) # the number of predictors fitted in final model 
# inner functions  
 div.2<-function(x, var.sub) {
    var.2 <- x
   var.2[x <= var.sub[1]] <- 1
   var.2[x > var.sub[1]] <- 2
   return(var.2)
   }
 div.3<-function(x, var.sub) {
 var.3 <- x
   var.3[x<=var.sub[1]] <- 1
   var.3[(x > var.sub[1]) & (x< var.sub[2])] <- 2
   var.3[(x >= var.sub[2])] <- 3
   return(var.3)
 } 
 div.4<-function(x, var.sub) {
 var.4 <- x
   var.4[x <= var.sub[1]] <- 1
   var.4[(x > var.sub[1]) & (x < var.sub[2])] <- 2
   var.4[(x >= var.sub[2]) & (x < var.sub[3])] <- 3
   var.4[(x >= var.sub[3])] <- 4
   return(var.4)
 }
div.4<-function(x, var.sub) {
var.4 <- x
   var.4[x <= var.sub[1]] <- 1
   var.4[(x > var.sub[1]) & (x < var.sub[2])] <- 2
   var.4[(x >= var.sub[2]) & (x < var.sub[3])] <- 3
   var.4[(x >= var.sub[3])] <- 4
   return(var.4)
 }
div.5<-function(x, var.sub) {
var.5 <- x
   var.5[x <= var.sub[1]] <- 1
   var.5[(x > var.sub[1]) & (x < var.sub[2])] <- 2
   var.5[(x >= var.sub[2]) & (x < var.sub[3])] <- 3
   var.5[(x >= var.sub[3]) & (x < var.sub[4])] <- 4
   var.5[(x >= var.sub[4])] <- 5
   return(var.5)
 }
# end inner functions
 
 types<-as.data.frame(matrix(NA,nrow=nrow(domain),ncol=justified))
 name.vars<-vector("character")
 name.levs<-vector("character")

 
      for (i in 1:justified) { 
      func <- get(paste("div",as.numeric(matcor.object$fit.pred.Nsubs[i]) , sep="." ))
      sel.pred <- as.vector(subset(domain, select = matcor.object$fit.preds[i]))  # select the ith fitted variables     
      types[,i] <- func(sel.pred, matcor.object$fit.pred.subs[[i]]) 
      name.vars <- c(name.vars, paste("var", i, matcor.object$fit.preds[i], sep="."))
      }
      names(types) <- name.vars 
      
      
      types.levs <-as.data.frame(matrix(NA,nrow=nrow(domain),ncol=0))
      for (k in 1:justified) {       # concatenate the types
      if (k==1) conc <- paste(types[ ,1], sep="")
      if (k==2) conc <- paste(types[ ,1], types[,2], sep="")
      if (k==3) conc <- paste(types[ ,1], types[,2], types[,3], sep="")
      if (k==4) conc <- paste(types[ ,1], types[,2],types[,3], types[,4],  sep="")
      if (k==5) conc <- paste(types[ ,1], types[,2],types[,3], types[,4], types[,5  ],  sep="") 
      types.levs<-cbind.data.frame(types.levs, conc)
      name.levs<-c(name.levs, paste("lev", k, sep="."))
      } 
      names(types.levs) <- name.levs      
         
sol$domain.cat <- types
sol$domain.class <- types.levs
sol$levs <- apply(types.levs, 2, function(x) levels(as.factor(x)))     # the levels
sol$NumType<- as.vector(apply(types.levs, 2, function(x) length(levels(as.factor(x))))) # the number of types in each level
sol
}    # end matcor predict function      
      
 lda.class <- function(resp, pred, method="canberra", level=seq(5,100,10),  toolong=NULL, clusters="pam", domain, divs=10, keep=F)  {
# fit linear discriminant function to traing data and predict for domain

if(nrow(resp)!=nrow(pred)) stop("Wrong data dimensions")

require(vegan)
require('MASS')
require(cluster)

N.cases <- nrow(resp)


sol <- c(call = match.call())
class(sol) <- "lda.class"

   resp.d <-  vegdist(resp, method=method) #nonextended RESPONSE distances
      if(!is.null(toolong)) {
      cat("........................................................", "\n")      
      cat("Extending distances.....", "\n")
      resp.dx <-  stepacross(resp.d, path = "extended", toolong) # compute extended RESPONSE distances with cutoff default = 0.9        
      } else {  
      resp.dx <- resp.d 
      }
      rm(resp.d) # tidy up

if(keep) {
lda.mods <- vector("list", length=length(level)) #  keep each lda model
class.assign <-matrix(nrow=nrow(resp), ncol=length(level))  # keep class assignmments
}

if(clusters!="pam") {
resp.clus <- hclust(resp.dx, method = clusters)   #
} #

pick <- rep(1:divs, length.out=nrow(domain))   #
class.frame <-  as.data.frame(matrix(NA, nrow=nrow(domain), ncol=length(level)))
row.names(class.frame) <- row.names(domain) 

for(j in 1:length(level)) {

 if(clusters=="pam") {
resp.grp.l <- as.factor(pam(resp.dx, k=level[j])$clustering)
 } else {
resp.grp.l <- as.factor(cutree(resp.clus, k = level[j]))
}

if(keep) {
lda.mod <- lda(x=pred, grouping=resp.grp.l)
lda.mods[[j]] <- assign(paste("lda.mod.", j, sep=""), lda.mod)   # store this model
class.assign[,j] <- resp.grp.l # store the site to class assignments
} else {
lda.mod <- lda(x=pred, grouping=resp.grp.l)
}

         for (d in 1:divs) {
                  class.frame[pick==d, j] <- predict(lda.mod, newdata=domain[pick==d, ])$class  # add prediction to frame
                  cat("preds made for subset", d, "level", j, "\n")
                  }
  }

 names(class.frame) <- paste("Lev_", as.character(level), sep="")

    sol$domain.class <-class.frame    # keep these for else where

    if (keep) {
    sol$lda.mods <- lda.mods # keep the trees
    sol$train.class <- class.assign
    }
    sol

 } # end lda.class


DescribeClass <- function(domain, classes=NULL, level = 1, coords = coords, RGB = c(1,2,3), sense = c(1,-1,1),  choices = c(1,2), pos="topleft", res = 30000) {
 # domain = the predictor data used to define the classification 
  # classes = matrix of sites by classes with columns representing diffrent levels of detail or classifications  
 # level = classification level to be plotted
 # coords = coordinates for sites if distance clusters and geographic regions are to be evaluated 
 # sense = to be used to control dirrections of colours                        
  # res = restricts number of points used in analysis and plotting (used when number of rows in classes > res 
   # choices PCA axes to plot
   
 nvars <- ncol(domain)  
   if(nrow(coords)>res) {
              pick <- sample(1:nrow(coords), size=res, replace = F)
              } else {
              pick<-c(1:nrow(coords))
              } 
   domain.s <- domain[pick,]            
   coords.s <- coords[pick,] 
               
if(!is.null(classes)) {     
 ob.name <- deparse(substitute(classes))   #name
 classes.s <- classes[pick,]              
 tclass <-  classes.s[,level]
 n.class <- length(which(as.vector((table(tclass) >0))))
 class.labs <- as.numeric(names(table(tclass) >0))
 class.names <-  names(table(tclass) >0)
 } else {
 ob.name <- "Classification"
 }

 # pca
   domain.pca <- princomp(domain.s, cor = FALSE, scores = TRUE)   # covariances NOT correlations
   rgb.control <- decostand(domain.pca$scores[,1:3], method="range", MARGIN=2) # range standradise the scores on axis 1:3
      
# extract loadings 
VarLoads <- domain.pca$loadings[names(domain.s), paste("Comp.", 1:nvars, sep="")]
# df with highest wiegthed variables on first three PCA components    
HiVars <-  as.data.frame(apply(VarLoads[,1:3], 2, function(x) row.names(VarLoads)[match(sort(abs(x), decreasing = T), abs(x))]))
 # get the sign of the correlations
HiVarSign <-  as.data.frame(apply(VarLoads[,1:3], 2, function(x) sign(x)[match(abs(x), sort(abs(x), decreasing = T))]))   
 HiVarSignC <- replace(HiVarSign, HiVarSign==1, "+")
 HiVarSignC <- replace(HiVarSignC, HiVarSign==-1, "-")
 
HiVarsSigned  <- as.data.frame(matrix(ncol=3, nrow=nvars))
row.names(HiVarsSigned) <-  paste("Variable ", 1:nvars, sep="") 
names(HiVarsSigned) <- paste("PCA Comp ", 1:3, sep="")
 for(j in 1:3)  HiVarsSigned[, j] <- paste(HiVarSignC[,j], HiVars[,j], sep="")
 
cat("Variable loadings on first 3 PCA axis", "\n")
cat("..........................................................................." , "\n")
print(VarLoads[,1:3]) 
cat("..........................................................................." , "\n")
cat("Variables ranked by loadings on first 3 PCA axes", "\n")
cat("..........................................................................." , "\n")
print(HiVarsSigned)  # print to screen
cat("..........................................................................." , "\n")

rgb.col <-c("RED", "GREEN", "BLUE")
col.ord <- rgb.col[RGB]
cat("..........................................................................." , "\n")
cat(ifelse(sense[1]<0, "Low", "High"), " values of ", as.character(HiVars[1,1]), "are ", col.ord[1], "\n")
cat(ifelse(sense[2]<0, "Low", "High"), " values of ", as.character(HiVars[1,2]), "are ", col.ord[2], "\n")
cat(ifelse(sense[3]<0, "Low", "High"), " values of ", as.character(HiVars[1,3]), "are ", col.ord[3], "\n")

rgb.sense <- rgb.control  # reverses the sense of the variable if negative
for (j in 1:3)  if(sense[j] < 0) rgb.sense[,j] <- abs(rgb.control[,j] - max(rgb.control[,j]))

point.col <- rgb(red=rgb.sense[,RGB], maxColorValue = 1) 

if(is.null(classes)) { class.pch <- rep(15, length(pick)) }
if(!is.null(classes)){ 
plot.pch <- 1:n.class 
class.pch <- vector("numeric", length(tclass))
for(i in 1:n.class) class.pch[tclass==class.labs[i]] <- plot.pch[i]
}
                                       
# plotting map
windows(width = 20, height = 20) # ; par(mfrow=c(2,2))
    
  plot(coords.s, col="white",xlab = "east", ylab = "north", main=paste("Map of PCA", ob.name))    # map of RGB
 points(coords.s, col=point.col, pch=class.pch)

 if(!is.null(classes)) {
 legend(pos, legend=class.names, title="Classes",
 text.col = "black",  pch = class.pch , bg = 'gray90',
            inset = .01) }
                                                     
# plotting PCA
windows(width = 4, height = 4) # ; par(mfrow=c(2,2)) 

 plot(domain.pca$scores[, choices], col="white", main=paste("PCA ordination of ", ob.name))    # plot of RGB coloured PCA
  points(domain.pca$scores[, choices], col=point.col, pch=class.pch)
 
 if(!is.null(classes)) {
 legend(pos, legend=class.names, title="Classes", cex=0.5,
 text.col = "black",  pch = class.pch , bg = 'gray90',
            inset = .01) } 

 } # end DescribeClass function
 
MapVar <- function(x=domain, divs=10, Var=1, coords, res = 30000, reverse = F, pos="topleft", MyCex = 0.5, ...) {
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
   for (i in 1:divs) points(coords.s[cut.var==OrdVar[i],],col=HeatCol[match(OrdVar[i], NamesCutVar)], pch=15, cex=MyCex) 
   
legend(pos, legend=div.mean, title=paste("Divisions of ", ob.name, sep=""),
 text.col = "black",  pch = 15 , bg = 'gray90', col=HeatCol,
            inset = .01) 
            
} # end MapVar

######################################################################################################################
 
taxa.class.fit <- function(resp, pred, Occ=20, Type = "ab", ...) {
# Occ = percent occupancy of individual taxa to be modelled.
# Type = pa will produce P/A models

if(nrow(resp)!=nrow(pred)) stop("Wrong data dimensions")

PerOccupy <- apply(resp, 2, function(x) sum(x>0)/nrow(resp)*100)
TaxaToModel <- PerOccupy >= Occ
NMods <- sum(TaxaToModel)
TaxaMods <- colnames(resp)[TaxaToModel]
RespToUse <- resp[,TaxaToModel]


cat("........................................................", "\n")
cat("Fitting Random Forest models for ", NMods, "taxa", "\n")

forest.mods <- vector("list", NMods) #  keep each forest model
names(forest.mods) <- TaxaMods

if(Type=="pa") {    # convert to PA data
      for (i in 1:NMods) {
    cat("........................................................", "\n")
    cat("Fitting presence/absence model to taxa ",i, ":", TaxaMods[i], "\n")
        forest.mods[[i]] <- randomForest(x=pred, y=as.factor(as.numeric(RespToUse[,i]>0)), ntree=500, ...) # store this model
         }
} else {
    for (i in 1:NMods) {
    cat("........................................................", "\n")
    cat("Fitting abundance model to taxa: ",i, ":", TaxaMods[i], "\n")
        forest.mods[[i]] <- randomForest(x=pred, y= RespToUse[,i], ntree=500, ...) # store this model
         }
}

out <- list(TaxaModels = forest.mods, Type = Type)
return(out)

} # end taxa.class.fit

######################################################################################################################

taxa.class.predict <- function(x, domain, divs=10, ...) {
  
Mods <- x$TaxaModels # get RF models for each taxa from the fit object
Type <-  x$Type 
  
pick <- rep(1:divs, length.out=nrow(domain))   # for subsetting the predictions
TaxaPredFrame <-  as.data.frame(matrix(NA, nrow=nrow(domain), ncol=length(Mods)))

row.names(TaxaPredFrame) <-  row.names(domain) 
names(TaxaPredFrame) <- paste("ModTaxa", 1:length(Mods)) 

if(Type=="pa") {     
    for (j in 1:length(Mods)) {
        for (d in 1:divs) { 
        TaxaPredFrame[pick==d, j] <- predict(Mods[[j]] , newdata=domain[pick==d, ], type="prob",...)[,2]  # add prediction to frame
        cat("preds made for subset", d, "taxa", names(Mods)[j], "\n")
        }
        } 
        } else {
    for (j in 1:length(Mods)) {
        for (d in 1:divs) { 
        TaxaPredFrame[pick==d, j] <- predict(Mods[[j]] , newdata=domain[pick==d, ], ...)  # add prediction to frame
        cat("preds made for subset", d, "taxa", names(Mods)[j], "\n")
        }
        } 
        }            
                           
return(TaxaPredFrame)                       
}   # end 

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



    

