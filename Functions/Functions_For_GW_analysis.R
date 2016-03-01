GetDaily <- function(Daily) {  # this has the natural median as a column
  # make factors for time-series
  
  Daily$chrono <- chron(as.character(Daily$Date), format = c(dates = "d/m/y"), out.format = c(dates = "d-m-y"))
  Daily$Dateformat <- as.Date(Daily$chrono, "%d/%m/%y")
  Daily$Julian <- julian(Daily$Dateformat)
  
  DayMonthYear <- unlist(sapply(as.character(Daily$Date), strsplit, split = "-"))
  Daily$day <- as.numeric(DayMonthYear[seq(1, length(DayMonthYear), by = 3)])
  Daily$month <- factor(DayMonthYear[seq(2, length(DayMonthYear), by = 3)], levels = MonthString)
  Daily$year <- as.numeric(DayMonthYear[seq(3, length(DayMonthYear), by = 3)])
  Daily$YYYY <- Daily$year
  
  Daily$Date1Jan <- as.Date(paste(as.character(Daily$YYYY), "01", "01", sep = "-") )
  Daily$Julian1Jan <- julian(Daily$Date1Jan)
  Daily$Julian <- julian(Daily$Dateformat)
  Daily <- Daily[order(Daily$Julian), ]
  Daily$DayOfYear <- Daily$Julian - Daily$Julian1Jan + 1
  
  names(Daily)[which(names(Daily) == "Data")] <- "Daily"
  
  # Get day to day changes in flow
  Daily$Daily <- as.numeric(Daily$Daily)
  Daily$DailyChange <- c(diff(Daily$Daily), NA )
  Daily$NaturalMedian <- median(Daily$Daily, na.rm=T)  # put the natural median into the frame for later use
  
  return(Daily)
}

#########
# function to rbind all the elements in a list of data frames
# works even when all NA's
Doug.rbind.list <- function (dfs)
{
  # dfs <- list(...)
  if (length(dfs) == 0)
    return(list())
  all.names <- unique(unlist(lapply(dfs, names)))
  data.frame(do.call("rbind", lapply(dfs, function(df) df)))
}
#########

#########
DougSplit <- function(mytext, myPart = 1, ...) {
  myOut <- strsplit(mytext, ...)[[1]][myPart]
  return(myOut)
}
#########

# function to cbind all the elements in a list of data frames
# works even when all NA's
Doug.cbind.list <- function (dfs)
{
  data.frame(do.call("cbind", lapply(dfs, function(df) df)))
}
#########
#########
