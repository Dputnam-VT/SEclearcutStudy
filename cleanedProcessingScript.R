setwd("C:/R_workspace/Collection2_data")
library(DBEST)
library(zoo)
library(lubridate)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(Kendall)
library(mblm)

dataDate1 = "2022_09_23" # using for the latest version of the time series csv

# Bringing in vegetation index time series
inputDF = read.csv(paste("NBR_timeSeries_",dataDate1,".csv",sep = ""),header = TRUE)
TS_DF = inputDF[order(inputDF$X),]

numSamples = length(inputDF$X)

## creating the dates for the time series index
compositeMonth = 2 # first month used for annual composite

startYear = colnames(TS_DF)[2]
startYear = substr(startYear, 2,5)
startDate = paste(startYear,'-',compositeMonth,'-01',sep = '')
startDate = as.Date.character(startDate,format = '%Y-%m-%d')

endYear = colnames(TS_DF)[length(TS_DF)]
endYear = substr(endYear, 2,5)
endDate = paste(endYear,'-',compositeMonth,'-01',sep = '')
endDate = as.Date.character(endDate,format = '%Y-%m-%d')

dates = seq(startDate,endDate, by = 'years')

# Dataframe to store recovery metrics for each stand
rMetrics <- data.frame(matrix(nrow = length(rownames(TS_DF)), ncol = 14,byrow = TRUE))
colnames(rMetrics) = c('DistY','Recovery_min','Recovery_max','Recovery_mag',
                       'Disturbance_mag','Y2R_seg','Recovery_slope','NBRy3','NBRy5',
                       'NBRy7','NBRy10','pre_distMean','post_recovVal','Y2R_80'
)
rownames(rMetrics) = TS_DF$X

# Dataframe to store auxiliary data for plotting (start dates for segments, etc.)
auxSegData <- data.frame(matrix(nrow = length(rownames(TS_DF)), ncol = 6,byrow = TRUE))
colnames(auxSegData) = c('DistStartIndex','DistEndIndex','RecoveryStartIndex','RecoveryEndIndex','gainLossGap','fitRMSE')
rownames(auxSegData) = TS_DF$X

## A function to calculate recovery metrics, will be looped over the output
#   objects from DBEST, outputs to be passed to plotting function
#   also outputs index values for the segments I'm selecting to an auxiliary dataframe
calcRecoveryMetrics = function(DBEST_object, index) {
  
  changeSeg_VIdif = DBEST_object$Change
  
  # Need to find the greatest magnitude decline in the time series first
  distIndex = which.min(changeSeg_VIdif) # index value for breakpoint arrays for disturbance
  distSegEndDate = dates[DBEST_object$End[distIndex]] # date object of end of disturbance segment
  distSegStartDate = dates[DBEST_object$Start[distIndex]] # date object of start of disturbance segment
  distStartIndex = DBEST_object$Start[distIndex]
  distEndIndex = DBEST_object$End[distIndex]
  # index value for breakpoint arrays for recovery, should be an increasing segment directly after disturbance
  gainSegIndices = which(changeSeg_VIdif > 0)
  gainSegStartDates = dates[DBEST_object$Start[gainSegIndices]]
  
  ## need to find the gain segment directly after the disturbance
  # making a list of the possible dates for possible gain segments
  potentialGainSegStartDates = gainSegStartDates[gainSegStartDates >= distSegEndDate]
  
  # adding the disturbance end date to the list
  potentialGainSegStartDates = sort(potentialGainSegStartDates,decreasing = FALSE)
  gainIndex = gainSegIndices[which(gainSegStartDates == potentialGainSegStartDates[1])]
  gainSegStartDate = dates[DBEST_object$Start[gainIndex]]
  gainSegEndDate = dates[DBEST_object$End[gainIndex]]
  gainSegEndIndex = DBEST_object$End[gainIndex]
  gainSegStartIndex = DBEST_object$Start[gainIndex]
  
  ## Adding a statement to remove data from bad time series
  # i.e time series with no gain, no loss, or only one segment
  if ((length(changeSeg_VIdif) <= 1) |
      (length(changeSeg_VIdif[changeSeg_VIdif < 0]) == 0) |
      (length(changeSeg_VIdif[changeSeg_VIdif > 0]) == 0) |
      (length(gainIndex) == 0)
  ) {
    rMetricsArray = rep(NA,14)
    auxSegArray = rep(NA,6)
  } else {
    # Calculating recovery metrics and segmentation outputs
    ts_len = length(DBEST_object$Data)
    
    DistStartY = as.double(substr(as.character.Date(distSegStartDate), start = 0, stop = 4))
    
    recovMin = as.double(DBEST_object$Data[distEndIndex]) # Recovery min
    
    recovMax = as.double(DBEST_object$Data[gainSegEndIndex]) # Recovery max
    
    recoveryMag = as.double(recovMax - recovMin) # Recovery magnitude
    
    # Using Max() and Min() in these calculations to enable the 5 year calculation period without going out of range of the dataset
    preDistMean = as.double(mean(DBEST_object$Data[max(1,(distStartIndex-5)):distStartIndex])) # Pre-disturbance mean VI value
    PostRecovMean = as.double(quantile(DBEST_object$Data[gainSegEndIndex:ts_len],c(0.85))) # Post-recovery mean VI value
    
    distMag = as.double(preDistMean - DBEST_object$Data[distEndIndex]) # Disturbance magnitude
    
    Y2R = as.double(gainSegEndIndex - distEndIndex) # Years to recovery
    
    Y2R80 = min(which(output$Data[distEndIndex:ts_len] >= (0.8*preDistMean))) - 1
    Y2R80 = min(Y2R80,(ts_len-distEndIndex))
    
    recovSlope = as.double(recoveryMag / Y2R) # Recovery slope
    
    NBR_year3 = as.double(DBEST_object$Data[(min((distEndIndex+3),ts_len))]) # absolute NBR value at 3,5,7,10 years post disturbance
    NBR_year5 = as.double(DBEST_object$Data[(min((distEndIndex+5),ts_len))])
    NBR_year7 = as.double(DBEST_object$Data[(min((distEndIndex+7),ts_len))])
    NBR_year10 = as.double(DBEST_object$Data[(min((distEndIndex+10),ts_len))])
    
    gapInGainLoss = as.double(DBEST_object$Start[gainIndex] - DBEST_object$End[distIndex]) # Gap between the end of the disturbance segment and the start of the gain segment
    
    RMSE = sqrt(mean((DBEST_object$Fit - DBEST_object$Data)^2))
    
    ## Adding recovery metrics to an array to be returned from the function
    rMetricsArray = c(DistStartY,recovMin,recovMax,recoveryMag,distMag,Y2R,recovSlope,NBR_year3,NBR_year5,NBR_year7,NBR_year10)
    
    ## Adding auxiliary data to dataframe
    auxSegArray = c(DBEST_object$Start[distIndex],DBEST_object$End[distIndex],DBEST_object$Start[gainIndex],gainSegEndIndex,gapInGainLoss,RMSE)
    
    rMetricsArray = append(rMetricsArray,c(preDistMean,PostRecovMean,Y2R80)) # These have to be calculated seperately, since their indexing data
  }                                                                   # if a time series doesn't have a gain or loss it'll crash without this 
  
  return(list(rMetricsArray,auxSegArray))
  
}

# Storing the max and min values of the input vegetation index / band, to be used to scale DBEST
VImax = max(TS_DF[2:length(TS_DF[1,])],na.rm = TRUE)
VImin = min(TS_DF[2:length(TS_DF[1,])],na.rm = TRUE)
span = VImax - VImin

## Function to plot the DBEST output correctly
plot_DBEST = function(DBEST_object,index,metricDF,auxDF,minVIval) {
  # getting output device set up
  svg(filename = paste('C:/R_workspace/Collection2_data/images/plot_',index,'_TimeSeriesGraph.svg',sep = ''),
      width = 4,
      height = 2,
      pointsize = 4,
      onefile = FALSE
  )
  
  stableColor = "#3B999B"
    # plotting the generalization lines furthest back
  plot.default(dates,as.vector(DBEST_object$Fit),
               type = 'l',
               xlab = 'YEAR',
               ylab = 'NBR',
               main = paste("DBEST Result : Stand #",index,sep = ""),
               lwd = 3,
               col = stableColor,
               ylim = c(VImin,1),
               axes = FALSE,
  )
  
  # getting axes customized
  box()
  axis.Date(1,at = seq(from = (startDate+365), to = endDate, by = '5 years'),labels = TRUE)
  axis(2)
  
  distColor = '#EE5873'
    # adding lines for the disturbance segments
  for (i in c(1,3)) {
    segments(x0 = dates[auxDF[index+1,i]],
             y0 = DBEST_object$Fit[auxDF[index+1,i]],
             x1 = dates[auxDF[index+1,i+1]],
             y1 = DBEST_object$Fit[auxDF[index+1,i+1]],
             col = distColor,
             lty = 1,
             lwd = 3
    )
  }
  
  # adding the original data last
  # raw time-series from original dataframe
  rawValues = as.numeric(as.vector(TS_DF[(index+1),2:length(TS_DF)]))
  # adding to plot
  lines(dates,as.vector(rawValues),
        col = 'black',
        lwd = 1,
        lty = 2
  )
  
  # adding bar to indicate disturbance duration
  segments(x0 = dates[auxDF[index+1,2]],
           y0 = minVIval,
           x1 = dates[auxDF[index+1,4]],
           y1 = minVIval,
           col = 'red',
           lty = 1,
           lwd = 4.5
  )
  
  # adding text to the graph for the disturbance duration
  text(x = dates[round((mean(c(auxDF[index+1,2],auxDF[index+1,4]))),digits = 0)],
       y = (minVIval + abs(minVIval*0.1)),
       labels = paste(as.character(metricDF$Y2R_seg[index+1]),'Year Recovery'),
       adj = 0.5
  )
  # recovery min / max indicators on axis
  points(x = c(as.Date(as.yearmon(dates[1]) -0.95, frac = 1),
               as.Date(as.yearmon(dates[1]) -0.99, frac = 1)
  ),
  y = c(metricDF[index+1,2],
        metricDF[index+1,3]
  ),
  col = c('orange','green'),
  pch = '-',
  cex = 4
  )
  
  
  # adding a legend
  legend(x = dates[(length(dates)-7.5)],
         y = 0.1,
         legend = c('Original Time-series','DBEST Generalization','Disturbance Event','Recovery Min','Recovery Max'),
         col = c('black',stableColor,distColor,'orange','green'),
         lty = c(2,1,1,NA,NA),
         pch = c(NA,NA,NA,'-','-'),
         pt.cex = c(NA,NA,NA,4,4),
         lwd = c(1,2,2,NA,NA),
         bg = 'light grey'
  )
  
  # closing the graphical device
  dev.off()
}

# implimenting DBEST, plotting results, and storing recovery metrics

for (i in TS_DF$X) {
  
  ## creating the time series object
  NBRvals = as.numeric(as.vector(TS_DF[(i+1),2:length(TS_DF)]))
  
  ts1 = zoo(x = NBRvals,order.by = dates)
  ts1 = na.locf(ts1, na.rm = FALSE, fromLast = TRUE) # filling in end values
  ts1 = na.locf(ts1, na.rm = FALSE, fromLast = FALSE)
  
  generalized = DBEST(data = ts1,
                      data.type = 'non-cyclical',
                      algorithm = "generalization",
                      breakpoints.no = 4,
                      first.level.shift = (span*0.15),
                      second.level.shift = (span*0.20),
                      duration = 5,
                      distance.threshold = 'default',
                      alpha = 0.05,
                      plot = 'off'
  )
  
  output = DBEST(data = generalized$Fit,
                 data.type = 'non-cyclical',
                 algorithm = "change detection",
                 breakpoints.no = 2,
                 first.level.shift = (span*0.30),
                 second.level.shift = (span*0.35),
                 duration = 5,
                 distance.threshold = (span*0.10),
                 alpha = 0.05,
                 plot = 'off'
  )
  
  # Calculate recovery metrics and auxilary data
  listOfArrays = calcRecoveryMetrics(output,i)
  rMetrics[(i+1),] = listOfArrays[[1]]
  auxSegData[(i+1),] = listOfArrays[[2]]
  
  # Plotting
  plot_DBEST(output,i,rMetrics,auxSegData,VImin) ##### rerun when needed #####
  
  cat("Stand #",i, "complete","\n")
}

#### Post Segmentation evaluation of the time series data / metrics ####
# need to add the identifier as a column for tracking things during filtering
rMetrics$'UniqueID' = row.names(rMetrics)
rMetrics$'K_shift' = rMetrics$post_recovVal - rMetrics$pre_distMean
auxSegData$'UniqueID' = row.names(rMetrics)
auxSegData$'DistDuration' = auxSegData$DistEndIndex - auxSegData$DistStartIndex

# Adding the NBR regrowth metric
rMetrics$NBRregrowth = rMetrics$NBRy5 - rMetrics$Recovery_min

# Joining table with information about sample location
locationData = read.csv("AllsamplesLocationData.csv")
rMetrics = merge.data.frame(rMetrics,locationData, 
                            by = intersect(rMetrics$UniqueID,locationData$UniqueID)
                            )

#### Cleaning data based on practical relationships ####

# Number of TS which did not have qualifying segments for metric calculation
areNAindicies = which(is.na(rMetrics$Y2R_seg) == TRUE)
areNAstands = rMetrics$UniqueID[areNAindicies]

# Number of TS which had disturbances lasting longer than three years
longDistStands = auxSegData$UniqueID[which(auxSegData$DistDuration > 4)]

# Removing low magnitude disturbance
smalldist = rMetrics$UniqueID[which(rMetrics$Disturbance_mag < 0.20)]

# Going to remove stands that were probably not pine prior to disturbace
hardwoodCNV = rMetrics$UniqueID[which(rMetrics$pre_distMean < 0.20)]

# Number of TS which had a recovery within 4 years or less
shortRecovStands = rMetrics$UniqueID[which(rMetrics$Y2R_seg < 3)]

# An array which will hold all UniqueID's of stands removed from statistics
badStands2 = sort(as.numeric(unique(c(areNAstands,longDistStands,shortRecovStands,hardwoodCNV,smalldist))),decreasing = FALSE)

## Creating a new dataframe with only good values
goodDataDF = rMetrics[which(is.na(rMetrics$DistY) == FALSE),]
goodDataDF$fit_RMSE = auxSegData$fitRMSE[which(is.na(auxSegData$DistStartIndex) == FALSE)]

for (index in goodDataDF$UniqueID){
  ID = as.numeric(index)
  IDindex = which(goodDataDF$UniqueID == ID)
  bool = ID %in% as.numeric(badStands2)
  if (bool == TRUE){
    goodDataDF = goodDataDF[-IDindex,]
  }
}

plottingDF = goodDataDF

# final printing of sources of removed stands from statistics
cat("Total Number of Bad time series :", length(badStands2),"\n", 
    "Unusable segmentation filter :", length(areNAindicies),"\n",
    "Disturbance segment > 4 years :", length(longDistStands),"\n",
    "Recovery segment <= 2 years :",length(shortRecovStands), "\n",
    "Disturbance magnitude < 0.20 :",length(smalldist),"\n",
    "Pre-disturbance mean < 0.20 :",length(hardwoodCNV),"\n",
    "Good time series fit RMSE :", mean(goodDataDF$fit_RMSE)
)

# simple function for randomly pulling up stand images for validation of samples
randomLook = function(Nsamples){
  for (i in seq(1,Nsamples)){
    randomStand = sample(goodDataDF$UniqueID,1,replace = FALSE)
    shell(paste("C:/R_workspace/Collection2_data/images/GoodStands/",
                'plot_',
                randomStand,
                '_TimeSeriesGraph.svg',
                sep = ''),
          wait = FALSE
    )
  }
}

# different but similar function for looking at plots of a given set of samples
lookAtSamples = function(sampleNameVector){
  if (length(sampleNameVector) > 30){
    input = readline(prompt = paste("Number of samples is large (",length(sampleNameVector),"), display anyway? (Y/N) : "))
    if (input != 'Y'){
      print("Cancelled")
    } else {
      for (sampleN in sampleNameVector){
        shell(paste("C:/R_workspace/Collection2_data/images/",
                    'plot_',
                    sampleN,
                    '_TimeSeriesGraph.svg',
                    sep = ''),
              wait = FALSE
        )
      }
    }
  } else {
    for (sampleN in sampleNameVector){
      shell(paste("C:/R_workspace/Collection2_data/images/",
                  'plot_',
                  sampleN,
                  '_TimeSeriesGraph.svg',
                  sep = ''),
            wait = FALSE
      )
    }
  }
}

# a function that will make it easier to look up stands in google earth
GE_coords = function(standID) {
  attributeEntry = ecoRegionInfo[which(ecoRegionInfo$UniqueID == standID),]
  cat("Stand #",standID," : ",attributeEntry$lat,", ",attributeEntry$long, sep = "")
}

write.csv(plottingDF,file = "RecoveryMetrics_outputCOL2.csv")
## the below line can be used to read in previously written stand calculations
# plottingDF = read.csv("RecoveryMetrics_outputCOL2.csv")

#### Temporal trend analysis ####

# a function to print a normal Q-Q plot, metric label
#     should be what the name should appear as in the thesis. Note that for
#     NBR regrowth you have to hardcode the label because of the subscript
ggqqplot <- function(resid, metricLabel) {
  resid_sd <- sd(resid)
  resid_mean <- mean(resid)
  
  qq_data <- data.frame(sample = qnorm(ppoints(length(resid))), 
                        residual = sort((resid - resid_mean) / resid_sd))
  
  QQ = ggplot(qq_data, aes(sample, residual)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, 
                color = "red", linetype = "dashed") +
    ggtitle(metricLabel) +
    # uncomment below line for the NBR regrowth label to show up with subscript
    #ggtitle(expression(paste(NBR[regrowth]))) +
    xlab("Theoretical Quantiles") +
    ylab("Standardized Residual Quantiles")+
    theme(axis.text=element_text(size=11),
          axis.title=element_text(size=11),
          title = element_text(size=12)
          )
  
  ggsave(filename = paste(metricLabel,"_QQ.png"), 
         QQ,
         units = "in",
         height = 3.05, 
         width = 3.97, 
         dpi = 300, 
         device='png'
         )
  print(QQ)
}

# Need a quick funciton to calculate RMSE
calcRMSE = function(residualVector){
  sqError = residualVector^2
  meanSqError = mean(sqError)
  return(sqrt(meanSqError))
}

# input the analysis metric for this run
trendMetric = 'Y2R_seg' #'NBRregrowth', 'K_shift', 'Y2R_seg'
metricColumn = match(trendMetric, colnames(plottingDF))

# plotting, running OLS, adding OLS regression line to plot
plot.default(plottingDF$DistY, plottingDF[,metricColumn])
temporalOLS = lm(plottingDF[,metricColumn] ~ plottingDF$DistY)
abline(reg = temporalOLS)

## might need to include the thein-sen slope estimator for Y2R
theilSenEst = mblm(Y2R_seg ~ DistY, # CHANGE HERE
                   dataframe = plottingDF,
                   repeated = TRUE
                   )

# Printing results and displaying Normal q-q plot
cat('Trend results for',trendMetric)
summary(temporalOLS)
cat('RMSE = ',calcRMSE(temporalOLS$residuals))
print(shapiro.test(temporalOLS$residuals))
summary(theilSenEst)
cat('RMSE = ',calcRMSE(theilSenEst$residuals))
ggqqplot(temporalOLS$residuals,'Y2R') # CHANGE HERE

#### Plotting ####
statsPlotsPath = "C:/R_workspace/Collection2_data/StatsPlots"

metrics = c('NBRregrowth', 'K_shift', 'Y2R_seg')

# Need to create a subset of the data that fall within the time window I want
distWindowDF = plottingDF[which(plottingDF$DistY >= 1989 & plottingDF$DistY <= 2011),]

## Another graph to show the distribution of disturbances during the time frame
png(file = 'disturbanceFrequency.png',
    width = 6.5,
    height = 3.97,
    unit = 'in',
    res = 300
    )
par(cex.axis = 1, cex.lab = 1.3)
distFreq = hist.default(as.numeric(distWindowDF$DistY),
                        main = NA,
                        breaks = 18,
                        xlab = NA
)
dev.off()

distWindowDF$DistY = as.factor(distWindowDF$DistY)
distWindowDF$'DistYNUM' = as.numeric(distWindowDF$DistY)

for (aMetric in metrics){
  
  yearsPlot = ggplot(data = distWindowDF, 
                     aes_string(x = 'DistY', y = aMetric) #, color = 'groupID' ,group = 1
  )+
    geom_bin2d()+
    scale_fill_gradient(low = , high = )+
    stat_summary(fun = mean, geom ="point", size = 5, color = 'red')+
    stat_smooth(mapping = aes_string(x = 'DistYNUM', y = aMetric),
                size = 1.5,
                method = 'lm',
                se = TRUE,
                n = 5,
                show.legend = TRUE,
                color = 'yellow',
                linetype = 'solid'
    )+
    labs(title = paste("Variation in",aMetric,"by Disturbance Year"),
         xlab = "Year of Disturbance",
         fill = "Observation Count"
    )+
    theme(axis.text.x = element_text(angle = 45,vjust = 0.95,hjust = 0.9), 
          aspect.ratio = 4/8,
          legend.position = 'right'
    )
  print(yearsPlot)
}

## Creating a simple plot to show relationship between K-Shift and pre-disturbance mean
K_shiftRelationship = ggplot(aes(pre_distMean, K_shift),
                             data = plottingDF
                             )+
  geom_bin2d()+
  geom_hline(yintercept = 0,
             color = 'red',
             size = 1.5
  )+
  xlab("Pre-Disturbance Mean (NBR)")+
  ylab("K-Shift (NBR)")+
  labs(fill = "Observation Count")

print(K_shiftRelationship)
