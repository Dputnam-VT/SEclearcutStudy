#---------------------------COLLECTION 2 INVESTIGATION -------------------------
setwd("C:/R_workspace")
library(DBEST)
library(zoo)
library(lubridate)
library(tidyverse)
library(RColorBrewer)

dataDate1 = "2022_06_30" # using for the latest version of the time series csv
dataDate2 = "2022_08_12" # using for the latest version of the stand attribute dataset
dataDate3 = "2022_08_04" # using for the latest version of the climate variables csv

# Bringing in vegetation index time series
inputDF = read.csv(paste("NBR_timeSeriesDF2000_",dataDate1,".csv",sep = ""),header = TRUE)
TS_DF = inputDF[order(inputDF$X),]

EVI2_DF = read.csv(paste("EVI2_timeSeriesDF2000_",dataDate1,".csv",sep = ""),header = TRUE)
EVI2_DF = EVI2_DF[order(EVI2_DF$X),]

B5_DF = read.csv(paste("B5_timeSeriesDF2000_",dataDate1,".csv",sep = ""),header = TRUE)
B5_DF = B5_DF[order(B5_DF$X),]

TCC_DF = read.csv("AllCleanedData_ONLYtcc.csv",header = TRUE)
TCC_DF = TCC_DF[order(TCC_DF$UniqueID),]

# sznLength_DF = read.csv("seasonLength_timeSeriesDF_2022_07_14.csv",header = TRUE)
# sznLength_DF = sznLength_DF[order(sznLength_DF$X),]
# 
# peakDay_DF = read.csv("PeakDay_timeSeriesDF_2022_07_14.csv",header = TRUE)
# peakDay_DF = peakDay_DF[order(peakDay_DF$X),]
# 
# midSznLength_DF = read.csv("midSeasonLength_timeSeriesDF_2022_07_14.csv",header = TRUE)
# midSznLength_DF = midSznLength_DF[order(midSznLength_DF$X),]
# 
# EVI2area_DF = read.csv("EVI2_Area_timeSeriesDF_2022_07_14.csv",header = TRUE)
# EVI2area_DF = EVI2area_DF[order(EVI2area_DF$X),]

numSamples = length(inputDF$X)

reformatClimateData = function(unformattedDataframe,Nsamples){
  newDF <- data.frame(matrix(nrow = Nsamples, ncol = 39,byrow = TRUE))
  colnames(newDF) = c('UniqueID',as.character(seq(1984,2021)))
  rownames(newDF) = unique(unformattedDataframe$'UniqueID')
  
  lastIncluded = 0
  for (multiplier in seq(1,38)){
    newDF[,(multiplier+1)] = unformattedDataframe$'mean'[(lastIncluded+1):(Nsamples*multiplier)]
    lastIncluded = Nsamples*multiplier
  }
  newDF$UniqueID = unique(unformattedDataframe$'UniqueID')
  return(newDF)
}

freezeDays_DF = read.csv(paste("gridmetMetrics/freeze_TS_",dataDate3,".csv", sep = ''),header = TRUE)
freezeDays_DF = freezeDays_DF[,4:5]
freezeDays_DF = reformatClimateData(freezeDays_DF,numSamples)

GDD_DF = read.csv(paste("gridmetMetrics/GDD_TS_",dataDate3,".csv", sep = ''),header = TRUE)
GDD_DF = GDD_DF[,4:5]
GDD_DF = reformatClimateData(GDD_DF,numSamples)

meanVPD_DF = read.csv(paste("gridmetMetrics/meanVPD_TS_",dataDate3,".csv", sep = ''),header = TRUE)
meanVPD_DF = meanVPD_DF[,4:5]
meanVPD_DF = reformatClimateData(meanVPD_DF,numSamples)

totalPrecip_DF = read.csv(paste("gridmetMetrics/totalPrecip_TS_",dataDate3,".csv", sep = ''),header = TRUE)
totalPrecip_DF = totalPrecip_DF[,4:5]
totalPrecip_DF = reformatClimateData(totalPrecip_DF,numSamples)

meanTemp_DF = read.csv(paste("gridmetMetrics/meanTemp_TS_",dataDate3,".csv", sep = ''),header = TRUE)
meanTemp_DF = reformatClimateData(meanTemp_DF,numSamples)

# Bringing in topographic, soil, locational, and ecotonal data
ecoRegionInfo = read.csv(paste("standAttribute_values_",dataDate2,".csv",sep = ""),header = TRUE)

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
colnames(rMetrics) = c('YearOfCut','Recovery_min','Recovery_max','Recovery_mag',
                       'Disturbance_mag','Y2R_seg','Recovery_slope','NBRy3','NBRy5',
                       'NBRy7','NBRy10','pre_CutMean','post_CutMean','Y2R_80'
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
  svg(filename = paste('C:/R_workspace/images/ClearCutPlotsOnly/plot_',index,'_TimeSeriesGraph.svg',sep = ''),
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
               ylim = c(VImin,VImax),
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
  # 
  # # adding lines to indicate pre-disturbance and post-recovery means
  # abline(h = metricDF$pre_CutMean[(index+1)],
  #        col = 'azure3',
  #        lty = 3,
  #        lwd = 1
  #        )
  # abline(h = metricDF$post_CutMean[(index+1)],
  #        col = 'azure4',
  #        lty = 3,
  #        lwd = 1
  #        )
  # 
  # # adding line to show NBRy5 value
  # segments(x0 = dates[(auxDF[[index+1,2]]+5)],
  #          y0 = minVIval,
  #          x1 = dates[(auxDF[[index+1,2]]+5)],
  #          y1 = DBEST_object$Data[(auxDF[[index+1,2]]+5)]
  #          )
  # 
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
  legend(x = dates[(length(dates)-8)],
         y = -0.1,
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

# a function that will make it easier to look up stands in google earth
GE_coords = function(standID) {
  attributeEntry = ecoRegionInfo[which(ecoRegionInfo$UniqueID == standID),]
  cat("Stand #",standID," : ",attributeEntry$lat,", ",attributeEntry$long, sep = "")
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
  
  # ### REMOVE WHEN DONE TESTING ###
  # plot.default(ts1,type = 'l', col = 'black', ylim = c(-0.2,0.8))
  # lines(generalized$Fit,type = 'l', col = 'blue')
  # lines(output$Fit,type = 'l', col = 'red')

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
rMetrics$'K_shift' = rMetrics$post_CutMean - rMetrics$pre_CutMean
auxSegData$'UniqueID' = row.names(rMetrics)
auxSegData$'DistDuration' = auxSegData$DistEndIndex - auxSegData$DistStartIndex

# Experimenting with other recovery

rMetrics$'recovPcent3y' = (rMetrics$NBRy3 - rMetrics$Recovery_min) / (abs(rMetrics$Disturbance_mag)) * 100
rMetrics$'recovPcent5y' = (rMetrics$NBRy5 - rMetrics$Recovery_min) / (abs(rMetrics$Disturbance_mag)) * 100
rMetrics$'recovPcent7y' = (rMetrics$NBRy7 - rMetrics$Recovery_min) / (abs(rMetrics$Disturbance_mag)) * 100
rMetrics$'recovPcent10y' = (rMetrics$NBRy10 - rMetrics$Recovery_min) / (abs(rMetrics$Disturbance_mag)) * 100

# A function that finds index values for 3, 5, 7, 10 years post disturbance for a single stand
returnSyncVIs = function(standID,VItext,ts_length,sumMetric = FALSE) {
  VI_values = get(paste(VItext,'_DF',sep = ''))[which(get(paste(VItext,'_DF',sep = ''))[,1] == standID),2:length(get(paste(VItext,'_DF',sep = ''))[1,])]
  segmentationInfo = auxSegData[which(auxSegData$UniqueID == standID),]
  if (sumMetric == FALSE){
  
    if (length(VI_values) == 38) {
      
      timeSeries = zoo(x = as.numeric(VI_values),order.by = dates)
      timeSeries = na.locf(timeSeries, na.rm = FALSE, fromLast = TRUE) # filling in end values
      timeSeries = na.locf(timeSeries, na.rm = FALSE, fromLast = FALSE)
  
      VI_year3 = as.numeric(timeSeries[(min((segmentationInfo$DistEndIndex+3),ts_length))])
      VI_year5 = as.numeric(timeSeries[(min((segmentationInfo$DistEndIndex+5),ts_length))])
      VI_year7 = as.numeric(timeSeries[(min((segmentationInfo$DistEndIndex+7),ts_length))])
      VI_year10 = as.numeric(timeSeries[(min((segmentationInfo$DistEndIndex+10),ts_length))])
      
    } else if (length(VI_values) < 38){
      
      startY = substr(colnames(VI_values)[1],start = 2, stop = 5)
      endY = substr(colnames(VI_values)[length(colnames(VI_values))],start = 2, stop = 5)
      holder = as.vector(rep(NA,ts_length))
      timeSeries = zoo(x = holder,order.by = dates)
      startIndex = match(startY,year(timeSeries))
      endIndex = match(endY,year(timeSeries))
      timeSeries[startIndex:endIndex] = as.numeric(VI_values)
  
      VI_year3 = as.numeric(timeSeries[(min((segmentationInfo$DistEndIndex+3),ts_length))])
      VI_year5 = as.numeric(timeSeries[(min((segmentationInfo$DistEndIndex+5),ts_length))])
      VI_year7 = as.numeric(timeSeries[(min((segmentationInfo$DistEndIndex+7),ts_length))])
      VI_year10 = as.numeric(timeSeries[(min((segmentationInfo$DistEndIndex+10),ts_length))])
      
    } else {
      cat('Please use time series less than or equal to 38 years')
    }
  } else if (is.na(segmentationInfo$fitRMSE) == FALSE) {
    
    startY = substr(colnames(VI_values)[1],start = 1, stop = 5)
    endY = substr(colnames(VI_values)[length(colnames(VI_values))],start = 1, stop = 5)
    holder = as.vector(rep(NA,ts_length))
    timeSeries = zoo(x = holder,order.by = dates)
    startIndex = match(startY,year(timeSeries))
    endIndex = match(endY,year(timeSeries))
    timeSeries[startIndex:endIndex] = as.numeric(VI_values)
    
    VI_year3 = as.numeric(sum(timeSeries[segmentationInfo$DistEndIndex:(min((segmentationInfo$DistEndIndex+3),ts_length))],na.rm = TRUE))
    VI_year5 = as.numeric(sum(timeSeries[segmentationInfo$DistEndIndex:(min((segmentationInfo$DistEndIndex+5),ts_length))],na.rm = TRUE))
    VI_year7 = as.numeric(sum(timeSeries[segmentationInfo$DistEndIndex:(min((segmentationInfo$DistEndIndex+7),ts_length))],na.rm = TRUE))
    VI_year10 = as.numeric(sum(timeSeries[segmentationInfo$DistEndIndex:(min((segmentationInfo$DistEndIndex+10),ts_length))],na.rm = TRUE))
  } else {
    VI_year3 = NA
    VI_year5 = NA
    VI_year7 = NA
    VI_year10 = NA
  }
  return(c(VI_year3,VI_year5,VI_year7,VI_year10))
}

rMetrics$'EVI2y3' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'EVI2y5' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'EVI2y7' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'EVI2y10' = rep(NA,length(rMetrics$UniqueID))

rMetrics$'B5y3' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'B5y5' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'B5y7' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'B5y10' = rep(NA,length(rMetrics$UniqueID))

rMetrics$'TCCy3' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'TCCy5' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'TCCy7' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'TCCy10' = rep(NA,length(rMetrics$UniqueID))

rMetrics$'freezeDaysy3' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'freezeDaysy5' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'freezeDaysy7' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'freezeDaysy10' = rep(NA,length(rMetrics$UniqueID))

rMetrics$'GDDy3' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'GDDy5' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'GDDy7' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'GDDy10' = rep(NA,length(rMetrics$UniqueID))

rMetrics$'meanVPDy3' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'meanVPDy5' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'meanVPDy7' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'meanVPDy10' = rep(NA,length(rMetrics$UniqueID))

rMetrics$'totalPrecipy3' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'totalPrecipy5' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'totalPrecipy7' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'totalPrecipy10' = rep(NA,length(rMetrics$UniqueID))

rMetrics$'meanTempy3' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'meanTempy5' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'meanTempy7' = rep(NA,length(rMetrics$UniqueID))
rMetrics$'meanTempy10' = rep(NA,length(rMetrics$UniqueID))

# rMetrics$'sznLengthy3' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'sznLengthy5' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'sznLengthy7' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'sznLengthy10' = rep(NA,length(rMetrics$UniqueID))
# 
# rMetrics$'midSznLengthy3' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'midSznLengthy5' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'midSznLengthy7' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'midSznLengthy10' = rep(NA,length(rMetrics$UniqueID))
# 
# rMetrics$'peakDayy3' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'peakDayy5' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'peakDayy7' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'peakDayy10' = rep(NA,length(rMetrics$UniqueID))
# 
# rMetrics$'EVI2areay3' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'EVI2areay5' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'EVI2areay7' = rep(NA,length(rMetrics$UniqueID))
# rMetrics$'EVI2areay10' = rep(NA,length(rMetrics$UniqueID))

for (standNum in rMetrics$UniqueID) {
    for (VI in c('EVI2','B5','TCC')) {
      output = returnSyncVIs(standNum, VI, 38)
      startCol = match(paste(VI,'y3',sep = ''),colnames(rMetrics))
      endCol = startCol + 3
      standIndex = which(rMetrics$UniqueID == standNum)
      rMetrics[standIndex,startCol:endCol] = output
  }
}

for (standNum in rMetrics$UniqueID) {
  for (VI in c('freezeDays','GDD','meanVPD','totalPrecip','meanTemp')) {
    output = returnSyncVIs(standNum, VI, 38, TRUE)
    startCol = match(paste(VI,'y3',sep = ''),colnames(rMetrics))
    endCol = startCol + 3
    standIndex = which(rMetrics$UniqueID == standNum)
    rMetrics[standIndex,startCol:endCol] = output
  }
}

### Adding in the temporal change in these metrics ###

rMetrics$'freezeDaysChange' = apply(freezeDays_DF[34:39], MARGIN = 1,FUN = mean) - apply(freezeDays_DF[2:7], MARGIN = 1,FUN = mean)
rMetrics$'GDDChange' = apply(GDD_DF[34:39], MARGIN = 1,FUN = mean) - apply(GDD_DF[2:7], MARGIN = 1,FUN = mean)
rMetrics$'meanVPDChange' = apply(meanVPD_DF[34:39], MARGIN = 1,FUN = mean) - apply(meanVPD_DF[2:7], MARGIN = 1,FUN = mean)
rMetrics$'totalPrecipChange' = apply(totalPrecip_DF[34:39], MARGIN = 1,FUN = mean) - apply(totalPrecip_DF[2:7], MARGIN = 1,FUN = mean)
rMetrics$'meanTempChange' = apply(meanTemp_DF[34:39], MARGIN = 1,FUN = mean) - apply(meanTemp_DF[2:7], MARGIN = 1,FUN = mean)

#plottingDF$'meanTempChange'

rMetrics$'NonPSAy3' = atan2(rMetrics$EVI2y3,rMetrics$B5y3)
rMetrics$'NonPSAy5' = atan2(rMetrics$EVI2y5,rMetrics$B5y5)
rMetrics$'NonPSAy7' = atan2(rMetrics$EVI2y7,rMetrics$B5y7)
rMetrics$'NonPSAy10' = atan2(rMetrics$EVI2y10,rMetrics$B5y10)

# export table I can add back into arcgis pro
exportDF = rMetrics[,c(14,32:length(colnames(rMetrics)))]
write.csv(exportDF,file = "additionalVIvals2.csv")

write.csv(rMetrics,file = "RecoveryMetrics_output.csv")

#### Cleaning data based on practical relationships ####

# Number of TS which did not have qualifying segments for metric calculation
areNAindicies = which(is.na(rMetrics$Y2R_seg) == TRUE)
areNAstands = rMetrics$UniqueID[areNAindicies]

# Number of TS which had large gaps between disturbance and recovery segments
bigGapIndices = which(auxSegData$gainLossGap > 5)
bigGapLengths = auxSegData$gainLossGap[bigGapIndices]
bigGapStands = auxSegData$UniqueID[bigGapIndices]

# Number of TS which had disturbances lasting longer than three years
longDistStands = auxSegData$UniqueID[which(auxSegData$DistDuration > 4)]

# Experimental, there are many stands which because of segmentation or noise have longer disturbances
#     seem to be able to keep the good ones by including a disturbance magnitude condition
expFilter = auxSegData$UniqueID[which((auxSegData$DistDuration > 4) & (rMetrics$Disturbance_mag < 0.3))]

# Removing low magnitude disturbance
smalldist = rMetrics$UniqueID[which(rMetrics$Disturbance_mag < 0.20)]
  
# Going to remove stands that were probably not pine prior to disturbace
hardwoodCNV = rMetrics$UniqueID[which(rMetrics$pre_CutMean < 0.20)]

# Number of TS which had a recovery within 4 years or less
shortRecovStands = rMetrics$UniqueID[which(rMetrics$Y2R_seg < 3)]

# An array which will hold all UniqueID's of stands removed from statistics
badStands = sort(as.numeric(unique(c(areNAstands,bigGapStands,longDistStands,shortRecovStands))),decreasing = FALSE)
badStands2 = sort(as.numeric(unique(c(areNAstands,longDistStands,shortRecovStands,hardwoodCNV,smalldist))),decreasing = FALSE)

# # trying out an outlier test
# findExtremeStands = function(metricValues,IDs,c){
#   stats = summary(metricValues)
#   anIQR = IQR(metricValues, na.rm = TRUE)
#   lower = stats[2] - (c*anIQR)
#   upper = stats[5] + (c*anIQR)
#   indices = which(metricValues > upper | metricValues < lower)
#   return(IDs[indices])
# } # can use either 1 or 3 for c (3 for extreme values)

## Creating a new dataframe with only good values
goodDataDF = rMetrics[which(is.na(rMetrics$YearOfCut) == FALSE),]
goodDataDF$fit_RMSE = auxSegData$fitRMSE[which(is.na(auxSegData$DistStartIndex) == FALSE)]
# 
# goodDataDF = goodDataDF[-which(goodDataDF$UniqueID == badStands2),]

for (index in goodDataDF$UniqueID){
  ID = as.numeric(index)
  IDindex = which(goodDataDF$UniqueID == ID)
  bool = ID %in% as.numeric(badStands2)
  if (bool == TRUE){
    goodDataDF = goodDataDF[-IDindex,]
  }
}

#### the above code is a much simplier and more transparent way of selecting the good stands
# for (ID in goodDataDF$UniqueID) {
#   gapLength = auxSegData$gainLossGap[which(auxSegData$UniqueID == ID)]
#   distDuration = as.numeric(auxSegData$DistDuration[which(auxSegData$UniqueID == ID)])
#   recovDuration = goodDataDF$Y2R[goodDataDF$UniqueID == ID]
#   if ((distDuration > 3) | (recovDuration < 5)) {
#     indexValue = which(goodDataDF$UniqueID == ID)
#     goodDataDF = goodDataDF[-indexValue,]
#   }
# } # removing stands with large gaps between disturbance end and recovery start

# # Copying only the good files into a seperate directory within all the plot images folder
# allPlotsDirPath = "C:/R_workspace/images/ClearCutPlotsOnly"
# goodDirPath = "C:/R_workspace/images/ClearCutPlotsOnly/CleanedData"
# for (standI in goodDataDF$UniqueID){
#   filePath = paste('plot_',standI,'_TimeSeriesGraph.svg',sep = '')
#   file.copy(from = paste(allPlotsDirPath,'/',filePath, sep = ''),
#             to = goodDirPath,
#             overwrite = TRUE,
#             recursive = FALSE,
#             copy.mode = TRUE,
#             copy.date = TRUE
#         )
# }
# 
# # Copying only the bad stands graphs into a seperate directory within all the plot images folder
# badDirPath = "C:/R_workspace/images/ClearCutPlotsOnly/BadStands"
# for (standI in badStands2){
#   filePath = paste('plot_',as.character(standI),'_TimeSeriesGraph.svg',sep = '')
#   file.copy(from = paste(allPlotsDirPath,'/',filePath, sep = ''),
#             to = badDirPath,
#             overwrite = TRUE,
#             recursive = FALSE,
#             copy.mode = TRUE,
#             copy.date = TRUE
#   )
# }

# merging the ecoregion info and good data for plotting #
ecoRegionInfo$intLat = floor(ecoRegionInfo$lat)
ecoRegionInfo$ownership_mode = replace_na(ecoRegionInfo$ownership_mode,999)
ecoRegionInfo$ownership_mode = as.factor(ecoRegionInfo$ownership_mode)
plottingDF = merge(x = goodDataDF,
                   y = ecoRegionInfo,
                   by.x = 'UniqueID',
                   by.y = 'UniqueID',
                   all.x = TRUE
                   )

SSURGOdata = read.csv(file = 'allStandsSSURGOdata.csv')
plottingDF = merge(x = plottingDF,
                   y = SSURGOdata,
                   by.x = 'UniqueID',
                   by.y = 'UniqueID',
                   all.x = TRUE
)

SCMCdata = read.csv(file = 'cleanedData_SCMC.csv')
colnames(SCMCdata) = c('UniqueID','groupID','memProb')
plottingDF = merge(x = plottingDF,
                   y = SCMCdata,
                   by.x = 'UniqueID',
                   by.y = 'UniqueID',
                   all.x = TRUE
)

# Setting some of the categorical data columns to be factors #
plottingDF$aspect_mode = as.factor(plottingDF$aspect_mode)
plottingDF$soil0_mode = as.factor(plottingDF$soil0_mode)
plottingDF$soil30_mode = as.factor(plottingDF$soil30_mode)
plottingDF$soil100_mode = as.factor(plottingDF$soil100_mode)
plottingDF$ecoR_code = as.character(plottingDF$Classified)
plottingDF$plantation_class = as.factor(plottingDF$plantation_class)
plottingDF$soilGG_mode = as.factor(plottingDF$soilGG_mode)
plottingDF$groupID = as.factor(plottingDF$groupID)
plottingDF$'rel_Kshift' = ((plottingDF$post_CutMean - plottingDF$pre_CutMean) / plottingDF$pre_CutMean) * 100 #(plottingDF$post_CutMean / plottingDF$pre_CutMean) * 100 #

write.csv(plottingDF,file = "RecoveryMetrics_output.csv")

#### Experimenting with CART analysis ####
library(rpart)
library(rpart.plot)
library(ranger)
library(lctools)
library(dplyr)
library(corrplot)

# Looking at correlation matrix for all variables
onlyNumeric = unlist(lapply(plottingDF, is.numeric))  
correlationMatrix = cor(plottingDF[,onlyNumeric], 
                        use = 'everything'
                        )
corrplot(correlationMatrix)

# Determining spatial autocorrelation
spatialWeights = w.matrix(Coords = data.frame(plottingDF$long,plottingDF$lat),
                          Bandwidth = 50,
                          WType = 'Binary',
                          family = 'adaptive'
                          )

morans = moransI.w(plottingDF$K_shift,spatialWeights)

localMorans = l.moransI(Coords = data.frame(plottingDF$long,plottingDF$lat),
                        Bandwidth = 500,
                        x = plottingDF$K_shift,
                        WType = 'Binary',
                        scatter.plot = TRUE,
                        family = 'adpative'
                        )

plot.default(plottingDF$pre_CutMean, plottingDF$rel_Kshift, col = plottingDF$YearOfCut)
testModel = lm(formula = 'rel_Kshift ~ pre_CutMean', data = log(plottingDF[,c('rel_Kshift', 'pre_CutMean')]))
plot.default(plottingDF$rel_Kshift, exp(testModel$fitted.values))
abline(0,1)

# Model formula
formulaString = "recovPcent5y ~ aspect_mode + elevation_mean + slope_mean + 
                  soil0_mode + soil30_mode + soil100_mode + pH30_mean + 
                  ecoR_code + ownership_mode + freezeDaysy3 + freezeDaysy5 + 
                  freezeDaysy7 + freezeDaysy10 + GDDy3 + GDDy5 + GDDy7 + 
                  GDDy10 + meanVPDy3 + meanVPDy5 + meanVPDy7 + meanVPDy10 + 
                  totalPrecipy3 + totalPrecipy5 + totalPrecipy7 + 
                  totalPrecipy10 + meanTempy3 + meanTempy5 + meanTempy7 +
                  meanTempy10 + freezeDaysChange + GDDChange + meanVPDChange + 
                  totalPrecipChange + meanTempChange + plantation_class + 
                  soilGG_mode + SI_25 + SI_50 + Recovery_min + YearOfCut + 
                  pre_CutMean + lat + long"

# # Seperating train from test data 70% train, 30% test
# trainIndices = sample(x = (1:length(plottingDF$UniqueID)), size = (0.7*length(plottingDF$UniqueID)))
# trainData = plottingDF[trainIndices,]
# testIndices = (1:length(plottingDF$UniqueID))[-trainIndices]
# testData = plottingDF[testIndices,]

# Fitting the regression CART model for K_shift
CART_fit = rpart(formula = formulaString,
              data = plottingDF,
              method = 'anova',
              cp = 0
              )
rsq.rpart(CART_fit)
# Finding the lowest cp value within 1SE of the lowest cross validated error
lowERRcpI = as.numeric(which.min(CART_fit$cptable[,4]))
lowERRcp = as.numeric(CART_fit$cptable[lowERRcpI,1])
optimumCP = which(CART_fit$cptable[,4] <= (CART_fit$cptable[lowERRcpI,4] + CART_fit$cptable[lowERRcpI,5]))
optimumCP = as.double(CART_fit$cptable[min(optimumCP),1])

# CART_test = rpart.predict(CART_fit, type = 'matrix')
# 
# # Going to look at the highest residuals from this model
# CARTResids = sort(abs(residuals(CART_fit)),decreasing = TRUE)
# residQs = quantile(CARTResids,c(0.80,0.85,0.90,0.95))
# Qresids = CARTResids[which(CARTResids >= residQs[4])]
# residsI = as.numeric(names(Qresids))
# residsStands = plottingDF$UniqueID[residsI]

# # Testing the regression tree
# CART_test = rpart.predict(object = CART_fit,
#                           type = 'matrix',
#                           newdata = testData
#                           )

# # RMSE of the test fit
# CART_rmse = sqrt(mean((CART_test - testData$K_shift)^2))
# 
# # Finding the number of splits with the lowest complexity
# lowestCP = CART_fit$cptable[which.min(CART_fit$cptable[,'xerror']),'CP']

# Pruning the tree to the lowest complexity parameter
CART_fit2 = prune.rpart(CART_fit, cp = optimumCP)
paste('Pruned CART R^2',round(1-min(CART_fit2$cptable[,4]),digits = 4))
rpart.plot(CART_fit2)
# CART_test2 = rpart.predict(object = CART_fit2,
#                           type = 'matrix',
#                           newdata = testData
# )
# CART2_rmse = sqrt(mean((CART_test2 - testData$K_shift)^2))

noSIstands = which(is.na(plottingDF$SI_25))

### Trying out random forest ###
RF_1 = ranger(formula = formulaString,
              data = plottingDF[c(-178,-3448,-2596,-noSIstands),],
              num.trees = 2000,
              importance = 'impurity',
              replace = TRUE,
              respect.unordered.factors = TRUE,
              oob.error = TRUE,
              num.threads = 8,
              write.forest = TRUE
)
# Checking out accuracy of Random Forest
paste("Random Forest R2",round(RF_1$r.squared, digits = 3))
plot.default(plottingDF$recovPcent5y [c(-178,-3448,-2596,-noSIstands)],RF_1$predictions) #xlim = c(-50,200),ylim = c(-50,200)
abline(0,1)

## Plotting variable importance c(RF_1$variable.importance,CART_fit$variable.importance)
varIdata = sqrt(sort(CART_fit2$variable.importance, decreasing = TRUE))[0:floor(0.6*length(RF_1$variable.importance))]
varGraph = ggplot(mapping = aes(x = reorder(names(varIdata),as.numeric(varIdata)), y = as.numeric(varIdata))) +
  geom_bar(width = 1, 
           stat="identity",
           colour = "black", 
           fill="lightblue",
           position = position_dodge()
           )+
  theme(aspect.ratio = 2/3)+
  coord_polar(theta = "x", start=0)+
  labs(x = NULL,y = 'Variable Importance Metric')
print(varGraph)

## Plotting variable importance c(RF_1$variable.importance,CART_fit$variable.importance)
varIdata = sqrt(sort(RF_1$variable.importance, decreasing = TRUE))[0:floor(0.6*length(RF_1$variable.importance))]
varGraph = ggplot(mapping = aes(x = reorder(names(varIdata),as.numeric(varIdata)), y = as.numeric(varIdata),fill = as.numeric(varIdata))) +
  geom_bar(width = 1, 
           stat="identity",
           colour = "black",
           position = position_dodge()
  )+
  coord_flip()+
  scale_fill_continuous(low="blue", high="red")+
  theme(aspect.ratio = 2/3,
        text = element_text(size = 16)
        )+
  guides(fill=guide_legend(title=""))+
  labs(x = NULL,y = 'Gini Variable Importance')
print(varGraph)

# Going to try the CART model with binned K_shift rather than a continuous regression
test = ecdf(plottingDF$K_shift)

calcConfusion = function(confusionMatrix){
  inputDim = dim(confusionMatrix)
  accuracyTable = matrix(ncol = 3,nrow = inputDim[1])
  colnames(accuracyTable) = c('Users','producers','overall')
  rownames(accuracyTable) = seq(1,inputDim[1])
  for (classN in seq(1,inputDim[1])){
    accuracyTable[classN,1] = confusionMatrix[classN,classN] / sum(confusionMatrix[classN,])
    accuracyTable[classN,2] = confusionMatrix[classN,classN] / sum(confusionMatrix[,classN])
  }
  total = sum(confusionMatrix[,c(seq(1,inputDim[1]))])
  accuracyTable[1,3] = sum(diag(confusionMatrix)) / total
  return(accuracyTable)
}

### trying to bin values with a given number of bins (equal number of samples per bin)
rowOrder = data.frame(sort(plottingDF$K_shift,decreasing = FALSE,index.return = TRUE))[,2]
plottingDF$'binKshift'[rowOrder] = ntile(sort(plottingDF$K_shift, decreasing = FALSE), n = 3)
#plottingDF$'binKshift'[c(91,2374,3393)] = 2 # the above code misclassifies 3 values as negative
plottingDF$'binKshift' = as.factor(plottingDF$'binKshift')

### this code can be used to bin values based on quantiles
# KshiftQs = quantile(plottingDF$K_shift,c(0.0,0.1241812,0.345,0.540,0.75,1.0), digits = 3) # this is just for < 0 OR > 0
# binnedKshift = data.frame(findInterval(plottingDF$K_shift,KshiftQs,all.inside = TRUE))
# plottingDF$'binKshift' = as.factor(binnedKshift$findInterval.plottingDF.K_shift..KshiftQs.)

# I'm trying out removing the other soil modes (soil0_mode + soil100_mode)
formulaString2 = "binKshift ~ aspect_mode + elevation_mean + slope_mean + 
                  soil0_mode + soil30_mode + soil100_mode + pH30_mean + 
                  ecoR_code + ownership_mode + freezeDaysy3 + freezeDaysy5 + 
                  freezeDaysy7 + freezeDaysy10 + GDDy3 + GDDy5 + GDDy7 + GDDy10 + 
                  meanVPDy3 + meanVPDy5 + meanVPDy7 + meanVPDy10 + 
                  totalPrecipy3 + totalPrecipy5 + totalPrecipy7 + totalPrecipy10 +
                  meanTempy3 + meanTempy5 + meanTempy7 + meanTempy10 + 
                  freezeDaysChange + GDDChange + meanVPDChange + 
                  totalPrecipChange + meanTempChange + plantation_class + 
                  soilGG_mode"
CART_fit3 = rpart(formula = formulaString2,
                  data = plottingDF,
                  method = 'class',
                  cp = 0
                  )
paste('Pruned binned K_shift CART R^2',round(1-min(CART_fit3$cptable[,4]),digits = 4))
CART_fit3 = prune.rpart(CART_fit3,cp = 0.00470925)
CART_confM = table(predict(CART_fit3, type = 'class'),plottingDF$binKshift)
print(CART_confM)
print(calcConfusion(CART_confM))

### Trying out random forest for binned K_shift ###
RF_2 = ranger(formula = formulaString2,
            data = plottingDF[c(-178,-3448,-2596,-noSIstands),],
            num.trees = 2000,
            importance = 'impurity_corrected',
            replace = TRUE,
            respect.unordered.factors = TRUE,
            oob.error = TRUE,
            num.threads = 8,
            write.forest = TRUE,
            classification = TRUE,
            verbose = TRUE
)
print(RF_2$confusion.matrix)
print(calcConfusion(RF_2$confusion.matrix))

#### Temporal trend analysis ####
library(Kendall)

yearsMedian = aggregate.data.frame(plottingDF$K_shift,
                                   by = list(plottingDF$YearOfCut),
                                   FUN = median,
                                   simplify = TRUE
                                   )

trendTest = Kendall(yearsMedian$Group.1,yearsMedian$x)
summary(trendTest)
plot.default(yearsMedian$Group.1,yearsMedian$x)
temporalTrend = lm('x ~ Group.1',data = yearsMedian)
summary(temporalTrend)
abline(temporalTrend$coefficients)

#### Looking at regional groups ####
group1 = c(8,6,4,5,3)
group2 = c(7,2,1)

temporalVariable = 'freezeDaysy3'

## Looking at temporal trend for coastal regions
rows1 = plottingDF$groupID %in% group1
G1data = plottingDF[rows1,]
yearsMedianG1 = aggregate.data.frame(G1data[temporalVariable],
                                   by = list(G1data$YearOfCut),
                                   FUN = median,
                                   simplify = TRUE
)

trendTestG1 = Kendall(yearsMedianG1[,1], yearsMedianG1[,2])
summary(trendTestG1)
g1zoo = zoo(x = yearsMedianG1[,2],
            order.by = yearsMedianG1$Group.1
            )
g1ACF = acf(x = g1zoo,
        type = 'correlation',
        plot = TRUE
        )
print(g1ACF)
plot.default(yearsMedianG1[,1], yearsMedianG1[,2], type = 'l', main = paste('Coastal Stands',temporalVariable,'Temporal Trend'))
temporalTrendG1 = lm(paste(temporalVariable,' ~ Group.1', sep = ''),data = yearsMedianG1)
abline(temporalTrendG1$coefficients, col = 'blue')

## Looking at temporal trend for inland regions
rows2 = plottingDF$groupID %in% group2
G2data = plottingDF[rows2,]
yearsMedianG2 = aggregate.data.frame(G2data[temporalVariable],
                                     by = list(G2data$YearOfCut),
                                     FUN = median,
                                     simplify = TRUE
)

trendTestG2 = Kendall(yearsMedianG2[,1], yearsMedianG2[,2])
summary(trendTestG2)
g2zoo = zoo(x = yearsMedianG2[,2],
            order.by = yearsMedianG1$Group.1
)
g2ACF = acf(x = g2zoo,
            type = 'correlation',
            plot = TRUE
)
print(g2ACF)
plot.default(yearsMedianG2[,1], yearsMedianG2[,2], type = 'l', main = paste('In-land Piedmont Stands',temporalVariable,'Temporal Trend'))
temporalTrendG2 = lm(paste(temporalVariable,' ~ Group.1', sep = ''),data = yearsMedianG2)
abline(temporalTrendG2$coefficients, col = 'blue')

## Plot to look at cluster temporal trends together
combinedTrends = ggplot(mapping = aes_string(x = 'YearOfCut', y = temporalVariable, color = 'groupID'),
                        data = plottingDF
                        )+
                        geom_point(shape = 4)+
                        stat_smooth(method = 'lm')
print(combinedTrends)

## Looking at TCC values
onlyTCCrows = which(is.na(plottingDF$TCCy10) == FALSE)
onlyTCC_DF = plottingDF[onlyTCCrows,]

onlyNumeric = unlist(lapply(onlyTCC_DF, is.numeric))  
correlationMatrix = cor(onlyTCC_DF[,onlyNumeric], 
                        use = 'everything'
)
corrplot(correlationMatrix)


## Boxplots to describe spatial cluster
arcmapSymb = c('#78AAFF','#FF6455','#7DDC55','#FFB400','#C864E1','#BEA064','#FABEC8','#AFAFAF')
clusterAttri = ggplot(data = plottingDF[-104,],
                      mapping = aes_string(y = 'pre_CutMean', fill = 'groupID')
                      )+
  geom_boxplot(notch = TRUE,
               varwidth = TRUE,
               coef = 1.0
               )+
  scale_fill_manual(values = arcmapSymb
                    )+
  theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank()
  )+
  labs(xlab = '', fill = 'Cluster ID')
print(clusterAttri)

## Plot for relationship between pre-disturbance mean and relative K_shift
Krelation = ggplot(data = plottingDF,
                   mapping = aes_string(y = 'pre_CutMean', x = 'rel_Kshift') # color = 'groupID'
                   )+
  #geom_point(shape = 4, stroke = 0.967)+
  geom_hex()+
  theme(legend.position = 'right',
        aspect.ratio = 2/3,
        text = element_text(size = 16)
  )+
  #guides(fill=guide_legend(title="Cluster ID"))+
  labs(x = "% Change in NBR Post-Disturbance",
       y = 'Pre-Disturbance Mean',
       fill = "Observation Count")
print(Krelation)
  
#### new data visualization plotting methods ####

# making a list of recovery metrics to plot
metrics = c('Recovery_max','Recovery_mag','Recovery_slope','K_shift','Y2R_seg','NonPSAy3','NonPSAy5','NonPSAy7','NonPSAy10')
metrics = append(metrics,colnames(plottingDF)[match('recovPcent3y',colnames(plottingDF)) : match('B5y10',colnames(plottingDF))])
metrics = append(metrics,colnames(plottingDF)[match('NBRy3',colnames(plottingDF)):match('NBRy10',colnames(plottingDF))])

# path for export of stats plots
statsPlotsPath = "C:/R_workspace/images/StatsPlots"

#---------------- Ungrouped Stats plots for each metric ----------------------#
# histograms
for (aMetric in metrics){
  summaryStats = ggplot(data = plottingDF, 
                         aes_string(x = aMetric))+ 
    geom_histogram(bins = 20)+
    #geom_vline(aes(xintercept = mean(aMetric)), color = 'red')
    labs(title = paste("Distribution of",aMetric, "Metric"))
    print(summaryStats)
    ggsave(filename = paste(aMetric,'_histogram.svg',sep = ''),
           plot = summaryStats,
           device = 'svg',
           path = statsPlotsPath
           )
}
# boxplots
# have to make a seperate DF for just the metrics for this one
metricsDF = plottingDF %>% select(all_of(metrics[c(-5,-10,-11,-12,-13)]))
summaryStats = ggplot(data = stack(metricsDF), 
                      aes(x = ind, y = values)
                      )+ 
  geom_boxplot()+
  #geom_vline(aes(xintercept = mean(aMetric)), color = 'red')
  labs(title = paste("Distribution of Selected Recovery Metrics"))+
  xlab(NULL)+
  ylab("NBR Value")+
  theme(axis.text.x = element_text(angle = 20,vjust = 0.95,hjust = 0.9)
  )
print(summaryStats)
ggsave(filename = 'selectMetricsBoxplots.svg',
       plot = summaryStats,
       device = 'svg',
       path = statsPlotsPath
      )

#---------------- Stratifying by disturbance year ----------------------------#
# Need to create a subset of the data that fall within the time window I want
distWindowDF = plottingDF[which(plottingDF$YearOfCut >= 1989 & plottingDF$YearOfCut <= 2011),]
distWindowDF$YearOfCut = as.factor(distWindowDF$YearOfCut)
distWindowDF$intLat =as.factor(distWindowDF$intLat)
distWindowDF$'YearOfCutNUM' = as.numeric(distWindowDF$YearOfCut)

# Plotting for poster
yearsMedian = aggregate.data.frame(distWindowDF$K_shift,
                                   by = list(as.numeric(distWindowDF$YearOfCut)),
                                   FUN = median,
                                   simplify = TRUE
)

trendTest = Kendall(yearsMedian$Group.1,yearsMedian$x)
summary(trendTest)
plot.default(yearsMedian$Group.1,yearsMedian$x)
temporalTrend = lm('x ~ Group.1',data = yearsMedian)
summary(temporalTrend)
abline(temporalTrend$coefficients)

distWindowDF$YearOfCut = as.factor(distWindowDF$YearOfCut)

for (aMetric in metrics){
  
  # columnIndex = which(colnames(distWindowDF) == aMetric)
  # linearModel = lm(paste(aMetric,"~ as.numeric.distWindowDF.YearOfCut.", sep = ""),data = data.frame(as.numeric(distWindowDF$YearOfCut),distWindowDF[aMetric]))
  
  yearsPlot = ggplot(data = distWindowDF, 
                     aes_string(x = 'YearOfCut', y = aMetric) #, color = 'groupID' ,group = 1
                     )+
    geom_bin2d()+
    scale_fill_gradient(low = , high = )+
    #geom_point(mapping = aes(fill = 'Stand Value', shape = "Stand Value"), size = 2.5, stroke = 0.5)+
    stat_summary(fun = median, geom ="point", size = 5, color = 'red')+
    stat_smooth(mapping = aes_string(x = 'YearOfCutNUM', y = aMetric),
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
    #scale_color_brewer(type = 'qual', palette = 'Set2')
    #geom_text(mapping = aes(x = 15,y = ),label = eq(distWindowDF$'YearOfCutNUM',distWindowDF[,columnIndex]), parse = TRUE)
    # scale_fill_manual(name = "", values = c("Stand Value" = "black", "Yearly_Median" = 'red'))+
    # scale_color_manual(name = "", values = c('Stand Value' = 'black', "Temporal_Trend" = 'blue'))+
    # scale_shape_manual(name = "", values = c("Stand Value" = 4, "Yearly_Median" = 21))+
    # scale_linetype_manual(name = "", values = c("Temporal_Trend" = 1))
    #guides(shape = guide_legend(override.aes = list(linetype = 1)))
    print(yearsPlot)
    ggsave(filename = paste(aMetric,'_temporalTrend.svg',sep = ''),
           plot = yearsPlot,
           device = 'svg',
           path = statsPlotsPath,
           width = 8,
           height = 4,
           limitsize = FALSE
           )
}

#---------------- Stratifying by whole lattitudes ----------------------------#
plottingDF$intLat =as.factor(plottingDF$intLat)
for (aMetric in metrics){
  lattitudePlot = ggplot(data = plottingDF, 
                         aes_string(x = 'intLat', y = aMetric)
                         ) + 
    geom_violin(trim = TRUE)+
    stat_summary(fun = mean,geom = 'crossbar', color = 'black',width = 0.7)+
    labs(title = paste("Variation in",aMetric,"Within Degrees of Latitude"))
    xlab(NULL)+
    theme(axis.text.x = element_text(angle = 20,vjust = 0.95,hjust = 0.9), 
          legend.position = 'none'
    )+
  print(lattitudePlot)
}

#-------------------- Stratifying by Ecoregions ------------------------------#
plottingDF$ecoR_code = as.character(plottingDF$Classified)
for (aMetric in metrics){
  ecoRegionPlot = ggplot(data = plottingDF, 
                         aes_string(x = 'ecoR_code', y = aMetric, fill = 'ecoR_code')
                         )+ 
    geom_violin(trim = TRUE)+
    stat_summary(fun = mean,geom = 'crossbar', color = 'black',width = 0.7)+
    labs(title = paste("Variation in",aMetric,"Within EcoRegions"))+
    xlab(NULL)+
    theme(axis.text.x = element_text(angle = 20,vjust = 0.95,hjust = 0.9), 
          legend.position = 'none'
          )+
    scale_fill_brewer(palette = 'Paired')
  print(ecoRegionPlot)
}

#---------------- Looking at K-Shift by toil texture class -------------------#
soilDF = plottingDF[-which(is.na(plottingDF$soil0_mode) == TRUE),]
soilDF$K_shift = as.numeric(soilDF$K_shift)
#plottingDF$soil0_mode = as.factor(plottingDF$soil0_mode)
kshiftSoil = ggplot(data = soilDF, 
                    mapping = aes(group = soilDF$soil30_mode, y = soilDF$K_shift)
                    )+ 
  geom_boxplot()+
  #geom_vline(aes(xintercept = mean(aMetric)), color = 'red')
  labs(title = paste("K-Shift Distribution by Soil Texture Class"))+
  xlab(NULL)+
  ylab("K-Shift Value (NBR)")+
  theme(axis.text.x = element_text(angle = 20,vjust = 0.95,hjust = 0.9)
  )
print(kshiftSoil)

write.csv(plottingDF,file = "AllCleanedData_output.csv")

# final printing of sources of removed stands from statistics
totalNumberOfBadDataPoints = length(bigGapIndices) + length(areNAindicies) + length(longDistStands) + length(shortRecovStands)
cat("Total Number of Bad time series :", length(badStands),"\n", 
    "Unusable segmentation filter :", length(areNAindicies),"\n",
    "Disturbance segment > 3 years :", length(bigGapIndices),"\n",
    "Recovery segment <= 2 years :",length(shortRecovStands), "\n",
    "Disturbance -> recovery gap > 3 years :",length(bigGapIndices),"\n",
    "Good time series fit RMSE :", mean(plottingDF$fit_RMSE)
    )
