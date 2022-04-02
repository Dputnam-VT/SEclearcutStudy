setwd("C:/R_workspace")
library(DBEST)
library(zoo)
library(lubridate)
TS_DF = read.csv("timeSeriesDF500_3.csv",header = TRUE)

## creating the dates for the time series index
startYear = colnames(TS_DF)[2]
startYear = substr(startYear, 2,5)
startDate = paste(startYear,'-11-01',sep = '')
startDate = as.Date.character(startDate,format = '%Y-%m-%d')

endYear = colnames(TS_DF)[length(TS_DF)]
endYear = substr(endYear, 2,5)
endDate = paste(endYear,'-11-01',sep = '')
endDate = as.Date.character(endDate,format = '%Y-%m-%d')

dates = seq(startDate,endDate, by = 'years')

# back filling missing values at beginning of time-series with first value
for (arrayNum in TS_DF[,1]) {
  aArray = TS_DF[(arrayNum+1),]
  NonNAindex = which(!is.na(aArray))
  if (length(NonNAindex) != length(aArray)){
    firstNonNA = NonNAindex[2]
    fillVal = as.double(aArray[firstNonNA])
    for (i in seq(from = 2, to = (firstNonNA-1))) {
      aArray[i] = fillVal
    }
    TS_DF[(arrayNum+1),] = aArray
  }
}

# Dataframe to store recovery metrics for each stand
rMetrics <- data.frame(matrix(nrow = length(rownames(TS_DF)), ncol = 8,byrow = TRUE))
colnames(rMetrics) = c('YearOfCut','Recovery_min','Recovery_max','Recovery_mag','Y2R','Recovery_slope','pre_CutMean','post_CutMean')
rownames(rMetrics) = TS_DF$X

# Dataframe to store auxiliary data for plotting (start dates for segments, etc.)
auxSegData <- data.frame(matrix(nrow = length(rownames(TS_DF)), ncol = 5,byrow = TRUE))
colnames(auxSegData) = c('DistStartIndex','DistEndIndex','RecoveryStartIndex','RecoveryEndIndex','gainLossGap')
rownames(auxSegData) = TS_DF$X

## Function to plot the DBEST output correctly
plot_DBEST = function(DBEST_object,index,metricDF,auxDF) {
  # getting output device set up
  svg(filename = paste('C:/R_workspace/images/ClearCutPlotsOnly/plot_',index,'_TimeSeriesGraph.svg',sep = ''),
      width = 4,
      height = 2,
      pointsize = 4,
      onefile = FALSE
      )
  
  # plotting the generalization lines furthest back
  plot.default(dates,as.vector(DBEST_object$Fit),
               type = 'l',
               xlab = 'YEAR',
               ylab = 'NBR',
               main = paste("DBEST Result : Stand #",index,sep = ""),
               lwd = 3,
               col = 'blue',
               ylim = c(-0.2,1),
               axes = FALSE,
  )
  
  # getting axes customized
  box()
  axis.Date(1,at = seq(from = (startDate+365), to = endDate, by = '5 years'),labels = TRUE)
  axis(2)
  
  # adding lines for the disturbance segments
  for (i in c(1,3)) {
    segments(x0 = dates[auxDF[index+1,i]],
             y0 = DBEST_object$Fit[auxDF[index+1,i]],
             x1 = dates[auxDF[index+1,i+1]],
             y1 = DBEST_object$Fit[auxDF[index+1,i+1]],
             col = 'yellow',
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
             y0 = -0.19,
             x1 = dates[auxDF[index+1,4]],
             y1 = -0.19,
             col = 'red',
             lty = 1,
             lwd = 4.5
    )
  
  # adding text to the graph for the disturbance duration
  text(x = dates[round((mean(c(auxDF[index+1,2],auxDF[index+1,4]))),digits = 0)],
       y = -0.145,
       labels = paste(as.character(metricDF[index+1,5]),'Year Recovery'),
       adj = 0.5
       )
  
  # recovery min / max indicators on axis
  points(x = c(as.Date(as.yearmon(dates[1]) -0.95, frac = 1),
               as.Date(as.yearmon(dates[1]) -0.99, frac = 1)
               ),
         y = c(metricDF[3,2],
               metricDF[3,3]
               ),
         col = c('orange','green'),
         pch = '-',
         cex = 4
  )
  
  
  # adding a legend
  legend(x = dates[(length(dates)-8)],
         y = 0.2,
         legend = c('Original Time-series','DBEST Generalization','Disturbance Event','Recovery Min','Recovery Max'),
         col = c('black','blue','yellow','orange','green'),
         lty = c(2,1,1,NA,NA),
         pch = c(NA,NA,NA,'-','-'),
         pt.cex = c(NA,NA,NA,4,4),
         lwd = c(1,2,2,NA,NA),
         bg = 'light grey'
  )
  
  # closing the graphical device
  dev.off()
}

## A function to calculate recovery metrics, will be looped over the output
#   objects from DBEST, outputs to be passed to plotting function
#   also outputs index values for the segments I'm selecting to an auxiliary dataframe
calcRecoveryMetrics = function(DBEST_object, index) {
  
  changeSeg_VIdif = DBEST_object$Change
  
  # Need to find the greatest magnitude decline in the time series first
  distIndex = which.min(changeSeg_VIdif) # index value for breakpoint arrays for disturbance
  distSegEndDate = dates[DBEST_object$End[distIndex]] # date object of end of disturbance segment
  distSegStartDate = dates[DBEST_object$Start[distIndex]] # date object of start of disturbance segment
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
  
  # This statement allows for an additional segment up to 2 years directly after the first gain segment to be identified as the official
  #     end of the disturbance, many time series have their full recovery period split into two segments
  if (length(gainSegIndices) > 1) {
    maxDateForGain2 = dates[DBEST_object$End[gainIndex]] %m+% years(2)
    gainIndex2 = gainSegIndices[which((gainSegStartDates >= dates[DBEST_object$End[gainIndex]]) & (gainSegStartDates <= maxDateForGain2))]
    if (length(gainIndex2) == 0){ 
      gainIndex2 = gainIndex # for time series which do not have multiple gains after loss
    } else if (length(gainIndex2) > 1) {
      # for time series which have two gain segments within 2 years of the end of the disturbance, 
      #   takes the gain segment with the newest start date
      gainIndex2 = gainIndex2[which.max(dates[DBEST_object$Start[gainIndex2]])] 
    }
    gainSegEndDate = dates[DBEST_object$End[gainIndex2]]
    gainSegEndIndex = DBEST_object$End[gainIndex2]
  }
  
  DistStartY = as.double(substr(as.character.Date(distSegStartDate), start = 0, stop = 4))
  
  distEndYindex = DBEST_object$End[distIndex] # Recovery min
  recovMin = as.double(DBEST_object$Fit[distEndYindex])

  recovMax = as.double(DBEST_object$Fit[gainSegEndIndex]) # recovery max
  
  recoveryMag = as.double(recovMax - recovMin) # Recvoery magnitude
  
  Y2R = as.double(gainSegEndIndex - distEndYindex) # Years to recovery
  
  recovSlope = as.double(recoveryMag / Y2R) # Recovery slope
  
  gapInGainLoss = as.double(DBEST_object$Start[gainIndex] - DBEST_object$End[distIndex]) # Gap between the end of the disturbance segment and the start of the gain segment
  
  ## Adding recovery metrics to an array to be returned from the function
  rMetricsArray = c(DistStartY,recovMin,recovMax,recoveryMag,Y2R,recovSlope)
  
  ## Adding auxiliary data to dataframe
  auxSegArray = c(DBEST_object$Start[distIndex],DBEST_object$End[distIndex],DBEST_object$Start[gainIndex],gainSegEndIndex,gapInGainLoss)
  
  ## Adding a statement to remove data from bad time series
  # i.e time series with no gain, no loss, or only one segment
  if ((length(changeSeg_VIdif) <= 1) |
      (length(changeSeg_VIdif[changeSeg_VIdif < 0]) == 0) |
      (length(changeSeg_VIdif[changeSeg_VIdif > 0]) == 0) |
      (length(gainIndex) == 0)
      ) {
    rMetricsArray = rep(NA,8)
    auxSegArray = rep(NA,5)
  } else {
    preDistMean = as.double(mean(DBEST_object$Data[0:DBEST_object$Start[distIndex]])) # Pre-disturbance mean VI value
    PostDistMean = as.double(mean(DBEST_object$Data[DBEST_object$End[distIndex]:length(dates)])) # Post-disturbance mean VI value
    rMetricsArray = append(rMetricsArray,c(preDistMean,PostDistMean)) # These have to be calculated seperately, since their indexing data
  }                                                                   # if a time series doesn't have a gain or loss it'll crash without this 
  
  return(list(rMetricsArray,auxSegArray))

}

# implimenting DBEST, plotting results, and storing recovery metrics

for (i in TS_DF$X) {
  
  ## creating the time series object
  NDVIvals1 = as.numeric(as.vector(TS_DF[(i+1),2:length(TS_DF)]))
  
  ts1 = zoo(x = NDVIvals1,order.by = dates)
  
  output = DBEST(data = ts1, # data = generalization$Fit,
                  data.type = 'non-cyclical',
                  algorithm = "change detection",
                  change.magnitude = 0.25,
                  first.level.shift = 0.25,
                  second.level.shift = 0.15,
                  duration = 3,
                  distance.threshold = 0.15,
                  alpha = 0.05,
                  plot = 'no'
  )
  
  # Calculate recovery metrics and auxilary data
  listOfArrays = calcRecoveryMetrics(output,i)
  rMetrics[(i+1),] = listOfArrays[[1]]
  auxSegData[(i+1),] = listOfArrays[[2]]
  
  # Plotting
  plot_DBEST(output,i,rMetrics,auxSegData)
  
  cat("Plot #",i,'\n',
      "Year of disturbance :",rMetrics$'YearOfCut'[i],'\n',
      'direction / magnitude :',output$Change,'\n',
      'Start Years :', format(dates[output$Start], '%Y'),'\n',
      "durations :", output$Duration,'\n',
      '(0) = no abrupt, (1) = abrupt :', output$ChangeType,'\n',
      sep = ' '
      )
}

### Post Segmentation evaluation of the time series data / metrics ###
goodDataY2R = rMetrics$Y2R[is.na(rMetrics$Y2R) == FALSE]

areNAindicies = which(is.na(rMetrics$Y2R) == TRUE)

bigGapData = auxSegData$gainLossGap[is.na(auxSegData$gainLossGap) == FALSE]
bigGapData = bigGapData[bigGapData > 3]

totalNumberOfBadDataPoints = length(bigGapData) + length(areNAindicies)
cat("Number of bad time series :",totalNumberOfBadDataPoints)

write.csv(rMetrics,file = "RecoveryMetrics_output.csv")
