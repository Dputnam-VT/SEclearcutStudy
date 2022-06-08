setwd("C:/R_workspace")
library(DBEST)
library(zoo)
library(lubridate)
library(tidyverse)
library(RColorBrewer)

dataDate1 = "2022_05_26" # using for the latest version of the time series csv
dataDate2 = "2022_06_03" # using for the latest version of the stand attribute csv
dataDate3 = "2022_06_08" # using for the latest version of the climate variables csv

# Bringing in vegetation index time series
inputDF = read.csv(paste("timeSeriesDF500_",dataDate1,".csv",sep = ""),header = TRUE)
TS_DF = inputDF[order(inputDF$X),]

# Bringing in topographic, soil, locational, and ecotonal data
ecoRegionInfo = read.csv(file = paste("standAttributes_",dataDate2,".csv",sep = ""))

# Bringing in the climate data for a single scenario first
scenario = 'rcp45'
precipTS_DF = read.csv(file = paste("precipDF_",scenario,"_",dataDate3,".csv",sep = ""))
precipTS_DF = precipTS_DF[order(precipTS_DF$X),]
minTempTS_DF = read.csv(file = paste("minTempDF_",scenario,"_",dataDate3,".csv",sep = ""))
minTempTS_DF = minTempTS_DF[order(minTempTS_DF$X),]
maxTempTS_DF = read.csv(file = paste("maxTempDF_",scenario,"_",dataDate3,".csv",sep = ""))
maxTempTS_DF = maxTempTS_DF[order(maxTempTS_DF$X),]

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

# filling in missing values not calculated by the pandas interpolation 
#   (start and end values)
for (arrayNum in TS_DF[,1]) {
  aArray = TS_DF[(arrayNum+1),]
  NonNAindex = which(!is.na(aArray))
  if (length(NonNAindex) != length(aArray)){
    firstNonNA = NonNAindex[2]
    fillVal = as.double(aArray[firstNonNA])
    for (i in seq(from = 2, to = (firstNonNA))) {
      aArray[i] = fillVal
    }
    
  NonNAindex2 = which(!is.na(aArray)) # check to see if NA's remain (end values)
  if (length(NonNAindex2) != length(aArray)){
    lastNonNA = NonNAindex2[length(NonNAindex2)]
    fillVal = as.double(aArray[lastNonNA])
    for (i in seq(from = (lastNonNA+1), to = length(aArray))) {
      aArray[i] = fillVal
    }
  }
    TS_DF[(arrayNum+1),] = aArray
  }
}

# Dataframe to store recovery metrics for each stand
rMetrics <- data.frame(matrix(nrow = length(rownames(TS_DF)), ncol = 8,byrow = TRUE))
colnames(rMetrics) = c('YearOfCut','Recovery_min','Recovery_max','Recovery_mag','Y2R','Recovery_slope','pre_CutMean','post_CutMean')
rownames(rMetrics) = TS_DF$X

# Dataframe to store auxiliary data for plotting (start dates for segments, etc.)
auxSegData <- data.frame(matrix(nrow = length(rownames(TS_DF)), ncol = 6,byrow = TRUE))
colnames(auxSegData) = c('DistStartIndex','DistEndIndex','RecoveryStartIndex','RecoveryEndIndex','gainLossGap','fitRMSE')
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
         y = c(metricDF[index+1,2],
               metricDF[index+1,3]
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
  
  # # This statement allows for an additional segment up to 2 years directly after the first gain segment to be identified as the official
  # #     end of the disturbance, many time series have their full recovery period split into two segments
  # if (length(gainSegIndices) > 1) {
  #   maxDateForGain2 = dates[DBEST_object$End[gainIndex]] %m+% years(2)
  #   gainIndex2 = gainSegIndices[which((gainSegStartDates >= dates[DBEST_object$End[gainIndex]]) & (gainSegStartDates <= maxDateForGain2))]
  #   if (length(gainIndex2) == 0){ 
  #     gainIndex2 = gainIndex # for time series which do not have multiple gains after loss
  #   } else if (length(gainIndex2) > 1) {
  #     # for time series which have two gain segments within 2 years of the end of the disturbance, 
  #     #   takes the gain segment with the newest start date
  #     gainIndex2 = gainIndex2[which.max(dates[DBEST_object$Start[gainIndex2]])] 
  #   }
  #   gainSegEndDate = dates[DBEST_object$End[gainIndex2]]
  #   gainSegEndIndex = DBEST_object$End[gainIndex2]
  # }
  
  DistStartY = as.double(substr(as.character.Date(distSegStartDate), start = 0, stop = 4))
  
  distEndYindex = DBEST_object$End[distIndex] # Recovery min
  recovMin = as.double(DBEST_object$Fit[distEndYindex])

  recovMax = as.double(DBEST_object$Fit[gainSegEndIndex]) # recovery max
  
  recoveryMag = as.double(recovMax - recovMin) # Recvoery magnitude
  
  Y2R = as.double(gainSegEndIndex - distEndYindex) # Years to recovery
  
  recovSlope = as.double(recoveryMag / Y2R) # Recovery slope
  
  gapInGainLoss = as.double(DBEST_object$Start[gainIndex] - DBEST_object$End[distIndex]) # Gap between the end of the disturbance segment and the start of the gain segment
  
  RMSE = sqrt(mean((DBEST_object$Fit - DBEST_object$Data)^2))
  
  ## Adding recovery metrics to an array to be returned from the function
  rMetricsArray = c(DistStartY,recovMin,recovMax,recoveryMag,Y2R,recovSlope)
  
  ## Adding auxiliary data to dataframe
  auxSegArray = c(DBEST_object$Start[distIndex],DBEST_object$End[distIndex],DBEST_object$Start[gainIndex],gainSegEndIndex,gapInGainLoss,RMSE)
  
  ## Adding a statement to remove data from bad time series
  # i.e time series with no gain, no loss, or only one segment
  if ((length(changeSeg_VIdif) <= 1) |
      (length(changeSeg_VIdif[changeSeg_VIdif < 0]) == 0) |
      (length(changeSeg_VIdif[changeSeg_VIdif > 0]) == 0) |
      (length(gainIndex) == 0)
      ) {
    rMetricsArray = rep(NA,8)
    auxSegArray = rep(NA,6)
  } else {
    preDistMean = as.double(mean(DBEST_object$Data[0:DBEST_object$Start[distIndex]])) # Pre-disturbance mean VI value
    PostDistMean = as.double(mean(DBEST_object$Data[DBEST_object$End[distIndex]:length(dates)])) # Post-disturbance mean VI value
    rMetricsArray = append(rMetricsArray,c(preDistMean,PostDistMean)) # These have to be calculated seperately, since their indexing data
    }                                                                   # if a time series doesn't have a gain or loss it'll crash without this 
  
  return(list(rMetricsArray,auxSegArray))

}

### Experimental function, finding the confidence of the segmentation ###

DBEST_rep = function(inputZooTS, numIterations) {
  
  segStarts = matrix(nrow = numIterations, ncol = 3) # segments start index
  
  segEnds = matrix(nrow = numIterations, ncol = 3) # segments end index
  
  fits = matrix(nrow = numIterations, ncol = 38) # fit data values
  
  parameterMatrix = matrix(nrow = numIterations, ncol = 5) # Recording the parameter values
  
  for (rep in seq(1,numIterations)) {
    
    firstLS = runif(1,0.10,0.50)
    secondLS = runif(1,0.15,0.55)
    dur = round(runif(1,2,10),digits = 0)
    dist = runif(1,0.01,0.20)
    
    iterationOut = DBEST(data = ts1,
                         data.type = 'non-cyclical',
                         algorithm = "change detection",
                         breakpoints.no = 3,
                         first.level.shift = firstLS,
                         second.level.shift = secondLS,
                         duration = dur,
                         distance.threshold = dist,
                         alpha = 0.05,
                         plot = 'off'
    )
    
    RMSE = sqrt(mean((iterationOut$Fit - iterationOut$Data)^2))
    parameterMatrix[rep,] = c(firstLS,secondLS,dur,dist,RMSE)
    
    startsArray = iterationOut$Start
    endsArray = iterationOut$End
    fitsArray = iterationOut$Fit
    while (length(startsArray) < 3) {
      endsArray = append(endsArray,NA)
      startsArray = append(startsArray,NA)
    }
    segStarts[rep,] = sort(startsArray, na.last = TRUE)
    segEnds[rep,] = sort(endsArray, na.last = TRUE)
    fits[rep,] = fitsArray
  }
  return(list(segStarts,segEnds,fits,parameterMatrix))
}

# Also a function to calculate the mode #
# Create the function.
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# implimenting DBEST, plotting results, and storing recovery metrics

for (i in TS_DF$X) {
  
  ## creating the time series object
  NDVIvals1 = as.numeric(as.vector(TS_DF[(i+1),2:length(TS_DF)]))
  
  ts1 = zoo(x = NDVIvals1,order.by = dates)
  
  #repList = DBEST_rep(ts1,50) # repeatly trying out randomly selected parameter values for segmentation
  # 
  # #### test plotting to visualize the variations in segmentation ####
  # 
  # 
  # plot.default(repList[[3]][1,],type = 'l', ylim = c(-0.2,1))
  # for (n in seq(2,50)) {
  #   lines(repList[[3]][n,])
  # }
  # points(x = repList[[1]], y = coredata(ts1)[repList[[1]]],
  #        col = 'green', cex = 2.0)
  # points(x = repList[[2]], y = coredata(ts1)[repList[[2]]],
  #        col = 'red')
  # lines(coredata(ts1),col = 'blue')
  # 
  # # Trying out a mean annual VI value
  # meanValues = c()
  # for (year in seq(1,50)) {
  #   yearVals = repList[[3]][,year]
  #   meanVal = mean(yearVals)
  #   meanValues = append(meanValues,meanVal)
  # }
  # lines(x = seq(1,37),y = meanValues, col = 'red')
  # 
  # # trying out a majority vote for segment starts and ends
  # majorityStarts = c()
  # for (segN in seq(1,3)) {
  #   start = getmode(repList[[1]][,segN])
  #   majorityStarts[segN] = start
  # }
  # majorityEnds = c()
  # for (segN in seq(1,3)) {
  #   end = getmode(repList[[2]][,segN])
  #   majorityEnds[segN] = end
  # }
  # plot.default(coredata(ts1),col = 'blue', type = 'l')
  # # adding lines for the disturbance segments
  # for (v in c(1,2,3)) {
  #   segments(x0 = majorityStarts[v],
  #            y0 = coredata(ts1)[majorityStarts[v]],
  #            x1 = majorityEnds[v],
  #            y1 = coredata(ts1)[majorityEnds[v]],
  #            col = 'yellow',
  #            lty = 1,
  #            lwd = 3
  #   )
  # }
  # 
  # # trying out a matching approach which finds the total number of
  # #   unique start and end combinations
  # uniqueStarts = unique(repList[[1]])
  # uniqueEnds = unique(repList[[2]])

  output = DBEST(data = ts1,
                  data.type = 'non-cyclical',
                  algorithm = "change detection",
                  breakpoints.no = 4,
                  first.level.shift = 0.20,
                  second.level.shift = 0.30,
                  duration = 5,
                  distance.threshold = 0.15,
                  alpha = 0.05,
                  plot = 'off'
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

write.csv(rMetrics,file = "RecoveryMetrics_output.csv", )

### Post Segmentation evaluation of the time series data / metrics ###
# need to add the identifier as a column for tracking things during filtering
rMetrics$'UniqueID' = row.names(rMetrics)
rMetrics$'K_shift' = rMetrics$post_CutMean - rMetrics$pre_CutMean
auxSegData$'UniqueID' = row.names(rMetrics)

areNAindicies = which(is.na(rMetrics$Y2R) == TRUE)

bigGapIndices = which(auxSegData$gainLossGap > 3)
bigGapLengths = auxSegData$gainLossGap[bigGapIndices]

totalNumberOfBadDataPoints = length(bigGapIndices) + length(areNAindicies)
cat("Number of bad time series :",totalNumberOfBadDataPoints)

## Creating a new dataframe with only good values
goodDataDF = rMetrics[which(is.na(rMetrics$YearOfCut) == FALSE),]
goodDataDF$fit_RMSE = auxSegData$fitRMSE[which(is.na(auxSegData$DistStartIndex) == FALSE)]

auxSegData$'DistDuration' = auxSegData$DistEndIndex - auxSegData$DistStartIndex

for (ID in goodDataDF$UniqueID) {
  gapLength = auxSegData$gainLossGap[which(auxSegData$UniqueID == ID)]
  distDuration = as.numeric(auxSegData$DistDuration[which(auxSegData$UniqueID == ID)])
  recovDuration = goodDataDF$Y2R[goodDataDF$UniqueID == ID]
  if ((gapLength > 3) | (distDuration > 3) | (recovDuration <=2)) {
    indexValue = which(goodDataDF$UniqueID == ID)
    goodDataDF = goodDataDF[-indexValue,]
  }
} # removing stands with large gaps between disturbance end and recovery start


# merging the ecoregion info and good data for plotting #
ecoRegionInfo$intLat = floor(ecoRegionInfo$lat)
plottingDF = merge(x = goodDataDF,
                   y = ecoRegionInfo,
                   by.x = 'UniqueID',
                   by.y = 'X'
                   )

#### new data visualization plotting methods ####

# making a list of recovery metrics to plot
metrics = c('Recovery_min','Recovery_max','Recovery_mag','Y2R',
            'Recovery_slope','K_shift')

#---------------- Ungrouped Stats plots for each metric ----------------------#
# histograms
for (aMetric in metrics){
  summaryStats = ggplot(data = plottingDF, 
                         aes_string(x = aMetric))+ 
    geom_histogram(bins = 20)+
    #geom_vline(aes(xintercept = mean(aMetric)), color = 'red')
    labs(title = paste("Distribution of",aMetric, "Metric"))
    print(summaryStats)
}
# boxplots
# have to make a seperate DF for just the metrics for this one
metricsDF = plottingDF %>% select(all_of(metrics[-4]))
summaryStats = ggplot(data = stack(metricsDF), 
                      aes(x = ind, y = values)
                      )+ 
  geom_boxplot()+
  #geom_vline(aes(xintercept = mean(aMetric)), color = 'red')
  labs(title = paste("Distribution of Selected Recovery Metrics"))+
  xlab(NULL)+
  ylab("NBR Value")
print(summaryStats)

#---------------- Stratifying by disturbance year ----------------------------#
# Need to create a subset of the data that fall within the time window I want
distWindowDF = plottingDF[which(plottingDF$YearOfCut >= 1994 & plottingDF$YearOfCut <= 2010),]
distWindowDF$YearOfCut = as.factor(distWindowDF$YearOfCut)

for (aMetric in metrics){
  yearsPlot = ggplot(data = distWindowDF, 
                     aes_string(x = 'YearOfCut', y = aMetric,group = 1)
                     )+
    #geom_bin_2d(drop = FALSE)+
    geom_point()+
    stat_smooth(size = 1.5,
                method = 'loess',
                span = 0.80,
                show.legend = TRUE)+
    labs(title = paste("Variation in",aMetric,"by Disturbance Year")
         )+
    xlab("Year of Disturbance")
  print(yearsPlot)
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
plottingDF$ecoR_code = as.character(plottingDF$ecoR_code)
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
    scale_fill_brewer(palette = 'Set1')
  print(ecoRegionPlot)
}
