setwd("C:/R_workspace")
library(DBEST)
library(zoo)
TS_DF = read.csv("timeSeriesDF500.csv",header = TRUE)

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

# filling missing values at beginning of time-series with first value
for (arrayNum in TS_DF[,1]) {
  aArray = TS_DF[(arrayNum+1),]
  NonNAindex = which(!is.na(aArray))
  if (length(NonNAindex) != length(aArray)){
    aResult = 'this worked'
    firstNonNA = NonNAindex[2]
    fillVal = as.double(aArray[firstNonNA])
    for (i in seq(from = 2, to = (firstNonNA-1))) {
      aArray[i] = fillVal
    }
    TS_DF[(arrayNum+1),] = aArray
  }
}

# function to plot the DBEST output correctly
plot_DBEST = function(DBEST_object,index) {
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
  for (i in seq(from = 1, to = DBEST_object$BreakpointNo)) {
    segments(x0 = dates[(DBEST_object$Start[i])],
             y0 = DBEST_object$Fit[(DBEST_object$Start[i])],
             x1 = dates[(DBEST_object$End[i])],
             y1 = DBEST_object$Fit[(DBEST_object$End[i])],
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
    segments(x0 = dates[min(DBEST_object$End)],
             y0 = -0.18,
             x1 = dates[max(DBEST_object$End)],
             col = 'red',
             lty = 1,
             lwd = 4.5
    )
  # calculating disturbance duration
  start = min(DBEST_object$End)
  end = max(DBEST_object$End)
  duration = end - start
  
  # adding text to the graph for the disturbance duration
  text(x = dates[round((mean(c(start,end))),digits = 0)],
       y = -0.145,
       labels = paste(as.character(duration),'Year Recovery'),
       adj = 0.5
       )
  
  recovMin = DBEST_object$Fit[min(DBEST_object$End)]
  recovMax = DBEST_object$Fit[max(DBEST_object$End)]
  
  # recovery min / max indicators on axis
  points(x = c(as.Date(as.yearmon(dates[1]) -0.95, frac = 1),
               as.Date(as.yearmon(dates[1]) -0.99, frac = 1)
               ),
         y = c(DBEST_object$Fit[min(DBEST_object$End)],
               DBEST_object$Fit[(max(DBEST_object$End))]
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

# Dataframe to store recovery metrics for each stand
rMetrics <- data.frame(matrix(nrow = length(rownames(TS_DF)), ncol = 5,byrow = TRUE))
colnames(rMetrics) = c('YearOfCut','Recovery_min','Recovery_max','Recovery_slope','Y2R')
rownames(rMetrics) = TS_DF$X

# implimenting DBEST, plotting results, and storing recovery metrics

for (i in TS_DF$X) {
  
  ## creating the time series object
  NDVIvals1 = as.numeric(as.vector(TS_DF[(i+1),2:length(TS_DF)]))
  
  ts1 = zoo(x = NDVIvals1,order.by = dates)

# # implimenting DBEST change detection / generalization
#   generalization = DBEST(data = ts1,
#                           data.type = 'non-cyclical',
#                           algorithm = "generalization",
#                           change.magnitude = 0.25,
#                           first.level.shift = 0.3,
#                           second.level.shift = 0.15,
#                           duration = 5,
#                           distance.threshold = 'default',
#                           alpha = 0.05,
#                           plot = 'no'
# )
  
  output = DBEST(data = ts1, #generalization$Fit
                  data.type = 'non-cyclical',
                  algorithm = "change detection",
                  breakpoints.no = 2,
                  first.level.shift = 0.20,
                  second.level.shift = 0.15,
                  duration = 5,
                  distance.threshold = 0.15,
                  alpha = 0.05,
                  plot = 'no'
  )
  
  ## Calculating and storing recovery metrics
  distIndex = which.min(output$Change) # index value for breakpoint arrays for disturbance
  recovIndex = which.max(output$Change) # index value for breakpoint arrays for recovery
  
  DistStartYOBJ = dates[output$Start[distIndex] + 1] # year of disturbance, have to add one since values are from November it's unlikely the disturbance occurred within the last month of the year
  DistStartY =substr(as.character.Date(DistStartYOBJ), start = 0, stop = 4)
  
  distEndY = output$End[distIndex] # Recovery min
  recovMin = output$Fit[distEndY]
  
  recovEndY = output$End[recovIndex] # recovery max
  recovMax = output$Fit[recovEndY]
  
  Y2R = (recovEndY - distEndY) # Years to recovery
  
  recovSlope = output$Change[recovIndex] / Y2R # Recovery slope
  
  
  ## Adding recovery metrics to dataframe
  rMetrics[(i+1),1] = DistStartY
  rMetrics[(i+1),2] = recovMin
  rMetrics[(i+1),3] = recovMax
  rMetrics[(i+1),5] = Y2R
  rMetrics[(i+1),4] = recovSlope
  
  ## Plotting
  plot_DBEST(output,i)
  
  cat("Plot #",i,'\n',
      "Year of disturbance :",format(dates[distEndY],'%Y'),'\n',
      'direction / magnitude :',output$Change,'\n',
      "durations :", output$Duration,'\n',
      '(0) = no abrupt, (1) = abrupt :', output$ChangeType,'\n',
      sep = ' '
      )
}

write.csv(rMetrics,file = "RecoveryMetrics_output.csv")
