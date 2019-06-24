#'@author Chris Dejong <dejonc [at] mcmaster.ca>, RCPP update by Jennifer Cabral <cabraj3 [at] mcmaster.ca>
#'@title getPeaksAboveCutoff
#'@description gets all the peaks for the mass spec, each peak is for continous masses (1 apart) over a retention time range as one peak 
#'@param xc the mass spec data loaded as an xcms object
#'@param mz.range that is being checked for
#'@param rt.min minimum retention time to consider in seconds
#'@param mz.tolerance is the mass value that the different peaks can differ by
#'@param rt.tolerance is the retention time value that the different peaks can differ by
#'@param min.intensity the minimum intensity to be considered
#'@return significantPeakMasses
getPeaksAboveCutoff <- function(xc, mz.range = c(100, 1500), rt.min = 90, mz.tolerance = 0.005, rt.tolerance = 60, minimum.intensity = 1000, window.size = 100, above.average.multiplier = 5){
  intens = xc@env$intensity
  mzs = xc@env$mz
  ions.table <- cbind(xc@env$intensity, xc@env$mz, 1:length(xc@env$mz))
  colnames(ions.table) <- c("intensity", "mz", "ionIndex")
  #Remove those - and + window.size/2 from what is wanted
  ions.table <- ions.table[ions.table[,"mz"] >= mz.range[1] - (window.size/2), ]
  ions.table <- ions.table[ions.table[,"mz"] <= mz.range[2] + (window.size/2), ]
  #Add the rt
  scans <- cbind(xc@scanindex, xc@scantime) #Positional information from the intensities and mass for each scan, and the RT for the scans
  colnames(scans) <- c("scanIndex", "scanRT") #move this to after remove background
  
  ions.table <- addScansCPP(ions.table,scans)
  
  ions.table <- ions.table[order(ions.table[,"rt"], ions.table[,"mz"]),]
  ions.table <- ions.table[ions.table[,"rt"] >= rt.min,] #remove ions below the rt cutoff
  ions.table <- removeBackgroundCPP(ions.table, window.size, above.average.multiplier, minimum.intensity, mz.range) #Do sliding window inside the range
  #Get average time between scans
  mean.time.between.scans <- time.between.scans(scans)
  #change the grouping methods to take a matrix instead
  peaks.data.frame <- group.ions(ions.table, mz.tolerance, rt.tolerance, mean.time.between.scans)
  return(peaks.data.frame)
}

#'@title group.ions
#'@description groups ions into peaks
#'@param ions.table
#'@param mz.tolerance is the mass value that the different peaks can differ by to be grouped
#'@param rt.tolerance is the retention time value that the different peaks can differ by to be grouped in seconds
#'@param mean.time.between.scans average time between scans in the file in seconds
#'@return ions.data.frame
group.ions <- function(ions.table, mz.tolerance, rt.tolerance, mean.time.between.scans){
  ions.matrix <- groupIonsSameScanCPP(ions.table, mz.tolerance)
  ions.data.frame <- as.data.frame(ions.matrix)
  ions.data.frame <- group.peaks.across.scans(ions.data.frame, mz.tolerance, rt.tolerance, mean.time.between.scans)
  ions.matrix <- as.matrix(ions.data.frame)
  ions.matrix <- collectFinalPeaksCPP(ions.matrix, mz.tolerance, rt.tolerance)
  ions.data.frame <-data.frame(ions.matrix)
  names(ions.data.frame)[2] <- "intensMZ"
  names(ions.data.frame)[3] <- "intensRT"
  ions.data.frame <- ions.data.frame[ions.data.frame[,"NumScanRange"] > 1,]
  return(ions.data.frame)
}

#'@title group.peaks.across.scans
#'@description groups collected mz's together that overlap in mz and are close in rt
#'@param mz.intensities.mzCollected sorted by smallest to highest rt
#'@param mz.tolerance is the mass value that the different masses can be off by
#'@param rt.tolerance in seconds that the RTs can be different by
#'@param mean.time.between.scans is the mean of the time between each scan for an entire run
#'@return mz.intensities.rtCollected
group.peaks.across.scans <- function(mz.intensities.mzCollected, mz.tolerance = 0.005, rt.tolerance = 60, mean.time.between.scans){
  mz.intensities.rtCollected = data.frame()
  while(nrow(mz.intensities.mzCollected) > 0){
    row.query <- mz.intensities.mzCollected[1,]
    mz.intensities.mzCollected <- mz.intensities.mzCollected[-1,]
    
    rtRange <- c(as.numeric(row.query[3]),as.numeric(row.query[3]))
    mz.range <- as.numeric(row.query[4:5])
    intens.mz.rt <- as.numeric(row.query[1:3])
    rtCount <- 1
    reduced.mzCollected <- subset(mz.intensities.mzCollected, mz.intensities.mzCollected$minMZ > mz.range[1] - mz.tolerance * 5 
                                  & mz.intensities.mzCollected$minMZ < mz.range[2] + mz.tolerance * 5
                                  & mz.intensities.mzCollected$rt < rtRange[2] + rt.tolerance)
    if(!is.na(reduced.mzCollected[1,1])){
      for(j in 1:nrow(reduced.mzCollected)){
        row.next <- reduced.mzCollected[j,]
        intensity.next <- as.numeric(row.next[1])
        mzIntens.next <- as.numeric(row.next[2])
        rt.next <- as.numeric(row.next[3])
        mz.range.next <- as.numeric(row.next[4:5])
        
        if(sameMZ(mz.range,mz.range.next, mz.tolerance)){ #Check if any of the masses match, if match they are the same compound
          rtRange[2] <- rt.next
          if(mz.range.next[1] < mz.range[1]){ #Check if lower mz is lower than the query to extend
            mz.range[1] <- mz.range.next[1]
          }
          if(mz.range.next[2] > mz.range[2]){ #Check if higher mz is higher than the query to extend
            mz.range[2] <- mz.range.next[2]
          }
          if(intensity.next > intens.mz.rt[1]){ #Check if is more intense than query, to replace the most intense peak
            intens.mz.rt <- c(intensity.next,mzIntens.next,rt.next)
          }
          rowName.to.remove <- rownames(row.next)
          mz.intensities.mzCollected <- mz.intensities.mzCollected[!rownames(mz.intensities.mzCollected) %in% rowName.to.remove, ]
          rtCount <- rtCount + 1
        }
      } 
    }
    #Get the number of scans the ions appeared in as a fraction of the total over the range
    numScans <- ((rtRange[2] - rtRange[1])/mean.time.between.scans)+1
    rtFraction <- round(rtCount/numScans,2)
    #Add to the data.frame
    mz.intensities.rtCollected <- rbind(mz.intensities.rtCollected, c(intens.mz.rt, mz.range, rtRange, rtFraction, rtCount))
  }
  
  #Adds headers
  names(mz.intensities.rtCollected)[1] <- "highestIntensity"
  names(mz.intensities.rtCollected)[2] <- "mz"
  names(mz.intensities.rtCollected)[3] <- "rt"
  names(mz.intensities.rtCollected)[4] <- "minMZ"
  names(mz.intensities.rtCollected)[5] <- "maxMZ"
  names(mz.intensities.rtCollected)[6] <- "minRT"
  names(mz.intensities.rtCollected)[7] <- "maxRT"
  names(mz.intensities.rtCollected)[8] <- "RTfraction"
  names(mz.intensities.rtCollected)[9] <- "NumScanRange"
  return(mz.intensities.rtCollected)
}

#'@title sameMZ 
#'@param mz.range.1 query
#'@param mz.range.2 subject
#'@param mz.tolerance is the mass value that the different masses can be off by
#'@return overlap boolean
sameMZ <- function(mz.range.1,mz.range.2, mz.tolerance){
  overlap <- FALSE
  if((mz.range.1[1]-mz.tolerance) < mz.range.2[1] & (mz.range.1[2]+mz.tolerance) > mz.range.2[1]
     |(mz.range.1[1]-mz.tolerance) < mz.range.2[2] & (mz.range.1[2]+mz.tolerance) > mz.range.2[2]
     |(mz.range.1[1]+mz.tolerance) > mz.range.2[1] & (mz.range.1[2]-mz.tolerance) < mz.range.2[2]){ #Last check here is in case mz.range1 is inside mz.range2
    decimal1 <- revtrunc(mz.range.1[1])
    decimal2 <- revtrunc(mz.range.2[1])
    if(abs(decimal1 - decimal2) <= mz.tolerance*2){
      overlap <- TRUE 
    }
  }
  return(overlap)
}

#Get the time between scans given a list of scans and their retentiontime
time.between.scans <- function(scans){
  diff <- vector()
  for(i in 1:nrow(scans)-1){
    diff <- c(diff, scans[i+1,2]-scans[i,2])
  }
  return(mean(diff))
}

#'@title revtrunc extracts decimal number from float eg: 0.15 -> 15 and -0.15 -> 15
#'@param x float to have decimal extracted 
#'@return decimal values only
revtrunc <- function(x) {
  x <- abs(x)
  x - floor(x)
}
