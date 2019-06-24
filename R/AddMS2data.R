#'@author cDejong, mCannon
#'@title addMS2
#'@description Adds cooresponding MS2 scans to ms1 CLAMS output
#'@param ms1Peaks LCMSunit$IsotopeOvervew or LCMSunit$SignificantPeaks uses "intensMZ' and 'intensRT' to find the cooresponding ms2 scan
#'@param ms2.summary a summary of ms2scans including the scan number, "scanNum", "precursorRT", and "precursorMZ" are necessary columns.
#'@return ms1Peaks with a new column containing cooresponding MS2 scan number, ms2 mz values, and ms2 intensity values
addMS2 <- function(xcms, ms1Peaks){
  mz <- openMSfile(xcms@filepath[[1]])
  mz.header <- header(mz)
  ms1.header <- mz.header[(which(mz.header$msLevel == 1)),]
  ms2.header <- mz.header[(which(mz.header$msLevel == 2)),]
  if(nrow(ms2.header) < 1){
    ms2mz <- rep(NA, nrow(ms1Peaks))
    ms2int <- rep(NA, nrow(ms1Peaks))
    ms1Peaks <- cbind(ms1Peaks, ms2mz, ms2int)
    return(ms1Peaks)
  }
  ms1.df <- scansToMatrix(mz, ms1.header)
  ms2.df <- scansToList(mz, ms2.header)
  ms2.summary <- getMs2Summary(mz.header)
  ms2.summary <- refinePrecursor(ms2.summary,ms1.df)
  ms1Peaks <- linkMS2SpectraToPeaks(ms1Peaks, ms2.summary)
  ms1Peaks <- addMS2SpectraDetails(ms1Peaks, ms2.df)
  ms1Peaks
}

addMS2SpectraDetails <- function(ms1Peaks, ms2.summary, backgroundIntensity = 5, ionsToSave = 50){
  ms1Peaks.ms2Data <- data.frame()
  for(i in 1:nrow(ms1Peaks)){
    row <- ms1Peaks[i,]
    ms2ScanNum <- row[["ms2Scan"]]
    if(!is.na(ms2ScanNum)){
      ms2.spectra <- ms2.summary[[toString(ms2ScanNum)]]
      ms2.spectra <- ms2.spectra[,1:2]
      if(is.null(nrow(ms2.spectra))){
        row$ms2mz <- NA
        row$ms2int <- NA
      }else{
        ms2.spectra <- ms2.spectra[ms2.spectra[,"intensity"] > backgroundIntensity,]
        if(is.null(nrow(ms2.spectra))){
          row$ms2mz <- NA
          row$ms2int <- NA
        }else{
          ms2.spectra <- ms2.spectra[order(-ms2.spectra[,"intensity"]), ]
          if(nrow(ms2.spectra) > ionsToSave){
            ms2.spectra <- ms2.spectra[1:50,]
          }
          mz <- ms2.spectra[,"mz"]
          mz <- round(mz, 3)
          int <- ms2.spectra[,"intensity"]
          int <- round(int, 0)
          row$ms2mz <- list(mz)
          row$ms2int <- list(int)
        } 
      }
    }else{
      row$ms2mz <- NA
      row$ms2int <- NA
    }
    ms1Peaks.ms2Data <- rbind(ms1Peaks.ms2Data, row)
  }
  ms1Peaks.ms2Data
}

#'@author cDejong, mCannon
#'@title linkMS2SpectraToPeaks
#'@description Adds cooresponding MS2 scans to ms1 CLAMS output
#'@param ms1Peaks LCMSunit$IsotopeOvervew or LCMSunit$SignificantPeaks uses "intensMZ' and 'intensRT' to find the cooresponding ms2 scan
#'@param ms2.summary a summary of ms2scans including the scan number, "scanNum", "precursorRT", and "precursorMZ" are necessary columns.
#'@return peaks.ms2 ms1Peaks with a new column containing cooresponding MS2 scan number
linkMS2SpectraToPeaks <- function(ms1Peaks, ms2.summary, precision = 0.02){ #UNTESTED
  
  peaks.ms2 <- do.call(rbind, lapply(1:nrow(ms1Peaks), function(i){
    peak <- ms1Peaks[i,]
    mz <- peak["intensMZ"][[1]]
    rt <- peak["intensRT"][[1]]
    matchCheck <- which((abs(ms2.summary[,'precursorMZ'] - mz) < precision) & (abs(ms2.summary[,'retentionTime'] - rt) < 15))
    if (length(matchCheck) == 0){ #If the precursor does not exist in ms1Peaks
      ms2Scan <- NA
    }
    else{
      matchIndex <- which((abs(ms2.summary[,'precursorMZ'] - mz) < precision) & (abs(ms2.summary[,'retentionTime'] - rt) < 15)) #There may be multiple here
      matchBest <- which.min(abs(ms2.summary[matchIndex,'precursorMZ'] - mz))
      matchIndex <- matchIndex[matchBest]
      ms2Scan <- ms2.summary[matchIndex, "ms2scanNum"]
    }
    return(cbind(peak,ms2Scan))    
  }))
  return(peaks.ms2)
}

refinePrecursor <- function(ms2.summary, ms1.df){ 
  ms1.df.scanNum <- ms1.df[,'scanNum']
  refinedPrecursors <- do.call(rbind, lapply(1:nrow(ms2.summary), function(i){
    row <- ms2.summary[i,]
    precursor.scan <- row["precursorScan"]
    roughPrecursor.MZ <- row["precursorMZ"]
    singleScan <- ms1.df[ms1.df.scanNum == precursor.scan,]
    singleScan <- rbind(singleScan, c(0,0,0,0))
    exactPrecursor.MZ <- which.min(abs(singleScan[,'mz'] - roughPrecursor.MZ))
    row["precursorMZ"] <- singleScan[exactPrecursor.MZ,'mz']
    row["precursorIntensity"] <- singleScan[exactPrecursor.MZ, 'intensity']
    row
  }))
  return(refinedPrecursors)
}


getMs2Summary <- function(mz.header){
  scan.num <- NULL
  ms2.summary <- data.frame()
  for(i in 1: nrow(mz.header)){ #Change this to lapply
    single <- mz.header[i,]
    if(single["msLevel"] == 1){
      scan.num <- single["seqNum"]
    }else if(single["msLevel"] == 2){
      entry <- c(single["seqNum"], scan.num, single["precursorMZ"], single["precursorIntensity"], single["retentionTime"])
      ms2.summary <- rbind(ms2.summary, entry)
    }else{
      stop(paste0("BAD MS LEVEL DETECTED AT: ", single))
    }
  }
  colnames(ms2.summary) <- c("ms2scanNum", "precursorScan", "precursorMZ", "precursorIntensity", "retentionTime")
  ms2.summary <- as.matrix(ms2.summary)
  return(ms2.summary)
}

#'@title scansToList
#'@description takes a scan's header (summary) and builds a list of data.frame for all the scans in a list named by the scanNum (seqNum)
#'@param mz is the msR.object which contains all the data for the lcms run
#'@param mz.header is the header for all scans (a summary of scans)
#'@return a list of data.frame containing all scans with proper headers and the list is named by the scanNum (seqNum)
scansToList <- function(mz, mz.header){
  scans.all <- lapply(1:nrow(mz.header), function(x){ #Can thread this as needed
    buildSingleScanMatrix(mz, mz.header, x)
  })
  names(scans.all) <- mz.header$seqNum
  return(scans.all)
}




#'@title scansToMatrix
#'@description takes a scan's header (summary from mzR) and builds the matrix for all the scans in a single data.frame
#'@param mz is the msR.object which contains all the data for the lcms run
#'@param mz.header is the header for all scans (a summary of scans)
#'@return a data.frame containing all scans with proper headers
scansToMatrix <- function(mz, mz.header){
  scans.all <- do.call(rbind, lapply(1:nrow(mz.header), function(x){ #Can thread this as needed
    buildSingleScanMatrix(mz, mz.header, x)
  }))
  scans.all <- scans.all[scans.all[,'intensity'] >= 20, ]
  return(scans.all)
}

#'@title buildSingleScanMatrix
#'@description takes mz data from a single header and builds the matrix object with proper column names
#'@param mz is the msR.object which contains all the data for the lcms run
#'@param mz.header is the header for all scans (a summary of scans)
#'@param scanNum is the scan number/index that is being extracted
#'@return returns a dataframe of the scan that is desired with labelled column names
buildSingleScanMatrix <- function(mz, mz.header, scan.num){
  scan.info <- mz.header[scan.num,]
  scan <- as.matrix(peaks(mz, scan.info$seqNum))
  if(nrow(scan) == 0){
    scan <- rbind(scan, c(0,0)) #For empty scans, keep it from being empty
  }
  colnames(scan) <- c("mz", "intensity")
  #   scan[,'retentionTime'] <- scan.info[,'retentionTime']
  #   scan[,'scanNum'] <- scan.info[,'seqNum']
  scan <- cbind(scan, retentionTime = scan.info[,'retentionTime'], scanNum = scan.info[,'seqNum'])
  return(scan)
}