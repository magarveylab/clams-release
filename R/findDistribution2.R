#####Let's try a 'scanblend' approach to finding distributions

###Input looks like:
# highestIntensity intensMZ intensRT    minMZ    maxMZ   minRT   maxRT RTfraction NumScanRange
# 1         1307.743 209.0598  311.649 209.0589 209.0598 309.050 313.263       1.18           14
# 2         1130.500 237.1208  310.668 237.1208 237.1215 310.359 311.649       1.16            5
# 3         4344.721 211.1442  320.624 211.1434 211.1443 316.176 327.979       0.99           31

#'@title getDistribution2
#'@description Different internal logic than findDistribution.  Focus on 'universal' logic that isn't hardcoded, and the ability to look left without headaches
#'@param xc xcmsRaw object of interest
#'@param sigPeaksTable The entire table from getSignificantPeaks
#'@param charges Vector of the charges to be considered. Ex/ c(1,2,3,4)
#'@param mzRadius Number of theoretical masses to generate on either side of intensMZ
#'@param weakPrecision Maximum distance between expected mz values and observed ones to be considered the same. This can be loose, and will be refined later in the algorithm.
#'@param nearbyThreshold When comparing relative distributions, how much percentage drift between scans is allowed while still being considered 'nearby'
#'@param precisionMultiplier When generating a 'strongPrecision', what multiple of the (mz difference + stdev) should we use?
#'@param minPrecision What is the smallest 'strongPrecision' that is allowed?
#'@param attendancePercentage How frequently must a peak be observed to 'replace' the leftmost peak as intensMZ? 
#'@param consistencyPercentage When multiple charges are plausible, what 'attendance record' is necessary in numOccurances and nearbyCount to be counted as a 'score'
getDistribution2 <- function(sigPeaksTable, xc, charges = c(1,2,3,4), mzRadius = 5,weakPrecision = 0.02,nearbyThreshold = 0.10, precisionMultiplier=2, minPrecision=0.001, attendancePercentage = 0.75, consistencyPercentage = 0.5){
  numPeaks <- nrow(sigPeaksTable)
  out <- c()
  for (i in 1:numPeaks){ #for each candidate
    #if (i%%100 == 0){print(i)}
    distributionData <- getSingleDistribution(xc, sigPeaksTable, i, charges = charges, mzRadius = mzRadius,weakPrecision = weakPrecision, nearbyThreshold = nearbyThreshold, precisionMultiplier = precisionMultiplier, minPrecision = minPrecision, attendancePercentage = attendancePercentage,consistencyPercentage = consistencyPercentage)
    out <- rbind(out,cbind(sigPeaksTable[i,1:2], oldIntensMZ = NA,sigPeaksTable[i,-c(1,2)], charge = NA, distribution = NA, distributionMasses = NA))
    if (!is.na(distributionData[[1]])){ #If a distribution was found
      out[i, "intensMZ"] <- as.numeric(unlist(distributionData[,"intensMZ"]))
      out[i, "oldIntensMZ"] <- as.numeric(unlist(sigPeaksTable[i,"intensMZ"]))
      out[i, "charge"] <- as.numeric(unlist(distributionData[,"charge"]))
      out[i, "distribution"][[1]] <- distributionData[,"distribution"]
      out[i, "distributionMasses"][[1]] <- distributionData[,"distributionMasses"]
    }
    else{  #Nothing found, but let's still fill the old CLAMS data
      out[i, "oldIntensMZ"] <- as.numeric(unlist(sigPeaksTable[i,"intensMZ"]))
      }
  }
  
  return(out)
}

#'@title getSingleDistribution
#'@description Finds the strongest distribution for a single row in getSigPeaks table.
#'@param xc xcmsRaw object of interest
#'@param sigPeaksTable The entire table from getSignificantPeaks
#'@param currentRowNumber The row of sigPeaksTable that we're currently interested in.
#'@param charges Vector of the charges to be considered. Ex/ c(1,2,3,4)
#'@param mzRadius Number of theoretical masses to generate on either side of intensMZ
#'@param weakPrecision Maximum distance between expected mz values and observed ones to be considered the same. This can be loose, and will be refined later in the algorithm.
#'@param nearbyThreshold When comparing relative distributions, how much percentage drift between scans is allowed while still being considered 'nearby'
#'@param precisionMultiplier When generating a 'strongPrecision', what multiple of the (mz difference + stdev) should we use?
#'@param minPrecision What is the smallest 'strongPrecision' that is allowed?
#'@param attendancePercentage How frequently must a peak be observed to 'replace' the leftmost peak as intensMZ? 
#'@param consistencyPercentage When multiple charges are plausible, what 'attendance record' is necessary in numOccurances and nearbyCount to be counted as a 'score'
getSingleDistribution <- function(xc, sigPeaksTable, currentRowNumber, charges = c(1,2,3,4), mzRadius = 5,weakPrecision = 0.02,nearbyThreshold = 0.10, precisionMultiplier=2, minPrecision=0.001, attendancePercentage = 0.75, consistencyPercentage = 0.5){
  a <- generatePredictedMasses(sigPeaksTable[currentRowNumber,"intensMZ"], mzRadius,charges)
  firstScan <- which.min(abs(xc@scantime - sigPeaksTable[currentRowNumber,"minRT"]))
  lastScan <- which.min(abs(xc@scantime - sigPeaksTable[currentRowNumber,"maxRT"]))
  intensScan <- which.min(abs(xc@scantime - sigPeaksTable[currentRowNumber,"intensRT"]))
  b <- getScansInWindow(xc, sigPeaksTable[currentRowNumber,"intensMZ"],  firstScan, lastScan)
  d <- getMatchingMasses(a,b,precision = weakPrecision)
  d2 <- getMassesInIntensScan(d, intensScan)
  e <- getSurvivingMasses(d2, mzRadius, nearbyThreshold)
  g <- recalibrateMZprecision(e, mzRadius,precisionMultiplier,minPrecision)
  h <- getStuffThatMightBeDistributions(g,mzRadius)
  j <- getTrimmedStuff(h, attendancePercentage)
  out <- reportDistribution(j,consistencyPercentage, scanMasses = d, intensScan)
  return(out)
}

#'@title getADecision
#'@description Chooses which distribution is strongest when multiple candidates exist. Scores based on 'consistency' in occurance and relative itensity.
#'@param trimmedStuff Output from getTrimmedStuff
#'@param consistencyPercentage When multiple charges are plausible, what 'attendance record' is necessary in numOccurances and nearbyCount to be counted as a 'score'
getADecision <- function(trimmedStuff, consistencyPercentage){
output <- sapply(trimmedStuff, function(RATEME){ #iterate over all charges
  if (is.na(RATEME[[1]])){#No distribution, score 0
    0}
  else{
    attendance <- as.numeric(RATEME[,"numOccurances"])
    neighbourhood <- as.numeric(RATEME[,"nearbyCount"])
    maxAttendance <- max(attendance)
    maxNeighbour <- max(neighbourhood)
    attendanceScore <- sum(attendance >= (consistencyPercentage*maxAttendance))
    neighbourhoodScore <- sum(neighbourhood >= (consistencyPercentage*maxNeighbour))
    score <- attendanceScore + neighbourhoodScore  ###The important part
    score
  }
  })
  winner <- which.max(output)
  if (length(winner) > 1){ #If we have a tie(rare)
    tiebreaker <- sapply(trimmedStuff[winner], nrow)
    tiebreaker <- which.max(tiebreaker)
    if (length(tiebreaker) > 1){ #If we're somehow STILL tied
      tiebreaker <- tiebreaker[1] #just take the first one I guess
      }
    winner <- winner[tiebreaker]
  }
  return(winner)
}

#'@title getReportFormat
#'@description Once a charge/distribution has been selected, this function converts the data to appropriate format for LCMSunit
#'@param distributionTable  A single element table from trimmedStuff
#'@param charge The charge.
#'@param scanMasses The exact masses observed in scans, used to get intensMZ 'right'. Part of gross fix and we are ashamed.
#'@param intensScan The original most intense scan. Also part of my embarassment.
getReportFormat <- function(distributionTable, charge, scanMasses, intensScan){
  if(is.na(distributionTable[1])) {
    return(NA)
  }
  masses <- distributionTable[,"averageMZ"]
  firstMZ <- masses[1]
  distribution <- distributionTable[, "averageProportion"]
  
  ##Gross fix so that reported intensMZ will be the exact value in intensScan
  initialIntensMZIndex <- distributionTable[1, "predictedMZ"] 
  realMasses <- scanMasses[[charge]][[initialIntensMZIndex]]
  intensMZ <- realMasses[which(realMasses[,"scan"] == intensScan), "mz"]
  
  
  
  return(cbind(intensMZ = intensMZ, distribution = list(as.numeric(as.vector(unlist(distribution)))), distributionMasses = list(as.numeric(as.vector(unlist(masses)))), charge = charge))
}

#'@title reportDistribution
#'@description Chooses the distribution to report back to LCMSunit. Chooses when multiple charges have plausible distributions
#'@param trimmedStuff Output from getTrimmedStuff
#'@param consistencyPercentage When multiple charges are plausible, what 'attendance record' is necessary in numOccurances and nearbyCount to be counted as a 'score'
#'@param scanMasses The exact masses observed in scans, used to get intensMZ 'right'. Part of gross fix and we are ashamed.
#'@param intensScan The original most intense scan. Also part of my embarassment.
reportDistribution <- function(trimmedStuff, consistencyPercentage = 0.5, scanMasses, intensScan){
  charges <- names(trimmedStuff)
  numCharges <- length(charges)
  empties <- unlist(which(is.na(trimmedStuff))) #unfound distributions
  numEmpties <- length(empties)
  if (numEmpties == numCharges){ #If there are no distributions at all
    return(NA)
  }
  else if (numEmpties == (numCharges - 1)){ #If there is only one distribution
    if (numEmpties == 0){#For case where only one charge is considered, and it has a distribution
      return(getReportFormat(trimmedStuff[[1]], charges, scanMasses, intensScan))
    }
    else{ #Normal scenario, multiple charges
      return(getReportFormat(trimmedStuff[-empties][[1]], charges[-empties], scanMasses, intensScan))
      
    }
  }
  else{ #If there's actually a decision to be made
    winnerIndex <- getADecision(trimmedStuff,0.5)
    return(getReportFormat(trimmedStuff[[winnerIndex]], charges[winnerIndex], scanMasses, intensScan))
  }
}

#'@title getTrimmedStuff
#'@description The leftmost peak is the most important, so removes it if it didn't show up enough relative to original intensMZ
#'@param stuffThatMight Output from getStuffThatMightBeADistribution
#'@param attendancePercentage How many times must the leftmost peak appeared relative to original intensMZ? 
getTrimmedStuff <- function(stuffThatMight, attendancePercentage = 0.75){
  output <- lapply(stuffThatMight, function(STUFFTHATMARGE){
    if (!is.na(STUFFTHATMARGE[[1]])){#If we have data
      attendance <- as.numeric(STUFFTHATMARGE[,"numOccurances"])
      perfectAttendance <- max(attendance, na.rm=TRUE)
      while(is.na(attendance[1]) | attendance[1] < (attendancePercentage*perfectAttendance)){ #While we have a first peak that should be trimmed
        attendance <- attendance[-1]
        STUFFTHATMARGE <- STUFFTHATMARGE[-1,]
      }
      if (is.null(nrow(STUFFTHATMARGE))){ #If we threw out everything but 1
        STUFFTHATMARGE <- NA
      }
    }
    STUFFTHATMARGE  
  })
  return(output)
}

#'@title getStuffThatMightBeDistributions
#'@description  Determines which charges should be further considered, based on whether adjacent isotopes are found (to the right OR to the left!)
#'@details Missing masses terminate a distribution search in all cases.
#'@param recalibratedSurvivingMasses Output from recalibrateMZTolerance
#'@param radius The number of masses that we checked on either side (used as coordinate to find original intensMZ)
getStuffThatMightBeDistributions <- function(recalibratedSurvivingMasses, radius){
  pivotIndex <- radius + 1 #index for original intensMZ
  out <- lapply(recalibratedSurvivingMasses, function(SURVIVORSFORACHARGE){ #Iterate over each charge
    unfoundIndex <- which(is.na(SURVIVORSFORACHARGE[,"averageMZ"]))
    leftStartIndex <- unfoundIndex[which((pivotIndex - unfoundIndex) > 0)]  #NA that is closest to pivotIndex without going over
    leftStartIndex <- ifelse(length(leftStartIndex) > 0, max(leftStartIndex) + 1, 1)
    rightEndIndex <- unfoundIndex[which((pivotIndex - unfoundIndex) < 0)] #NA that is closest to pivotIndex without going under
    rightEndIndex <- ifelse(length(rightEndIndex) > 0, min(rightEndIndex) - 1, pivotIndex + radius)
    if (leftStartIndex == rightEndIndex){ #If there are no possible adjacent isotopes
      NA
    }
    else{ #If we found some good stuff
      SURVIVORSFORACHARGE[leftStartIndex:rightEndIndex,]
      }
    })
}

#'@title recalibrateMZprecision
#'@description Based on the observed variance in the original intensMZ, we can narrow our mzPrecision tolerance, and remove other peaks.
#'@param survivingMasses Output from getSurvivingMasses
#'@param radius The number of masses that we checked on either side (used as coordinate to find original intensMZ)
#'@param precisionMultiplier The new precision will be this multiple of the intensMZ mass deviation (plus difference between averageMZ and predictedMZ).
#'@param minPrecision What's the lowest that precision can go? To avoid missing legitimate stuff when intensMZ matches perfectly.
recalibrateMZprecision <- function(survivingMasses, radius, precisionMultiplier, minPrecision = 0.001){
  out <- lapply(survivingMasses, function(WILLISURVIVETHISCHARGE){ #Iterate over each charge
      targetIndex <- radius + 1
      targetDeviation <- as.numeric(WILLISURVIVETHISCHARGE[targetIndex, "mzDeviation"])
      predicteds <- as.vector(as.numeric(WILLISURVIVETHISCHARGE[,"predictedMZ"]))
      averages <- as.vector(as.numeric(WILLISURVIVETHISCHARGE[,"averageMZ"]))
      newPrecision <- precisionMultiplier*(abs(predicteds[targetIndex] - averages[targetIndex]) + targetDeviation)   ####The important part
      newPrecision <- max(minPrecision, newPrecision)
      differences <- abs(predicteds - averages)
      test <- (differences <= newPrecision)
      WILLISURVIVETHISCHARGE[!test,] <- NA
      WILLISURVIVETHISCHARGE
  })
  return(out)
}

#'@title compareRelativeIntensities
#'@description Compares relative intensities to the original intensMZ, returns several measures of consistency over time. Important for determining if the peak is related.
#'@param queryIntensities Vector of the intensities for the mass currently being investigated (NAs removed)
#'@param referenceIntensities Vector for the intensities of the original intensMZ, with appropriate values removed to align with query Intensities
#'@param nearbyThreshold If the average proportion is 0.32, how close to 0.32 is considered 'about the same'? Default 0.05.
compareRelativeIntensities <- function(queryIntensities, referenceIntensities, nearbyThreshold){
  relativeIntensities <- round(queryIntensities/referenceIntensities, digits = 3)
  avrg <- round(mean(relativeIntensities), digits = 3)
  inHood <- length(which(abs(relativeIntensities - avrg) <= nearbyThreshold))
  return(c(averageProportion = avrg, stdev = round(sd(relativeIntensities), digits = 3), nearbyCount = inHood, lowestProportion = min(relativeIntensities), highestProportion = max(relativeIntensities)    ))
}

#'@title getSurvivingMasses
#'@description Returns masses that were found at least once in getMatchingMasses, along with some basic frequency data, and intensity consistency data
#'@param matchingMasses Output from getMatchingMasses
#'@param radius The number of masses that we checked on either side (used as coordinate to find original intensMZ)
#'@param nearbyThreshold When monitoring relative intensity changes, if the average proportion is 0.32, how close to 0.32 is considered 'about the same'? Default 0.05.
getSurvivingMasses <- function(matchingMasses, radius, nearbyThreshold = 0.05){
  targetIntensities <- as.vector(as.numeric(matchingMasses[[1]][[radius+1]][,"intensity"]))
  missingTargets <- which(is.na(targetIntensities))  #In some weird cases the mass won't have appeared in all scans
  missedCheck <- (length(missingTargets) > 0) #If something was missing, remove the missing stuff.
  if (missedCheck == TRUE){targetIntensities <- targetIntensities[-missingTargets]}
  out <- lapply(matchingMasses, function(MATCHESFORACHARGE) { #iterate over each charge
    singleResult <- lapply(MATCHESFORACHARGE, function(TABLE){ #Iterate over every predicted mass's results
      if (missedCheck == TRUE){  #remove the scans where the intensMZ didn't show up (for whatever reason). Should be used very rarely.
        TABLE <- matrix(TABLE[-missingTargets,], ncol = 3)  #ncol helps prevent one row matrix issue (which should be crazy rare)
        }
      matchIndex <- which(!is.na(TABLE[,1]))
      if (length(matchIndex) > 0){ #If anything was found
        matches <- matrix(TABLE[matchIndex,], ncol  = 3)  #To deal with stupid one-match matrix issue
        matches <- apply(matches, c(1,2), as.numeric)
        relatives <- compareRelativeIntensities(matches[,2], targetIntensities[matchIndex], nearbyThreshold)
        c(averageMZ = round(mean(matches[,1]), digits = 4), mzDeviation = round(sd(matches[,1]), digits = 4),averageIntensity = round(mean(matches[,2])) ,lowIntensity = min(matches[,2]), highIntensity = max(matches[,2]), highIntensityScan =  matches[which.max(matches[,2]),3],numOccurances = length(matchIndex), firstScan = min(matches[,3]), lastScan = max(matches[,3]), relatives)
      }
      else{ #Nothing found
        c(averageMZ = NA, mzDeviation = NA, averageIntensity = NA, lowIntensity = NA, highIntensity = NA,highIntensityScan = NA, numOccurances = NA, firstScan = NA, lastScan = NA, averageProportion = NA, stdev = NA, nearbyCount = NA, lowestProportion = NA, highestProportion = NA)
      }
    })
    singleResult <- do.call(rbind, singleResult)
    singleResult <- cbind(predictedMZ = names(MATCHESFORACHARGE), singleResult)
  })
  # out <- do.call(rbind, out)
  return(out)
}

#'@title getMassesInIntensScan
#'@description Removes any masses that did not appear in the most intense scan for intensMZ
#'@param matchingMasses Output from getMatchingMasses
#'@param intensScan Scan number of intensRT
getMassesInIntensScan <- function(matchingMasses, intensScan){
  out <- lapply(matchingMasses, function(MATCHESFORACHARGE){ #Iterate over every charge
    singleResult <- lapply(MATCHESFORACHARGE, function(TABLE){ #Iterate over every mass
      keyIndex <- which(TABLE[,"scan"] == intensScan)
      if (is.na(TABLE[keyIndex, "mz"])){ #If we didn't find the mass in the key scan
        TABLE[,"mz"] <- NA  #remove this mass from contention
        TABLE[,"intensity"] <- NA
      }
      TABLE
      })
    singleResult
    })
  return(out)
}

#'@title getMatchingMasses
#'@description Locates any observed peaks that match to predicted masses, for a single query mass from sigPeaks
#'@param predictedMasses Output from generatePredictedMasses
#'@param scanList Output from getScansInWindow
#'@param precision How close to the predicted mass must an observed mass be?
getMatchingMasses <- function(predictedMasses, scanList, precision = 0.02){
  scanNums <- names(scanList)
  out <- lapply(predictedMasses, function(MASSESFORACHARGE) { #Iterate over each charge
    singleChargeResults <- lapply(MASSESFORACHARGE,  function(PREDICTEDMASS) { #Iterate over all predicted masses, each gets it's own list element
      matches <- sapply(scanList, function(SCAN){ #Iterate over every scan
        test <- which(abs(SCAN[,"mz"] - PREDICTEDMASS) <= precision)
        if (length(test) == 0){ #If there was no match
          c(mz = NA, intensity = NA)
        }
        else { #If there were multiple matches, take most intense
          SCAN[test[which.max(SCAN[test,2])],]
        }
      })
      cbind(round(t(matches), digits = 4), scan = scanNums)
    })
    names(singleChargeResults) <- round(MASSESFORACHARGE, digits = 4)
    singleChargeResults
  })
  return(out)
}

#'@title getScansInWindow
#'@description Retrieves the scans indicated as relevant by sigPeaks. Returns as list of tables.
#'@param xcmsRaw The object.
#'@param intensMZ Single mass from sigPeaks
#'@param startScan First scan in window.
#'@param endScan Last scan in window
#'@param scanMzRadius How much of the scan to retrieve
getScansInWindow <- function(xcmsRaw, intensMZ, startScan, endScan, scanMzRadius = 15){
  out <- lapply(seq(startScan, endScan), function(SCAN){
    partialOut <- getScan(xcmsRaw, SCAN, c(intensMZ - scanMzRadius, intensMZ + scanMzRadius))
    partialOut
  })
  names(out) <- seq(startScan,endScan)
  return(out)
}

#'@title generatePredictedMasses
#'@description  Returns list of vectors with all masses that would be predicted based on charge and starting mass.  Looks left and right
#'@param intensMZ  Single mass from sigPeaks table
#'@param radius Number of masses to generate on either side of intensMZ
#'@param charges Which charges to consider i.e. 1,2,3,4,...
#'@return List of vectors of all possible mz values, in ascending order, each list element named by its charge.
generatePredictedMasses <- function(intensMZ, radius, charges = c(1)){
  out <- lapply(charges,function(CHARGE) {
    difference <- round(1.00335/CHARGE, digits = 5)
    partialOut <- seq(intensMZ - (difference*radius), intensMZ + (difference*radius), by = difference)  
    partialOut
  })
  names(out) <- charges
  return(out)
}