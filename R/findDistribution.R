#'findDistribution
#'
#'Finds the charge, rough distribution, then filters this distribution to find candidate isotope peaks
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>

#'@title prepareForOverview
#'@description Rearranges getDistributionOutput into IsotopeOverview format.
#'@param getDistributionOutput Output from a 12C or 13C run of getDistribution
#'@param significantPeaksTable The significantPeaks table of interest from LCMSunit (single table). For 13C, use IsotopeOverview. 
#'@return isotopeOverview
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
prepareForOverview <- function(getDistributionOutput ,significantPeaksTable) {
  allDistributions <- lapply(getDistributionOutput, function(X) X$filteredDistribution)
  nothingFoundIndex <- which(is.na(allDistributions))
  if (length(nothingFoundIndex) == 0){#If there were no NAs. Everything had a distribution.
    significantPeaksFound <- significantPeaksTable
    relevantGetDistributionOutput <- getDistributionOutput
    foundDistributions <- allDistributions
  } else {
    significantPeaksFound <- significantPeaksTable[-(nothingFoundIndex),]
    relevantGetDistributionOutput <- getDistributionOutput[-(nothingFoundIndex)]
    foundDistributions <- allDistributions[-(nothingFoundIndex)]
  }
  overviewData <- NULL
  #From here on, there should only be peaks that have at least one candidate isotope. 
  overviewData$Charge <- sapply(relevantGetDistributionOutput, function(X) X$Charge[[1]])
  overviewData$halogenCheck <- sapply(relevantGetDistributionOutput, function(X) X$halogenCheck[[1]])
  overviewData$startMass <- round(sapply(foundDistributions, function(X) X[1,'mz']), digits = 4)
  overviewData$endMass <- round(sapply(foundDistributions, function(X) X[nrow(X),'mz']), digits = 4)
  overviewData$numIsotopes <- sapply(foundDistributions, function(X) nrow(X))
  overviewData <- as.data.frame(overviewData, stringsAsFactors = FALSE)
  overviewData$distribution <- (lapply(foundDistributions, function(X) (round(as.numeric(X[,'relativeIntensity']), digits = 3))))
  overviewData$isotopeMasses <- (lapply(foundDistributions, function(X) (round(as.numeric(X[,'mz']), digits = 4))))
  isotopeOverview <- cbind(significantPeaksFound, overviewData)
  isotopeOverview$maxMZ <- isotopeOverview$endMass #Update maxMZ value (yes, this is redundant. Will eventually be removed.)
  return(isotopeOverview)
}

#'@title getDistribution 
#'@description The top-level function that finds isotope distributions for a single file within a LCMSunit. Works for both labelled and unlabelled files.
#'@details The user function. Identifies most likely isotope distributions for M/Zs of interest, assuming that they correspond to real compounds.
#'@param tableOfCandidates The table of interest. A single SignificantPeaks table from LCMSunit
#'@param xcmsRaw The xcmsRaw output for the file of interest.
#'@param precision When looking for a particular mz, how close must we be? Default 0.02
#'@return bestCandidates List of lists. Each top-level list corresponds to candidate from tableOfCandidates, each bottom level list features information about the corresponding isotope distribution search. Final infered distribution is $filteredDistribution.
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
getDistribution <- function(tableOfCandidates, xcmsRaw, precision = 0.02){
  simpleSP <- tableOfCandidates[,1:3]
  bestCandidates <- getCandidateInfo(simpleSP, xcmsRaw)
  bestCandidates <- getCharge(bestCandidates, precision) #Find charge (from single to quadruple)
  bestCandidates <- getRoughDistribution(bestCandidates, precision)
  bestCandidates <- refineDistribution(bestCandidates)
  return(refineDistribution(bestCandidates))
}

#'@title getCandidateInfo
#'@description Generates list containing relevant information for each peak in significantPeaks table.
#'@details Creates the list format with $bestMZ, $bestScan, etc. Essentially a rearranging of significantPeaks to the desired format, with the relevant scan added. 
#'@param simpleSP Data-frame. First three columns of significantPeaks table. Will also work with the full table.
#'@param xcmsRaw The xcmsRaw output for the 12C file of interest.
#'@return bestCandidates List of lists. Each element in top-level corresponds to a row in significantPeaks. Each bottom-level element is one of $bestMass, $bestRT, or $bestScan. 
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
getCandidateInfo <- function(simpleSP, xcmsRaw){
  numCandidates <- nrow(simpleSP)
  bestCandidates <- list()
  for (i in (1:numCandidates)){#Look at each peak in SP individually
    bestCandidates[[i]] <- list() #List associated with a single row in SignificantPeaks
    bestCandidates[[i]]$bestMass <- simpleSP$intensMZ[i]
    bestCandidates[[i]]$bestRT <- simpleSP$intensRT[i]
    massRange <- c(bestCandidates[[i]]$bestMass - 10, bestCandidates[[i]]$bestMass + 25)
    bestCandidates[[i]]$bestScan <- getScan(xcmsRaw, which.min(abs((slot(xcmsRaw, "scantime") - bestCandidates[[i]]$bestRT))), massRange)                                        
  }
  return(bestCandidates)
}

#'@title getCharge
#'@description Assuming that candidates masses correspond to real compounds, this function identifies the most probable charge of each. Searches up to quadruple charge.
#'@param bestCandidates Output from either getCandidateInfo or findMatchingPeaks, depending on whether a 12C or 13C file.
#'@param precision As set in getDistribution, how close to expected M/Z must a peak be to be considered a possible isotope?
#'@return bestCandidates Adds $Charge to each list within the list. Will be either NULL (no possible isotopes found), 1, 2, 3, or 4.
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
getCharge <- function(bestCandidates, precision){
  numCandidates <- length(bestCandidates)
  for (i in 1:numCandidates){
    singleData <- bestCandidates[[i]]
    if (length(singleData) <= 1){ 
      next
    } #If something strange has happened
    else{
      startMass <- singleData$bestMass
      RT <- singleData$bestRT
      scan <- singleData$bestScan
      
      startIndex <- which(scan[,'mz'] == startMass)
      charges <- list(1 * 1.00335, 0.5 * 1.00335, 0.33 * 1.00335, 0.25 * 1.00335)
      chargeCandidates <- lapply(X = charges
                                 , FUN = function(X) scan[which((0 < (scan[,'mz'] - startMass)) 
                                                                & (abs(scan[,'mz'] - startMass) <= (X + precision)) 
                                                                & (abs(scan[,'mz'] - startMass) >= (X - precision))),])
      chargeCandidates <- lapply(X = chargeCandidates, FUN = function(X) if (length(X) > 2) {X <- X[which.max(X[,'intensity']),]} else{X <- X}) #for situation where multiple hits are found for a single charge, take highest intensity hit.
      names(chargeCandidates) <- c(1,2,3,4)
      chargeCandidatesTable <- do.call(rbind, chargeCandidates)
      
      if(length(chargeCandidatesTable) != 0){ #Safety against situation where nothing was found at any charge
        #Now determine charge based on the candidate with highest intensity.
        chargeFinalIndex <- which.max(chargeCandidatesTable[,'intensity'])
        if (length(chargeFinalIndex) > 1) {
          chargeFinalIndex <- chargeFinalIndex[1]
        }
        
        chargeFinal <- rownames(chargeCandidatesTable)[chargeFinalIndex]
        bestCandidates[[i]]$Charge <- as.numeric(chargeFinal)
      }
      
    }
  }
  return(bestCandidates)
}

#'@title getRoughDistribution
#'@description Using the original starting mass, charge, and precision, looks for the isotope peaks at the next 10 predicted M/Zs.
#'@details Note that even if isotopes are not found at some masses, subsequant isotope masses will still be examined.
#'@param bestCandidates Output from getCharge
#'@param precision See above
#'@return bestCandidates Adds $roughIsotopes to bestCandidates. Is a matrix with mz/intensity for identified isotope candidates. Unfound masses will show NA/NA. 
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
getRoughDistribution <- function(bestCandidates, precision){
  numCandidates <- length(bestCandidates)
  for (i in (1:numCandidates)){
    singleData <- bestCandidates[[i]]
    if (length(singleData$Charge) == 0){ #If I did not identify a charge, there is no distribution to search for.
      bestCandidates[[i]]$roughIsotopes <- NA
    }
    else{
      startMass <- singleData$bestMass
      difference <- 1/(as.numeric(singleData$Charge))
      scan <- singleData$bestScan

      startMassIndex <- which(scan[,'mz'] == startMass)
      isotopesToSearchFor <- seq(startMass, startMass + (10*(difference*1.00335)), by = difference*1.00335) #STRATEGY: Look for up to 10 isotopes individually rather than sequentially (might find third isotope is a hit, but second is mysteriously missing). 1.00335 is 13C mass - 12C mass 
      isotopeCandidates <- lapply(isotopesToSearchFor, function(X) scan[which(abs(scan[,'mz'] - X) <= precision),])
      isotopeCandidates <- lapply(isotopeCandidates, function(X) if(length(X) > 2){X <- X[which.max(X[,'intensity']),]} else{X <- X}) #If multiple candidates found (ridiculous mass precision/sampling rate on machine), take highest intensity peak.
      isotopeCandidates <- lapply(isotopeCandidates, function(X) if(length(X) == 0){X <- c(0,0)} else{X <- X}) #Put place-fillers for the unfound masses
      names(isotopeCandidates) <- isotopesToSearchFor
      isotopeCandidates <- do.call(rbind, isotopeCandidates)
      bestCandidates[[i]]$roughIsotopes <- isotopeCandidates
    }
  }
  return(bestCandidates)
}

#'@title refineDistribution
#'@description Takes the 'roughDistributions' and applies several filters to them. Also searches for halogens.
#'@param bestCandidates Output from getRoughDistribution
#'@return bestCandidates Adds $filteredDistribution, the inferred 'real' isotope distribution for the candidate M/Z. Also adds $halogenCheck, is '1' if distribution indicates that there may be a halogen in the compound, '0' if not. 
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
refineDistribution <- function(bestCandidates){
  numCandidates <- length(bestCandidates)
  distributions <- lapply(bestCandidates, function(X) X$roughIsotopes)
  for (i in (1:numCandidates)){#In this loop, apply filters to individual distributions
    singleDistribution <- distributions[[i]]
    if (is.na(singleDistribution[1])) { #Case where no roughIsotopes found
      bestCandidates[[i]]$filteredDistribution <- NA
      next}
    filteredDistribution <- filter.UntilNoIsotopesFound(singleDistribution) #Basic filter, only includes isotopes until interruption.
    filteredDistribution <- filter.SufficientChangeInRelativeIntensity(filteredDistribution, requiredIntensityChange = 0.10)
    filteredDistribution <- filter.SufficientChangeInAbsoluteIntensity(filteredDistribution, requiredIntensityChange = 50) #Tries to avoid scenario where large relative intensity changes are a consequence of sub-100 absolute intensities
    filteredDistribution <- filter.GreaterThanOneRelativeIntensity(filteredDistribution)
    filteredDistribution <- filter.FallingIntensities(filteredDistribution)
    if ((ncol(as.data.frame(filteredDistribution))) == 1) { #If the startMass is the only peak remaining after filters (used this funny looking condition because the matrix transposes when there's only one row)
      bestCandidates[[i]]$filteredDistribution <- NA
    }
    else{
      bestCandidates[[i]]$halogenCheck <- chlorine.bromineCheck(filteredDistribution)
      bestCandidates[[i]]$filteredDistribution <- filteredDistribution
    }
  }
  return(bestCandidates)
}

#'@title filter.UntilNoIsotopesFound
#'@description Simple filter used in refineDistribution. Only stores the peaks that were found immediately, without interruption, to the right of the startMass. 
#'@details (there may be some cases where this filter is not ideal, but these will be very rare. I have yet to see any examples.)
#'@param singleDistribution A single distribution from bestCandidates that has begun to filtered.
#'@return filteredDistribution Input that may have had some isotope candidates removed. 
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
filter.UntilNoIsotopesFound <- function(singleDistribution){
  indexOfUnfoundIsotopes <- which(singleDistribution[,'mz'] == 0)
  if(length(indexOfUnfoundIsotopes) == 0){
    filteredDistribution.1 <- singleDistribution
    return(filteredDistribution.1)
  }
  else{
    endIsotopeMassIndex <- indexOfUnfoundIsotopes[[1]] - 1
    filteredDistribution.1 <- singleDistribution[1:endIsotopeMassIndex,]
    return(filteredDistribution.1)
  }
}



#'@title filter.SufficientChangeInRelativeIntensity
#'@description Medium intensity filter used in refineDistribution. Meant to identify the point at which relative intensities stagnate (prevents overly long distributions)
#'@details Strategy: 'moving window', three isotopes wide, checks whether there is a sufficient change in subsequent intensities.
#'@param filteredDistribution Output from the function that adds relative intensities to the distribution. 
#'@param requiredIntensityChange How much average movement must there be in relativeIntensity inside a three peak window before it is considered 'stagnate'? (Default = 0.10, peaks should move by more than 10% of their neighbours, on average)
#'@return filteredDistribution Same format, further refined isotope candidates.
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
filter.SufficientChangeInRelativeIntensity <- function(filteredDistribution
                                                       ,requiredIntensityChange = 0.10){
  lengthDistribution <- nrow(filteredDistribution)
  windowSize <- 3
  
  if(length(lengthDistribution) == 0) {lengthDistribution <- 1} #Matrix gets adjusted strangely when it only has one row.
  if (lengthDistribution == 1){ #Case where all isotopes have been filtered out.
    filteredDistribution.2 <- filteredDistribution
    return(filteredDistribution.2)
  }
  filteredDistribution <- cbind(filteredDistribution, relativeIntensity = (filteredDistribution[,'intensity']/filteredDistribution[1,'intensity']))  #Add relative intensity column
  relativeIntensities <- filteredDistribution[,'relativeIntensity'] 
  if (lengthDistribution == 2){ #Case with only one isotope peak identified.
    
    if (abs(1 - relativeIntensities[2]/relativeIntensities[1]) > requiredIntensityChange){
      filteredDistribution.2 <- filteredDistribution
    }
    else {
      filteredDistribution.2 <- filteredDistribution[1,]
    }
    return(filteredDistribution.2)
  }
  else{
    
    windowStarts <- as.matrix(seq(1, lengthDistribution- windowSize + 1, by = 1))
    windowSearchIndexes <- as.matrix(seq(1, windowSize-1, by = 1))
    requiredChangeCheck <- apply(windowStarts,1, function(X) ((mean(apply(windowSearchIndexes,1, function(Y) abs(1 - relativeIntensities[X+Y-1]/relativeIntensities[X+Y])))) > requiredIntensityChange)) #If FALSE, indicates that there was not an average requiredChange inside the window starting at that element in windowStarts
    firstFailedWindowIndex <- which(requiredChangeCheck == FALSE)[1]
    if (is.na(firstFailedWindowIndex)){ #If nothing failed
      filteredDistribution.2 <- filteredDistribution
      return(filteredDistribution.2)
    }
    else{ #If something failed
      #       browser()
      filteredDistribution.2 <- filteredDistribution[1:(windowStarts[firstFailedWindowIndex] + (windowSize - 2)),] #Only eliminate the final isotope in the failed window
      return(filteredDistribution.2)
    }
  }
}


#'@title filterSufficientChangeInAbsoluteIntensity
#'@description Medium intensity filter used in refineDistribution. Intended for scenario where relativeIntensity appears to be changing a lot simply because absolute intensities are very, very small (sub-100 generally)
#'@details As soon as absolute intensity stops fluctuating the right amount, stops. This should not be a problem in most cases, since at high intensities absolute intensity still fluctuates more than, say, 50. 
#'@param filteredDistribution Output from filter.SufficientChangeInRelativeIntensity
#'@param requiredIntensityChange How much must intensity change from one peak to the next? Default = 50
#'@return filteredDistribution Further filtered individual distribution.
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
filter.SufficientChangeInAbsoluteIntensity <- function(filteredDistribution
                                                       ,requiredIntensityChange = 50){
  lengthDistribution <- nrow(filteredDistribution)
  if(length(lengthDistribution) == 0) {lengthDistribution <- 1} #Matrix gets adjusted strangely when it only has one row.
  if (lengthDistribution == 1){ #Case where all isotopes have been filtered out.
    filteredDistribution.3 <- filteredDistribution
    return(filteredDistribution.3)
  }
  else{
    absoluteIntensities <- filteredDistribution[,'intensity']
    sufficientChangeCheck <- apply(as.matrix(seq(2,lengthDistribution)), 1, function(X) (abs(absoluteIntensities[X] - absoluteIntensities[X-1]) > requiredIntensityChange)) #If FALSE, intensity did not change enough
    firstFailedIndex <- which(sufficientChangeCheck == FALSE)[1]
    if (is.na(firstFailedIndex)){ #If nothing failed
      filteredDistribution.3 <- filteredDistribution
      return(filteredDistribution.3)
    }
    else{ #If something failed
      filteredDistribution.3 <- filteredDistribution[1:(firstFailedIndex),]
      return(filteredDistribution.3)
    }
  }
}

#'@title filter.GreaterThanOneRelativeIntensity
#'@description 12C filter used in refineDistribution. Makes sure that no peaks have relative intensity > 1, unless there is evidence that there is a chlorine/bromine causing m+2, m+4 peaks to increase.
#'@details Logic: If m+2 is >1, that's evidence of bromine/chlorine. If m+4 is >1, m+2 must also be >1.
#'@param filteredDistribution Single distribution
#'@return filteredDistribution.4 A futher filtered distribution.
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
filter.GreaterThanOneRelativeIntensity <- function(filteredDistribution){
  lengthDistribution <- nrow(filteredDistribution)
  if(length(lengthDistribution) == 0) {lengthDistribution <- 1} #Matrix gets adjusted strangely when it only has one row.
  if (lengthDistribution == 1){ #Case where all isotopes have been filtered out.
    filteredDistribution.4 <- filteredDistribution
    return(filteredDistribution.4)
  }
  else{
    relativeIntensities <- filteredDistribution[,'relativeIntensity']
    failedPeaksIndex <- which(relativeIntensities > 1) #The test
    #     browser()
    if (length(failedPeaksIndex) == 0){#If nothing failed
      filteredDistribution.4 <- filteredDistribution
      return(filteredDistribution.4)
    }
    else{ #If there are peaks >1
      singleFailedPeak <- failedPeaksIndex[1] #Start with the earliest failed peak
      if(singleFailedPeak == 3){ #3 is m+2 peak. At this stage in analysis, this is forgiveable as a chlorine/bromine
        singleFailedPeak <- failedPeaksIndex[2]
        if (is.na(singleFailedPeak)){#If there was no second failure to find
          filteredDistribution.4 <- filteredDistribution
          return(filteredDistribution.4)
        }
        else if (singleFailedPeak == 5){#5 is m+4 peak. Should only be >1 if m+2 is also.}
          singleFailedPeak <- failedPeaksIndex[3]
        }
      }
      if (is.na(singleFailedPeak)){#If nothing failed that wasn't a potential chlorine/bromine
        filteredDistribution.4 <- filteredDistribution
        return(filteredDistribution.4)
      }
      else{#If there is a legitimate failure
        filteredDistribution.4 <- filteredDistribution[1:(singleFailedPeak-1),]
        return(filteredDistribution.4)
      }
    }
  }
}

#'@title filter.FallingIntensities
#'@description filter used in refineDistribution. Checks that intensities are falling as expected. Exceptions made for chlorine/bromine candidates (m+2, m+4)
#'@param filteredDistribution A single distribution.
#'@return filteredDistribution Distribution passed through filter.
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
filter.FallingIntensities <- function(filteredDistribution){
  lengthDistribution <- nrow(filteredDistribution)
  if(length(lengthDistribution) == 0) {lengthDistribution <- 1} #Matrix gets adjusted strangely when it only has one row.
  if (lengthDistribution == 1){ #Case where all isotopes have been filtered out.
    filteredDistribution.5 <- filteredDistribution
    return(filteredDistribution.5)
  }
  else{ #If a distribution exists
    intensities <- filteredDistribution[,'intensity']
    differences <- apply(as.matrix(seq(1,(lengthDistribution-1))), 2, function(X) (intensities[X] - intensities[X+1])) #Element i is how much greater row i is than i+1
    failedDifferenceIndex <- which(differences < 0) #The test. Did intensity increase?
    #     browser()
    if (length(failedDifferenceIndex) == 0){#If nothing failed
      filteredDistribution.5 <- filteredDistribution
      return(filteredDistribution.5)
    }
    else{ #If there are non-declining peaks
      singleFailedDifference <- failedDifferenceIndex[1] #Start with the earliest failed peak
      if(singleFailedDifference == 2){ #2 is m+2 peak. At this stage in analysis, this is forgiveable as a chlorine/bromine
        singleFailedDifference <- failedDifferenceIndex[2]
        if (is.na(singleFailedDifference)){#If there was not a second failure to find
          filteredDistribution.5 <- filteredDistribution
          return(filteredDistribution.5)
        }
        else if (singleFailedDifference == 4){#4 is m+4 peak. Should only be increased intensity if m+2 is also.}
          singleFailedDifference <- failedDifferenceIndex[3]
        }
      }
      if (is.na(singleFailedDifference)){#If nothing failed that wasn't a potential chlorine/bromine
        filteredDistribution.5 <- filteredDistribution
        return(filteredDistribution.5)
      }
      else{#If there is a legitimate failure
        filteredDistribution.5 <- filteredDistribution[1:(singleFailedDifference),] #Remove starting at the failed isotope
        return(filteredDistribution.5)
      }
    }
  }
}

#'@title chlorine.bromineCheck
#'@description A very rough check of m+2 peak for evidence of bromine/chlorine. Used in refineDistribution.
#'@details Logic: Might be halogen if m+2 has relIntensity >1 OR if third peak is atleast 90% of the second.
#'@param filteredDistribution A single, fully filtered distribution
#'@return 1 if candidate found, 0 if not. 
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>
chlorine.bromineCheck <- function(filteredDistribution){
  #Everything that enters here should be a distribution with atleast 2 peaks
  lengthDistribution <- nrow(filteredDistribution)
  if (lengthDistribution == 2){
    return(0)
  }
  else{
    relativeIntensities <- filteredDistribution[,'relativeIntensity']
    mPlusTwo <- relativeIntensities[3]
    if (mPlusTwo > 1){ #Very conservative. If this doesn't work well, explore exact ratios as I do below with 0.33
      return(1)
    }
    else if ((mPlusTwo/relativeIntensities[2]) >= (.90)){ #Test. If third peak is greater than or within 10% of second peak. 10% Arbitrarily chosen.
      return(1)
    }
    else(
      return(0))
  }
}