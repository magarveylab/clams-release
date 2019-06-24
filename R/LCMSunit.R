#Author: cDejong, mCannon
#Object class for LCMS analysis

#'@author Chris Dejong <dejonc [at] mcmaster.ca>
#'@author Michael Cannon <cannonmj [at] mcmaster.ca>

#'@title LCMSunit object
#'Creates the LCMSunit object
#'@param MZXMLs list of mzmls for input
#'@return unit
LCMSunit <- function(MZMLs){
  unit <- list("xcmsRaw" = sapply(MZMLs, xcmsRaw),
               "SignificantPeaks" = NA, "Unique" <- NA, "IsotopeOverview" = NA, "IsotopeDetails" = NA, "FeedingDetails" = NA, "Adducts" = NA, "MiscData" = NA)
  class(unit) <- "LCMSunit"
  return(unit)
}

print.LCMSunit <- function(x){
  cat("LCMSunit files: ")
  print(names(x$xcmsRaw))
}

#'@title LCMSunitMaxIntens
#'Finds the maximum intensity value in the entire LCMS run of the First input file
#'@param LCMSunit being analyzed
#'@return LCMSunit with Max intensity added in MiscData
LCMSunitMiscData <- function(LCMSunit){ # Can thread this
  #Gets the maximum intensity from a run
  maxIntensity <- sapply(LCMSunit$xcmsRaw,function(xcms) max(xcms@env$intensity))
  minIntensity <- sapply(LCMSunit$xcmsRaw,function(xcms) min(xcms@env$intensity))
  medianIntensity <- sapply(LCMSunit$xcmsRaw,function(xcms) median(xcms@env$intensity))
  #Gets the time difference between the scans for each run and averages them
  meanRTdiff <- sapply(LCMSunit$xcmsRaw, function(xcms){
    rts <- data.frame(scanIndex = xcms@scanindex, scanRT <- xcms@scantime)[,2]
    diff <- vector()
    for(i in 1:length(rts)-1){ diff <- c(diff, rts[i+1]-rts[i]) }
    mean(diff)
  })
  miscData <- data.frame(MaxIntensity = maxIntensity, MinIntensity = minIntensity, MedianIntensity = medianIntensity, MeanRTdiff = meanRTdiff)
  LCMSunit$MiscData <- miscData
  return(LCMSunit)
}

#'@title LCMSunitAnalyze
#'@param LCMSunit being analyzed through the entire CLAMs analysis
#'@param masses window to analyze
#'@param rtMin minimum retention time to consider in second
#'@param MZtolerance is the mass value that the different masses can be off by
#'@param RTtolerance in seconds that the RTs can be different by
#'@param backgroundIntensity the minimum intensity to consider
#'@param feeding If doing shift analysis, vector of 13C filepaths
#'@param correctedMZCutoff correctedMZ value used to merge the 'High' and 'Low' peak refinement methods
#'@return LCMSunit analyzed via the below methods
LCMSunitAnalyze <- function(LCMSunit,
                            masses = c(100,1500),
                            rtMin = 300,
                            MZtolerance = 0.005,
                            RTtolerance = 60,
                            backgroundIntensity = 1000,
                            correctedMZCutoff = 1500) {
  LCMSunit <- LCMSunitSignificantPeaks(LCMSunit, masses = masses, rtMin = rtMin, MZtolerance = MZtolerance, RTtolerance = RTtolerance, backgroundIntensity = backgroundIntensity)
  LCMSunit <- LCMSunitRefineIsotope(LCMSunit, correctedMZCutoff) #Adds IsotopeOverview and IsotopeDetails
  return(LCMSunit)
}

#'@title LCMSunitSignificantPeaks
#'Runs LCMSunit through getSignificantPeaks.R
#'See getSignificantPeaks.R for full description of the analysis
#'@param LCMSunit being analyzed through Ranges
#'@param masses window to analyze
#'@param rtMin minimum retention time to consider in seconds
#'@param MZtolerance is the mass value that the different masses can be off by
#'@param RTtolerance in seconds that the RTs can be different by
#'@param backgroundIntensity the minimum intensity
#'@return LCMSunit with ranges data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAkAAAAJCAYAAADgkQYQAAAAMElEQVR42mNgIAXY2Nj8x8cHC8AwMl9XVxe3QqwKcJmIVwFWhehW4LQSXQCnm3ABAHD6MDrmRgfrAAAAAElFTkSuQmCCadded
LCMSunitSignificantPeaks <- function(LCMSunit, masses = c(100,1500), rtMin = 300, MZtolerance = 0.005, RTtolerance = 60, backgroundIntensity = 1000){
  LCMSunit$SignificantPeaks  <- lapply(LCMSunit$xcmsRaw,function(x) getPeaksAboveCutoff(x, masses, rtMin, MZtolerance, RTtolerance, backgroundIntensity))
  return(LCMSunit)
}

#'@title LCMSunitRefineIsotope
#'Runs LCMSunit through findDistribution.R
#'Builds $IsotopeOverview, $IsotopeDetails
#'@param LCMSunit Following LCMSunitSignificantPeaks
#'@param correctedMZCutoff correctedMZ value used to merge the 'High' and 'Low' peak refinement methods
#'@return LCMSunit with IsotopeOverview and IsotopeDetails added.
LCMSunitRefineIsotope <- function(LCMSunit, correctedMZCutoff = 1500){
  numFiles <- length(LCMSunit$xcmsRaw)
  #new
  LCMSunit$IsotopeOverviewHighMass <- lapply(as.matrix(seq(1,numFiles)), function(X) getDistribution2(LCMSunit$SignificantPeaks[[X]], LCMSunit$xcmsRaw[[X]], charges = c(1,2,3,4,5,6,7), minPrecision = 0.005, precisionMultiplier = 3))
  names(LCMSunit$IsotopeOverviewHighMass) <- paste0(names(LCMSunit$SignificantPeaks))
  #old
  LCMSunit$IsotopeDetailsLowMass <- lapply(as.matrix(seq(1,numFiles)), function(X) getDistribution(LCMSunit$SignificantPeaks[[X]], LCMSunit$xcmsRaw[[X]]))                      
  LCMSunit$IsotopeOverviewLowMass <- lapply(as.matrix(seq(1,numFiles)), function(X) prepareForOverview(LCMSunit$IsotopeDetailsLowMass[[X]], LCMSunit$SignificantPeaks[[X]]))
  names(LCMSunit$IsotopeOverviewLowMass) <- paste0(names(LCMSunit$SignificantPeaks))
  LCMSunit <- LCMSunitMergeMethods(LCMSunit, correctedMZCutoff)
  return(LCMSunit)
}

#'@title LCMSunitMergeMethods
#'@param LCMSunit After running HighMass and LowMass
#'@param correctedMZCutoff correctedMZ value used to merge the 'High' and 'Low' peak refinement methods
#'@return LCMSunit with two methods merged.
LCMSunitMergeMethods <- function(LCMSunit, correctedMZCutoff = 1500) {
  numFiles <- length(LCMSunit$xcmsRaw)
  LCMSunit$IsotopeOverview <- lapply(as.matrix(seq(1,numFiles)), function(X) mergeOverviews(LCMSunit$IsotopeOverviewLowMass[[X]],
                                                                                            LCMSunit$IsotopeOverviewHighMass[[X]],
                                                                                            correctedMZCutoff))
  return(LCMSunit)
}

#'@title LCMSunitAddMS2Spectra
#'@param LCMSunit After running addMS2 
#'@return LCMSunit with MS2 data
LCMSunitAddMS2Spectra <- function(LCMSunit){ #Need to add more options here from the defaults
  LCMSunit$IsotopeOverview <- lapply(1:length(LCMSunit$xcmsRaw), function(x){
    xcms <- LCMSunit$xcmsRaw[[x]]
    overview <- LCMSunit$IsotopeOverview[[x]]
    addMS2(xcms, overview)
  })
  LCMSunit
}
