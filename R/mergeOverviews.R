#'@title mergeOverviews
#'@param Isotopeoverview done with the original findDistribution method
#'@param Isotopeoverview done with the new findDistribution2 method
#'@param correctedMZCutoff correctedMZ value used to merge the 'High' and 'Low' peak refinement methods
#'@return merged dataframe.
mergeOverviews <- function(isoLow, isoHigh, correctedMZCutoff = 1500) {
  dropCols <- c("halogenCheck","startMass","endMass","numIsotopes","oldIntensMZ")
  isoLow <- isoLow[ , !(names(isoLow) %in% dropCols)]
  colnames(isoLow)[10] <- "charge"
  colnames(isoLow)[12] <- "distributionMasses"
  isoLow <- isoLow[isoLow[ ,"intensMZ"] * isoLow[ ,"charge"] <= correctedMZCutoff, ]
  isoHigh <- isoHigh[ , !(names(isoHigh) %in% dropCols)]
  isoHigh <- isoHigh[isoHigh[ ,"intensMZ"] * isoHigh[ ,"charge"] > correctedMZCutoff, ]
  isoHigh <- na.omit(isoHigh)
  return(rbind(isoLow, isoHigh))
}
