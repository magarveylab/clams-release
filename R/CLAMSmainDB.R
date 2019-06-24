#'@author cDejong
#'@title CLAMSmainDB
#'@description runs CLAMS as a peak finder and outputs the settings, TIC and peaks as a JSON string
#'@param mz.file the location of the mzML file
#'@param minMZ minMZ
#'@param maxMZ maxMZ
#'@param rtMin rtMin
#'@param MZtolerance is the mass value that the different peaks can differ by
#'@param RTtolerance is the retention time value that the different peaks can differ by
#'@param backgroundIntensity minimum intensity
#'@param correctedMZCutoff correctedMZ value used to merge the 'High' and 'Low' peak refinement methods
#'@return return results as JSON as a string
CLAMSmainDB <- function(mz.file, minMZ = 100, maxMZ = 2500, rtMin = 300, MZtolerance = 0.005, RTtolerance = 60, backgroundIntensity = 2000, correctedMZCutoff = 1500){
  lc <- LCMSunit(mz.file)
  lc <- LCMSunitAnalyze(lc, c(minMZ, maxMZ), rtMin , MZtolerance , RTtolerance , backgroundIntensity, correctedMZCutoff)
  overview <- LCMSunitAddMS2Spectra(lc)$IsotopeOverview[[1]]
  overview <- overview[,c("intensMZ", "intensRT", "highestIntensity", "distribution", "distributionMasses", "charge", "maxRT", "minRT", "RTfraction", "ms2mz", "ms2int")]
  mz.tic <- list(RT = lc$xcmsRaw[[1]]@scantime, int = lc$xcmsRaw[[1]]@tic)
  settings <- list(mzML = mz.file, mz_min = minMZ, mz_max = maxMZ, rt_min = rtMin, mz_tolerance = MZtolerance, rt_tolerance = RTtolerance, background_intensity = backgroundIntensity, corrected_mz_cutoff = correctedMZCutoff)
  json <- lcToJSON(overview, mz.tic, settings)
  return(json)
}

#'@title jsonToFile
#'@description writes a JSON text object to a file
#'@param json the json string to write
#'@param file.name the name/path of the file to write to
jsonToFile <- function(json, file.name){
  f <- file(file.name)
  writeLines(json, f)
  close(f)
}
