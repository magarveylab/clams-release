#'@title lcToJson
#'@author cDejong
#'@description Converts the below input into a json object for pasing
#'@param isotopeOverview is the data with the peaks
#'@param tic is the total ion chromatogram for the lc run
#'@param settings are the clams settings used to get the results, as a list of named strings
#'@return all the data in JSON format
lcToJSON <- function(isotopeOverview, tic, settings){

	lc.jsonParse <- list(settings = settings,
                      TIC =  tic,
                      data = lapply(seq(nrow(isotopeOverview)), function(x) list(MS1peak = isotopeOverview[x,]))
                      )
  clams_results <- list(clams_results = lc.jsonParse)
	toJSON(clams_results)
}
