mz.file <- "/mnt/storage/mass_spec_samples/data/2017_9/AMK_17-09-15_15613_1006_5_PA-2117.mzml"
minMZ <- 150
maxMZ <- 2000
rtMin <- 50
MZtolerance <- 0.005
RTtolerance <- 30
backgroundIntensity <- 1000


lc <- LCMSunit(mz.file)
lc <- LCMSunitAnalyze(lc, c(minMZ, maxMZ), rtMin , MZtolerance , RTtolerance , backgroundIntensity)
lc <- LCMSunitPDAData(lc)
