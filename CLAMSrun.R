#Runs CLAMS from commandline
suppressMessages(library(clams))

args <- commandArgs(trailingOnly = TRUE)
mz.file <- args[1]
minMZ <- as.numeric(args[2])
maxMZ <- as.numeric(args[3])
rtMin <- as.numeric(args[4])
MZtolerance <- as.numeric(args[5])
RTtolerance <- as.numeric(args[6])
backgroundIntensity <- args[7]


CLAMSmainDB(mz.file, minMZ, maxMZ, rtMin, MZtolerance, RTtolerance, backgroundIntensity)

