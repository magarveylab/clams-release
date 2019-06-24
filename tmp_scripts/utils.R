#very simple TIC
tic <- cbind(mz.xc@scantime, mz.xc@tic)
colnames(tic) <- c("rt", "int")
plot(tic, type = "o")

#to minutes
tic[,1] <- tic[,1] / 60
plot(tic, type = "o")

########plotting scans
plotScan(mz.xc, 1300) #based on scan number
plotScan.simp(mz.xc, 570, c(870,880)) #based on time in seconds

#get scan data
scan <- getScan(mz.xc, 100)
scan <- getScan.simp(mz.xc, 100)

#A convienience method for plotting a scan based on time rather than index
plotScan.simp <- function(mz.xcms, time, ...) {
  scanIndex = which(abs(mz.xcms@scantime - time) == min(abs(mz.xcms@scantime - time)))
  print(scanIndex)
  plotScan(mz.xcms, scanIndex, ...)
}

#A convienience method for getting a scan's values based on time rather than index
getScan.simp <- function(mz.xcms, time, ...) {
  scanIndex = which(abs(mz.xcms@scantime - time) == min(abs(mz.xcms@scantime - time)))
  getScan(mz.xcms, scanIndex, ...)
}

plotEIC.simp <- function(mz.xcms, mz, time_range, tol = 0.005, ...) {
  mzrange <- c(mz - tol, mz + tol)
  plotEIC(mz.xcms, mzrange, time_range, ...)
}


