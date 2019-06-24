library(xcms)
library(ggplot2)
library(ggrepel)
#METHODSSSSSSSSSS

parse_peaks <- function(peaks) {
  lapply(peaks, function(x) {
    c(x$intens_rt, x$rel.clams_results$rel.mzml$file_path, x$id, x$rt_min, x$rt_max)
  })
}
plotScan.cust <- function(object, time, mzrange = numeric()) {
  scan = which(abs(object@scantime - time) == min(abs(object@scantime - time)))

  if (scan<1 || scan>length(object@scanindex) ) {
    warning("scan out of range")
    return()
  }
  
  ## handle last spectrum
  if (scan == length(object@scanindex)) {
    followingScanIndex <- length(object@env$mz)
  } else {
    followingScanIndex <- object@scanindex[scan+1]
  }
  
  ## hendle empty spectra
  if (object@scanindex[scan] == length(object@env$mz) ||
      object@scanindex[scan] == followingScanIndex) {
    warning("empty scan")
    return()
  }
  
  idx <- (object@scanindex[scan]+1):min(followingScanIndex,
                                        length(object@env$mz), na.rm=TRUE)

  points <- cbind(object@env$mz[idx], object@env$intensity[idx])
  title = paste("Mass Spectrum: ", round(object@scantime[scan], 1),
                " seconds (scan ", scan, ")", sep = "")
  longData <- data.frame(points)
  top_points_a <- data.frame(top_points(points))
  colnames(top_points_a) <- c("X1","X2")
  return(ggplot(longData, aes(x = X1, y = X2, label=round(X1, digits = 2))) +
          geom_col(width = 1) + ylab("Intensity") + xlab("m/z") + ggtitle(title) +
          geom_text_repel(data=top_points_a, aes(X1, X2, label = round(X1, digits = 3)), segment.size = 0.1, size = 2, direction = 'y') +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x = element_text(), axis.text.y = element_text()))
  
  #plot(points, type="h", main = title, xlab="m/z", ylab="Intensity")
  #text(top_points[,1],top_points[,2],round(top_points[,1],3), cex=2)
  #invisible(points)
}

top_points <- function(points, max_shown = 10, num_consider = 50) {
  top_ints <- tail(sort(points[,2]), num_consider)
  top_indx <- sapply(top_ints, function(x) match(x, points[,2]))
  to_draw <- top_indx
  for(q_indx in top_indx) {
    if (q_indx %in% to_draw) {
      for(s_indx in to_draw) {
        if(q_indx != s_indx) {
          if(points[q_indx,2] <= points[s_indx,2]) {
            if(abs(points[q_indx,1] - points[s_indx,1]) < 4) {
              to_draw <- to_draw[!to_draw==q_indx]
            }
          }
        }
      }
    }
  }
  
  out <- points[to_draw,]
  top_ints <- tail(sort(out[,2]), max_shown)
  top_indx <- sapply(top_ints, function(x) match(x, out[,2]))
  out[top_indx,]
}

#CODEEEEEEE

peaks_f <- "peaks_for_scan.json"
peaks <- rjson::fromJSON(file = peaks_f)
peaks_parsed <- parse_peaks(peaks)
width <- 2327
height <- 899
for(pp in peaks_parsed) {
  rt <- as.numeric(pp[1])
  xc <- xcmsRaw(paste0('/mnt/storage/mass_spec_samples/data/',pp[2]))
  bn <- paste0('scan_figs/', pp[3], '_')
  p <- plotScan.cust(xc, rt - 10)
  ggsave(filename=paste0(bn,'_', rt - 10, '.png'), plot = p, width = width/300, height = height/300, device = "png")
  p <- plotScan.cust(xc, rt)
  ggsave(filename=paste0(bn,'_', rt, '.png'), plot = p, width = width/300, height = height/300, device = "png")
  p <- plotScan.cust(xc, rt + 10)
  ggsave(filename=paste0(bn,'_', rt + 10, '.png'), plot = p, width = width/300, height = height/300, device = "png")
}
