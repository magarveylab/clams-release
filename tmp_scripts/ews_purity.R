pda_rt_offset = 9.45
auc_window = 2.5 * 60

#1,2,4,10 have no blanks
peaks_f <- "peaks_for_scan.json"
peaks <- rjson::fromJSON(file = peaks_f)
peaks_parsed <- parse_peaks(peaks)
width <- 2327
height <- 899
for(pp in peaks_parsed) {
  mz.file <- paste0("/mnt/storage/mass_spec_samples/data/", pp[2])
  min_rt <- as.numeric(pp[4])
  max_rt <- as.numeric(pp[5])
  pda_data_blank <- read.blank.pda(mz.file)
  #pda_data_blank <- read.pda.file(mz.file.blank)
  if (!isempty(pda_data_blank)) {
    pda_data <- read.pda.file(mz.file)
    pda_220 <- standardize.ews(pda_data, pda_data_blank)
    pda_rts <- get.pda.rts(pda_data, pda_rt_offset)
    scans_per_second <- 1 / (pda_rts[2] - pda_rts[1])
    ews_rt = cbind(pda_rts, pda_220)
    AUC.total <- get.auc(ews_rt, get.rt.scan.index(pda_rts, min_rt - auc_window), get.rt.scan.index(pda_rts, max_rt + auc_window))
    AUC.peak <- get.auc(ews_rt, get.rt.scan.index(pda_rts, min_rt), get.rt.scan.index(pda_rts, max_rt))
    purity <- AUC.peak / AUC.total
    png(paste0('pda_figs/', pp[3], '_purity-', round(purity, 2) ,'.png'), width = width, height = height)
    plot(ews_rt)
    dev.off()
  }
  else{
    print(mz.file)
  }
}

get.auc <- function(xy, min_x = NULL, max_x = NULL) {
  if (is.null(min_x)) {
    min_x = xy[1,1]
  }
  if (is.null(max_x)) {
    max_x = xy[nrow(xy), 1]
  }
  trapz(xy[min_x:max_x,1], xy[min_x:max_x,2])
}

standardize.ews <- function(pda_data, pda_data_blank, wavelength = 220) {
  #indexes are the same between the query and the blank
  wavelength_indexs <- pda_data$pda_wavelengths
  wavelength_index <- match(wavelength, wavelength_indexs)
  ews <- get.extracted.wavelength.spectra(pda_data, wavelength_index)
  ews_blank <- get.extracted.wavelength.spectra(pda_data_blank, wavelength_index)
  ews = remove.blank(ews, ews_blank)
  ews = ews - runmed(ews,751)
  ews[ews < 0] <- 0
  ews
}

#removes blank from query based on normalized to max intensities
remove.blank <- function(ews, ews_blank) {
  ews_max <- max(abs(ews))
  ews <- ews / ews_max
  ews_blank <- ews_blank / max(abs(ews_blank))
  ews <- ews - ews_blank
  ews <- ews * ews_max
  ews
}

read.blank.pda <- function(mz.file) {
  split_name <- strsplit(mz.file, '/')[[1]]
  dir <- paste(split_name[1:length(split_name) -1 ], collapse = '/')
  name <- split_name[length(split_name)]
  date <- strsplit(name, '_')[[1]][2]
  date = gsub('17-10-16', '17-10-18', date)
  dir_files <- list.files(dir)
  mzml_files <- dir_files[grepl('mzml', dir_files)]
  blank_files <- mzml_files[grepl('_blank', mzml_files)]
  same_date <- blank_files[grepl(date, blank_files)]
  blank_file <- paste(c(dir,same_date[1]), collapse = '/') #can just choose the first blank, any will work
  read.pda.file(blank_file)
}
