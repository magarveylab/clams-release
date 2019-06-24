#mz.file <- "/mnt/storage/mass_spec_samples/data/2017_5/RM_17-05-12_11648_NE-042944-SC.mzml"
mz.file <- "/mnt/storage/mass_spec_samples/data/2017_10/AMK_17-10-16_17081_BTA-0514_NE-076270_590mz_MeOH_bin2_fractionation_first_round.mzml"
split_name <- strsplit(mz.file, '/')[[1]]
dir <- paste(split_name[1:length(split_name) -1 ], collapse = '/')
name <- split_name[length(split_name)]
date <- strsplit(name, '_')[[1]][2]
date = gsub('17-10-16', '17-10-18', date)
dir_files <- list.files(dir)
mzml_files <- dir_files[grepl('mzml', dir_files)]
blank_files <- mzml_files[grepl('_blank', mzml_files)]
same_date <- blank_files[grepl(date, blank_files)]
mz.blank.file <- paste(c(dir,same_date[1]), collapse = '/') #can just choose the first blank, any will work
maxMZ <- 2000
minMZ <- 100
rtMin <- 30
MZtolerance <- 0.005
RTtolerance <- 30
backgroundIntensity <- 1000
correctedMZCutoff <- 1500

#test standards
mz_files <- c("pda_tmp/HL_2017-09-14_polymyxin_standard_02mM.mzml", "pda_tmp/HL_2017-09-14_rifamycin_standard_02mM.mzml", "pda_tmp/HL_2017-09-14_tetracycline_standard_02mM.mzml")
pda_files <- c("pda_tmp/HL_2017-09-14_polymyxin_standard_02mM_pda.json", "pda_tmp/HL_2017-09-14_rifamycin_standard_02mM_pda.json", "pda_tmp/HL_2017-09-14_tetracycline_standard_02mM_pda.json")
mz_files2 <- c("pda_tmp/HL_17-09-18_PDA_standard_mixture.mzml", "pda_tmp/HL_17-09-18_PDA_standard_novobiocin.mzml", "pda_tmp/HL_17-09-18_PDA_standard_podophylotoxin.mzml")
pda_files2 <- c("pda_tmp/HL_17-09-18_PDA_standard_mixture_pda.json", "pda_tmp/HL_17-09-18_PDA_standard_novobiocin_pda.json", "pda_tmp/HL_17-09-18_PDA_standard_podophylotoxin_pda.json")
mz_files <- c(mz_files, mz_files2)
pda_files <- c(pda_files, pda_files2)
lc <- LCMSunit(mz.file)
lc <- LCMSunitAnalyze(lc, list(mz.blank.file), c(minMZ, maxMZ), rtMin , MZtolerance , RTtolerance , backgroundIntensity, correctedMZCutoff)

ios <- lcs$IsotopeOverview
drops <- c("distributionMasses", "distribution")
ios <- lapply(ios, function (x) x[ , !(names(x) %in% drops)])

pda_poly <- analyse.pda(ios[[1]], pda_files[1])
pda_poly_out <- as.matrix(pda_poly)
write.table(pda_poly_out, file="pda_tmp/HL_2017-09-14_polymyxin_standard_02mM_pda_peaks_4.tsv", quote=FALSE, sep="\t")

pda_rif <- analyse.pda(ios[[2]], pda_files[2])
pda_rif_out <- as.matrix(pda_rif)
write.table(pda_rif_out, file="pda_tmp/HL_2017-09-14_rifamycin_standard_02mM_pda_peaks_4.tsv", quote=FALSE, sep="\t")

pda_tetra <- analyse.pda(ios[[3]], pda_files[3])
pda_tetra_out <- as.matrix(pda_tetra)
write.table(pda_tetra_out, file="pda_tmp/HL_2017-09-14_tetracycline_standard_02mM_pda_peaks_4.tsv", quote=FALSE, sep="\t")

pda_stand <- analyse.pda(ios[[4]], pda_files[4])
pda_stand_out <- as.matrix(pda_stand)
write.table(pda_stand_out, file="pda_tmp/HL_17-09-18_PDA_standard_mixture_pda_peaks_4.tsv", quote=FALSE, sep = '\t')

pda_novo <- analyse.pda(ios[[5]], pda_files[[5]])
pda_novo_out <- as.matrix(pda_novo)
write.table(pda_novo_out, file="pda_tmp/HL_17-09-18_PDA_standard_novobiocin_pda_peaks_4.tsv", quote=FALSE, sep = '\t')

pda_podo <- analyse.pda(ios[[6]], pda_files[[6]])
pda_podo_out <- as.matrix(pda_podo)
write.table(pda_podo_out, file="pda_tmp/HL_17-09-18_PDA_standard_podophylotoxin_pda_peaks_4.tsv", quote=FALSE, sep = '\t')

