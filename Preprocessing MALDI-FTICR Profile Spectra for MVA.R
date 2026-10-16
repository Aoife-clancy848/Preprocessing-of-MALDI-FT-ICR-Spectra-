#List Spectra files 

Spectra <- list.files("#path to spectral files", 
                      pattern = "\\.mzml$", 
                      full.names = TRUE)


library(MALDIquant)
library(MALDIquantForeign)  # for mzML, mzXML, Bruker file

#Import and Preprocess Spectra

profile.spectraNP <- importMzMl(Spectra, centroided = FALSE)


spectra.smooth <- smoothIntensity(
  profile.spectraNP,
  method = "SavitzkyGolay",
  halfWindowSize = 5,
  polynomialOrder = 3
)


spectra.transform <- transformIntensity(spectra.smooth, method = "sqrt")

spectra.base <- removeBaseline(
  spectra.transform,
  method = "SNIP",
  iterations = 100
)

spectra.cal <- calibrateIntensity(spectra.base, method = "TIC")

#Generating Average Spectrum

spectra.centroid <- detectPeaks(
  spectra.cal,
  method = "MAD",
  halfWindowSize = 25,
  SNR = 7
)

ref <- averageMassSpectra(spectra.cal, method = "mean")

ref.peaks <- detectPeaks(
  ref,
  method = "MAD",
  halfWindowSize = 25,
  SNR = 7
)

#align spectra to reference spectrum

align_tol <- 2e-6 * 800  # ~2 ppm at m/z 800

spectra.aligned <- alignSpectra(
  spectra.cal,
  peaks = spectra.centroid,
  reference = ref.peaks,
  tolerance = align_tol   
)

sapply(spectra.centroid[1:5], length)

# Re-detect peaks
peaks.aligned <- detectPeaks(
  spectra.aligned,
  method = "MAD",
  halfWindowSize = 25,
  SNR = 7
)

#Checking Alignment and Deviations
idx <- c(1, 10, 21)
plot(
  mass(spectra[[idx[1]]]),
  intensity(spectra[[idx[1]]]),
  type = "l",
  col = "black",
  xlab = "m/z",
  ylab = "Intensity",
  main = "Unaligned spectra"
)

for (i in idx[-1]) {
  lines(mass(spectra[[i]]), intensity(spectra[[i]]), col = "grey50")
}


plot(
  mass(spectra.aligned[[idx[1]]]),
  intensity(spectra.aligned[[idx[1]]]),
  type = "l",
  col = "black",
  xlab = "m/z",
  ylab = "Intensity",
  main = "Aligned spectra"
)

for (i in idx[-1]) {
  lines(mass(spectra.aligned[[i]]),
        intensity(spectra.aligned[[i]]),
        col = "grey50")
}


# pick a strong reference peak
mz663.384 <- mass(ref.peaks)[1]

# deviation before alignment
dev.before <- sapply(spectra, function(s)
  min(abs(mass(s) - mz663.384)))

# deviation after alignment
dev.after <- sapply(spectra.aligned, function(s)
  min(abs(mass(s) - mz663.384)))

summary(dev.before)
summary(dev.after)

#Binning

#Check optimal binning ppm
ppm_vals <- c(3, 5, 7)

results <- lapply(ppm_vals, function(ppm) {
  tol <- ppm * 1e-6 * 1000
  
  binned <- binPeaks(peaks.aligned, tolerance = tol)
  mat <- intensityMatrix(binned, spectra.aligned)
  
  data.frame(
    ppm = ppm,
    features = ncol(mat),
    NA_fraction = mean(is.na(mat))
  )
})

do.call(rbind, results)

sapply(peaks[1:5], length)

bin_tol <- 3e-6 * 1000   # 3 ppm

binned.peaks <- binPeaks(
  peaks.aligned,
  tolerance = bin_tol   
)

#Generate a feature matrix

feature.matrix <- intensityMatrix(binned.peaks, spectra.aligned)
mean(is.na(feature.matrix))

#Export Matrix as CSV
sapply(profile.spectraNP, metaData)
sample_names <- basename(Spectra)
sample_names <- tools::file_path_sans_ext(sample_names)
length(sample_names)
nrow(feature.matrix)
length(profile.spectraNP)
rownames(feature.matrix) <- sample_names

write.table(feature.matrix,
            file = "FTICR_MALDI_Nanoparticles.csv",
            sep = ";",
            row.names = TRUE,
            col.names = NA,
            quote = FALSE)

getwd()