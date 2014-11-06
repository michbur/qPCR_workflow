if(!grepl("figures", getwd()))
  setwd(paste0(getwd(), "/figures/"))
# Supplement to 'R as Platform for the Analysis of qPCR experiments' for the R
# Journal by Rödiger et al. 2014
#################################
# Example one
#################################
# Load the required packages for the data import and analysis.
# Load the chipPCR package for the pre-processing and curve data quality
# analysis.
require(chipPCR)

# Collect information about the R session used for the analysis of the qPCR
# experiment.
current.session <- sessionInfo()

# Next we load the 'guescini1' dataset from the qpcR package the to
# workspace and assign it to the object gue.
require(qpcR)
gue <- guescini1

# Define the diltuion of the sample DNA quantity for
# the calibration curve.

dil <- 10^(2:-4)

# Preprocess the amplification curve data with the CPP function from the chipPCR
# package.
res.CPP <- cbind(gue[, 1], apply(gue[, -1], 2, function(x) {
    CPP(gue[, 1], x, trans = TRUE, method.norm = "minm", bg.range = c(1,7))[["y.norm"]]
}))

# Use the th.cyc function from the chipPCR package to calculate the Cq values
# by the cycle threshold method. The threshold level r was set to 0.05.

Cq.Ct <- apply(gue[, -1], 2, function(x) 
  th.cyc(res.CPP[, 1], x, r = 0.05)[1])
Cq.SDM <- apply(gue[, -1], 2, function(x) 
  summary(inder(res.CPP[, 1], x))[2])

res.Cq <- lm(Cq.Ct ~ Cq.SDM)

pdf("dilution_Cq.pdf", width = 9.5, height = 14)
layout(matrix(c(1,2,3,3,4,5), 3, 2, byrow = TRUE))

matplot(gue[, -1], type = "l", lty = 1, col = 1, xlab = "Cycle", 
	    ylab = "RFU", main = "Raw data")
legend("topleft", "A", cex = 3, bty = "n")

matplot(res.CPP[, -1], type = "l", lty = 1, col = 1, xlab = "Cycle", 
	ylab = "RFU", main = "Pre-processed data")
legend("topleft", "B", cex = 3, bty = "n")
abline(h = 0.05, col = "red", lwd = 2)

plot(Cq.SDM, Cq.Ct, xlab = "Ct method", ylab = "SDM method", 
     main = "Comparison of Cq methods")
abline(res.Cq)
legend("topleft", "C", cex = 3, bty = "n")

plot(effcalc(dil, t(matrix(Cq.Ct, nrow = 12, ncol = 7))), CI = TRUE)
legend("topright", "D", cex = 3, bty = "n")

plot(effcalc(dil, t(matrix(Cq.SDM, nrow = 12, ncol = 7))), CI = TRUE)
legend("topright", "E", cex = 3, bty = "n")

dev.off()
#################################
# Example two
#################################
# Load the required packages for the data import and analysis.

# Import the qPCR and melting curve data via the RDML package
require(RDML)

# Load the chipPCR package for the pre-processing and curve data quality
# analysis.
require(chipPCR)

# Load the MBmca package for the melting curve analysis.
require(MBmca)

# Collect information about the R session used for the analysis of the qPCR
# experiment.
current.session <- sessionInfo()


# Load the BioRad_qPCR_melt.rdml file form RDML package and assign the data to the
# object BioRad.

filename <- paste(path.package("RDML"), "/extdata/", "BioRad_qPCR_melt.rdml", sep = "")
BioRad <- RDML(filename, name.pattern = "%TUBE%_%NAME%_%TYPE%_%TARGET%")

# Fetch cycle dependent fluorescence for the EvaGreen channel of the 
# Mycobacterium tuberculosis katG gene and aggregate the data in the 
# object qPCR. 
qPCR <- cbind(BioRad[["qPCR"]][["EvaGreen"]][["pos"]], 
	      BioRad[["qPCR"]][["EvaGreen"]][["unkn"]][, -1], 
	      BioRad[["qPCR"]][["EvaGreen"]][["ntc"]][, -1])
# Leave data only from row 'D' that contains target 'Cy5-2' at channel 'Cy5'
qPCR <- cbind(qPCR[,1], qPCR[, grep("^D", names(qPCR))])

# Use plotCurves function from the chipPCR package to get an overview of the
# amplification curve samples.

pdf("plotCurves.pdf", width = 6, height = 4)

plotCurves(qPCR[, 1], qPCR[, -1], type = "l")

dev.off()
# Fetch temperature dependent fluorescence for the Cy5 channel of the 
# probe that can hybridize with Mycobacterium tuberculosis katG gene (codon 315)
# and aggregate the data in the object melt.
melt <- cbind(BioRad[["Melt"]][["EvaGreen"]][["pos"]], 
	      BioRad[["Melt"]][["EvaGreen"]][["unkn"]][, -1], 
	      BioRad[["Melt"]][["EvaGreen"]][["ntc"]][, -1])

# Calculate the melting temperature with the diffQ function
# from the MBmca package. Use simple logical conditions to find out
# if a wild-type sample with the expected Tm of circa 54.5 degree 
# Celsius is found.
res.Tm <- apply(melt[, -1], 2, function(x) {
  res <- mcaSmoother(melt[, 1], x, Trange = c(70, 95))
  res.Tm <- diffQ(res, fct = max, inder = TRUE)
  Decission <- ifelse(res.Tm[1] > 82 & res.Tm[1] < 87 & res.Tm[2] > 10, 1, 0)
  out <- data.frame(res.Tm[c(1,2)], Decission)
})     
# Present the results in a tabular output as matrix "results".	      
(results.Tm <- matrix(unlist(res.Tm), nrow = length(res.Tm), byrow = TRUE, 
       dimnames = list(colnames(melt[, -1]),
       c("Tm", "Height", "Decission"))))
       
pdf("amp_melt.pdf", width = 8, height = 6)

# Convert the Decission from the "relsults" object in a color code:
# Negative, black; Positive, red.

color <- c(results.Tm[, 3] + 1)

# Arrange the results of the calculations in plot.
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))

# Use the CPP function to preporcess the amplification curve data.
plot(NA, NA, xlim = c(1, 40), ylim = c(0,200), xlab = "Cycle", ylab = "RFU")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)
lapply(2L:ncol(qPCR), function(i) 
  lines(qPCR[, 1], 
        CPP(qPCR[, 1], qPCR[, i], trans = TRUE, 
            bg.range = c(1,9))[["y.norm"]],
        col = color[i - 1]))
matplot(melt[, 1], melt[, -1], type = "l", col = color, 
	lty = 1, xlab = "Temperature [°C]", ylab = "RFU")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)
	
plot(NA, NA, xlim = c(35, 95), ylim = c(-15,30), xlab = "Temperature [°C]", 
     ylab = "-d(RFU)/dT")
mtext("C", cex = 2, side = 3, adj = 0, font = 2)
lapply(2L:ncol(melt), function(i)
  lines(diffQ(cbind(melt[, 1], melt[, i]), verbose = TRUE, 
              fct = max, inder = TRUE)[["xy"]], col = color[i - 1]))
dev.off()

res.Cq <- lapply(2L:ncol(qPCR), function(i) {
	      res <- CPP(qPCR[, 1], qPCR[, i], trans = TRUE, bg.range = c(10, 20))[["y.norm"]]
	      th.cyc <- th.cyc(qPCR[, 1], res, r = 5)[1]
	      })
	      
result.Cq <- matrix(unlist(res.Cq), nrow = length(res.Cq), byrow = TRUE, 
       dimnames = list(colnames(qPCR[, -1]),
       "Cq"))
       
result.Cq
#################################
# Example three
#################################
require(dpcR)

# Analysis of a digital PCR experiment. The density estimation
# analysis of results of biorad experiment; of 16800 counted droplets (n) 4601 were positive (k).
pdf("dpcR.pdf")

(dens <- dpcr_density(k = 4601, n = 16800, average = TRUE, methods = "wilson"))

dev.off()
# Let us assume, that every droplet has roughly 5 nl 
# total concentration (and its confidence intervals) in molecules/ml
dens[4:6] / 5 * 1e-6


#################################
# Example four
#################################

pdf("qIA.pdf")

# Drawn in an 2-by-1 array on the device by two columns and one row.
par(mfrow = c(2, 1))

# Plot the raw data from the C81 dataset to the first array and add
# a legend.
plot(NA, NA, xlim = c(0, 120), ylim = c(0.4, 1.2), xlab = "Time (min)", ylab = "RFU")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)
lapply(c(2, 4), function(i) {
    lines(C81[, i]/60, C81[, i + 1], type = "b", pch = 20, col = i - 1)
})
legend(10, 0.8, c("D1: 1x", "D2: 1:10 diluted sample"), pch = 19, col = c(1, 3), 
    bty = "n")

# Prepare a plot on the second array for the pre-proccessed data.
plot(NA, NA, xlim = c(0, 120), ylim = c(0, 0.8), xlab = "Time (min)", ylab = "RFU")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

# Apply the CPP functions to pro-process the raw data.
# 1) Basline data to zero, 2) Smooth data with spline,
# 3) Remove outliers in background range between 
# entry 1 and 190.
res <- lapply(c(2, 4), function(i) {
  y.s <- CPP(C81[, i]/60, C81[, i + 1],
             trans = TRUE, 
             method = "spline",
             bg.outliers = TRUE,
             bg.range = c(1, 190))
  lines(C81[, i]/60, y.s[["y.norm"]], type = "b", pch = 20, col = i - 1)
  # Use the th.cyc function to calculate the cycle threshold time. 
  # The threshold level r was set to 0.05.
  paste(round(th.cyc(C81[, i]/60, y.s[["y.norm"]], r = 0.05)[1], 2), "min")
})

# Add the cycle threshold time and the threshold level to plot.

abline(h = 0.05, lty = 2)
text(10, 0.55, "Cq:")
legend(10, 0.5, paste(c("D1: ", "D2: "), res), pch = 19, col = c(1, 3), bty = "n")
dev.off()
