setwd("/home/tux/Work/paper/R_Journal/qPCR_workflow/figures/")
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
# workspace and assign it to the object tmp.
require(qpcR)
tmp <- guescini1

# Define the threshold value for the th.cyc function
Ct <- 0.05

# Define the diltuion of the sample DNA quantity for
# the calibration curve.

dil <- sapply((2:-4), function(i) {10^i})

# Preporcess the amplification curve data with the CPP function from the chipPCR
# package.
res.CPP <- cbind(tmp[, 1], apply(tmp[, -1], 2, function(x) {
    CPP(tmp[, 1], x, trans = TRUE, method.norm = "minm", bg.range = c(1,7))$y.norm
}))

# Use the th.cyc function from the chipPCR package to calculate the Cq values
# by the cycle threshold method.

Cq.Ct <- apply(tmp[, -1], 2, function(x) {th.cyc(res.CPP[, 1], x, r = Ct)[1]})
Cq.SDM <- apply(tmp[, -1], 2, function(x) {summary(inder(res.CPP[, 1], x))[2]})

res.Cq <- lm(Cq.Ct ~ Cq.SDM)

pdf("dilution_Cq.pdf", width = 9.5, height = 14)
layout(matrix(c(1,2,3,3,4,5), 3, 2, byrow = TRUE))

matplot(tmp[, -1], type = "l", lty = 1, col = 1, xlab = "Cycle", 
	    ylab = "RFU", main = "Raw data")
legend("topleft", "A", cex = 3, bty = "n")

matplot(res.CPP[, -1], type = "l", lty = 1, col = 1, xlab = "Cycle", 
	ylab = "RFU", main = "Pre-processed data")
legend("topleft", "B", cex = 3, bty = "n")
abline(h = Ct, col = "red", lwd = 2)

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


# Load lc96_bACTXY.rdml dataset form RDML package and assign the data to the
# object LC96.dat. The data were measured with CFX96 (Bio-Rad). The data set
# contains qPCR data with four targets and two types.

path <- path.package("RDML")
filename <- paste(path, "/extdata/", "BioRad_qPCR_melt.rdml", sep = "")
BioRad <- RDML(filename, name.pattern = "%TUBE%_%NAME%_%TYPE%_%TARGET%")
# 
# 
# # Fetch cycle dependent fluorescence for EvaGreen chanel of NAME_GENE_HERE.
qPCR <- cbind(BioRad$qPCR$EvaGreen$pos, BioRad$qPCR$EvaGreen$ntc[, -1])
melt <- cbind(BioRad$Melt$EvaGreen$pos, BioRad$Melt$EvaGreen$ntc[, -1])

res.diffQ <- lapply(2:ncol(melt), function(x) {
						res <- mcaSmoother(melt[, 1], melt[, x], Trange = c(70, 95))
						diffQ(res, verbose = TRUE, inder = TRUE)
						}
	     )

# # Use plotCurves function from the chipPCR package to get an overview of the
# # amplification curve samples.
# 
pdf("plotCurves.pdf", width = 6, height = 4)

plotCurves(qPCR[, 1], qPCR[, -1], type = "l")

dev.off()

pdf("amp_melt.pdf", width = 8, height = 6)

layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
matplot(qPCR[, 1], qPCR[, -1], type = "l", col = c(rep(1,12), rep(2,12)), lty = 1, xlab = "Cycle", 
	    ylab = "RFU")

matplot(melt[, 1], melt[, -1], type = "l", col = c(rep(1,12), rep(2,12)), lty = 1, xlab = "Temperature", 
	    ylab = "RFU")
plot(NA, NA, xlim = c(70, 93), ylim = c(0,25), xlab = "Temperature", ylab = "-d(RFU)/dT")
color <- c(rep(1,12), rep(2,12))
lapply(1L:24, function(i) {lines(res.diffQ[[i]]$xy, col = color[i])})

dev.off()
# 
# dil <- as.vector(lc96[["Dilutions"]][["FAM"]])


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
Ct <- 0.05

pdf("qIA.pdf")
par(mfrow = c(2,1))
plot(NA, NA, xlim = c(0, 120), ylim = c(0.4,1.2), xlab = "Time (min)", ylab = "RFU")
legend("topleft", "A", cex = 3, bty = "n")
lapply(c(2,4), function(i) {lines(C81[, i]/60, C81[, i + 1], type = "b", pch = 20)})

plot(NA, NA, xlim = c(0, 120), ylim = c(0,0.8), xlab = "Time (min)", ylab = "RFU")
legend("topleft", "B", cex = 3, bty = "n")
res <- lapply(c(2,4), function(i) {
			    y.s <- CPP(C81[, i]/60, C81[, i + 1], trans = TRUE, method = "spline", bg.outliers = TRUE, bg.range = c(1, 190))
			    lines(C81[, i]/60, y.s$y.norm, type = "b", pch = 20)
			    th.cyc(C81[, i]/60, y.s$y.norm, r = Ct)
			    })
abline(h = Ct, lty = 2)
dev.off()
