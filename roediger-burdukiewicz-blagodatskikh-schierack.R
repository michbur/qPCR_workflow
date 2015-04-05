if(!grepl("figures", getwd()))
  setwd(paste0(getwd(), "/figures/"))
# Supplement to 'R as an Environment for the Analysis of dPCR and qPCR Experiments'
# Journal by Rödiger et al. 2014
#################################
# Case study one
#################################
# Load the required packages for the data import and analysis.
# Load the chipPCR package for the pre-processing and curve data quality
# analysis and load the qpcR package as data resource.
require(chipPCR)
require(qpcR)

# Collect information about the R session used for the analysis of the
# experiment.
current.session <- sessionInfo()

# Next load the 'guescini1' dataset from the qpcR package to the
# workspace and assign it to the object 'gue'.
gue <- guescini1

# Define the dilution of the sample DNA quantity for the calibration curve 
# and assign it to the object 'dil'.
dil <- 10^(2:-4)

# Pre-process the amplification curve data with the CPP function from the 
# chipPCR package. The trans parameter was set TRUE to perform a baselining and 
# the method.norm parameter was set to minm for a min-maximum normalization. All
# amplification curves were smoothed by Savitzky-Golay smoothing.

res.CPP <- cbind(gue[, 1], apply(gue[, -1], 2, function(x) {
  CPP(gue[, 1], x, trans = TRUE, method.norm = "minm", method.reg = "least", 
      bg.range = c(1,7))[["y.norm"]]
}))

# Use the th.cyc function from the chipPCR package to calculate the Cq values
# by the cycle threshold method at a threshold signal level "r" of 0.05.
Cq.Ct <- apply(gue[, -1], 2, function(x) 
  th.cyc(res.CPP[, 1], x, r = 0.05)[1])

# Use the inder function from the chipPCR package to calculate the Cq values
# by the SDM method.
Cq.SDM <- apply(gue[, -1], 2, function(x)
  summary(inder(res.CPP[, 1], x))[2])

# Fit a linear model to carry out a regression analysis.
res.Cq <- lm(Cq.Ct ~ Cq.SDM)

summary(res.Cq)

#pdf("dilution_Cq.pdf", width = 9.5, height = 14)

# Arrange and plot the results in a convenient way.
layout(matrix(c(1,2,3,3,4,5), 3, 2, byrow = TRUE))
# Set bigger top margin.
par(mar = c(5.1, 4.1, 6.1, 2.1))

# Plot the raw amplification curve data.
matplot(gue[, -1], type = "l", lty = 1, col = 1, xlab = "Cycle", 
        ylab = "RFU", main = "Raw data")
legend("topleft", "A", cex = 3, bty = "n")

# Plot the pre-processed amplification curve data.
matplot(res.CPP[, -1], type = "l", lty = 1, col = 1, xlab = "Cycle", 
        ylab = "RFU", main = "Pre-processed data")
legend("topleft", "B", cex = 3, bty = "n")
abline(h = 0.05, col = "red", lwd = 2)

# Plot Cq.SDM against Cq.Ct and add the trendline from the linear regression
# analysis.

plot(Cq.SDM, Cq.Ct, xlab = "Ct method", ylab = "SDM method", 
     main = "Comparison of Cq methods")
abline(res.Cq)

legend("topleft", "C", cex = 3, bty = "n")

# Use the effcalc function from the chipPCR package to calculate the
# amplification efficiency.
plot(effcalc(dil, t(matrix(Cq.Ct, nrow = 12, ncol = 7))), CI = TRUE)
legend("topright", "D", cex = 3, bty = "n")

plot(effcalc(dil, t(matrix(Cq.SDM, nrow = 12, ncol = 7))), CI = TRUE)
legend("topright", "E", cex = 3, bty = "n")

# Set top margin to default value.
par(mar = c(5.1, 4.1, 4.1, 2.1))
#dev.off()
#################################
# Case study two
#################################
# Import the qPCR and melting curve data via the RDML package.
# Load the chipPCR package for the pre-processing and curve data quality
# analysis and the MBmca package for the melting curve analysis.
require(RDML)
require(chipPCR)
require(MBmca)
require(dplyr)

# Collect information about the R session used for the analysis of the qPCR
# experiment.
current.session <- sessionInfo()

# Load the BioRad_qPCR_melt.rdml file form RDML package and assign the data to the
# object BioRad.
filename <- paste(path.package("RDML"), "/extdata/", "BioRad_qPCR_melt.rdml", sep = "")
BioRad <- RDML$new(filename)

# Structure of experiment can be overviewed by AsDendrogram() function
# We can see that our experiment contains two detection channels
# ('EvaGreen' and 'Cy5' at 'Run ID'). 'EvaGreen' channel has one
# probe (target) - 'EvaGreen'. 'Cy5' has: 'Cy5', 'Cy5-2' and 'Cy5-2_rr'.
# each target has three sample types (positive, unknown, negative).
# And each sample type has qPCR ('adp') and melting ('mdp') data.
# Last column shows how many samples of this type at this experiment.
BioRad$AsDendrogram()

# Fetch cycle dependent fluorescence for the EvaGreen channel and row 'D'
# (that contains target 'Cy5-2' at channel 'Cy5') of the 
# katG gene and aggregate the data in the object qPCR. 

qPCR <- BioRad$AsTable() %>%
  filter(target == "EvaGreen",
         grepl("^D", position))  %>% 
  BioRad$GetFData(.)

# Use plotCurves function to get an overview of the amplification curve samples.
#pdf("plotCurves.pdf", width = 6, height = 4)

plotCurves(qPCR[, 1], qPCR[, -1], type = "l")

#dev.off()
# Detect positive samples - calculate Cq values by the cycle threshold method. 
# The threshold signal level r was set to 10.
Cq.Positive <- t(apply(qPCR[, -1], 2, function(x)
{
  res <- CPP(qPCR[, 1], x, trans = TRUE, bg.range = c(2, 8),
             method.reg = "least")[["y.norm"]]
  # The th.cyc fails when the threshold exceeds maximum 
  # observed fluorescence values, so it must be used with try()
  th.cycle <- try(th.cyc(qPCR[, 1], res, r = 10)[1], silent = TRUE)
  cq <- ifelse(class(th.cycle) != "try-error", as.numeric(th.cycle), NA)
  pos <- !is.na(cq)
  c(Cq=cq, M.Tub_positive = pos)
}
))

# Fetch temperature dependent fluorescence for the Cy5 channel of the 
# probe 'Cy5-2' that can hybridize with Mycobacterium tuberculosis 
# katG gene (codon 315) and aggregate the data in the object 'melt'.
melt <- BioRad$AsTable() %>%
  filter(target == "Cy5-2")  %>% 
  BioRad$GetFData(., data.type = "mdp")

# Calculate the melting temperature with the diffQ function from the MBmca 
# package. Use simple logical conditions to find out if a positive sample with 
# the expected Tm of circa 54.5 degree Celsius is found. The result of the test
# is assigned to the object 'positive'.
Tm.Positive <- matrix(nrow = ncol(melt) - 1,
                      byrow = TRUE,
                      dimnames = list(colnames(melt)[-1]),
                      unlist(apply(melt[, -1], 2, function(x) {
                        res.Tm <- diffQ(cbind(melt[, 1], x), 
                                        fct = max, inder = TRUE)
                        positive <- ifelse(res.Tm[1] > 54 & 
                                             res.Tm[1] < 55 & 
                                             res.Tm[2] > 80, 1, 0)
                        c(res.Tm[1], res.Tm[2], positive)
                      })))

# Present the results in a tabular output as data.frame 'results.tab'.
# Result of analysis logic is:
# Cq.Positive && Tm.Positive = Wild-type
# Cq.Positive && !Tm.Positive = Mutant
# !Cq.Positive && !Tm.Positive = NTC
# !Cq.Positive && Tm.Positive = Error
results <- sapply(1:length(Cq.Positive[,1]), function(i) {
  if(Cq.Positive[i, 2] == 1 && Tm.Positive[i, 3] == 1)
    return("Wild-type")
  if(Cq.Positive[i, 2] == 1 && Tm.Positive[i, 3] == 0)
    return("Mutant")
  if(Cq.Positive[i, 2] == 0 && Tm.Positive[i, 3] == 0)
    return("NTC")
  if(Cq.Positive[i, 2] == 0 && Tm.Positive[i, 3] == 1)
    return("Error")
})

results.tab <- data.frame(cbind(Cq.Positive, Tm.Positive, results))
names(results.tab) <- c("Cq", "M.Tub DNA", "Tm", "Height", 
                        "Tm positive", "Result")

results.tab[["M.Tub DNA"]] <- factor(results.tab[["M.Tub DNA"]], 
                                     labels=c("Not Detected", "Detected"))

results.tab[["Tm positive"]] <- factor(results.tab[["Tm positive"]], 
                                       labels=c(TRUE, FALSE))
results.tab

#pdf("amp_melt.pdf", width = 8, height = 6)

# Convert the decision from the "results" object in a color code:
# Negative, black; Positive, red.

color <- c(Tm.Positive[, 3] + 1)

# Arrange the results of the calculations in plot.
layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))

# Use the CPP function to preporcess the amplification curve data.
plot(NA, NA, xlim = c(1, 41), ylim = c(0,200), xlab = "Cycle", ylab = "RFU")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)
lapply(2L:ncol(qPCR), function(i) 
  lines(qPCR[, 1], 
        CPP(qPCR[, 1], qPCR[, i], trans = TRUE, 
            bg.range = c(1,9))[["y.norm"]],
        col = color[i - 1]))

matplot(melt[, 1], melt[, -1], type = "l", col = color, 
        lty = 1, xlab = "Temperature [degree Celsius]", ylab = "RFU")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

plot(NA, NA, xlim = c(35, 95), ylim = c(-15, 120), xlab = "Temperature [degree Celsius]", 
     ylab = "-d(RFU)/dT")
mtext("C", cex = 2, side = 3, adj = 0, font = 2)

lapply(2L:ncol(melt), function(i)
  lines(diffQ(cbind(melt[, 1], melt[, i]), verbose = TRUE, 
              fct = max, inder = TRUE)[["xy"]], col = color[i - 1]))

#dev.off()
#################################
# Case study three
#################################

#pdf("qIA.pdf")

# Drawn in an 2-by-1 array on the device by two columns and one row.
par(mfrow = c(2, 1))

# Plot the raw data from the C81 dataset to the first array and add
# a legend. Note: The abcsissa values (time in seconds) was divided 
# by 60 (C81[, i] / 60) to convert to minutes.
plot(NA, NA, xlim = c(0, 120), ylim = c(0, 1.2), xlab = "Time (min)", ylab = "RFU")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)
lapply(c(2, 4), function(i) {
  lines(C81[, i] / 60, C81[, i + 1], type = "b", pch = 20, col = i - 1)
})
legend(0, 0.35, c("D1: 1x", "D2: 1:10 diluted sample"), pch = 19, col = c(1, 3), 
       bty = "n")

# Prepare a plot on the second array for the pre-processed data.
plot(NA, NA, xlim = c(0, 120), ylim = c(0, 1.2), xlab = "Time (min)", ylab = "RFU")
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

# Apply the CPP functions to pro-process the raw data.1) Baseline data to zero, 
# 2) Smooth data with a spline, 3) Remove outliers in background range between 
# entry 1 and 190. Assign the results of the analysis to the object 'res'.
res <- lapply(c(2, 4), function(i) {
  y.s <- CPP(C81[, i] / 60, C81[, i + 1],
             trans = TRUE, 
             method = "spline",
             bg.outliers = TRUE,
             bg.range = c(1, 190))
  lines(C81[, i] / 60, y.s[["y.norm"]], type = "b", pch = 20, col = i - 1)
  # Use the th.cyc function to calculate the cycle threshold time (Cq.t). 
  # The threshold signal level r was set to 0.05.
  paste(round(th.cyc(C81[, i] / 60, y.s[["y.norm"]], r = 0.05)[1], 2), "min")
})

# Add the cycle threshold time from the object 'res' to the plot.

abline(h = 0.05, lty = 2)
text(10, 0.55, "Cq.t:")
legend(10, 0.5, paste(c("D1: ", "D2: "), res), pch = 19, col = c(1, 3), 
       bty = "n")
#dev.off()

#################################
# Case study four
#################################
# Load the dpcR package for the analysis of the digital PCR experiment.
require(dpcR)

# Analysis of a digital PCR experiment. The density estimation.
# In our in-silico experiment we counted in total 16800 droplets (n). 
# Thereof, 4601 were positive (k).
#pdf("dpcR.pdf")

(dens <- dpcr_density(k = 4601, n = 16800, average = TRUE, methods = "wilson"))

#dev.off()
# Let us assume, that every droplet has roughly a volume of 5 nL.
# The total concentration (and its confidence intervals) in molecules/ml is:
dens[4:6] / 5 * 1e-6

##################
# dPCR demo
##################
# Generate an amplitude plot for the first fluorescence channel (e.g., FAM)
# fluos1 <- sim_ddpcr(m = 7, n = 20, times = 100, pos_sums = FALSE, n_exp = 1,
#   fluo = list(0.1, 0))
# 
# # Generate an amplitude plot for the second fluorescence channel (e.g., VIC)
# fluos2 <- sim_ddpcr(m = 10, n = 20, times = 100, pos_sums = FALSE, n_exp = 1,
#   fluo = list(0.1, 0))
# #pdf("dpcR_sim.pdf", width = 12, height = 6.5)
# # Plot the amplitudes of both fluorescence channel in an aligned fashion
# plot_vic_fam(fam = fluos1, vic = fluos2, col_vic = "green", col_fam = "pink")
# #dev.off()

# NEW CASE STUDY
# Load the dpcR package for the analysis of the digital PCR experiment.

pdf("dpcR_bioamp.pdf", width = 8, height = 12)

require(dpcR)
# Load the dpcR package for the analysis of the digital PCR experiment.
# Analysis of a droplet dPCR experiment. Data were taken from the pds_raw dataset.

# Select the wells for the analysis.
# A01 to D01 are four replicate dPCR reactions and G04 is the 
# no template control (NTC).
wells <- c("A01", "B01", "C01", "D01", "G04")

# Set the arrangement for the plots. The first column contains the amplitude 
# plots, column two the density functions and column three the concentration
# calculated on according to the droplet volume as defined in the QX100 system,
# or the method proposed by Corbisier et al. (2015).

par(mfrow = c(5,3))

for (i in 1L:length(wells)) {
  cluster.info <- unique(pds_raw[wells[i]][[1]]["Cluster"])
  res <- bioamp(data = pds_raw[wells[i]][[1]], amp_x = 2, amp_y = 1, 
		main = paste("Well", wells[i]), xlab = "Amplitude of ileS (FAM)",
		ylab = "Amplitude of styA (HEX)", xlim = c(500,5500), 
		ylim = c(0,3000), pch = 19)
  # Draw threshold line to visualize between positive and negative droplts.
  abline(h = max(with(pds_raw[wells[i]][[1]], 
		 subset(Assay1.Amplitude, Cluster == 4))), lty = 2)
  legend("topleft", as.character(cluster.info[, 1]), col = cluster.info[, 1], 
	 ncol = 4, pch = 19)
  
  k.tmp <- res[1, "Cluster.3"]
  n.tmp <- sum(res[1, ])
  dens <- dpcr_density(k = k.tmp, n = n.tmp, 
			average = TRUE, methods = "wilson")
  legend("topleft", paste("k:", k.tmp,"\nn:", n.tmp))
  res.conc <- rbind(original = dens[4:6] /  0.90072 * 1e-6, 
		    corrected = dens[4:6] / 0.834 * 1e-6)
  barplot(res.conc[, 1], col = c("white","grey"), 
	  names = c("Bio-Rad", "Corbisier"), 
	  main = "Influence of\nDroplet size", ylab = "molecules/ml", ylim = c(0,8*10E-9))
    arrows(c(0.7,1.9), res.conc[, 2], c(0.7,1.9), res.conc[, 3], angle = 90, 
	   code = 3, lwd = 2)
}
dev.off()

