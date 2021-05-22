# Isotope Dendrochronology
# Calculation and interpretation of physiological indices derived from d13C
# COURSE: DISC, LTRR, UNIVERSITY OF ARIZONA
# INSTRUCTORS: S. BELMECHERI & P. SZEJNER

# The following code and comments should further your experience to use tree ring and atmospheric d13C data to calculate WUE and other
# physiological indices using R. This exercise relies upon raw tree-ring d13C measurements, atmospheric d13C measurements, and uses
# the Farqhuar 1982 model.

# There are a few questions throughout to catalyze understanding and critical thinking.
#
# The analysis is broken into 4 main parts - at the end you should have
# gone through many of the steps commonly used for reconstructing and interpreting iWUE derived from tree ring d13C
# You will also be familiar with uncertainties related to such calculations and best practices to perform them.

###################################################################################################################
#
# PART 0: Setting up directories and files.
#
#
###################################################################################################################

# For the following exercises, you will need  the following files in one directory.
# 1. A file containing tree-ring d13C measurements. The filename is "LIL.txt"
# 2. A file containing atmospheric 13C data. The filename is "ATM.txt"
# 3. A file containing regional temperature data. The filename is "CRU_tt.txt"

# You will also need to install and load the dplR package, and the rpmodel Correction Package
library(dplR)
library(rpmodel)


# You will need to change the directory in R so that this is the working directory.
# The command "setwd" is used for this as in the following line. 
setwd("~/Documents/DISC/NOTEBOOK3")


# You should see the above file names (and anything else in this directory) when you 
# type (or copy and paste) the following line:
list.files()

#If you see these files listed, you're ready to roll!

###################################################################################################################
#
# PART 1: Importing and plotting raw tree-ring isotope data [d13C] & Atmospheric Data
#
#
###################################################################################################################

# import the raw tree ring d13C dataset
TR13 <- ts((read.table("LIL.txt", header=T)), start=1901, frequency=1)

#plot the raw data
ts.plot(TR13.ts,ylab="d13C (permil.VPDB)",xlab="Years", col="#9C964A")

# import the atmospheric d13C and CO2 dataset
ATM <- ts((read.table("ATM.txt", header=T)), start=)

# subset d13C
ATM.13C <- ts(ATM[,"atm13C"],start=1901, frequency=1)
# subset CO2
ATM.CO2 <- ts(ATM[,"atmCO2"],start=1901, frequency=1)


# plot the raw tree ring d13C dataset and the atmospheric d13C dataset
layout(matrix(1:2,nrow = 1, ncol = 2,byrow = T), widths = c(3,3),heights = c(3,3),respect = TRUE)
ts.plot(ATM.13C,ylab="d13C (permil.VPDB)",xlab="Years", col="#9C964A",main="d13C")
ts.plot(ATM.CO2,ylab="CO2 (ppm)",xlab="Years", main="CO2")

# Question: compare both trends. Describe and interpret the relationship between
#           the decline in d13C and the increase in CO2 concentration. 

###################################################################################################################
#
# PART 2: Calculating iWUE- SIMPLE APPROACH
#
#
###################################################################################################################

# the Farquhar 1982 biochemical model  describes the isotopic discrimination against 13C (Δ13C) 
# during carbon diffusion and fixation by plants.

# The Farquhar equation can be written as follows:
# Δ13C = a+(b-a)*(ci/ca)

a=4.4 # the fractionation due to CO2 diffusion in air through the stomata.
b= 28 # the apparent net fractionation by RuBisCO during carboxylation.

# ci and ca are are the leaf intercellular and ambient partial pressure of CO2 (Pa), respectively
# ci <-  ca*(Δ13C-a)/(b-a)

# We know ca from measurements.ca is in ppm. you use use it as is 
ca <- ATM.CO2 # in ppm
# or convert to Pa. For this you need a r package "p model"
library("rpmodel")
# and the function calc_patm {rpmodel} to calculate atmospheric pressure at sea level. The variable needed is elevation.
elv= 1600# Elevation above sea-level (m.a.s.l.)
patm<- calc_patm(elv, patm0 = 101325)
ca <- ( 1.e-6 ) * ca * patm  # in Pa
# convert to a TS
ca <- ts(ca, start=1901, frequency=1)

# WARNING: comment the ca (ppm or Pa) you will not use. Preferably, use ca in Pa

# We can calculate Δ13C as follows:
d13C.disc <- (ATM.13C-TR13)/(1+TR13/1000)

# now calculate ci as follows:
ci <-  ca*(d13C.disc-a)/(b-a)

# now calculate iWUE (umol.mol)
iwue <- (ca-ci)/1.6

# plot the physiological indices

layout(matrix(1:6,nrow = 3, ncol = 2,byrow = F), widths = c(3,3,3),heights = c(3,3,3),respect = TRUE)

ts.plot(TR13,ylab="d13C (permil.VPDB)",xlab="Years", main="OBSERVATIONS")
ts.plot(d13C.disc,ylab="D13C (permil.VPDB)",xlab="Years", main="")
ts.plot(ca,ylab="CO2 (pa)",xlab="Years", main="")
ts.plot(ci,ylab="ci (pa)",xlab="Years", main="FARQUHAR-MODEL")
ts.plot(iwue,ylab=" iwue (umol/mol)",xlab="Years", main="")
ts.plot(ci/ca,ylab="ci/ca",xlab="Years", main="")


# Question: Describe the trends of the various physiological indices
# Question: Based on what you have learned in the Isotope theory lectures,
#           describe the trees'physiological response over time at this site.

###################################################################################################################
#
# PART 3: Calculating iWUE- COMPREHENSIVE APPROACH
#
#
###################################################################################################################


# In the Farquhar model, the d13C is the isotopic ratio for sugars fixed in the leaf. 
# For tree-ring cellulose, a correction factor is used to account for the offset of d13C between leaf sugars and tree-ring cellulose. 
# The offset results from post-photosynthetic fractionation processes (Gessler et al., 2014). 
# The offset between whole wood and leaves is ~1.3 ± 0.2‰ for oak and conifer species 
# Additional isotopic offsets include the difference between tree-ring cellulose and bulk wood with an average value of ~1.3 ± 0.2‰, 
# and isotopic depletion between primary assimilates and bulk leaf with values of -0.5 ± 1‰. 


#  Upscaling d13C measurements from tree rings to the leaf level can significantly improve
#  estimates of Δ13C and hence reduce uncertainties in the determination of the ratio of ci to ca.

# we cam scale tree ring d13 to leaf level using a factor d.
# d represents the sum of post-photosynthetic isotope fractionations between
# the leaf organic matter and the plant material considered.
d=1.9 # for d13C measured in wood
d=2.1 # for d13C measured in cellulose

# This correction can be made as follows:

d13C.disc.leaf <- (ATM.13C-(TR13-d))/(1+(TR13-d)/1000)
              
# now calculate ci:
ci.leaf <-  ca*(d13C.disc.leaf-a)/(b-a)
# and iWUE (umol.mol)
iwue.leaf <- (ca-ci.leaf)/1.6

# now add the leaf level physiological indices to the previous plot

layout(matrix(1:6,nrow = 3, ncol = 2,byrow = F), widths = c(3,3,3),heights = c(3,3,3),respect = TRUE)
# Note you might need to adjust the ylim values.
ts.plot(TR13,ylab="d13C (permil.VPDB)",xlab="Years", main="OBSERVATIONS")
ts.plot(d13C.disc,ylab="D13C (permil.VPDB)",xlab="Years", main=""); lines(d13C.disc.leaf, col="#0B775E")
ts.plot(ca,ylab="CO2 (ppm)",xlab="Years", main="")
ts.plot(ci,ylab="ci (pa)",xlab="Years", main="FARQUHAR-MODEL"); lines(ci.leaf, col="#0B775E")
ts.plot(iwue,ylab=" iwue (umol/mol)",xlab="Years", main=""); lines(iwue.leaf, col="#0B775E")
ts.plot(ci/ca,ylab="ci/ca",xlab="Years", main=""); lines(ci.leaf/ca, col="#0B775E")

# Question: What differences do you observe between physiological indices
#           at stem level and leaf level?


# The version of the Farquhar model we used above is a simplified description of the isotopic discrimination.
# It does not include the fractionation effects during the transfer of CO2 from substomatal cavities
# to the site of fixation via the mesophyll, and during mitochondrial respiration and photorespiration
# While the mesophyll fractionation and their values remain highly unconstrained, Recent studies have recommended 
# the inclusion of the photorespiratory effect in the discrimination model as this term contributes to an increase of Δ13C 
# with ca rise by 0.004‰ ppm-1 (Keeling et al., 2017; Lavergne et al., 2019).
# 
# A Farquhar model that includes the photorespiratory term is as follows:
# Δ13C = a+(b-a)*(ci/ca) -f*(Gst/ca)
# The photorespiratory term is the: f*(Gst/ca)

f=12 # the fractionation due to photorespiration.

# Gst (or Γ*, Gamma Star) is the CO2 compensation point in the absence of mitochondrial respiration (Pa), 
# it is calculated from the temperature (T) and atmospheric pressure response.

# you can use the calc_gammastar(tc, patm) function from the rpmodel to calculate Gamma Star
# It depends on temperature (tc) and atmospheric pressure (patm)

# Based on Notebook 1 and the temperature file, define the appropriate temperature window (e.g.average growing season).
# This has to be a Time series of site level temperature variability during the summer or growing season 
# or relevant season for trees at this site.
# below is an example for summer (June-July-August)
gridded.tt.data <- ts(read.table("CRU_TT.txt",header=T),start=1901,frequency=1)
tt.JJA <- rowMeans(gridded.tt.data[,c(6,7,8)])
tt.JJA <- ts(tt.Mar_Jun,start=1901,frequency =1)
tc <- tt.JJA

# now you can calculate gamma start as follows:
Gst <- calc_gammastar(tc, patm)
# in Help, you can look at the details of this calculation (Full equation)
# and the photorespiratory term as follows:
photoresp_term <- 12*(Gst/ca) # ca is  in Pa

 
# now calculate ci to include photorespiration:
ci.leaf.photo <-  ca*(d13C.disc.leaf-a+photoresp_term)/(b-a)
# and iWUE (umol.mol)
iwue.leaf.photo <- (ca-ci.leaf.photo)/1.6

# now plot the leaf level physiological parameters using the simple and photorespiration models:
layout(matrix(1:6,nrow = 3, ncol = 2,byrow = F), widths = c(3,3,3),heights = c(3,3,3),respect = TRUE)
# Note you might need to adjust the ylim values.
ts.plot(TR13,ylab="d13C (permil.VPDB)",xlab="Years", main="OBSERVATIONS")
ts.plot(d13C.disc.leaf,ylab="D13C (permil.VPDB)",xlab="Years", main="")
ts.plot(ca,ylab="CO2 (pa)",xlab="Years", main="")
ts.plot(ci.leaf,ylab="ci (pa)",xlab="Years", main="FARQUHAR_MODEL"); lines(ci.leaf.photo, col="#0B775E")
ts.plot(iwue.leaf,ylab=" iwue (umol/mol)",xlab="Years", main=""); lines(iwue.leaf.photo, col="#0B775E")
ts.plot(ci.leaf/ca,ylab="ci/ca",xlab="Years", main=""); lines(ci.leaf.photo/ca, col="#0B775E")
legend("bottomleft",c("SIMPLE", "PHOTO"),lwd= c(1,1), bty = "n", col=c("black","#0B775E"),
       text.col=c("black","#0B775E"),ncol=1)

# Question: Do the trends and amplitude of interannual variability vary between the simple model
#           and the one including photorespiration?
# 

###################################################################################################################
#
# PART 4: interpreting iWUE and leaf gas exchange strategies
#
#
###################################################################################################################

# From the parts above, you have now calculated physiological parameters using best practices
# to account for source data, scaling to leaf level, and most fractionation factors.

# In the following, we are going to interpret variations of this physiological parameters 
# We will focus on those physiological parameters using the model with photorespiration.


# STEP1: compute the trends of iWUE, ci and ci/ca over the period of record.
# you can use the lm function and the summary to look at and report the statistics (R2, p value).

# Here is a first example with iWUE
iwue.linear <- lm(iwue.leaf.photo ~ seq(1901,2002,by=1))
summary(iwue.linear)
#repeat these steps for ci and ci/ca

# STEP2: Plot time-series of iWUE, ci and ci/ca and add the linear trend
# Here is a first example with iWUE. You can use the layout to have 3 plots 

layout(matrix(1:3,nrow = 1, ncol = 3,byrow = T), widths = c(3,3,3),heights = c(3,3,3),respect = TRUE)
ts.plot(iwue.leaf.photo,ylab=" iwue (umol/mol)",xlab="Years", main="") 
abline(iwue.linear, col="#0B775E",lwd=2,lty=3)


# Question: Using computed statistics above, Describe the trends 
#           and their significance for each of the physiological indices?
# Question: How do you interpret the ci/ca trends in term of stomatal conductance 
# Question: Was the iWUE increase steady over time?
# Question: What was the rate of the ci (pa.year-1) increase over time? 
# Question: Was this rate steady/constant over time?
# Question: compared to the rate of ca, was the rate of ci:
#           a) proportional   b) similar
#           To answer this question, you need to estimate the rate of ca using lm and summary functions


# STEP3: Estimate the ci increase relative to ca.
# you can do this by regressing ci by ca
ci.linear <- lm(ci ~ ca)
summary(ci.linear)

# Question: Report the rate and describe it in terms of a) proportional   b) similar
# Note that this rate translate the ci increase for each ppm/pa of ca increase
# Question: Looking at all physiological indices, their trends and rates, how do you interpret
#           the temporal variations of stomatal conductance and photosynthesis of these trees?



# With rising ca, variations in plant Δ13C have been grouped into three leaf gas-exchange strategies 
# (1) constant ci, (2) constant ci/ca, and (3) constant ca − ci. 
# The first two strategies are considered active, and the third passive. 
# These theoretical considerations serve as the basis for the interpretation of
# the physiological mechanisms underlying iWUE trends.

# calculate ci following scenario 1
# first, estimate the average ci from tree ring for the first decade of the record
ci.average_1decade <- mean(ci.leaf.photo[1:10])
ci.scenario1 <- ts(rep(ci.average_1decade, 102), start=1901, frequency=1)


# calculate ci following scenario 2
# first, estimate the average ci/ca from tree ring for the first decade of the record
cica.average_1decade <- mean((ci.leaf.photo/ca)[1:10])
ci.scenario2 <-ca*cica.average_1decade

# calculate ci following scenario 3
# first, estimate the average ca-ci from tree ring for the first decade of the record
caci.average_1decade <- mean((ca-ci.leaf.photo)[1:10])
ci.scenario3 <-ca-caci.average_1decade

# now, plot the observed ci with the predicted ci following the 3 scenarios

ts.plot(ci.leaf.photo,ylab="ci (pa)",xlab="Years", main="", ylim=c(15,25))
lines(ci.scenario1, lty=2, col="#F2AD00",lwd=2)
lines(ci.scenario2, lty=6, col="#F2AD00",lwd=2)
lines(ci.scenario3, lty=3, col="#F2AD00",lwd=2)
legend("topleft",c("ci-observed", "ci constant","ci/ca constant","ca-ci constant"),
       lty= c(1,2,6,3), lwd=c(1,2,2,2),bty = "n", col=c("black","#F2AD00","#F2AD00","#F2AD00"),
       text.col=c("black","#F2AD00","#F2AD00","#F2AD00"), ncol=1)


# Question: Does the observed ci (derived from tree rings) follow any of the 3 scenarios?
#           Does the observed ci follow a scenario consistently throughout the record?

# Now, compute iwue following the three scenarios (using ci for 3 scenarios from above) and plot 
# with observed iwue.

# From Notebook 1, you assessed climate correlations of 13C and climate trends
# With the climate trends and tree ring d13C sensitivity to these climate factors
# Answer the following:

# Question: How do you interpret the change (increase) in iwue
#           is it driven by an increase is photosynthetic assimilation?
#           is it driven by a decrease is stomatal conductance?
#           are iwue, photosynthetic assimilation,stomatal conductance driven by
#           a) rising ca, b) trends in climate, c) both?
# Question: discuss how trends in ca and/or climate are affecting tree physiology
#           

