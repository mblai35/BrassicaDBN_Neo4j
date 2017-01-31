# RNAinterventionalDBN.R
# R version 3.3.1 (2016-06-21)
# January 22, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Create a Dynamic Bayesian Network using the R package 
# interventionalDBN and RNA seq data from well-watered and drought
# conditions. 

#-----------------------------------------------------------------------
library(interventionalDBN)
library(stringr)
library(data.table)
#-----------------------------------------------------------------------

# Read in the data. 
genes <- read.csv(file = "BrassicaDEgenes.csv", row.names = 1)

# Create directory to store output. 
dir.create(paste(getwd(), "/RinterDBN", sep = ""))
setwd(paste(getwd(), "/RinterDBN", sep = ""))

# Remove cluster column. 
genes <- genes[, -49]

# Transform genes to get in correct format for interventionalDBN.
genes <- t(genes)

# Turn into a data.table for ease of manipulation. 
genes <- as.data.table(genes, keep.rownames = T)

# Split rownames to extract Timepoint and Treatment. 
genes <- genes[, 'rn'  := str_split_fixed(genes$rn, "_", 2)[, 1]]
genes <- genes[, 'Treatment'  := str_split_fixed(genes$rn, "", 2)[, 1]]
genes <- genes[, 'Timepoint'  := str_split_fixed(genes$rn, "", 2)[, 2]]
genes <- genes[, 'rn' := NULL]

# Average replicates. 
genes <- genes[, lapply(.SD, mean), by= .(Timepoint, Treatment) ]
# NOTE: this changes data.table so that Timepoint and Treatment are the 
# first two columns instead of the last. 

# Make Timepoint an integer and Treatment a factor. 
genes$Timepoint <- as.integer(genes$Timepoint)
genes$Treatment <- as.factor(genes$Treatment)

# Create data.frame in appropriate format: Cell.line, Inhibitor, 
# Stimuli, Timepoint, and expression values. 
rnaNet <- data.frame(Cell.line = rep(1, 24), 
                       Inhibitor = genes$Treatment, 
                       Stimuli = rep("W", 24), 
                       Timepoint = genes$Timepoint, 
                       genes[, 3:(dim(genes)[2])])

# Define baseline as WW and inhibited as Dry.
droughtEffects <- interventionEffects(rnaNet, 1, "W", "D")

# Update phenoNet to the correct format using the interventionalDBN
# formatData function. 
Net <- formatData(rnaNet)

# Create an empty n by P (# of obs by # of cols of expression data) 
# Z matrix.
Z <- matrix(0, 24, 230)

# Perform inference for DBN.
Network <- interventionalInference(Net$y,Net$X0, 
                                    Net$X1, Z, max.indeg = 3,
                                    perfectOut = T, fixedEffectOut = T)

# Write output to csv file. 
write.csv(Network$pep, file = "pep.csv")
write.csv(Network$MAP, file = "MAP.csv")
write.csv(Network$MAPprob, file = "MAPprob.csv")
write.csv(Network$marginal.likelihood, 
          file = "marginalLikelihood.csv")

