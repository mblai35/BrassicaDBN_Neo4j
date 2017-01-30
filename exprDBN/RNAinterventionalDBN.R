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

#avgRNA <- genes[, mean(genes[, 1:4583]), by = .(Timepoint, Treatment)]

# Create data.frame in appropriate format: Cell.line, Inhibitor, 
# Stimuli, Timepoint, and expression values. 
rnaNet <- data.frame(Cell.line = rep(1, 48), 
                       Inhibitor = genes$Treatment, 
                       Stimuli = rep("WW", 48), 
                       Timepoint = genes$Timepoint, 
                       genes[, 1:230])

# Define baseline as WW and inhibited as Dry.
droughtEffects <- interventionEffects(rnaNet, 1, "WW", "Dry")

# Update phenoNet to the correct format using the interventionalDBN
# formatData function. 
Net <- formatData(rnaNet)

# Create an empty n by P (# of obs by # of cols of expression data) 
# Z matrix.
Z <- matrix(0, 48, 230)


# Perform inference for DBN.
Network <- interventionalInference(Net$y,Net$X0, 
                                    Net$X1, Z, max.indeg = 3,
                                    perfectOut = T, fixedEffectOut = T)

