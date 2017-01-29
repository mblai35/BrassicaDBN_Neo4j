# DEexprGRENITS.R
# R version 3.3.1 (2016-06-21)
# January 27, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating transcriptomic network for Brassica data using 
# simone package. Expression data taken from Brassica under well-watered
# and droughted conditions. 

#-----------------------------------------------------------------------
library(GRENITS)
library(dplyr)
#-----------------------------------------------------------------------

#### Pre-processing: 

# Read in the data. 
Brassica <- read.csv(file = "BrassicaDEgenes.csv", row.names = 1)

# Manipulate data so rep2 set follows rep1 set.
Brassica <- data.frame(select(Brassica, -ends_with("rep2")), 
                      select(Brassica, -ends_with("rep1")))

#### Network Inference:

# Create tempdir folder to hold MCMC files in output.folder. 
resultsFolder <- paste(getwd(), "/DEexpGRENTIS_gauss", sep = "")

# Run MCMC function; 2 chains, default parameters
ReplicatesNet_gauss(resultsFolder, Brassica, numReps = 2)

# Analys raw results, place analysis plots and files in resultsFolder
analyse.output(resultsFolder)

# View contents of output.folder
dir(resultsFolder)

# Load inferred network probabilities.
# Create prob.file to store full path name for 
# NetworkProbability_Matrix.txt.
prob.file <- paste(resultsFolder, "/NetworkProbability_Matrix.txt", 
                   sep = "")

# Read in the NetworkProbability_Matrix.txt to R as a data.frame.
prob.mat <- read.table(prob.file)

# Apply a threshold of 0.8 on the link probability for a network
# prediction. 
inferred.net <- 1 * (prob.mat > 0.8)

# View inferred network. 
write.csv(inferred.net, file = 
            paste(getwd(), "/DEexpGRENTIS_gauss/inferred.csv", sep = ""))

# For information on direction of the interaction, see 
# NetworkProbability_List.txt
prob.list.file <- paste(resultsFolder, "/NetworkProbability_List.txt",
                        sep = "")

# Read in NetworkProbability_List.txt.
prob.list <- read.table(prob.list.file, header = T)

# Find genes with a probability greater than .8 from the Probability
# column (3rd column).
above.08 <- (prob.list[, 3] > 0.8)

# Print out genes where the probability is greater than .8.
write.csv((prob.list[above.08, ]), file = 
            paste(getwd(), "/DEexpGRENTIS_gauss/prob_above08.csv", sep = ""))

