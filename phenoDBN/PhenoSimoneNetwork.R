# PhenoSimoneNetwork.R
# R version 3.3.1 (2016-06-21)
# January 15, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating phenotypic network for Brassica data using simone package. 

#-----------------------------------------------------------------------
library(simone)
library(data.table)
#-----------------------------------------------------------------------

# Read in phenotype file. 
Pheno <- read.csv(file = "/Users/mblai/Documents/PhenoBrassicaImp.csv")

# Convert to data.table for ease of manipulation. 
Pheno <- as.data.table(Pheno)

# Average replicates across treatment and timepoints. 
avgPheno <- Pheno[, .(Photo = mean(Photo), gs = mean(gs), 
             FvFm = mean(Fv.Fm.), Starch = mean(Starch),
             NSC = mean(NSC), SM = mean(SM...)), 
         by = .(Timepoint, Treatment)]

# Save Treatment column as a vector of factors for simone.
Tmnt <- as.factor(avgPheno$Treatment)

# Remove columns that aren't phenotypic observations. 
avgPheno <- avgPheno[, 3:8]

# Combine phenotypic observations with Treatment vector as a list for 
# simone. 
PhenoList <- list(expr = avgPheno, status = Tmnt)

# Run simone network analysis. 
PhenoDBN <- simone(PhenoList$expr, type = "time-course", 
               tasks = Tmnt)

# Plot PhenoDBN.
plot(PhenoDBN)

# Use AIC to get network. 
PhenoNet <- getNetwork(PhenoDBN, "AIC")

# Plot AIC networks. 
plot(PhenoNet[[1]])
plot(PhenoNet[[2]])

# Use BIC to get network. 
PhenoNet <- getNetwork(PhenoDBN, "BIC")

# Plot BIC networks. 
plot(PhenoNet[[1]])
plot(PhenoNet[[2]])

# Use default settings to get network. 
PhenoNet <- getNetwork(PhenoDBN)

# Plot networks. 
plot(PhenoNet[[1]])
plot(PhenoNet[[2]])

