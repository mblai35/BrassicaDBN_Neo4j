# BrassicaBN.R
# R version 3.3.1 (2016-06-21)
# February 3, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating combined pheno & RNA seq BN for Brassica data timepoint 8   
# using bnlearn package. Data taken from Brassica control
# and droughted conditions. 

#-----------------------------------------------------------------------
library(bnlearn)
library(stringr)
#-----------------------------------------------------------------------

#### Preprocessing: 

# Read in phenotype file. 
Pheno <- read.csv(file.choose(), row.names = 1)

# Rename SM... to get rid of periods. 
colnames(Pheno)[8] <- "SM"

# Subset Timepoint 8. 
TP8 <- Pheno[Pheno$Timepoint == 8, ]

# Remove Timepoint column. 
TP8 <- TP8[, -2]

#### Discretize data. 

# Discretize the data. 
discTP8 <- discretize(TP8[, 2:7], method = "interval",
                      breaks = c(5, 5, 3, 5, 3, 5))

# Add INT column. 
discTP8$INT <- as.factor(TP8$Treatment)





