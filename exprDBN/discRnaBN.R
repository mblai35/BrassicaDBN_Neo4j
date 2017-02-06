# discRnaBN.R
# R version 3.3.1 (2016-06-21)
# February 5, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating genetic network for Brassica data using 
# bnlearn package. Data taken from Brassica control
# and droughted conditions. 

#-----------------------------------------------------------------------
library(bnlearn)
library(stringr)
#-----------------------------------------------------------------------

#### Preprocessing: 

# Read in phenotype file. 
rna <- read.csv(file.choose(), row.names = 1)

# Remove cluster column. 
rna <- rna[, -49]

# Transpose data. 
rna <- t(rna)

# Convert to data.frame. 
rna <- as.data.frame(rna)
          
# Discretize rna data. 
          
  # Discretize the data. 
  discRNA <- discretize(rna, method = "interval", breaks = 5)
          
# Split rownames into Treatment/Timepoint and Replicate. 
rnaNames <- str_split_fixed(rownames(Rna), '_', 2)

# Create a Timepoint column. 
discRNA$Timepoint <- as.numeric(str_split_fixed(rnaNames[, 1], '', 2)[, 2])
  
# Create a treatment column named INT. 
discRNA$INT <- as.factor(str_split_fixed(rnaNames[, 1], '', 2)[, 1])
  
# Subset only Timepoint 8. 
TP8 <- discRNA[discRNA$Timepoint == 8, ]

  # Remove Timepoint column from TP8. 
  TP8 <- TP8[, -231]

  # Create list for blacklist.
  tiers <- list("INT", names(TP8)[1:230])
  
  # Create blacklist.
  bl <- tiers2blacklist(nodes = tiers)
  
  # Search for network.
  bn.t8 <- tabu(TP8, blacklist = bl, score = "bde",
                iss = 10, tabu = 50)
  
  # Plot network. 
  plot(bn.t8, main = "TP8 Learned Interventions CO v DR")
  
  
  
  
  
  
  
  
  
