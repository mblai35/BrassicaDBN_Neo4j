# discPhenoBN.R
# R version 3.3.1 (2016-06-21)
# February 3, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating pheno network for Brassica data using 
# bnlearn package. Data taken from Brassica control
# and droughted conditions. 

#-----------------------------------------------------------------------
library(bnlearn)
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
  discTP8$INT <- TP8$Treatment
  
  # Let structure learning algorithm decide which arcs 
  # connect to INT.  
  
  # Create list for blacklist.
  tiers <- list("INT", names(discTP8)[1:6])
  
  # Create blacklist.
  bl <- tiers2blacklist(nodes = tiers)
  
  # Search for network.
  bn.t8 <- tabu(discTP8, blacklist = bl, score = "bde",
                iss = 10, tabu = 50)
  
  # Plot network. 
  plot(bn.t8, main = "TP8 Learned Interventions CO v DR")

# Subset Timepoint 9. 
TP9 <- Pheno[Pheno$Timepoint == 9, ]

  # Remove Timepoint column. 
  TP9 <- TP9[, -2]
  
  # Soil Moisture level of 0 seems like an error. 
  # Replace with lowest value. 
  TP9[which(TP9$SM == 0), "SM"] <- NA
  TP9[which(is.na(TP9$SM)), "SM"] <- min(TP9$SM, na.rm = T)
  
  
  #### Discretize data. 
  
  # Discretize the data. 
  discTP9 <- discretize(TP9[, 2:7], method = "interval",
                             breaks = c(5, 5, 3, 5, 3, 5))
  
  # Add INT column. 
  discTP9$INT <- TP9$Treatment
  
  # Let structure learning algorithm decide which arcs 
  # connect to INT.  

  # Create list for blacklist.
  tiers <- list("INT", names(discTP9)[1:6])
  
  # Create blacklist.
  bl <- tiers2blacklist(nodes = tiers)
  
  # Search for network.
  bn.t9 <- tabu(discTP9, blacklist = bl, score = "bde",
                   iss = 10, tabu = 50)
  
  # Plot network. 
  plot(bn.t9, main = "TP9 Learned Interventions CO v DR")


# Subset Timepoint 10. 
TP10 <- Pheno[Pheno$Timepoint == 10, ]
  
  # Remove Timepoint column. 
  TP10 <- TP10[, -2]

  #### Discretize data. 
  
  # Discretize the data. 
  discTP10 <- discretize(TP10[, 2:7], method = "interval",
                        breaks = c(5, 5, 3, 5, 3, 5))
  
  # Add INT column. 
  discTP10$INT <- TP10$Treatment
  
  # Let structure learning algorithm decide which arcs 
  # connect to INT.  
  
  # Create list for blacklist.
  tiers <- list("INT", names(discTP10)[1:6])
  
  # Create blacklist.
  bl <- tiers2blacklist(nodes = tiers)
  
  # Search for network.
  bn.t10 <- tabu(discTP10, blacklist = bl, score = "bde",
                   iss = 10, tabu = 50)
  
  # Plot network. 
  plot(bn.t10, main = "TP10 Learned Interventions CO v DR")
  
# Subset Timepoint 11. 
TP11 <- Pheno[Pheno$Timepoint == 11, ]
  
  # Remove Timepoint and Fv.Fm. columns. 
  TP11 <- TP11[, -c(2, 5)]
    
  #### Discretize data. 
  
  # Discretize the data. 
  discTP11 <- discretize(TP11[, 2:6], method = "interval",
                         breaks = c(5, 5, 5, 3, 5))
  
  # Add INT column. 
  discTP11$INT <- TP11$Treatment
  
  # Let structure learning algorithm decide which arcs 
  # connect to INT.  
  
  # Create list for blacklist.
  tiers <- list("INT", names(discTP11)[1:5])
  
  # Create blacklist.
  bl <- tiers2blacklist(nodes = tiers)
  
  # Search for network.
  bn.t11 <- tabu(discTP11, blacklist = bl, score = "bde",
                 iss = 10, tabu = 50)
  
  # Plot network. 
  plot(bn.t11, main = "TP11 Learned Interventions CO v DR")
  
# Subset Timepoint 12. 
TP12 <- Pheno[Pheno$Timepoint == 12, ]
  
  # Remove Timepoint and Fv.Fm. columns. 
  TP12 <- TP12[, -c(2, 5)]
  
  #### Discretize data. 
  
  # Discretize the data. 
  discTP12 <- discretize(TP12[, 2:6], method = "interval",
                         breaks = c(5, 5, 5, 3, 5))
  
  # Add INT column. 
  discTP12$INT <- TP12$Treatment
  
  # Let structure learning algorithm decide which arcs 
  # connect to INT.  
  
  # Create list for blacklist.
  tiers <- list("INT", names(discTP12)[1:5])
  
  # Create blacklist.
  bl <- tiers2blacklist(nodes = tiers)
  
  # Search for network.
  bn.t12 <- tabu(discTP12, blacklist = bl, score = "bde",
                 iss = 10, tabu = 50)
  
  # Plot network. 
  plot(bn.t12, main = "TP12 Learned Interventions CO v DR")
  
# Create pdf of plots. 
  pdf("TP8_12BN.pdf")  
  plot(bn.t8, main = "TP8 Learned Interventions CO v DR")
  plot(bn.t9, main = "TP9 Learned Interventions CO v DR")
  plot(bn.t10, main = "TP10 Learned Interventions CO v DR")
  plot(bn.t11, main = "TP11 Learned Interventions CO v DR")
  plot(bn.t12, main = "TP12 Learned Interventions CO v DR")
  dev.off()

# Write csv files for arcs.
  write.csv(bn.t8.arcs, file = "PhenoT8arcs.csv")
  write.csv(bn.t9$arcs, file = "PhenoT9arcs.csv")
  write.csv(bn.t10$arcs, file = "PhenoT10arcs.csv")
  write.csv(bn.t11$arcs, file = "PhenoT11arcs.csv")
  write.csv(bn.t12$arcs, file = "PhenoT12arcs.csv")

#### Try making a prediction. 
  
  # Fit parameters of Bayesian Network conditional on structure. 
  fit <- bn.fit(bn.t8, discTP8, method = "bayes")
  
  # Fit parameters of Bayesian Network for TP11, 
  # conditioned on structure.
  fit2 <- bn.fit(bn.t11, discTP11, method = "bayes")
  