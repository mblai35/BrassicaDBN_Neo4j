# DEexprSimoneNetwork.R
# R version 3.3.1 (2016-06-21)
# January 27, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating transcriptomic network for Brassica data using 
# simone package. Expression data taken from Brassica under well-watered
# and droughted conditions. 

#-----------------------------------------------------------------------
library(simone)
#-----------------------------------------------------------------------

#### Pre-processing: 

# Read in the data. 
Brassica <- read.csv(file = "BrassicaDEgenes.csv", row.names = 1)

# Remove cluster from Brassica. 
Brassica <- Brassica[, 1:48]

# Create vector of Treatment values. 
Tmnt <- as.factor(c(rep("Dry", 24), rep("WW", 24)))

# Combine expression data and treatment vector in a list. 
BrassicaExpr <- list(expr = data.frame(Brassica), status = Tmnt)

#### Network analysis. 

# Run multiple iterations. 

for (i in 5:6)
{
  # Create folder name using iteration number. 
  folder <- paste("/simoneDBN", i, sep = "")
  
  # Create folder for storing output. 
  dir.create(paste(getwd(), folder, sep = ""))
  setwd(paste(getwd(), folder, sep = ""))
  
  # Set options to cluster using BIC. 
  ctrl <- setOptions(clusters.crit = "BIC")
  DBNsimoneBICcl <- simone(BrassicaExpr$expr, type = "time-course", 
                           tasks = Tmnt,
                           clustering = T, control = ctrl)
  
  # Use BIC to get network. 
  PhenoNet <- getNetwork(DBNsimoneBICcl, "BIC")
  
  # Write A and Theta to csv. 
  write.csv(PhenoNet[[1]]$A, file = 
              paste(getwd(), "/PhenoNetBIC1A.csv", 
                    sep = ""))
  write.csv(PhenoNet[[2]]$A, file = 
              paste(getwd(), "/PhenoNetBIC2A.csv", 
                    sep = ""))
  write.csv(PhenoNet[[1]]$Theta, file = 
              paste(getwd(), "/PhenoNetBIC1Theta.csv", 
                    sep = ""))
  write.csv(PhenoNet[[2]]$Theta, file = 
              paste(getwd(), "/PhenoNetBIC2Theta.csv", 
                    sep = ""))
  
  # Plot BIC networks. 
  pdf(file = paste(getwd(), "/BICplots.pdf", sep = ""))
  plot(PhenoNet[[1]])
  plot(PhenoNet[[2]])
  plot(PhenoNet[[1]], PhenoNet[[2]], type = "overlap")
  dev.off()
  
  # Set options to cluster using AIC. 
  ctrl <- setOptions(clusters.crit = "AIC")
  DBNsimoneAICcl <- simone(BrassicaExpr$expr, type = "time-course", 
                           tasks = Tmnt,
                           clustering = T, control = ctrl)
  
  # Use AIC to get network. 
  PhenoNet <- getNetwork(DBNsimoneAICcl, "AIC")
  
  # Write A and Theta to csv. 
  write.csv(PhenoNet[[1]]$A, file = 
              paste(getwd(), "/PhenoNetAIC1A.csv", 
                    sep = ""))
  write.csv(PhenoNet[[2]]$A, file = 
              paste(getwd(), "/PhenoNetAIC2A.csv", 
                    sep = ""))
  write.csv(PhenoNet[[1]]$Theta, file = 
              paste(getwd(), "/PhenoNetAIC1Theta.csv", 
                    sep = ""))
  write.csv(PhenoNet[[2]]$Theta, file = 
              paste(getwd(), "/PhenoNetAIC2Theta.csv", 
                    sep = ""))
  
  # Plot AIC networks. 
  pdf(file = paste(getwd(), "/AICplots.pdf", sep = ""))
  plot(PhenoNet[[1]])
  plot(PhenoNet[[2]])
  plot(PhenoNet[[1]], PhenoNet[[2]], type = "overlap")
  dev.off()
  
  # Use default settings. 
  DBNsimone <- simone(BrassicaExpr$expr, type = "time-course", 
                           tasks = Tmnt)
  
  # Use default settings to get network. 
  PhenoNet <- getNetwork(DBNsimone)
  
  # Write A and Theta to csv. 
  write.csv(PhenoNet[[1]]$A, file = 
              paste(getwd(), "/PhenoNet1A.csv", 
                    sep = ""))
  write.csv(PhenoNet[[2]]$A, file = 
              paste(getwd(), "/PhenoNet2A.csv", 
                    sep = ""))
  write.csv(PhenoNet[[1]]$Theta, file = 
              paste(getwd(), "/PhenoNet1Theta.csv", 
                    sep = ""))
  write.csv(PhenoNet[[2]]$Theta, file = 
              paste(getwd(), "/PhenoNet2Theta.csv", 
                    sep = ""))
  
  # Plot networks. 
  pdf(file = paste(getwd(), "/plots.pdf", sep = ""))
  plot(PhenoNet[[1]])
  plot(PhenoNet[[2]])
  plot(PhenoNet[[1]], PhenoNet[[2]], type = "overlap")
  dev.off()
  
  # Write clusters to csv. 
  #write.csv(DBNsimoneBICcl$clusters, file = 
  #            paste(getwd(), "/clusters.csv", sep = ""))
  
  # Plot DBNsimoneBICcl.
  #pdf("DBNsimoneplots.pdf")
  #plot(DBNsimone)
  #dev.off()
  
  setwd("/Users/mblai/Documents")

}


