# Gene_Neo4j_visNetwork.R
# R version 3.3.1 (2016-06-21)
# February 5, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Create interactive visualization for  pheno network for 
# Brassica data using bnlearn package. Data taken from Brassica control
# and droughted conditions. 

#-----------------------------------------------------------------------
library(igraph)
library(visNetwork)
library(RNeo4j)
#-----------------------------------------------------------------------

# Convert to the proper format for importing: 

  # Read in relationship csv from simone output. 
  aicNet1 <- read.csv(file.choose(), row.names = 1)
  
  # Read in theta values for relationships from csv from simone output.
  aicNet1Theta <- read.csv(file.choose(), row.names = 1)
  
  # Find indices for 1's in sparse matrix that indicate relationships. 
  indexRel <- which(aicNet1 != 0, arr.ind = T)
  
  # Create data.frame for Genes for importation. 
  relGenes <- data.frame(From = rownames(indexRel), 
                         To = colnames(aicNet1)[indexRel[, 2]])
  
  # Add the R package used to find these relationships as a column. 
  relGenes$Method <- rep("simone", dim(relGenes)[1])
  
  # Add theta values to data.frame. 
  indexTheta <- which(aicNet1Theta != 0, arr.ind = T)
  
  # Grab values of theta using the index values. 
  Theta <- as.matrix(aicNet1Theta[indexTheta[, 1], indexTheta[, 2]])
  
  # Append to relGenes data.frame. 
  relGenes$Theta <- diag(Theta)
  
  # Write csv to file. 
  write.csv(relGenes, file = "AICsimoneDBN.csv")
  
  

# Start Neo4j graph. 
graph = startGraph("http://localhost:7474/LinaDB/data/", username = "neo4j",
                   password = "plantanalytics")




