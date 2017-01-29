# PhenoCombinedDataImputation.R
# R version 3.3.1 (2016-06-21)
# January 15, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Investigate accuracy of current values, remove extreme outliers, and
# impute missing Brassica phenotype data for analysis.

#-----------------------------------------------------------------------
library(CoImp)
#-----------------------------------------------------------------------

# Read in combined phenotype file. 
Pheno <- read.csv(file = "/Users/mblai/Documents/PhenotypeBrassica.csv")

# Investigate the values for Photo. 
boxplot(Pheno$Photo ~ Pheno$Timepoint)

# Photosynthesis should be negative at night, meaning that respiration
# is occurring. It seems appropriate to impute Timepoint 6 by sampling
# values from Timepoint 5 and 12 by Treatment. 
# Sample for WW Treatment. 
w <- sample(c(Pheno[Pheno$Timepoint == 5 & Pheno$Treatment == "WW", 
                    "Photo"],
              Pheno[Pheno$Timepoint == 12 & Pheno$Treatment == "WW", 
                    "Photo"]), 12, replace = T)

# Update data.frame with sampled values.    
Pheno[Pheno$Timepoint == 6 & Pheno$Treatment == "WW", "Photo"] <- w

# Repeat for Dry Treatment.
d <- sample(c(Pheno[Pheno$Timepoint == 5 & Pheno$Treatment == "Dry", 
                    "Photo"],
              Pheno[Pheno$Timepoint == 12 & Pheno$Treatment == "Dry", 
                    "Photo"]), 12, replace = T)
Pheno[Pheno$Timepoint == 6 & Pheno$Treatment == "Dry", "Photo"] <- d

# Re-investigate the values for Photo. 
boxplot(Pheno$Photo ~ Pheno$Timepoint)
# Seems appropriate.

# Investigate the values for gs. 
boxplot(Pheno$gs ~ Pheno$Timepoint)

# Repeat the same process as above for gs. 
# Sample for WW Treatment. 
w <- sample(c(Pheno[Pheno$Timepoint == 5 & Pheno$Treatment == "WW", 
                    "gs"],
              Pheno[Pheno$Timepoint == 12 & Pheno$Treatment == "WW", 
                    "gs"]), 12, replace = T)

# Update data.frame with sampled values.    
Pheno[Pheno$Timepoint == 6 & Pheno$Treatment == "WW", "gs"] <- w

# Repeat for Dry Treatment.
d <- sample(c(Pheno[Pheno$Timepoint == 5 & Pheno$Treatment == "Dry", 
                    "gs"],
              Pheno[Pheno$Timepoint == 12 & Pheno$Treatment == "Dry", 
                    "gs"]), 12, replace = T)
Pheno[Pheno$Timepoint == 6 & Pheno$Treatment == "Dry", "gs"] <- d

# Re-investigate the values for Photo. 
boxplot(Pheno$gs ~ Pheno$Timepoint)
# Seems appropriate.

# Investigate the values for Fv'Fm'. 
boxplot(Pheno$Fv.Fm. ~ Pheno$Timepoint)
# Nothing looks too out of the ordinary here. 

# Investigate the values for Starch. 
boxplot(Pheno$Starch ~ Pheno$Timepoint)

# Looks like there's an extreme outlier at Timepoint 6. Remove this 
# value and replace with NA. 
Pheno[which(Pheno[, "Starch"] == max(Pheno[Pheno$Timepoint == 6, 
                                           "Starch"])), "Starch"] <- NA

# Investigate the values for NSC. 
boxplot(Pheno$NSC ~ Pheno$Timepoint)
# Nothing looks too out of the ordinary here. 

# Convert to matrix for CoImp function. 
mPheno <- as.matrix(Pheno[, c(4:9)])

# Store CoImp function results in missingPheno. 
missingPheno <- CoImp(mPheno, type.data = "continuous")

# Create data.frame containing imputed values. 
PhenoCombined <- data.frame(Treatment = Pheno$Treatment, 
                            Timepoint = Pheno$Timepoint,
                            missingPheno@Imputed.data.matrix)

# Investigate the appropriateness of imputed values for Starch. 
boxplot(PhenoCombined$Starch ~ PhenoCombined$Timepoint)
# Seems reasonable. 

# Investigate the appropriateness of imputed values for NSC. 
boxplot(PhenoCombined$NSC ~ PhenoCombined$Timepoint)
# Seems reasonable. 

# Investigate the appropriateness of imputed values for Soil Moisture. 
boxplot(PhenoCombined$SM... ~ PhenoCombined$Timepoint)
# Seems reasonable. 

# Write to csv file. 
write.csv(PhenoCombined, file = "PhenoBrassicaImp.csv")

