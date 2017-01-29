# PhenoCombinedBrassicaDec13data.R
# R version 3.3.1 (2016-06-21)
# January 12, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Combining Brassica phenotype data for analysis

# This "data clean-up" is meant to combine separate files and get rid
# of columns and data that are not needed for analysis. 

#-----------------------------------------------------------------------
library(data.table)
#-----------------------------------------------------------------------

# Read in file 1.
Starch <-fread(file = "4Mall_NSC_Starchdec2013_edit.csv") 
Starch$NSC <- as.numeric(Starch$NSC)
Starch$Starch <- as.numeric(Starch$Starch)

# Read in file 2.
Soil <-fread(file = "4Mall_Soil_Moisture_dec_2013_edit.csv") 
Soil$Timepoint <- as.numeric(Soil$Timepoint)

# Read in file 3.
Photo <-fread(file = "4Mall_PhotoFv'Fm'gs_dec_2013_edit.csv") 

# File1 doesn't have replicate values. Since we're only interested 
# in leaf tissue so we'll want to remove any readings taken from the 
# roots. 
Starch <- Starch[Tissue == "leaves"]

# For each file, remove columns that aren't "Treatment",  
# "Timepoint", and output values. 
Starch <- Starch[, c("Treatment", "Timepoint", "NSC", "Starch")] 
Soil <- Soil[, c("Treatment", "Timepoint", "SM(%)")] 
Photo <- Photo[, c("Treatment", "Timepoint", "Photo", "gs", "Fv'Fm'")] 

# Examine the number of replicates at each Timepoint.
replicates <- data.table(Photo = Photo[ , .N, by = Timepoint], 
                Soil = Soil[ , .N, by = Timepoint], 
                Starch = Starch[ , .N, by = Timepoint])
# Note, the mismatched values have been recycled. 
# We can see that Photo is missing one timepoint and Soil and Starch
# are missing two. 
# Photo is missing eight replicates for the first three timepoints.
# Starch has an extra replicate for timepoint 4 and a missing replicate
# for timepoint 12.

# Build dataframe for combined data. 
# Start with timepoint 1, treatment "WW".
Pheno <- Photo[Timepoint == 1 & Treatment == "WW"]
Pheno <- cbind(Pheno, 
               Starch[Timepoint == 1 & Treatment == "WW", NSC, Starch], 
               Soil[Timepoint == 1 & Treatment == "WW", "SM(%)"])
# For now, missing replicates for Photo will be recycled. 

# Timepoint 1, Dry.
tp1dry <- Photo[Timepoint == 1 & Treatment == "Dry"]
tp1dry <- cbind(tp1dry,
               Starch[Timepoint == 1 & Treatment == "Dry", NSC, Starch], 
               Soil[Timepoint == 1 & Treatment == "Dry", "SM(%)"])
Pheno <- rbind(Pheno, tp1dry)
rm(tp1dry)

# Timepoint 2, WW.
tp2ww <- Photo[Timepoint == 2 & Treatment == "WW"]
tp2ww <- cbind(tp2ww, 
               Starch[Timepoint == 2 & Treatment == "WW", NSC, Starch], 
               Soil[Timepoint == 2 & Treatment == "WW", "SM(%)"])
Pheno <- rbind(Pheno, tp2ww)
rm(tp2ww)

# Timepoint 2, Dry.
tp2dry <- Photo[Timepoint == 2 & Treatment == "Dry"]
tp2dry <- cbind(tp2dry,
                Starch[Timepoint == 2 & Treatment == "Dry", NSC, Starch], 
                Soil[Timepoint == 2 & Treatment == "Dry", "SM(%)"])
Pheno <- rbind(Pheno, tp2dry)
rm(tp2dry)

# Timepoint 3, WW. Fill in missing starch timepoint with NA's.
tp3ww <- Photo[Timepoint == 3 & Treatment == "WW"]
tp3ww <- cbind(tp3ww, 
               Starch = NA, NSC = NA,
               Soil[Timepoint == 3 & Treatment == "WW", "SM(%)"])
Pheno <- rbind(Pheno, tp3ww)
rm(tp3ww)

# Timepoint 3, Dry.
tp3dry <- Photo[Timepoint == 3 & Treatment == "Dry"]
tp3dry <- cbind(tp3dry,
                Starch = NA, NSC = NA, 
                Soil[Timepoint == 3 & Treatment == "Dry", "SM(%)"])
Pheno <- rbind(Pheno, tp3dry)
rm(tp3dry)

# Timepoint 4, WW. Remove extra starch value.
# The repeated value in the first row of "WW" looks like a human error.
starch24 <- Starch[Timepoint == 4 & Treatment == "WW", NSC, Starch]
starch24 <- starch24[-1]
tp4ww <- Photo[Timepoint == 4 & Treatment == "WW"]
tp4ww <- cbind(tp4ww, 
               starch24, 
               Soil[Timepoint == 4 & Treatment == "WW", "SM(%)"])
Pheno <- rbind(Pheno, tp4ww)
rm(tp4ww)
rm(starch24)

# Timepoint 4, Dry.
tp4dry <- Photo[Timepoint == 4 & Treatment == "Dry"]
tp4dry <- cbind(tp4dry,
                Starch[Timepoint == 4 & Treatment == "Dry", NSC, Starch], 
                Soil[Timepoint == 4 & Treatment == "Dry", "SM(%)"])
Pheno <- rbind(Pheno, tp4dry)
rm(tp4dry)

# Timepoint 5. Fill in missing soil values with NA's. 
tp5 <- Photo[Timepoint == 5]
tp5 <- cbind(tp5, 
               Starch[Timepoint == 5, NSC, Starch], 
               "SM(%)" = NA)
Pheno <- rbind(Pheno, tp5)
rm(tp5)

# Timepoint 6. Fill in missing soil and Photo values with NA's. 
tp6 <- Starch[Timepoint == 6, Timepoint, Treatment]
tp6 <- cbind(tp6, Photo = NA, gs = NA, "Fv'Fm'" = 0,
               Starch[Timepoint == 6, NSC, Starch], 
               "SM(%)" = NA)
Pheno <- rbind(Pheno, tp6)
rm(tp6)

# Timepoint 7.
tp7 <- Photo[Timepoint == 7]
tp7 <- cbind(tp7,
                Starch[Timepoint == 7, NSC, Starch], 
                Soil[Timepoint == 7, "SM(%)"])
Pheno <- rbind(Pheno, tp7)
rm(tp7)

# Timepoint 8.
tp8 <- Photo[Timepoint == 8]
tp8 <- cbind(tp8,
             Starch[Timepoint == 8, NSC, Starch], 
             Soil[Timepoint == 8, "SM(%)"])
Pheno <- rbind(Pheno, tp8)
rm(tp8)

# Timepoint 9. Fill in missing starch timepoint with NA's.
tp9 <- Photo[Timepoint == 9]
tp9 <- cbind(tp9,
             Starch = NA, NSC = NA, 
             Soil[Timepoint == 9, "SM(%)"])
Pheno <- rbind(Pheno, tp9)
rm(tp9)

# Timepoint 10.
tp10 <- Photo[Timepoint == 10]
tp10 <- cbind(tp10,
             Starch[Timepoint == 10, NSC, Starch], 
             Soil[Timepoint == 10, "SM(%)"])
Pheno <- rbind(Pheno, tp10)
rm(tp10)

# Timepoint 11.
tp11 <- Photo[Timepoint == 11]
tp11 <- cbind(tp11,
             Starch[Timepoint == 11, NSC, Starch], 
             Soil[Timepoint == 11, "SM(%)"])
Pheno <- rbind(Pheno, tp11)
rm(tp11)

# Timepoint 12.
tp12 <- Photo[Timepoint == 12]
tp12 <- cbind(tp12,
             Starch[Timepoint == 12, NSC, Starch], 
             Soil[Timepoint == 12, "SM(%)"])
Pheno <- rbind(Pheno, tp12)
rm(tp12)

# Convert to data.frame. 
Pheno <- as.data.frame(Pheno)

# Remove data points containing errors. 
# Starch should never have a negative value. 
sum(Pheno$Starch < 0, na.rm = T)

# Only one value is negative. Replace negative value with NA.
Pheno[which(Pheno$Starch < 0), "Starch"] <- NA

# NSC should also never be negative. 
sum(Pheno$NSC < 0, na.rm = T)

# Only one value is negative. Replace negative value with NA.
Pheno[which(Pheno$NSC < 0), "NSC"] <- NA

# The Soil Moisture value of zero also seems to be an error. 
Pheno[which(Pheno$`SM(%)` == 0), "SM(%)"] <- NA

# Check that Photo is negative at timepoints 5, 6, 11, & 12. 
Pheno[which(Pheno$Photo > 0), "Timepoint" == 5]
Pheno[which(Pheno$Photo > 0), "Timepoint" == 6]
Pheno[which(Pheno$Photo > 0), "Timepoint" == 11]
Pheno[which(Pheno$Photo > 0), "Timepoint" == 12]

# Check that Fv'Fm' is zero at timepoints 5, 6, 11, 12. 
Pheno[which(Pheno$`Fv'Fm'` != 0), "Timepoint" == 5]
Pheno[which(Pheno$`Fv'Fm'` != 0), "Timepoint" == 6]
Pheno[which(Pheno$`Fv'Fm'` != 0), "Timepoint" == 11]
Pheno[which(Pheno$`Fv'Fm'` != 0), "Timepoint" == 12]

# Check that Fv'Fm' is positive at all timepoints. 
sum(Pheno$`Fv'Fm'` < 0)

# Write file to csv. 
write.csv(Pheno, file = "PhenotypeBrassica.csv")




