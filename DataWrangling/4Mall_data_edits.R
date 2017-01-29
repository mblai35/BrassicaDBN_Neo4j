# 4Mall_data_edits.R
# R version 3.2.2 (2015-08-14)
# August 23, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# 4Mall* data consistency fixes.

# This "data clean-up" is meant to add consistency to our data files. 
# eg. '4Mall_Dried_Biomass_dec_2013.csv' uses "W" and "D" for treatment
# conditions whereas other files use "WW" (for "well-watered") and "Dry."
# Change "W's" and "D's" were to be consistent with the other files. The
# "Line" column should also be changed to R500 for consistency. Similar 
# edits will be made to other "4Mall" files. 

library(data.table)
library(stringr)

#-----------------------------------------------------------------------

# Read in 4Mall_Dried_Biomass_dec_2013 file
DriedBiomass <- read.csv(file = "4Mall_Dried_Biomass_dec_2013.csv") 

# Change column names to get rid of periods. 
colnames(DriedBiomass) <- sub("Date.after.Oven", "DateAfterOven", 
                              colnames(DriedBiomass))
colnames(DriedBiomass) <- sub("Dry.Weight.Shoots", "DryWeightShoots", 
                              colnames(DriedBiomass))

# Change to data.table for following edits. 
DriedBiomass <- data.table(DriedBiomass)

# Remove Line column.
DriedBiomass <- DriedBiomass[, Line := NULL]

# Change "W" in Treatment to "WW" and "D" to "Dry".
DriedBiomass <- DriedBiomass[Treatment == "W", Treatment := "WW"]
DriedBiomass <- DriedBiomass[Treatment == "D", Treatment := "Dry"]

# Write to a csv file. 
write.csv(DriedBiomass, file = "4Mall_Dried_Biomass_dec_2013_edit.csv")

########################################################################

# Read in 4Mall_NSC_Starchdec2013 file
Starch <- read.csv(file = "4Mall_NSC_Starchdec2013.csv")

# Turn into data.table.
Starch <- data.table(Starch)

# Remove rows AVG, SD, and SE
Starch <- Starch[ Tissue != 'AVG']
Starch <- Starch[ Tissue != 'SD']
Starch <- Starch[ Tissue != 'SE']

# Remove column X and X.1.
Starch <- Starch[, grep("X", colnames(Starch)) := NULL]

# Rename Time.Point Timepoint
colnames(Starch) <- sub("Time.Point", "Timepoint", colnames(Starch))

# Move CARBOHYDRATES from row one up to column name.
colnames(Starch) <- sub("TOTAL", "TotalCarbohydrates", colnames(Starch))

# Change STARCH to Starch.
colnames(Starch) <- sub("STARCH", "Starch", colnames(Starch))

# Remove the first row. Make sure to keep mg/g as a property of NSC, 
# Starch, and Total Carbohydrates.
Starch <- Starch[-1, ]

# Write to a csv file. 
write.csv(Starch, file = "4Mall_NSC_Starchdec2013_edit.csv")

########################################################################

# Read in 4Mall_SoilMoisture_dec_2013 file
SoilMoisture <- read.csv(file = "4Mall_Soil Moisture_dec_2013.csv")

# Turn into data.table.
SoilMoisture <- data.table(SoilMoisture)

# Rename column "SM...." as "SM(%)".
colnames(SoilMoisture)[5] <- "SM(%)"

# Remove Growth.Chamber column.
SoilMoisture <- SoilMoisture[, "Growth.Chamber" :=NULL]

# Rename "Replicate.Number" as "Replicate".
colnames(SoilMoisture) <- sub("Replicate.Number", "Replicate", 
                              colnames(SoilMoisture))

# Split the timepoint from the Time column into its own column. 
SoilMoisture <- SoilMoisture[, 'Timepoint' := str_split_fixed(SoilMoisture$Time, " ", 4)[, 2]]

# Remove the TP from the Timepoint column 
SoilMoisture <- SoilMoisture[, 'Timepoint' := str_split_fixed(SoilMoisture$Timepoint, "TP", 2)[, 2]]

# Update the Time column to contain only the time.
SoilMoisture <- SoilMoisture[, 'Time' := str_split_fixed(SoilMoisture$Time, " ", 4)[, 4]]

# Change a.m. to AM
SoilMoisture$Time <- sub(" a.m.", ":00 AM", SoilMoisture$Time)

# Change p.m. to PM
SoilMoisture$Time <- sub(" p.m.", ":00 PM", SoilMoisture$Time)

# Write to a csv file. 
write.csv(SoilMoisture, file = "4Mall_Soil_Moisture_dec_2013_edit.csv")

########################################################################

# Read in 4Mall_PhotoFv'Fm'gs_dec_2013 file
Photo <- read.csv(file = "4Mall_PhotoFv'Fm'gs_dec_2013.csv")

# Turn into data.table.
Photo <- data.table(Photo)

# Remove rows AVG, SD, and SE
Photo <- Photo[ X != 'AVG']
Photo <- Photo[ X != 'SD']
Photo <- Photo[ X != 'SE']
Photo <- Photo[ X != 'tttest']

# Remove column X.
Photo <- Photo[, "X" := NULL]

# Remove empty rows.
Photo <- Photo[1:240]

# Rename Time.Point as Timepoint
colnames(Photo) <- sub("Time.Point", "Timepoint", colnames(Photo))

# Rename column "Fv.Fm." as "Fv'Fm'".
colnames(Photo)[10] <- "Fv'Fm'"

# Rename column "X..Rep" as "Replicate".
colnames(Photo)[7] <- "Replicate"

# Change column name "Time.of.day" to "Time". 
colnames(Photo) <- sub("Time.of.day", "Time", 
                       colnames(Photo))

# Change a.m. to AM
Photo$Time <- sub("a.m.", "AM", Photo$Time)

# Change p.m. to PM
Photo$Time <- sub("p.m.", "PM", Photo$Time)

# Remove TimeinhrsfromZT0 column.
Photo <- Photo[, "Time.in.hrs.from.ZT0" := NULL]

# Remove Growth.Chamber column.
Photo <- Photo[, "Growth.Chamber" := NULL]

# Remove Line column.
Photo <- Photo[, "Line" := NULL]

# Write to a csv file. 
write.csv(Photo, file = "4Mall_PhotoFv'Fm'gs_dec_2013_edit.csv")
