# GeneDEmaSigPro2.R
# R version 3.3.1 (2016-06-21)
# March 11, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating transcriptomic network for Brassica data using 
# simone package. Expression data taken from Brassica under well-watered
# and droughted conditions. 

#-----------------------------------------------------------------------
library(maSigPro)
library(stringr)
#-----------------------------------------------------------------------

#### Pre-processing: 

# Read in the data. 
dBrassicaFPKM <- read.table(file = "/Users/mblai/Documents/RNAseqData/R500_D_FPKM.txt", header = T)
wwBrassicaFPKM <- read.table(file = "/Users/mblai/Documents/RNAseqData/R500_WW_FPKM.txt", header = T)

#dBrassicaFPKM <- read.table(file = "R500_D_FPKM.txt", header = T)
#wwBrassicaFPKM <- read.table(file = "R500_WW_FPKM.txt", header = T)


# Remove non-unique genes. 
# TODO: Look into why there are repeated rows.
dBrassicaFPKM <- dBrassicaFPKM[duplicated(dBrassicaFPKM$gene_short_name) == 0, ]
wwBrassicaFPKM <- wwBrassicaFPKM[duplicated(wwBrassicaFPKM$gene_short_name) == 0, ]

# Merge data.frames.
BrassicaFPKM <- merge(dBrassicaFPKM, wwBrassicaFPKM, by = "gene_short_name")

# Remove unnecessary data.frames. 
rm(dBrassicaFPKM)
rm(wwBrassicaFPKM)

# Assign gene names to rownames. 
rownames(BrassicaFPKM) <- BrassicaFPKM$gene_short_name

# Remove columns without expression values. 
BrassicaFPKM <- BrassicaFPKM[, -c(1, 2, 27)]

# Remove first 6 timepoints. 
BrassicaFPKM <- BrassicaFPKM[, -c(1:12, 25:36)]

# Reduce data set to only genes with some expression. 
BrassicaFPKM <- BrassicaFPKM[rowSums(BrassicaFPKM > 5) >= 10,]

# Log2 transform wwBrassicaFPKM data set. 
log2FPKM <- log2(BrassicaFPKM +.1)

# Remove BrassicaFPKM data set.
rm(BrassicaFPKM)

# Coerce log2FPKM to a matrix. 
ExprMatrix <- as.matrix(log2FPKM)

# Create design dataframe. 
Time <- rep(c(rep(c(1:6), each = 2)), 2)
Replicates <- rep(c(1:12), each = 2)
D <- rep(c(1,0), each = 12)
W <- rep(c(0,1), each = 12)
ExprDesign <- cbind(Time, Replicates, D, W)

# Remove expression matrix vectors. 
rm(Time)
rm(Replicates)
rm(D)
rm(W)

# Make rownames the column names of the expression matrix. 
rownames(ExprDesign) <- colnames(ExprMatrix)

# Format the ExprDesign matrix. 
formatExprDesign <- make.design.matrix(ExprDesign, degree = 5)

# Perform regression fit for time series. 
fit <- p.vector(ExprMatrix, design = formatExprDesign, 
                Q = 0.01, counts = F)

# Select regression model by stepwise regression. 
model <- T.fit(fit, step.method = "two.ways.forward", alfa = .01)

# Get the significantly expressed genes. 
sigs <- get.siggenes(model, vars="each")

# Print the significantly expressed genes. 
summ <- sigs$summary

# Note: this will write multiple png files to the working directory.
# Make sure to set the directory to a directory you would like to keep
# the files. 

# Create png plots for each gene in the WvsD column of summary.
tiff(filename = "WvsD%03d.png", width = 8, height = 11, units = 'in', 
     res = 300)
par(mfrow = c(5,3))
  for (i in 1:(length(unique(summ$WvsD))-1)){
  k <- ExprMatrix[rownames(ExprMatrix) == as.character(summ$WvsD[i]), ]
  PlotGroups(k, edesign = ExprDesign, main = as.character(summ$WvsD[i]))
  }
dev.off()


# Create png plots for each gene in the Time5xW column of summary.
tiff(filename = "Time5xW%03d.png", width = 8, height = 11, units = 'in', 
     res = 300)
par(mfrow = c(5,3))
for (i in 1:(length(unique(summ$Time5xW))-1)){
  k <- ExprMatrix[rownames(ExprMatrix) == as.character(summ$Time5xW[i]), ]
  PlotGroups(k, edesign = ExprDesign, main = as.character(summ$Time5xW[i]))
}
dev.off()


# Create png plots for first 200 genes in the Ind column of summary.
tiff(filename = "IND%03d.png", width = 8, height = 11, 
     units = 'in', res = 300)
par(mfrow = c(5,3))
for (i in 1:(length(unique(summ$independ))-858)){
  k <- ExprMatrix[rownames(ExprMatrix) == 
                    as.character(summ$independ[i]), ]
  PlotGroups(k, edesign = ExprDesign, 
             main = as.character(summ$independ[i]))
}
dev.off()


# Write to csv file. 
#write.csv(DEgenes, file = "BrassicaDEgenes.csv")

  m <- sigs$sig.genes$WvsD$sig.pvalues

top <- rep("NA", dim(m)[1])  
  
for (i in 1:(dim(m)[1])){  
  if (sum(is.na(m[i, seq(4, 14, by = 2)])) < 3){
    top[i] <- rownames(m[i, ])
  }
}

i <- "Bra005337"
k <- ExprMatrix[rownames(ExprMatrix) == i, ]
PlotGroups(k, edesign = ExprDesign, 
             main = i)
  
# Don't particularly care about summ$Time or summ$independ

trmt <- summ[, seq(2, 12, by = 2)]
trmt <- trmt[-c(128:1058), ]

total <- c(as.character(trmt$WvsD), as.character(trmt$TimexW),
           as.character(trmt$Time2xW), as.character(trmt$Time3xW),
           as.character(trmt$Time4xW), as.character(trmt$Time5xW))

total <- unique(total)

tot <- rbind(sigs$sig.genes$WvsD$sig.profiles, 
             sigs$sig.genes$TimexW$sig.profiles,
             sigs$sig.genes$Time2xW$sig.profiles,
             sigs$sig.genes$Time3xW$sig.profiles,
             sigs$sig.genes$Time4xW$sig.profiles,
             sigs$sig.genes$Time5xW$sig.profiles)


# Run differential expression analysis. 
DEanalysis <- maSigPro(ExprMatrix, ExprDesign, 
                       Q = .025, counts = F, 
                       step.method = "two.ways.backward", 
                       pdf = T, degree = 5)

# Extract differentially expressed genes from output. 
DEgenes <- rbind(DEanalysis$sig.genes$D$sig.profiles, 
                 DEanalysis$sig.genes$WvsD$sig.profiles)


inf <- DEanalysis$influ.info

topInf <- inf[, which(colSums(inf) > 26)]

# Create png plots for each gene in the inf data slot of summary.
tiff(filename = "inf27%03d.png", width = 8, height = 11, units = 'in', 
     res = 300)
par(mfrow = c(5,3))
for (i in 1:(length(colnames(topInf)))){
  k <- ExprMatrix[rownames(ExprMatrix) == colnames(topInf)[i], ]
  PlotGroups(k, edesign = ExprDesign, main = colnames(topInf)[i])
}
dev.off()

genes <- DEanalysis$sig.genes$WvsD$sig.pvalues

# Subset genes with more than one significant p-value for treatment.
g <- genes[which(is.na(genes$p.valor_WvsD) == F 
                 & is.na(genes$p.valor_TimexW) == F
                 & is.na(genes$p.valor_Time2xW) == F
                 & is.na(genes$p.valor_Time3xW) == F
                 & is.na(genes$p.valor_Time4xW) == F
                 & is.na(genes$p.valor_Time5xW) == F), ]

# Create png plots for each gene in the inf data slot of summary.
tiff(filename = "allTP%03d.png", width = 8, height = 11, units = 'in', 
     res = 300)
par(mfrow = c(5,3))
for (i in 1:(length(rownames(g)))){
  k <- ExprMatrix[rownames(ExprMatrix) == rownames(g)[i], ]
  PlotGroups(k, edesign = ExprDesign, main = rownames(g)[i])
}
dev.off()


