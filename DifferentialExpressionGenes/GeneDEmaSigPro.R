# GeneDEmaSigPro.R
# R version 3.3.1 (2016-06-21)
# January 26, 2017. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating transcriptomic network for Brassica data using 
# simone package. Expression data taken from Brassica under well-watered
# and droughted conditions. 

#-----------------------------------------------------------------------
library(maSigPro)
library(MASS)
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

# Reduce data set to only genes with some expression. 
BrassicaFPKM <- BrassicaFPKM[rowSums(BrassicaFPKM > 10) >= 48,]

# Log2 transform wwBrassicaFPKM data set. 
log2FPKM <- log2(BrassicaFPKM)

# Remove BrassicaFPKM data set.
rm(BrassicaFPKM)

# Coerce log2FPKM to a matrix. 
ExprMatrix <- as.matrix(log2FPKM)
ExprMatrix <- as.matrix(BrassicaFPKM)

# Create design dataframe. 
Time <- rep(c(rep(c(1:12), each = 2)), 2)
Replicates <- rep(c(1:24), each = 2)
D <- rep(c(1,0), each = 24)
W <- rep(c(0,1), each = 24)
ExprDesign <- cbind(Time, Replicates, D, W)

# Make rownames the column names of the expression matrix. 
rownames(ExprDesign) <- colnames(ExprMatrix)

# Note: doing this will write a results pdf to the working directory.
# Make sure to set the directory to a directory you would like to keep
# the file. 
# Run differential expression analysis. 
DEanalysis <- maSigPro(ExprMatrix, ExprDesign, degree = 4)

# Extract differentially expressed genes from output. 
DEgenes <- rbind(DEanalysis$sig.genes$D$sig.profiles, 
                DEanalysis$sig.genes$WvsD$sig.profiles)

# Write to csv file. 
write.csv(DEgenes, file = "BrassicaDEgenes.csv")







# Alternatively...
########################################################################
# Format the ExprDesign matrix. 
formatExprDesign <- make.design.matrix(ExprDesign, degree = 4)

# Perform regression fit for time series. 
fit <- p.vector(ExprMatrix, design = formatExprDesign, 
                Q = 0.05, counts = F)

# Select regression model by stepwise regression. 
model <- T.fit(fit)

# Get the significantly expressed genes. 
get <- get.siggenes(model, vars="all")

# Print the significantly expressed genes. 
get$summary
########################################################################


