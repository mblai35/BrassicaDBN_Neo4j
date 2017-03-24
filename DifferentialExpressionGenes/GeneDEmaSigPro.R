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

# Format the ExprDesign matrix. 
formatExprDesign <- make.design.matrix(ExprDesign, degree = 11)

# Perform regression fit for time series. 
fit <- p.vector(ExprMatrix, design = formatExprDesign, 
                Q = 0.05, counts = F)

# Select regression model by stepwise regression. 
model <- T.fit(fit, step.method = "two.ways.forward", alfa = .05)

# Investigation of influential data was performed. 

# Get the significantly expressed genes. 
sigs <- get.siggenes(model, vars="groups", rsq = 0)

# Save p-values for genes with treatment differences to a vector g. 
g <- sigs$sig.genes$WvsD$sig.pvalues

# For each timepoint, grab genes with a p-value less than .0001.
g0 <- rownames(g[which(g$p.valor_WvsD < .0001), ])
g1 <- rownames(g[which(g$p.valor_TimexW < .0001), ])
g2 <- rownames(g[which(g$p.valor_Time2xW < .0001), ])
g3 <- rownames(g[which(g$p.valor_Time3xW < .0001), ])
g4 <- rownames(g[which(g$p.valor_Time4xW < .0001), ])
g5 <- rownames(g[which(g$p.valor_Time5xW < .0001), ])
g6 <- rownames(g[which(g$p.valor_Time6xW < .0001), ])
g7 <- rownames(g[which(g$p.valor_Time7xW < .0001), ])
g8 <- rownames(g[which(g$p.valor_Time8xW < .0001), ])
g9 <- rownames(g[which(g$p.valor_Time9xW < .0001), ])
g10 <- rownames(g[which(g$p.valor_Time10xW < .0001), ])
g11 <- rownames(g[which(g$p.valor_Time11xW < .0001), ])

highSigG <- rownames(g[which(g$`p-value` < .000000001), ])
DEgenes <- ExprMatrix[highSigG, ]

# Extract all unique gene names from each timepoint.
DEgenes <- unique(c(g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11))

# Extract differentially expressed genes from ExprMatrix.
DEgenes <- ExprMatrix[DEgenes, ]
########################################################################

# Plot final DEgenes.
tiff(filename = "highSigDEgenes%03d.png", width = 8, height = 11, 
     units = 'in', res = 300)
par(mfrow = c(5,3))
for (i in 1:(dim(DEgenes)[1])){
  k <- DEgenes[i, ]
  PlotGroups(k, edesign = ExprDesign, main = rownames(DEgenes)[i])
}
dev.off()  

p <- function(i){
g[which(rownames(g) == i),]
}

g[which(rownames(g) == "Bra002221"),]

p("Bra002221")
