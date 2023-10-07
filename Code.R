# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Differential expression analysis with limma
library(GEOquery)
#Loading required package: Biobase
#Loading required package: BiocGenerics
#Attaching package: ‘BiocGenerics’
library(limma)
#Attaching package: ‘limma’
#The following object is masked from ‘package:BiocGenerics’:
#plotMA
library(umap)

# load series and platform data from GEO
gset <- getGEO("GSE210064", GSEMatrix =TRUE, AnnotGPL=FALSE)
#retriving the data  by using the GEOquery package
if (length(gset) > 1) idx <- grep("GPL20844", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#finding the length by if condition and retrieving the gene names by comparing the GPL20844 
#checking the condition and storing the data to the variable gset

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "0000000000011111111111"
#defining the control and test sample sets
sml <- strsplit(gsms, split="")[[1]]
#converting into character array

# log2 transformation
ex <- exprs(gset)
#assign
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
#quality
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
#converting to log values based on the below condition
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }



# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("hai","bye"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
#comparing the expression levels

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SPOT_ID","GB_ACC","SEQUENCE","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")
#66 line making the schema of the table
#visualisation of the table

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
#histogram which diplay the p value is the probability of observing a test statistic 
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)
#adjusting the p value that is greater than the alpha value (0.05 ) 

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
#theoretical quantiles:ata set follows a particular theoretical distribution.
#It is used to compare the quantiles of the observed data to the quantiles of a theoretical distribution
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
#fold value is the change in the expression

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)
#Mean-Difference plot systematic differences between two measurement methods or instruments.

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE210064", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
#block plot visualize the spread, skewness, and central tendency of the data.of the control and test groups

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE210064", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")
#draw the graph by comparing the above plots

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 9, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=9", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)
#dimensionality reduction technique

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE210064")
#relationship between the expected return and the risk