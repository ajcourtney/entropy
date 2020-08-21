#Author Abigail J. Courtney, University of Georgia, 2020
library(BioQC)
library(scater)
library(ggplot2)

#set working directory to where your files are located
setwd("~/Dropbox/LewisLab/Research/EntropyAnalysis/EntropyAnalysis/")

#read in length of each gene to use in the TPM normalization below
#this can be collected directly from a single featureCounts file using the "length" column (column 4) or put into its own file like this example
lengthTable = read.table(file="./NewCountsNoName/FeatureCountsGeneLengthOutput.txt", header=TRUE, stringsAsFactors = FALSE)
#grab the gene names and use as row names
rnames_length=lengthTable[,1]
#convert to data matrix format (using all columns besides column with gene names)
matLength = data.matrix(lengthTable[2:ncol(lengthTable)])
#assign row names (gene names) to matrix
rownames(matLength)<-rnames_length

#read in table of counts to be used
countTable = read.table(file="./NewCountsNoName/combinedCountsGeneIDfinal.2.txt", header=TRUE, stringsAsFactors=FALSE)
#grab the gene names and use as row names
rnames_count=countTable[,1]
#convert to data matrix format (using all columns besides column with gene names)
matCount=data.matrix(countTable[,2:ncol(countTable)])
#assign row names (gene names) to matrix
rownames(matCount)<- rnames_count

#calculate TPM using the the counts matrix and the length matrix previously read in with scater
TPMtable <- calculateTPM(matCount, effective_length = matLength[,1])

#calculate entropy using entropySpecificity function from BioQC
tpmE <- data.frame(entropySpecificity(TPMtable))
#output gene IDs and entropy values to a file
write.table(tpmE, file="EntropyValuesGeneIDs.txt", quote = FALSE, sep="\t")

#Create Density plot of all entropy values
allDensity <- ggplot(tpmE, aes(x=tpmE$entropySpecificity.TPMtable.)) + 
  geom_density(fill="blue", alpha=0.6) +
  xlim((floor(min(tpmE$entropySpecificity.TPMtable.))),ceiling(max(tpmE$entropySpecificity.TPMtable.)))+
  labs(x = "Entropy", y = "Density") +
  theme(axis.title.x = element_text("Arial", size = 14))+
  theme(axis.title.y = element_text("Arial", size = 14))+
  theme(axis.text = element_text("Arial", size = 12))+
  geom_rug()

ggsave("./R_Plots/allDensity.png", plot=allDensity, dpi=300, width=5, height=6.5, device="png")

#create and save histogram of all entropy values
#determine bin width
all_bw <- ((max(na.omit(tpmE$entropySpecificity.TPMtable.))) - (min(na.omit(tpmE$entropySpecificity.TPMtable.)))) / ceiling(sqrt(length(tpmE$entropySpecificity.TPMtable.)))
#plot histogram
allHistogram <- ggplot(tpmE, aes(x=tpmE$entropySpecificity.TPMtable.)) +
  geom_histogram(alpha=0.1, color="blue", fill="white", binwidth = all_bw) +
  labs(x="Entropy", y="Count")+
  xlim(0,7)+
  theme(axis.title.x = element_text("Arial", size = 14))+
  theme(axis.title.y = element_text("Arial", size = 14))+
  theme(axis.text = element_text("Arial", size = 12))

ggsave("./R_Plots/allHistogram.png", plot=allHistogram, dpi=300, width=5, height=6.5, device="png")