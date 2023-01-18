#updated 2020_11_11
#adapted from chip_seq tutorial @physalia course
#Aarathy
#To find differential binding of Transcription factor IRF9 in 
#treatment vs untreated condition (treatment with IFN-I (IFN beta))
# a consensus bed file and reads from respective conditions are used as inputs
#Packages
library("csaw")
library("edgeR")
library("statmod")

#output
outFolder <- "./ChIP_seq/"
#files 
meta <- c()

meta[2] <- "WT_IRF9_IFNB"
meta[3] <- "WT_IRF9_UT"
#path
path <- "/home/decker/data/Nextflow_projects/Katrin/Katrin_ChIP_BMDM_without_dups/results/bwa/mergedLibrary/"

#read files
bamFiles <- paste(path,"/", meta[2:3], "_R", rep(1:2, each=2), ".mLb.clN.sorted.bam", sep="")

# peak file specifying regions where read counts are estimated
regions <- read.table("macs/narrowPeak/consensus/IRF9/IRF9.consensus_peaks_katrin.bed", sep="\t", header=FALSE)

# design matrix  
treatment <- rep(meta[2:3], 2)
treatment[treatment == meta[2]] <- 1
treatment[treatment == meta[3]] <- 0
treatment <- factor(treatment)
?model.matrix
design <- model.matrix(~treatment)
colnames(design) <- c("intercept", "treatment")

#creating genomic ranges
reg_id <- paste(regions[,1], ":", regions[,2], "-", regions[,3], sep="")
d <- GRanges(regions[,1], IRanges(regions[,2], regions[,3]), id=reg_id)

#assessing counts in given region
dCounts <- regionCounts(bamFiles, regions=d)

#normalize
#TMM
data <- normFactors(dCounts, se.out=TRUE)

#QC plots
#just counts per million
cpms_norm_N <- cpm(asDGEList(dCounts))
#Cpm+TMM
cpms_norm_Y <- cpm(asDGEList(data), normalized.lib.size=TRUE)
#replacing zero with leastvalue
cpms_norm_N[cpms_norm_N == 0] <- min(cpms_norm_N[cpms_norm_N > 0])
cpms_norm_Y[cpms_norm_Y == 0] <- min(cpms_norm_Y[cpms_norm_Y > 0])

ylab <- "Normalized counts"
nms <- rep("", 4)

#edgeR

y <- asDGEList(data)
?estimateDisp
y <- estimateDisp(y, design)
y
#fitting a quasi-likelihood negative binomial generalized log-linear model to count data
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)
#adjusteed pvalue (Benjamini & Hochberg)

Padj <- p.adjust(results$table$PValue, method = "BH")

#output data
out <- data.frame(chrom=regions[,1],
                  start=regions[,2],
                  end=regions[,3],
                  logCPM=results$table$logCPM,
                  logFC=results$table$logFC,
                  PValue=results$table$PValue,
                  Padj)
out["group"]<-"ns"
out[which(out['Padj'] <= 0.05 & out['logFC'] >= 1.0 ),"group"] <- "up"
out[which(out['Padj'] <= 0.05 & out['logFC'] <= 1.0 ),"group"] <- "down"

#write table
write.table(out, paste(outFolder,"results.chip-seq_diffbinding.txt",sep = ""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

