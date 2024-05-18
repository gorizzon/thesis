library(readr)
library(data.table)

load("raw_data/BRCAData.rda")

###### Separate the omics platforms ######

mRNA  <- as.data.frame(t(BRCAData[[1]]))
Meth  <- as.data.frame(t(BRCAData[[2]]))
RPPA  <- as.data.frame(t(BRCAData[[4]]))

setDT(mRNA, keep.rownames = "Sample")
setDT(Meth, keep.rownames = "Sample")
setDT(RPPA, keep.rownames = "Sample")

###### Shorten ID names - collection time not needed for analysis ######

Meth$Sample <- gsub("\\.01A", "", Meth$Sample)
Meth$Sample <- gsub("\\.01B", "", Meth$Sample)
Meth$Sample <- gsub("\\.", "", Meth$Sample)

mRNA$Sample <- gsub("\\.01A", "", mRNA$Sample)
mRNA$Sample <- gsub("\\.01B", "", mRNA$Sample)
mRNA$Sample <- gsub("\\.", "", mRNA$Sample)

RPPA$Sample <- gsub("\\.01A", "", RPPA$Sample)
RPPA$Sample <- gsub("\\.01B", "", RPPA$Sample)
RPPA$Sample <- gsub("\\.", "", RPPA$Sample)

###### Remove case with missing data (???) ######

Meth <- Meth[Meth$Sample != "TCGAB6A0RM", ]
mRNA <- mRNA[mRNA$Sample != "TCGAB6A0RM", ]
RPPA <- RPPA[RPPA$Sample != "TCGAB6A0RM", ]


###### Norm to 0-1 ######

Meth_norm <- apply(Meth[, 2:575], 2, function(x) (x - (min(x)))/max((x) - min(x)))
Meth <- cbind(Meth[, 1], Meth_norm)

mRNA_norm <- apply(mRNA[, 2:ncol(mRNA)], 2, function(x) (x - (min(x)))/max((x) - min(x)))
mRNA <- cbind(Meth[, 1], mRNA_norm)

RPPA_norm <- apply(RPPA[, 2:ncol(RPPA)], 2, function(x) (x - (min(x)))/max((x) - min(x)))
RPPA <- cbind(Meth[, 1], RPPA_norm)

###### Write csv files ######

write_csv(mRNA, file = "data/mRNA.csv")
write_csv(Meth, file = "data/Meth.csv")
write_csv(miRNA, file = "data/miRNA.csv")
write_csv(RPPA, file = "data/RPPA.csv")

