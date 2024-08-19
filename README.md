# Analysis of Microbial Communities Using DADA2

This document guides the analysis of microbial communities using the DADA2 pipeline in R, including quality filtering, error learning, dereplication, and merging of paired reads.

## Load Required Packages

```r
library(dada2)
```

## Define File Paths

```r
path <- "~/citl_fastq_data_2/"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
```

## Sample names

```r 
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
cat("The sequence files are:\n")
print(fnFs)
cat("\nThe sample names are:\n")
print(sample.names)
```

## Plot the quality of R1 and R2 files 

```r 
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

## Quality Filtering

```r
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                     compress=TRUE, multithread=TRUE)
```

## Learn the Error Rates

```r
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```

## Sample Inference

```r
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
```

## Merging Paired Reads

```r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
```

## Constructing the Sequence Table

```r
seqtab <- makeSequenceTable(mergers)
```

## Remove Chimeras

```r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
```

## Track Reads Through the Pipeline

```r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample.names
```

## Save the Sequence Table

```r
saveRDS(seqtab.nochim, file="path/to/seqtab_nochim.rds")
```

## Session Info

```r
sessionInfo()
```
