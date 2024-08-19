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

# Taxonomic Assignment

```r
taxa <- assignTaxonomy(seqtab.nochim, 
                       "greengenes2_trainset.fa.gz", 
                       multithread = TRUE)
```

## Removing sequence row names for visualization 

```r
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```


# Creating a Phyloseq Object

## Load necessary packages

```r
library(tidyverse)
library(vegan)
library(ggplot2)
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
library(phyloseq); packageVersion("phyloseq")
```

## Load metadata

```r
metadata <- read_csv("00.Data/Metadata.csv") %>%
  mutate_at(vars(depth), as.character) %>%
  column_to_rownames("seqR1")
```

## Create phylogeny

```r
ASVs.nochim = DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) = paste0("ASV", 1:ncol(seqtab.nochim))

alignment = AlignSeqs(ASVs.nochim, anchor = NA, processors = 30)

phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

fit = pml(treeNJ, data = phang.align)
fitGTR_green <- update(fit, k = 4, inv = 0.2)
```

## Create the object

```r
load("00.Data/tree.green.RData")

ASVs.nochim = DNAStringSet(colnames(seqtab.nochim))
names(ASVs.nochim) = paste0("ASV", 1:ncol(seqtab.nochim))

tmp.seqtab = seqtab.nochim
colnames(tmp.seqtab) = names(ASVs.nochim)
tmp.taxa = taxa.print
rownames(tmp.taxa) = names(ASVs.nochim)

ps.green.nochim = phyloseq(
  otu_table(tmp.seqtab, taxa_are_rows = FALSE),
  sample_data(metadata),
  tax_table(tmp.taxa),
  refseq(ASVs.nochim),
  phy_tree(fitGTR_green$tree)
)
```

## Create phylogenetic tree nodes

```r
ps = ps.green.nochim
set.seed(1)
is.rooted(phy_tree(ps))
phy_tree(ps) <- root(phy_tree(ps), sample(taxa_names(ps), 1), 
                     resolve.root = TRUE)
is.rooted(phy_tree(ps))
```

# Exploratory NMDS

```r
library(microbiome)
library(patchwork)
```

## Extract metadata

```r
metadata_green <- meta(ps) %>%
  rownames_to_column("seqR1") %>%
  mutate(season = str_replace(season, "flood", "Flood")) %>%
  mutate(season = str_replace(season, "dry", "Dry")) 

metadata_green %>%
  select(month, season)
```

## Extract OTU matrix from phyloseq object

```r
OTU1 = as(otu_table(ps), "matrix")
# Transpose if necessary 
if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Convert to data frame
OTUdf = as.data.frame(OTU1)
```

## MDS 

```r
OTUdf_prop <- OTUdf / rowSums(OTUdf)
mds <- metaMDS(OTUdf_prop, distance = "bray", trymax = 9000)
stressplot(mds) # checks that the fit is good

en = envfit(mds, metadata_green, permutations = 9999, na.rm = TRUE)
```

## Convert to long format for ggplot

```r
my.nmds <- data.frame("seqR1" = rownames(mds$points), "NMDS1" = mds$points[ ,1],
                      "NMDS2" = mds$points[ ,2])
my.nmds <- merge(my.nmds, metadata_green, by = "seqR1")
```

## Exploratory

```r
plot_nmds <- ggplot(data = my.nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = seqR1, shape = season), size = 3, alpha = 0.5) + 
  # scale_colour_manual(values = c("orange", "steelblue", "red")) +
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) + 
  labs(colour = "Season")
plotly::ggplotly(plot_nmds)

plot2 <- ggplot(data = my.nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = season), size = 5, alpha = 0.5) + 
  scale_colour_manual(values = c("darkgreen", "purple")) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
        axis.title.y = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30")) + 
  labs(colour = "Season")

plot1 <- ggplot(data = my.nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = zone), size = 5, alpha = 0.5) + 
  scale_colour_manual(values = c("darkorange", "steelblue", "red")) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30")) + 
  labs(colour = "Zone") +
  theme_bw()

plot3 <- ggplot(data = my.nmds, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = depth), size = 5, alpha = 0.5) + 
  scale_colour_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
  theme(axis.title = element_text(size = 12, face = "bold", colour = "grey30"),
        axis.title.y = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "grey30"), 
        legend.key = element_blank(), 
        legend.title = element_text(size = 12, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 12, colour = "grey30")) + 
  labs(colour = "Depth")

plot1 + plot2 + plot3 +
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A')

```

```r
my.stats <- anosim(x = as.data.frame(otu_table(ps)), 
                   grouping = metadata_green$depth, 
                   permutations = 9999,
                   distance = "bray")

my.stats
```
