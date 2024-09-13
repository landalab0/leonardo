# Analysis of Microbial Communities Using DADA2

This document guides the analysis of microbial communities using the DADA2 pipeline in R, including quality filtering, error learning, dereplication, and merging of paired reads.

## Load Required Packages

Load the required packages

```r
library(dada2)
```

## Define File Paths

Define the following path variable so that it points to the extracted directory and we read in the names of the fastq files, and perform some string manipulation to get lists of the forward and reverse fastq files in matched order

```r
path <- "~/citl_fastq_data_2/"
list.files(path)

fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
```

## Sample names

Extract the sample names of the fasq files

```r 
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)
cat("The sequence files are:\n")
print(fnFs)
cat("\nThe sample names are:\n")
print(sample.names)
```

## Plot the quality of R1 and R2 files 

View visualizing the quality profiles of the forward and reverse reads

```r 
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
```

## Quality Filtering

We proceed with filtering and trimming the reads according to what is observed in the quality plots. Specifically, we use the arguments truncLen=(indicate the average trimming length for each read), maxN=(the maximum number of ambiguous bases), maxEE= (the maximum number of errors) , and rm.phix=(whether we want to remove sequences belonging to the internal Illumina control).

```r
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), 
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, 
                     compress=TRUE, multithread=TRUE)
```

## Learn the Error Rates

Generating a probabilistic error model, it can filter erroneous reads and thus use the remaining ones directly for the taxonomic classification stage. This part of the method allows us to achieve higher resolution compared to OTU-based analyses

```r
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
```

## Sample Inference

The inference of Amplicon Sequence Variants (ASVs).

```r
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)
```

## Merging Paired Reads

We merge the R1 and R2 reads

```r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
```

## Constructing the Sequence Table

Generate a sequence table that we will use later

```r
seqtab <- makeSequenceTable(mergers)
```

## Remove Chimeras

Remove artifact sequences resulting from PCR amplification

```r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=TRUE, verbose=TRUE)
```

## Track Reads Through the Pipeline

Final check of our progress

```r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), 
               sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", 
                     "merged", "nonchim")
rownames(track) <- sample.names
```

# Taxonomic Assignment

The sequence table without chimeras is the table we use for taxonomic classification. We use Greengenes databes to sequence classification

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
  asv_table(tmp.seqtab, taxa_are_rows = FALSE),
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

## Extract ASV matrix from phyloseq object

```r
ASV1 = as(asv_table(ps), "matrix")
# Transpose if necessary 
if(taxa_are_rows(ps)){ASV1 <- t(ASV1)}
# Convert to data frame
ASVdf = as.data.frame(ASV1)
```

## MDS 

```r
ASVdf_prop <- ASVdf / rowSums(ASVdf)
mds <- metaMDS(ASVdf_prop, distance = "bray", trymax = 9000)
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
my.stats <- anosim(x = as.data.frame(asv_table(ps)), 
                   grouping = metadata_green$depth, 
                   permutations = 9999,
                   distance = "bray")

my.stats
```

# Prevalence 

We can look at is the prevalence of the taxonomic features

## Load packages

```r
library(tidyverse)
library(phyloseq)
library(ggplot2)

```

## Filtering Taxa

We generate a vector with all the taxa that we want to filter

```r
ps.pruned.1 <- subset_taxa(ps.pruned, !is.na(Phylum) & !Phylum %in% 
                             c("", "uncharacterized"))
```

## Compute prevalence of each feature, store as data.frame

```r
prevdf = apply(X = asv_table(ps.pruned.1),
               MARGIN = ifelse(taxa_are_rows(ps.pruned.1), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
```

## Add taxonomy and total read counts to this data.frame

```r
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps.pruned.1),
                    tax_table(ps.pruned.1))
```

## Viewing
```r
head(prevdf)
```

## Computing average and total prevalence

```r
temp_df <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
```

## Vector to be removed

```r
filterPhyla <- temp_df[temp_df$`2` < 2,]$Phylum
```

## Filter entries with unidentified Phylum.
```r
ps1 = subset_taxa(ps.pruned, !Phylum %in% filterPhyla)
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps1),color=Phylum)) +
## Include a guess for parameter
geom_hline(yintercept = 0.02, alpha = 0.5, linetype = 2) + 
geom_point(size = 1, alpha = 0.7) +
scale_x_log10() +  xlab("Total Abundance") + 
ylab("Prevalence [Frac. Samples]") +
facet_wrap(~Phylum) + theme(legend.position="none")
```

##  Define prevalence threshold as 2% of total samples

```r
prevalenceThreshold = 0.02 * nsamples(ps1)
prevalenceThreshold
```

## Execute prevalence filter, using `prune_taxa()` function

```r
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps1)

ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps2),color=Phylum)) +
  ## Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

ASV1 = as(asv_table(ps2), "matrix")
```

## transpose if necessary

```r
if(taxa_are_rows(ps2)){ASV1 <- t(ASV1)}
```

## Coerce to data.frame

```r
ASVdf = as.data.frame(ASV1)
ASVdf_prop <- ASVdf/rowSums(ASVdf)


which(is.na(ASVdf_prop))

ASVdf_prop[!complete.cases(ASVdf_prop), ]
ps3 = subset_samples(ps2, sample_names(ps2) != 'zr2502_16_R1')
```

# Alpha Diversity

## Load packages

```r
library("ggpubr")
library("phyloseq")
library("ggplot2")
library(patchwork)
```

## Manipulating phyloseq object
This code takes a phyloseq object (which has already been filtered and is stored
in ps.pruned) and further filters it by removing taxa where the phylum
information is missing, empty, or considered "uncharacterized"

```r
a_my_comparisons <- list( c("Fringe", "Basin"), c("Fringe", "Impaired"), 
                          c("Basin", "Impaired"))
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                   symbols = c("****", "***", "**", "*", "ns"))
```

## Plot reachnes per zone,season and depht at 5 cm

```r
ps3_5<-subset_samples(ps3, depth == "5")
plot1_5<-plot_richness(ps3_5, x = "zone", 
                     measures = c("Simpson", "Shannon")) +
  geom_boxplot(aes(fill = season), alpha=.7) + 
  scale_color_manual(values = c("steelblue", "darkorange", "red")) +
  scale_fill_manual(values = c("#606060FF","#D6ED17FF")) + 
  theme_bw()

newSTorder = c("Fringe", "Basin", "Impaired")

plot1_5$data$zone <- as.character(plot1_5$data$zone)
plot1_5$data$zone <- factor(plot1_5$data$zone, levels=newSTorder)
print(plot1_5)
```

## Plot reachnes per zone,season and depht at 20 cm

```r
ps3_20<-subset_samples(ps3, depth == "20")

plot1_20<-plot_richness(ps3_20, x = "zone", 
                        measures = c("Simpson", "Shannon")) +
  geom_boxplot(aes(fill = season), alpha=.7) + 
  scale_color_manual(values = c("steelblue", "darkorange", "red")) +
  scale_fill_manual(values = c("#606060FF","#D6ED17FF")) + 
  theme_bw() 

newSTorder = c("Fringe", "Basin", "Impaired")

plot1_20$data$zone <- as.character(plot1_20$data$zone)
plot1_20$data$zone <- factor(plot1_20$data$zone, levels=newSTorder)
print(plot1_20)
```

## Plot reachnes per zone,season and depht at 20 cm

```r
ps3_40<-subset_samples(ps3, depth == "40")

plot1_40<-plot_richness(ps3_40, x = "zone", 
                        measures = c("Simpson", "Shannon")) +
  geom_boxplot(aes(fill = season), alpha=.7) + 
  s0cale_color_manual(values = c("steelblue", "darkorange", "red")) +
  scale_fill_manual(values = c("#606060FF","#D6ED17FF")) + 
  theme_bw() 

newSTorder = c("Fringe", "Basin", "Impaired")

plot1_40$data$zone <- as.character(plot1_40$data$zone)
plot1_40$data$zone <- factor(plot1_40$data$zone, levels=newSTorder)
print(plot1_40)
```

## Anova per zone,season and depht at 5 cm

```r
ps3_5<-subset_samples(ps3, depth == "5")
ps3_5<-subset_samples(ps3_5, season == "dry")
richness <- estimate_richness(ps3_5)
head(richness)

write.table(richness, file = "Alpha_diversity_ps3_5.tsv", sep = "\t", row.names = T,
            quote = F)

anova.sh = aov(richness$Simpson ~ sample_data(ps3_5)$zone)
summary(anova.sh)

TukeyHSD(anova.sh)
pairwise.wilcox.test(richness$Simpson, sample_data(ps3_5)$zone, p.adj = "bonf")

ps3_5<-subset_samples(ps3, depth == "5")
ps3_5<-subset_samples(ps3_5, season == "flood")
richness <- estimate_richness(ps3_5)
head(richness)

write.table(richness, file = "Alpha_diversity_ps3_5.tsv", sep = "\t", row.names = T,
            quote = F)

anova.sh = aov(richness$Simpson ~ sample_data(ps3_5)$zone)
summary(anova.sh)

TukeyHSD(anova.sh)
pairwise.wilcox.test(richness$Simpson, sample_data(ps3_5)$zone, p.adj = "bonf")
```

## Anova per zone,season and depht at 20 cm

```r
ps3_20<-subset_samples(ps3, depth == "20")
ps3_20<-subset_samples(ps3_20, season == "dry")
richness <- estimate_richness(ps3_20)
head(richness)

write.table(richness, file = "Alpha_diversity_ps3_20.tsv", sep = "\t", row.names = T,
            quote = F)

anova.sh = aov(richness$Shannon ~ sample_data(ps3_20)$zone)
summary(anova.sh)

TukeyHSD(anova.sh)
pairwise.wilcox.test(richness$Shannon, sample_data(ps3_20)$zone, p.adj = "bonf")

ps3_20<-subset_samples(ps3, depth == "20")
ps3_20<-subset_samples(ps3_20, season == "flood")
richness <- estimate_richness(ps3_20)
head(richness)

write.table(richness, file = "Alpha_diversity_ps3_20.tsv", sep = "\t", row.names = T,
            quote = F)

anova.sh = aov(richness$Simpson ~ sample_data(ps3_20)$zone)
summary(anova.sh)

TukeyHSD(anova.sh)
pairwise.wilcox.test(richness$Simpson, sample_data(ps3_20)$zone, p.adj = "bonf")
```

## Anova per zone,season and depht at 40 cm

```r
ps3_40<-subset_samples(ps3, depth == "40")
ps3_40<-subset_samples(ps3_40, season == "dry")
richness <- estimate_richness(ps3_40)
head(richness)

write.table(richness, file = "Alpha_diversity_ps3_40.tsv", sep = "\t", row.names = T,
            quote = F)

anova.sh = aov(richness$Simpson ~ sample_data(ps3_40)$zone)
summary(anova.sh)

TukeyHSD(anova.sh)
pairwise.wilcox.test(richness$Simpson, sample_data(ps3_40)$zone, p.adj = "bonf")

ps3_40<-subset_samples(ps3, depth == "40")
ps3_40<-subset_samples(ps3_40, season == "flood")
richness <- estimate_richness(ps3_40)
head(richness)

write.table(richness, file = "Alpha_diversity_ps3_40.tsv", sep = "\t", row.names = T,
            quote = F)

anova.sh = aov(richness$Simpson ~ sample_data(ps3_40)$zone)
summary(anova.sh)

TukeyHSD(anova.sh)
pairwise.wilcox.test(richness$Simpson, sample_data(ps3_40)$zone, p.adj = "bonf")
```

# Beta diversity

## Load Library

```r
library(phyloseq)
library(microViz)
library(patchwork)
library(tidyverse)
```

## Principal component analysis (PCA) with transformation clr per season and zone

```r
plot3<-ps3 %>%
  tax_fix() %>%
  tax_transform(trans = "clr")  %>%
  ord_calc(method = "PCA") %>%
  ord_plot(color = "zone", shape = "season", size = 3)+ 
  scale_colour_manual(values = c("darkorange", "steelblue", "red"))  +
  theme_bw()  + 
  labs(colour = "Zone", shape = "Season") + 
  scale_shape_manual(values = c(15,17))
```

## PERMANOVA with Bray-Curtis distance

```r
ps3 %>%
  tax_transform("identity") %>% # don't transform!
  dist_calc("bray") %>%
  dist_permanova(
    seed = 99,
    variables = c("zone", "season", "depth"),
    n_processes = 1,
    n_perms = 9999 # only 99 perms used in examples for speed (use 9999+!)
  )
```

## Principal coordinate analysis (PCoA) with distance Bray-Curtis per Zone and Season with MDS1 vs MDS2

```r
plot1<-ps3 %>%
  tax_transform("identity") %>% # don't transform!
  dist_calc("bray") %>%
  ord_calc("PCoA") %>%
  ord_plot(axes = c(1, 2), color = "zone", shape ="season", size = 3) + 
  scale_colour_manual(values = c("darkorange", "steelblue", "red"))  +
  theme_bw() + 
  labs(colour = "Zone", shape = "Season") + 
  scale_shape_manual(values = c(15,17))
```

## Principal coordinate analysis (PCoA) with distance Bray-Curtis per Zone and Season with MDS1 vs MDS3

```r
plot4<-ps3 %>%
  tax_transform("identity") %>% # don't transform!
  dist_calc("bray") %>%
  ord_calc("PCoA") %>%
  ord_plot(axes = c(1, 3), color = "zone", shape ="season",
           size = 3) + 
  scale_colour_manual(values = c("darkorange", "steelblue", "red"))  +
  theme_bw() + 
  labs(colour = "Zone", shape = "Season") + 
  scale_shape_manual(values = c(15,17))
```

## Principal coordinate analysis (PCoA) with Bray-Curtis distance per Depth and Season with MDS1 vs MDS3

```r
plot2<-ps3 %>%
  tax_transform("identity") %>% # don't transform!
  dist_calc("bray") %>%
  ord_calc("PCoA") %>%
  ord_plot(axes = c(1, 3), color = "depth", shape ="season",
           size = 3) + 
  scale_colour_manual(values = c("#5e239d", "#00f0b5", "#f61067")) +
  theme_bw() + 
  labs(colour = "Depth", shape = "Season") + 
  scale_shape_manual(values = c(15,17))
```

## Principal coordinate analysis (PCoA) with Unifrac distance per Zone with MDS1 vs MDS2

```r
plot5<-ps3 %>%  
  phyloseq_validate(verbose = FALSE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("gunifrac", gunifrac_alpha = 0.5)  %>%
  ord_calc("auto")  %>%
  ord_plot(color = "zone",
                             size = 3) + 
  scale_colour_manual(values = c("darkorange", "steelblue", "red"))  +
  theme_bw() + 
  labs(colour = "Zone") + 
  scale_shape_manual(values = c(15,17))
```

## Principal coordinate analysis (PCoA) with Unifrac distance per Depth with MDS1 vs MDS2

```r

plot6<-ps3 %>%  
  phyloseq_validate(verbose = FALSE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("gunifrac", gunifrac_alpha = 0.5)  %>%
  ord_calc("auto")  %>%
  ord_plot(color = "depth",
           size = 3) + 
  scale_colour_manual(values = c("#5e239d", "#00f0b5", "#f61067")) +
  theme_bw() + 
  labs(colour = "Depth") + 
  scale_shape_manual(values = c(15,17))
```

## Principal coordinate analysis (PCoA) with Unifrac distance per Season with MDS1 vs MDS2

```r

plot7<-ps3 %>%  
  phyloseq_validate(verbose = FALSE) %>%
  tax_transform("identity", rank = "unique") %>%
  dist_calc("gunifrac", gunifrac_alpha = 0.5)  %>%
  ord_calc("auto")  %>%
  ord_plot(color = "season",
           size = 3) + 
  scale_color_manual(values = c("#606060FF","#D6ED17FF")) +
  theme_bw() + 
  labs(colour = "Season") + 
  scale_shape_manual(values = c(15,17))

```

## PERMANOVA with Unifrac distance

```r
ps3 %>%
  tax_transform("identity") %>% # don't transform!
  dist_calc("gunifrac") %>%
  dist_permanova(
    seed = 99,
    variables = c("zone", "season", "depth"),
    n_processes = 1,
    n_perms = 9999 # only 99 perms used in examples for speed (use 9999+!)
  )
```

# A differential abundance analysis for the comparison

## Load Library 

```r
library("microbiomeMarker")
library("phyloseq")
library(ggplot2)
library(patchwork)
library(dplyr)
library(tibble)
```
## Analysis ALDEx2 temporality 

Differential analysis using ALDEx2 of temporality with the statistical method of glm_anova

```r
ps4<-run_aldex(ps3, group = "season", taxa_rank = "all", 
               method = "glm_anova")

tabla_ps4<-marker_table(ps4) 

tabla_ps4_1<-tabla_ps4$feature
tabla_ps4_2<-tabla_ps4$enrich_group
tabla_ps4_3<-tabla_ps4$ef_F_statistic
tabla_ps4_4<-tabla_ps4$pvalue
tabla_ps4_5<-tabla_ps4$padj

tabla_ps4<-cbind(tabla_ps4_1, tabla_ps4_2, tabla_ps4_3, tabla_ps4_4, tabla_ps4_5) %>%
  as_tibble() %>%
  add_column(comparacion="Dry-Flood") %>%
  rename(feature=tabla_ps4_1) %>%
  rename(enrich_group=tabla_ps4_2) %>%
  rename(ef_F_statistic=tabla_ps4_3) %>%
  rename(pvalue=tabla_ps4_4) %>%
  rename(padj=tabla_ps4_5)

p_abd_ps4 <- plot_abundance(ps4, group = "season", label_level=1)
season<-p_abd_ps4 + 
  scale_fill_manual(values = c("#606060FF","#D6ED17FF"))
```

## Analysis ALDEx2 zone (basin - fringe)

```r
ps5<-ps3 %>%
  subset_samples(zone %in% c("Basin", "Fringe"))
ps6<-run_aldex(ps5, group = "zone", taxa_rank = "all", 
               method = "glm_anova")

tabla_ps6<-marker_table(ps6)

tabla_ps6_1<-tabla_ps6$feature
tabla_ps6_2<-tabla_ps6$enrich_group
tabla_ps6_3<-tabla_ps6$ef_F_statistic
tabla_ps6_4<-tabla_ps6$pvalue
tabla_ps6_5<-tabla_ps6$padj

tabla_ps6<-cbind(tabla_ps6_1, tabla_ps6_2, tabla_ps6_3, tabla_ps6_4, tabla_ps6_5) %>%
  as_tibble() %>%
  add_column(comparacion="MD-ND") %>%
  rename(feature=tabla_ps6_1) %>%
  rename(enrich_group=tabla_ps6_2) %>%
  rename(ef_F_statistic=tabla_ps6_3) %>%
  rename(pvalue=tabla_ps6_4) %>%
  rename(padj=tabla_ps6_5)

p_abd_ps6 <- plot_abundance(ps6, group = "zone",  label_level=1)
basin_fringe<-p_abd_ps6 + 
  scale_fill_manual(values = c("#ff8c00","#FF3333"))
```

## Analysis ALDEx2 depth (40 - 5)

```r
ps7<-ps3 %>%
  subset_samples(depth %in% c("40", "5"))

ps8<-run_aldex(ps7, group = "depth", taxa_rank = "all", 
               method = "glm_anova")

tabla_ps8<-marker_table(ps8)

tabla_ps8_1<-tabla_ps8$feature
tabla_ps8_2<-tabla_ps8$enrich_group
tabla_ps8_3<-tabla_ps8$ef_F_statistic
tabla_ps8_4<-tabla_ps8$pvalue
tabla_ps8_5<-tabla_ps8$padj

tabla_ps8<-cbind(tabla_ps8_1, tabla_ps8_2, tabla_ps8_3, tabla_ps8_4, tabla_ps8_5) %>%
  as_tibble() %>%
  add_column(comparacion="5-40") %>%
  rename(feature=tabla_ps8_1) %>%
  rename(enrich_group=tabla_ps8_2) %>%
  rename(ef_F_statistic=tabla_ps8_3) %>%
  rename(pvalue=tabla_ps8_4) %>%
  rename(padj=tabla_ps8_5)


p_abd_ps8 <- plot_abundance(ps8, group = "depth", label_level=1)
cuarenta_cinco<-p_abd_ps8 + 
  scale_fill_manual(values = c("#00f0b5","#f61067"))
```

## Analysis ALDEx2 depth (5 - 20)

```r
ps9<-ps3 %>%
  subset_samples(depth %in% c("5", "20"))

ps10<-run_aldex(ps9, group = "depth", taxa_rank = "all", 
                method = "glm_anova")

tabla_ps10<-marker_table(ps10)

tabla_ps10_1<-tabla_ps10$feature
tabla_ps10_2<-tabla_ps10$enrich_group
tabla_ps10_3<-tabla_ps10$ef_F_statistic
tabla_ps10_4<-tabla_ps10$pvalue
tabla_ps10_5<-tabla_ps10$padj

tabla_ps10<-cbind(tabla_ps10_1, tabla_ps10_2, tabla_ps10_3, tabla_ps10_4, tabla_ps10_5) %>%
  as_tibble() %>%
  add_column(comparacion="5-20") %>%
  rename(feature=tabla_ps10_1) %>%
  rename(enrich_group=tabla_ps10_2) %>%
  rename(ef_F_statistic=tabla_ps10_3) %>%
  rename(pvalue=tabla_ps10_4) %>%
  rename(padj=tabla_ps10_5)

p_abd_ps10 <- plot_abundance(ps10, group = "depth", label_level=1)
cinco_veinte<-p_abd_ps10 + 
  scale_fill_manual(values = c("#5e239d", "#f61067"))
```

## Analysis ALDEx2 zone (basin - impaired)

```r
ps11<-ps3 %>%
  subset_samples(zone %in% c("Basin", "Impaired"))

ps12<-run_aldex(ps11, group = "zone", taxa_rank = "all", 
                method = "glm_anova")

tabla_ps12<-marker_table(ps12)

tabla_ps12_1<-tabla_ps12$feature
tabla_ps12_2<-tabla_ps12$enrich_group
tabla_ps12_3<-tabla_ps12$ef_F_statistic
tabla_ps12_4<-tabla_ps12$pvalue
tabla_ps12_5<-tabla_ps12$padj

tabla_ps12<-cbind(tabla_ps12_1, tabla_ps12_2, tabla_ps12_3, tabla_ps12_4, tabla_ps12_5) %>%
  as_tibble() %>%
  add_column(comparacion="MD-D") %>%
  rename(feature=tabla_ps12_1) %>%
  rename(enrich_group=tabla_ps12_2) %>%
  rename(ef_F_statistic=tabla_ps12_3) %>%
  rename(pvalue=tabla_ps12_4) %>%
  rename(padj=tabla_ps12_5)

p_abd_ps12 <- plot_abundance(ps12, group = "zone",  label_level=1)
basin_impaired<-p_abd_ps12 + 
  scale_fill_manual(values = c("#ff8c00","#4682b4"))
```

## Analysis ALDEx2 zone (fringe - impaired)

```r
ps13<-ps3 %>%
  subset_samples(zone %in% c("Fringe", "Impaired"))

ps14<-run_aldex(ps13, group = "zone", taxa_rank = "all", 
                method = "glm_anova")

tabla_ps14<-marker_table(ps14)

tabla_ps14_1<-tabla_ps14$feature
tabla_ps14_2<-tabla_ps14$enrich_group
tabla_ps14_3<-tabla_ps14$ef_F_statistic
tabla_ps14_4<-tabla_ps14$pvalue
tabla_ps14_5<-tabla_ps14$padj

tabla_ps14<-cbind(tabla_ps14_1, tabla_ps14_2, tabla_ps14_3, tabla_ps14_4, tabla_ps14_5)%>%
  as_tibble() %>%
  add_column(comparacion="ND-D") %>%
  rename(feature=tabla_ps14_1) %>%
  rename(enrich_group=tabla_ps14_2) %>%
  rename(ef_F_statistic=tabla_ps14_3) %>%
  rename(pvalue=tabla_ps14_4) %>%
  rename(padj=tabla_ps14_5)


p_abd_ps14 <- plot_abundance(ps14, group = "zone",  label_level=1)
fringe_impaired<-p_abd_ps14 + 
  scale_fill_manual(values = c("#FF3333","#4682b4"))
```

## Combine multiple plots

```r
 (basin_fringe + basin_impaired ) / (fringe_impaired + season) /
  (cinco_veinte + cuarenta_cinco) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = 'collect') +
  plot_layout(axis_titles = "collect")
```
