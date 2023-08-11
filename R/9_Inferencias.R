file_path <- "~/Leonardo/metadata_citlalli.csv"
csv_content <- readLines(file_path)
csv_content_clean <- gsub("\uFFFD", "", csv_content)
samdf <- read.csv(text = csv_content_clean, header = TRUE, row.names = 4, sep = ",")
samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps)
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps