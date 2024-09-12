# Create a database from green genes to use in R 

## Load libraries

´´´r

library(Biostrings)
library(ShortRead)
library(dada2)

´´´

## Path where the created file will appear

´´´r

setwd("~/Desktop/GG2/final")

´´´

## Greengenes database download

´´´r

download.file("http://ftp.microbio.me/greengenes_release/current/2022.10.backbone.full-length.fna.qza", 
              "2022.10.backbone.full-length.fna.qza")
download.file("http://ftp.microbio.me/greengenes_release/current/2022.10.backbone.tax.qza",
              "2022.10.backbone.tax.qza")

´´´

## Decompress the downloaded database

´´´r

unzip("2022.10.backbone.full-length.fna.qza")
unzip("2022.10.backbone.tax.qza")
fn <- "a53d9300-5c5c-4774-a2e8-a5e23904f1ae/data/dna-sequences.fasta"
txfn <- "c16a953c-f24d-4d14-927c-40d90ced395e/data/taxonomy.tsv"

´´´

## Taxonomy is in the second column: ID and taxonomy. Taxonomy is at 7 levels: domain -- species

´´´r

sq <- getSequences(fn)
tdf <- read.csv(txfn, sep="\t", header=TRUE)
tax <- tdf[,2]
names(tax) <- tdf[,1]

identical(names(sq), names(tax)) # TRUE

´´´

## Analysis of identification chain taxonomies

All taxonomies are within 7 levels.
Note: A significant number of unassigned taxonomic levels here are coded as, e.g., 'g__'.
Note: The species level has the genus name duplicated, i.e., the species level has genus [SPACE] binomial species, instead of just species.

´´´r

taxes <- strsplit(tax, "; ")
tax.depth <- sapply(taxes, length)
table(tax.depth) 

´´´

## Removing the genus from the species assignment

´´´r

for(i in seq(length(taxes))) {
  gen <- taxes[[i]][[6]]
  gen <- substr(gen, 4, nchar(gen))
  taxes[[i]][[7]] <- gsub(gen, "", taxes[[i]][[7]])
  taxes[[i]][[7]] <- gsub("__ ", "__", taxes[[i]][[7]])
}

´´´
## Identification of unassigned taxonomic levels

Note: A relative number of entries have lower taxonomic levels assigned, although higher levels are not assigned,
e.g., MJ030-2-barcode58-umi83452bins-ubs-6, which is assigned as s__Spirochaeta aurantia, but with o__;f__;g__ for the order/family/genus designation. 
These assignments will be removed with the assignTaxonomy test

´´´r

tax_pre <- c("d__", "p__", "c__", "o__", "f__", "g__", "s__")
is.unassigned <- sapply(taxes, function(tx) {
  tx == tax_pre
}) |> t()
tax.depth <- apply(is.unassigned, 1, function(isu) { min(which(isu)-1L, 7L) })
tax.ids <- sapply(seq_along(taxes), function(i) {
  td <- tax.depth[[i]]
  id.str <- paste(taxes[[i]][1:td], collapse=";")
  id.str <- paste0(id.str, ";") # Add terminal semicolon
  id.str
})
names(tax.ids) <- names(taxes)

´´´

## Creating the FASTA file

´´´r

sq.out <- sq
names(sq.out) <- tax.ids
writeFasta(sq.out, "greengenes2_trainset.fa.gz", compress=TRUE)

´´´

## Testing the results of assignTaxonomy on this recently created file

´´´r

dada2:::tax.check("greengenes2_trainset.fa.gz", fn.test=system.file("extdata", "ten_16s.100.fa.gz", package="dada2"))

´´´
