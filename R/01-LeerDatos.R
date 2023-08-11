library(dada2)
packageVersion("dada2")

## Ruta de los archivos ---------------------------------------------------####
path <- "~/citl_fastq_data_2/"
list.files(path)

## Archivos de las secuencias ---------------------------------------------####
fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

## Nombres de las muestras ------------------------------------------------####
sample.names <- sapply(strsplit(basename(fnFs), "_R"), `[`, 1)

## Ver los archivos y nombres de las muestras -----------------------------####
cat("Los archivos de las secuencias son:\n")
print(fnFs)
cat("\nLos nombres de las muestras son:\n")
print(sample.names)

## Graficar la calidad de los archivos R1 y R2 ----------------------------####
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])
