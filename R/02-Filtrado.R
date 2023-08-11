## Asignar nombres a los archivos para los fastq filtrados ----------------####
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Filtrado ---------------------------------------------------------------####
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(275, 235),
                     maxN = 0,
                     maxEE = c(5, 4),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)

## Mostrar los primeros resultados ----------------------------------------####
head(out)
