## Inferencia de variantes ------------------------------------------------#### 
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

## Inspeccionando el objeto dada ------------------------------------------####
head(dadaFs[[1]])

## Fusionando lecturas pareadas
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

## Inspeccionar los datos de fusión de la primera muestra
head(mergers[[1]])

