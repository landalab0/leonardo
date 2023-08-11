## Aprender de las tasas de error -----------------------------------------####
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

## Graficar los errores ---------------------------------------------------####
plotErrors(errF, nominalQ = TRUE)
