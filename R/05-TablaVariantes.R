## Construir la tabla de secuencia (Una version nueva de la OTU table) ----####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

## Inspeccionando los largos de la secuencia de distribución --------------####
table(nchar(getSequences(seqtab)))

