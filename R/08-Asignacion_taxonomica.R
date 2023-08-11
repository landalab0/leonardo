## Asignando Taxonomia ----------------------------------------------------####
taxa <- assignTaxonomy(seqtab.nochim, 
        "~/Leonardo/PIPELINE_DADA_2/MiSeq_SOP/silva_nr99_v138.1_train_set.fa", 
                       multithread = TRUE)
taxa <- assignTaxonomy(seqtab.nochim,
"~/Leonardo/PIPELINE_DADA_2/MiSeq_SOP/silva_nr99_v138.1_wSpecies_train_set.fa",
                       multithread = TRUE)
taxa <- addSpecies(taxa, 
  "~/Leonardo/PIPELINE_DADA_2/MiSeq_SOP/silva_species_assignment_v138.1.fa")

## Eliminando nombre de las filas de secuencias para visualización --------####
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

