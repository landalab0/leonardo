# Remover Chimeras --------------------------------------------------------####
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method ="consensus", 
                                    multithread=TRUE, 
                                    verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
