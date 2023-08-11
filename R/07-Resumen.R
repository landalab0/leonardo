# Resumen del analisis ----------------------------------------------------####
getN <- function(x) sum(getUniques(x))

track<- cbind (out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers,getN), 
               rowSums(seqtab.nochim))


colnames(track)<- c("input", 
                    "filtered",
                    "denoisedF", 
                    "denoisedR", 
                    "merged",
                    "nochim")

rownames(track) <- sample.names

head(track)