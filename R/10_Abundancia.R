###---Filtrado por probables falsos positivos--------------------------########

# Definimos taxa a filtrar
filterPhyla = c("Altiarchaeota", "Cloacimonadota", "FCPU426", "Iainarchaeota",
                "Sumerlacota", "Thermotogota", "WOR-1", NA)

# Procedemos a filtrar
(psd1 = subset_taxa(ps, !Phylum %in% filterPhyla))

# Removiendo taxa no correspondiente como cloroplastos, mitocondrias y eukariota

filterPhyla2 <- c("Chloroplast", "Mitochondria", "Eukaryota")
psd1 <- subset_taxa(psd1, !Kingdom %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Phylum %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Class %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Order %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Family %in% filterPhyla2)
psd1 <- subset_taxa(psd1, !Genus %in% filterPhyla2)

# Seleccionamos las taxa de interés
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(psd1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Agregamos una línea para nuestro umbral
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

####-------UMBRAL DE PREVALENCIA-----------------------------------------######
(prevalenceThreshold = 0.05 * nsamples(psd1))

# Executando el filtro de prevalencia ---------------------------------########
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
(psd2 = prune_taxa(keepTaxa, psd1))

# Reemplazamos las secuencias por un nombre genérico-------------------------###
taxa_names(psd2) <- paste0("ASV", seq(ntaxa(psd2)))

sample_sum_df <- data.frame(sum = sample_sums(psd2))

ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "grey", binwidth = 15) +
  ggtitle("Distribution de largo de secuenciacion") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) 