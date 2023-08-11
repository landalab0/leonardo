# Esta es la tabla de cuentas o read counts
head(otu_table(psd2)) 
  kable(format = "html", col.names = colnames(otu_table(psd2))) 
  kable_styling() 
  kableExtra::scroll_box(width = "100%", height = "350px")

# Esta es la tabla de taxonomía
head(tax_table(psd2)) 
  kable(format = "html", col.names = colnames(tax_table(psd2)))
  kable_styling() 
  kableExtra::scroll_box(width = "100%", height = "320px")

# Esta es la metadata
as(sample_data(psd2), "data.frame") -> metad

metad 
  kable(format = "html", col.names = colnames(metad)) 
  kable_styling() 
  kableExtra::scroll_box(width = "100%", height = "400px")

# Esta es la filogenia asociada a las taxa en nuestro objeto phyloseq
summarize_phyloseq(psd2)

###----metadata_vs_taxonomia_vs_abundacia_de taxon mas abundante###############
df <- psmelt(psd2)

head(df) 
  kable(format = "html", col.names = colnames(df)) 
  kable_styling() 
  kableExtra::scroll_box(width = "100%", height = "400px")
