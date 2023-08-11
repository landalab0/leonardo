##############--------de acuerdo a Zona--------------------------#######
pro <- plot_frequencies(sample_data(psd2), "zone", "seqR1")
print(pro$plot)

# En formato tabla
pro$data 
kable(x = pro$data, format = "html", col.names = colnames(pro$data), digits = 2)
kable_styling() 
  kableExtra::scroll_box(width = "100%", height = "300px")
# Transformamos las cuentas en porcentaje
psd2r  = transform_sample_counts(psd2, function(x) x / sum(x) )

tax(psd2r)[,6] 
