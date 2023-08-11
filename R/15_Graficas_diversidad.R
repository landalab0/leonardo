#######-------------------Riqueza zona----------------------#############

plot_richness(psd2, color = "zone", x = "zone", measures = c("Observed", "Chao1", "Shannon")) +
  geom_boxplot(aes(fill = zone), alpha = 0.7) +
  scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) +
  scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f"))

#######-------------------Riqueza temporada----------------------#############

plot_richness(psd2, color = "season", x = "season", measures = c("Observed", "Chao1", "Shannon")) +
  geom_boxplot(aes(fill = season), alpha = 0.7) +
  scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) +
  scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f"))

#######-------------------Riqueza profundidad----------------------#############
plot_richness(psd2, color = "depht", x = "depth", measures = c("Observed", "Chao1", "Shannon")) +
  geom_boxplot(aes(fill = depth), alpha = 0.7) +
  scale_color_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f")) +
  scale_fill_manual(values = c("#a6cee3", "#b2df8a", "#fdbf6f"))

