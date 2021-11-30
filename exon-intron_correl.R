tiff("exon-intron-correl-plot.tiff", width =10, height = 6.5, units = "in",res = 400,compression = "lzw")
exon_libraries <- c(paste0("Exon_Total_S",seq(1,3)),paste0("Exon_PolyA_S",seq(1,3)))
intron_libraries <- c(paste0("Intron_Total_S",seq(1,3)),paste0("Intron_PolyA_S",seq(1,3)))
 
p <- list()
for (i in 2:7){
plt <- ggplot(data.frame(all_datasets), aes_string(x = log2(all_datasets[,i] + 0.01), y = log2(all_datasets[,i+6] + 0.01))) +
  geom_point(size=0.5) +  xlab(exon_libraries[i-1]) + ylab(intron_libraries[i-1]) +
  scale_y_continuous(breaks = seq(0, 20, by = 5)) +
  scale_x_continuous(breaks = seq(0, 50, by = 5)) +
  theme_bw() + labs(fill = "Level") +
  theme(axis.text=element_text(size=10),panel.grid.major = element_line(colour = "grey70", size = 0.2),
       panel.grid.minor = element_blank(),
       legend.title = element_text(color = "black", size = 6),
       legend.text = element_text(color = "black", size = 6),
       legend.position = c(0.91,0.15),
       legend.key.size = unit(0.2, "cm"),
       legend.background = element_rect(fill='transparent')) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + scale_fill_viridis(option="plasma",alpha = 0.5) +
  stat_cor(method = "spearman",size=4)
  print(plt)
p[[i-1]] <- plt
}
do.call(gridExtra::grid.arrange,c(p, ncol=3))
dev.off()
