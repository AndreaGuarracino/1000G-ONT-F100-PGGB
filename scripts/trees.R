library(tidyverse) # install.packages("tidyverse")
library(dendextend) # for tanglegram. install.packages('dendextend')
library(ape) # for phylogenetic analysis. install.packages('ape')
library(ggtree) # for phylogenetic tree plotting and annotation. BiocManager::install("ggtree")
library(ggrepel) # install.packages("ggrepel")

# Assign arguments
path_dist_tsv <- '/home/guarracino/Desktop/Garrison/1000G-ONT-F100-PGGB/1000G-30kbp+HPRC.chr20.haplotype.dist.tsv'
plot_title <- '"1000G-ONT + HPRC" chr20 pangenome'
prefix <- '1000G-30kbp+HPRC.chr20'
path_dist_tsv <- '/home/guarracino/Desktop/Garrison/1000G-ONT-F100-PGGB/1000G-30kbp+HPRC.chr20.clean.haplotype.dist.tsv'
plot_title <- '"1000G-ONT + HPRC" chr20 pangenome without euchromatic, non-centromeric regions'
prefix <- '1000G-30kbp+HPRC.chr20.clean'

path_dist_tsv <- '/home/guarracino/Desktop/Garrison/1000G-ONT-F100-PGGB/1000G-30kbp+HPRC.chr6.MHC.haplotype.dist.tsv'
plot_title <- '"1000G-ONT + HPRC" MHC pangenome'
prefix <- '1000G-30kbp+HPRC.chr6.MHC'

path_dist_tsv <- '/home/guarracino/Desktop/Garrison/1000G-ONT-F100-PGGB/1000G-ONT+HPRC.chr16.fa.gz.445f03b.c2fac19.432148a.smooth.final.dist.haplotype.tsv'
plot_title <- '"1000G-ONT + HPRC" chr16 pangenome'
prefix <- '1000G-ONT+HPRC.chr16'
path_dist_tsv <- '/home/guarracino/Desktop/Garrison/1000G-ONT-F100-PGGB/1000G-ONT+HPRC.chr16.p98.s20k.k79.clean.haplotype.dist.tsv'
plot_title <- '"1000G-ONT + HPRC" chr16 pangenome without low-/high-depth regions'
prefix <- '1000G-ONT+HPRC.chr16.clean'

path_dist_tsv <- '/home/guarracino/Desktop/Garrison/1000G-ONT-F100-PGGB/1000G-ONT+HPRC.chr20.fa.gz.445f03b.c2fac19.432148a.smooth.final.dist.haplotype.tsv'
plot_title <- '"1000G-ONT + HPRC" chr20 pangenome'
prefix <- '1000G-ONT+HPRC.chr20'
path_dist_tsv <- '/home/guarracino/Desktop/Garrison/1000G-ONT-F100-PGGB/1000G-ONT+HPRC.chr20.p98.s20k.k79.clean.haplotype.dist.tsv'
plot_title <- '"1000G-ONT + HPRC" chr20 pangenome without euchromatic, non-centromeric fraction'
prefix <- '1000G-ONT+HPRC.chr20.clean'

# Read matrices
meta_df <- read.table("/home/guarracino/Dropbox/1000G-ONT-F100-PGGB/data/graphs_1000g+hprc.populations.groups.tsv", header = F)
colnames(meta_df) <- c('Sample', 'Group')
meta2_df <- read.table("/home/guarracino/Dropbox/1000G-ONT-F100-PGGB/data/graphs_1000g+hprc.populations.tsv", header = T)
colnames(meta2_df) <- c('Sample', 'Sex', 'Superpopulation')
meta_df <- merge(meta_df, meta2_df)
rm(meta2_df)

jaccard_dist_df <- read_tsv(path_dist_tsv) %>%
  filter(!group.a %in% c('grch38', 'chm13') & !group.b %in% c('grch38', 'chm13')) %>%
  arrange(group.a, group.b) %>%
  select(group.a, group.b, jaccard.distance) %>%
  pivot_wider(names_from = group.b, values_from = jaccard.distance) %>%
  column_to_rownames(var = "group.a")

jaccard_hc <- as.dist(jaccard_dist_df) %>% hclust()

# Fit
fit <- cmdscale(as.dist(jaccard_dist_df), eig=TRUE, k = 5)

fit_df <- as.data.frame(fit$points) %>%
  rownames_to_column(var = "Name")

# Create a new column in fit_df for the prefix
fit_df$SamplePrefix <- sapply(strsplit(as.character(fit_df$Name), "#"), `[`, 1)
fit_and_meta_df <- merge(fit_df, meta_df, by.x = "SamplePrefix", by.y = "Sample")
#fit_and_meta_df <- fit_and_meta_df %>%
#  filter(Group %in% c('REFERENCE'))

fit_and_meta_df$Label <- fit_and_meta_df$Name
fit_and_meta_df$Label[fit_and_meta_df$Group != 'REFERENCE'] <- ''


plotD1D2 <- ggplot(
  data = fit_and_meta_df,
  aes(x = V1, y = V2, label = Label, shape = Group, color = Superpopulation)
) +
  #scale_color_manual(values=haplotype_colors) +
  geom_point(size=2, alpha=0.7) + 
  xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
  ylab(paste0("Dimension 2 (", trunc(fit$eig[2]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
  theme_bw() +
  ggtitle(plot_title)
plotD1D2 <- plotD1D2 + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position="right"
  ) + geom_text_repel(
    size=3,
    max.iter=1000,
    max.time=2,
    show.legend=FALSE, # to hide the `a` from the legend
    max.overlaps=10#Inf
  )
plotD1D2
ggsave(paste0(prefix, ".PCA.group+superpop.pdf"), plotD1D2, width = 10, height = 6, dpi = 300)


plotD1D2 <- ggplot(
  data = fit_and_meta_df,
  aes(x = V1, y = V2, label = Label, shape = Group, color = Group)
) +
  #scale_color_manual(values=haplotype_colors) +
  geom_point(size=2, alpha=0.7) + 
  xlab(paste0("Dimension 1 (", trunc(fit$eig[1]/sum(fit$eig) * 100 * 10^2)/10^2, "%)"))+
  ylab(paste0("Dimension 2 (", trunc(fit$eig[2]/sum(fit$eig) * 100 * 10^2)/10^2, "%)")) +
  theme_bw() +
  ggtitle(plot_title)
plotD1D2 <- plotD1D2 + 
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position="right"
  ) + geom_text_repel(
    size=3,
    max.iter=1000,
    max.time=2,
    show.legend=FALSE, # to hide the `a` from the legend
    max.overlaps=Inf
  )
plotD1D2
ggsave(paste0(prefix, ".PCA.only-ref-labels.png"), plotD1D2, width = 10, height = 6, dpi = 300)



if (FALSE){
  # Plot
  png("path/to/save/your_plot.png", width = 800, height = 600)
  plot(
    jaccard_hc,
    # Label at same height
    hang = -1,
    main = plot_title,
    xlab = 'Haplotype',
    ylab = 'Jaccard distance',
    sub = '',
    cex = 1.8,       # Adjusts the size of points/text in the plot
    cex.lab = 1.6,   # Adjusts the size of x and y labels
    cex.axis = 1.6,  # Adjusts the size of axis text
    cex.main = 1.6,  # Adjusts the size of the main title
    cex.sub = 1.6,    # Adjusts the size of the subtitle
    lwd = 2  # Increases the width of the branch lines
  )
  dev.off()
  
  png("path/to/save/your_plot.png", width = 800, height = 600)
  # Plotting the dendrogram
  plot(jaccard_hc, hang = -1, main = plot_title, xlab = 'Haplotype', ylab = 'Jaccard distance',
       sub = '', cex = 1.8, cex.lab = 1.6, cex.axis = 1.6, cex.main = 1.6, cex.sub = 1.6, lwd = 2)
  # Identifying clusters (e.g., cutting the tree into 3 clusters)
  clusters <- cutree(jaccard_hc, k = 3)
  # Adding colored rectangles
  rect.hclust(jaccard_hc, k = 3, border = 2:4)  # Here, 'border' specifies the colors
  dev.off()
}
