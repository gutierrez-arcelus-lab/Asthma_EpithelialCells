
meta_all <- readRDS("~/Documents/AADCRC/EpithelialRV/allDonor/results/D1-2-3-24h_qc_pca_umap_clust_harmony.rds")
#read scDRS results 
filenames_Scores <- Sys.glob("~/Documents/AADCRC/EpithelialRV/allDonor/results/scDRS_result/*.score")

Files_Scores <- lapply(filenames_Scores, read.table, sep="\t", header = T, row.names = 1)

#substract specific columns 
sub_meta_all <- metadata[,c("cell", "nUMI", "nGene", "percent_mito", "percent.ribo", "condition", "UMAP1", "UMAP2", "cell_type", "harmonized_UMAP1", "harmonized_UMAP2")]
head(sub_meta_all)

Files_Umap <- lapply(names(Files_Scores), function(x){
  files <- Files_Scores[[x]][meta_all$cell,]
  return(files)
})

names(Files_Umap) <- names_scores

Files_Umap <- lapply(names(Files_Umap), function(x){
  combine <- cbind(Files_Umap[[x]], sub_meta_all[,c("UMAP1", "UMAP2", "harmonized_UMAP1", "harmonized_UMAP2", "condition", "cell_type")])
  return(combine)
})

names(Files_Umap) <- names_scores

plot_list_trait_pval <- lapply(names(Files_Umap),function(x){
  plot <-ggplot(filter(Files_Umap[[x]], pval >= 0.05),
                aes(UMAP1, UMAP2)) +
    geom_point(size = 0.5, color = "white") +
    geom_point(data = filter(Files_Umap[[x]], pval < 0.05),
               aes(color = norm_score), size = 0.5) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    theme_minimal() +
    ggtitle(x) + 
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey"))
  return(plot)
})

names(plot_list_trait_pval) <- names_scores

plot_list_pvalue <- lapply(names(Files_Umap),function(var){
  plot <- ggplot(data=Files_Umap[[var]], aes_string(x='pval')) + 
    geom_histogram(color="black", fill="white") +
    ggtitle(var) +
    labs(x="pvalue") +
    theme_classic()
  return(plot)
})

grid.arrange(grobs=plot_list_pvalue)

library(qvalue)

#Calculate qvalues

pval <- lapply(Files_Umap, function(x) x$pval)
qobj <- lapply(pval, function(x) qvalue(p=x))
qvalues <- lapply(qobj, function(x) x$qvalues)

print("0.1")
lapply(qvalues, function(x) table(x < 0.1))
print("0.2")
lapply(qvalues, function(x) table(x < 0.2))


Files_Umap_qval <- lapply(names(Files_Umap), function(x) {
  combine <- cbind(Files_Umap[[x]], qval = qvalues[[x]])
  return(combine)
})

names(Files_Umap_qval) <- names(Files_Umap)


plot_list_trait_qval <- lapply(names(Files_Umap_qval),function(x){
  plot <-ggplot(filter(Files_Umap_qval[[x]], qval >= 0.1),
                aes(harmonized_UMAP1, harmonized_UMAP2)) +
    geom_point(size = 0.5, color = "white") +
    geom_point(data = filter(Files_Umap_qval[[x]], qval < 0.2),
               aes(color = norm_score), size = 0.5) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    theme_minimal() +
    ggtitle(x) + 
    coord_fixed() +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey")) +
    xlab("Harmonized UMAP1") + ylab("Harmonized UMAP2") 
  return(plot)
})

names(plot_list_trait_qval) <- names(Files_Umap)

grid.arrange(grobs=plot_list_trait_qval)

plot <- grid.arrange(plot_list_trait_qval$`UKB-asthma`,
                     plot_list_trait_qval$COA, 
                     plot_list_trait_qval$`Allergy-Eczema`, 
                     plot_list_trait_qval$AOA)

plot_list_trait_qval <- lapply(names(Files_Umap_qval),function(x){
  Files_Umap_qval[[x]] <- Files_Umap_qval[[x]][order(Files_Umap_qval[[x]]$norm_score),]
  plot <-ggplot(filter(Files_Umap_qval[[x]], qval >= 0.1),
                aes(harmonized_UMAP1, harmonized_UMAP2)) +
    geom_point(alpha= 0.5, size = 0.5, color = "grey") +
    geom_point(data = filter(Files_Umap_qval[[x]], qval < 0.1),
               aes(color = norm_score), size = 1) +
    scale_color_gradient(low = "#8F8CEC", high = "#9a3b98", limits=c(2.492536, 7.629132), space = "Lab", guide="colourbar") +
    labs(color = "Disease\nRelevant\nScore")+
    theme_minimal() +
    ggtitle(x) + 
    theme_void() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size=15),
      legend.text = element_text(size=20),
      legend.title = element_text(size=20)
    ) + 
    coord_fixed()
  return(plot)
})

filter_cells <- lapply(Files_Umap_qval, function(x){
  filter(x, qval < 0.1)
})

filter_cells_summary <- lapply(filter_cells, function(x) summary(x$norm_score))
min(unlist(lapply(filter_cells_summary, function(x) x[["Min."]])), na.rm = T)
max(unlist(lapply(filter_cells_summary, function(x) x[["Max."]])), na.rm = T)

melt_filter_cells_asthma <- list(filter_cells$`Adult-Onset Asthma`, 
                                 filter_cells$`Childhood-Onset Asthma`,
                                 filter_cells$`Allergy/Eczema`,
                                 filter_cells$`All asthma`)
names(melt_filter_cells_asthma) <- c("Adult-Onset Asthma", "Childhood-Onset Asthma", "Allergy/Eczema", "All asthma")


melt_filter_cells <- melt(melt_filter_cells_asthma)

melt_filter_cells_per <- melt_filter_cells %>% filter(variable == "norm_score") %>% group_by(L1, cell_type) %>%
  summarise(Nb = n()) %>% 
  mutate(C = sum(Nb)) %>% 
  mutate(percent = Nb/C*100)

#melt_filter_cells_per$cell_type <- factor(melt_filter_cells_per$cell_type, levels=c("Ciliated", "Deuterosomal", "Rare", "Basal","Secretory", "Transitional"))

ggplot(melt_filter_cells_per, aes(x=L1, y= percent, fill=cell_type)) +
  geom_bar(color="white", stat = "identity", width = 0.7)+
  geom_text_repel(aes(label = paste(round(percent),"%")), 
                  position = position_stack(vjust = 0.5), 
                  direction = "y",
                  box.padding = unit(0.01, "lines"),
                  size=7, 
                  max.overlaps = Inf) +
  theme_classic()+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=22), 
        legend.title = element_blank(), 
        legend.text = element_text(size=30),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22)) +
  scale_fill_manual(breaks = c("Ciliated", "Deuterosomal", "Rare", "Basal","Secretory", "Transitional"),
                    values = c("#009593", "#7AAF97", "#8a817c", "#D35C79", "#D69481", "#DEC0A3")) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  xlab("") +
  ylab("% of cells")

ggsave("~/Documents/AADCRC/EpithelialRV/allDonor/results/barplot-asthma_per.tiff", dpi = 500, width = 14, height = 10)

barplot_allcell <- meta_all %>% group_by(cell_type) %>%
  summarise(Nb=n()) %>% 
  mutate(C = sum(Nb)) %>% 
  mutate(percent = Nb/C*100)

barplot_allcell$L1 <- "All cells"

barplot_allcell$cell_type <- factor(barplot_allcell$cell_type, levels=c("Ciliated", "Deuterosomal", "Rare", "Basal","Secretory", "Transitional"))

ggplot(barplot_allcell, aes(x = "", y= percent, fill=cell_type)) +
  geom_bar(color="white", stat = "identity", width = 0.5)+
  geom_text_repel(aes(label = paste(round(percent),"%")), 
                  position = position_stack(vjust = 0.5),
                  box.padding = unit(0.01, "lines"),
                  size=7, 
                  max.overlaps = Inf) +
  theme_classic()+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=22), 
        legend.title = element_blank(), 
        legend.text = element_text(size=30),
        axis.text.x = element_text(size=22),
        axis.text.y = element_text(size=22),
        legend.position = "none") +
  scale_fill_manual(breaks = c("Ciliated", "Deuterosomal", "Rare", "Basal","Secretory", "Transitional"),
                    values = c("#009593", "#7AAF97", "#8a817c", "#D35C79", "#D69481", "#DEC0A3")) +
  xlab("All cells") +
  ylab("% of cells")

ggsave("~/Documents/AADCRC/EpithelialRV/allDonor/results/barplot-all_cells_asthma_per.tiff", dpi = 500, width = 11, height = 10)

library(ggpp)

all_data_asthma <- rbind(barplot_allcell, melt_filter_cells_per)
all_data_asthma$L1 <- ifelse(all_data_asthma$L1 == "Childhood-Onset Asthma", "Childhood-Onset\nAsthma", all_data_asthma$L1)

all_data_asthma$L1 <- factor(all_data_asthma$L1, levels=c("All cells", "Childhood-Onset\nAsthma", "Allergy/Eczema", "All asthma"))
all_data_asthma$color <- ifelse(all_data_asthma$L1 == "All cells", "white", "black")

ggplot(all_data_asthma, aes(x=L1, y= percent, fill=cell_type, label= paste0(round(percent), "%"))) +
  geom_bar(color= all_data_asthma$color, stat = "identity", width = 0.45, alpha=1, size=1.5)+
  theme_classic()+
  theme(axis.title = element_text(size=30),
        legend.title = element_blank(), 
        legend.position = "none",
        #legend.text = element_blank(),
        #legend.text = element_text(size=18),
        axis.text.x = element_text(size=30),
        axis.text.y = element_text(size=30)) +
  scale_fill_manual(breaks = c("Ciliated", "Deuterosomal", "Rare", "Basal","Secretory", "Transitional"),
                    values = c("#009593", "#7AAF97", "#8a817c", "#D35C79", "#D69481", "#DEC0A3")) +
  scale_x_discrete(position = "bottom")+
  coord_flip() +
  xlab("") +
  ylab("% of cells") +
  ggrepel::geom_text_repel(
    position = position_stacknudge(vjust=0.5, y=-0.5,
                                   x = ifelse(all_data_asthma$percent < 3, 0.5, 0), 
                                   kept.origin = "stacked"),
    #hjust=0.5,
    size = 8, 
    box.padding = unit(0.01, "lines"))

ggsave("~/Documents/AADCRC/EpithelialRV/allDonor/results/barplot-all_cells_asthma_per.tiff", dpi = 500, width = 15, height = 8)


ggplot() +
  geom_point(
    data = meta_all[sample(nrow(meta_all)),],
    mapping = aes_string(x = "harmonized_UMAP1", y = "harmonized_UMAP2", fill = "cell_type"),
    size = 1, stroke = 0.0001, shape = 21, alpha=0.5) +
  guides(fill = guide_legend(override.aes = list(size = 7, alpha=1))) +
  scale_fill_manual(breaks = c("Ciliated", "Deuterosomal", "Rare", "Basal","Secretory", "Transitional"),
                    values = c("#009593", "#7AAF97", "#8a817c", "#D35C79", "#D69481", "#DEC0A3"), name = "") +
  labs(
    x = "UMAP1",
    y = "UMAP2",
  ) +
  theme_bw(base_size =15) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom",
  ) + 
  coord_fixed()


barplot <- lapply(names(Files_Umap_qval), function(x) {
  Files_Umap_qval[[x]]$condition <- factor(Files_Umap_qval[[x]]$condition, levels=c("no-virus","24h"))
  filter_cells <- Files_Umap_qval[[x]][which(Files_Umap_qval[[x]]$qval < 0.1),]
  ggplot(filter_cells, aes(x=factor(condition), fill=cell_type)) +
    geom_bar(stat="count")+
    xlab("") + ylab("") + 
    ggtitle(x) +
    theme_bw () + 
    theme(text = element_text(size = 20),
          axis.text = element_text(size = 20))      
})

names(barplot) <- names_scores

grid.arrange(barplot$`UKB-asthma`, barplot$COA, barplot$`Allergy-Eczema`, ncol=3)



