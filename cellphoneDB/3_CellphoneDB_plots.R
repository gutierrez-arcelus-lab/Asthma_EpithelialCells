library(Seurat)

markers <- read.table("markers_ciliatedInfected_vs_CiliatedNonInfected_04292024.txt")

markers <- markers[order(markers$avg_log2FC, decreasing = T),]
markers$gene <- rownames(markers)
markers_up <- markers[which(markers$avg_log2FC > 0.5), ]
markers_up <- markers_up[which(markers_up$p_val_adj < 0.05), ]
markers_up <- markers_up[which(markers_up$pct.1 > 0.1), ]
markers_up <- markers_up[which(markers_up$pct.2 > 0.1), ]

markers_up$cluster <- "ciliatedInfected"

data_cpdb <- read.table("degs_analysis_relevant_interactions_05_02_2024_131803.txt", sep = "\t", header = TRUE)

means_cpdb <- read.table("degs_analysis_means_05_02_2024_131803.txt", sep = "\t", header = TRUE)

library(ktplots)

colnames(means_cpdb) <- gsub(pattern = "\\.", "|", colnames(means_cpdb))
colnames(data_cpdb) <- gsub(pattern = "\\.", "|", colnames(data_cpdb))

library(ggplot2)

plot_cpdb(
    scdata=seu_obj,
    cell_type1=".",
    cell_type2=".",
    celltype_key="cat",
    means=means_cpdb,
    pvals=data_cpdb,
    degs_analysis = TRUE,
    keep_significant_only = FALSE) +
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))

#ggsave("cellphoneDBres_l2fc0.5_allSig.tiff",  width = 10, height = 14)

data_cpdb_filter <- data_cpdb[data_cpdb$directionality=="Ligand-Receptor",]
data_cpdb_filter <- data_cpdb_filter[data_cpdb_filter$gene_a %in% markers_up$gene,]
data_cpdb_filter <- data_cpdb_filter[which(data_cpdb_filter$`ciliatedInfected|othersnonInfected` == 1),]
#data_cpdb_filter <- data_cpdb_filter[,c(1:13, 17)]

mean_filter <- means_cpdb[which(means_cpdb$id_cp_interaction %in% data_cpdb_filter$id_cp_interaction),]
#mean_filter <- mean_filter[,c(1:13,17)]

library(patchwork)
ligands=unique(data_cpdb_filter$gene_a)
ligands
DotPlot(seu_obj, ligands) + RotatedAxis() + plot_annotation(title = "Ligands") 

ggsave("cellphoneDBres_ligand_Sig.tiff",  width = 9, height = 7)

receptors=unique(data_cpdb_filter$gene_b)
receptors
DotPlot(seu_obj, receptors) + RotatedAxis() + plot_annotation("Receptors")
ggsave("cellphoneDBres_receptor_Sig.tiff",  width = 9, height = 7)

plot_cpdb(
    scdata=seu_obj,
    cell_type1=".",
    cell_type2=".",
    celltype_key="cat",
    means=mean_filter,
    pvals=data_cpdb_filter,
    degs_analysis = TRUE,
    keep_significant_only = TRUE,
    genes = ligands)+
  theme(legend.position="bottom",
        axis.text.x = element_text(angle = 90, vjust = 1))
  

ggsave("cellphoneDBres_l2fc0.5_filterSig.tiff",  width = 8, height = 10)