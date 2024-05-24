library(Seurat)
library(anndata)
meta_all <- readRDS("D1-2-3-24h_qc_pca_umap_clust_harmony_annotated.rds")
exprs_norm <- readRDS("exprs_norm.rds")
meta_all <- meta_all[which(meta_all$cell %in% colnames(exprs_norm)),]

exprs_norm <- exprs_norm[,meta_all$cell]

meta_all$cil <- ifelse(meta_all$cell_type == "Ciliated", "ciliated", "others")
meta_all$infection_status <- ifelse(meta_all$sample_id == "no_virus", "nonInfected", "Infected")

meta_all$cat <- paste0(paste0(meta_all$cil, meta_all$infection_status))

#write.table(meta_all[,c("cell", "cat")], "05022024_metaSubset_cilInfected_OthersnonInf_cpdb_anndataVersion.txt", sep="\t", row.names = F, quote = F)

seu_obj <- CreateSeuratObject(exprs_norm, project = "scRNA_24h", assay = "RNA")
seu_obj <- AddMetaData(seu_obj, meta_all)
seu_obj <- NormalizeData(seu_obj)
Idents(seu_obj) <- seu_obj$cat

seu_obj_subset <- subset(seu_obj, idents = c("ciliatedInfected", "othersnonInfected"))
#write.table(seu_obj_subset@meta.data[,c("cell", "cat")], "05022024_metaSubset_cilInfected_OthersnonInf_cpdb_anndataVersion.txt", sep="\t", row.names = F, quote = F)

reticulate::use_condaenv("/opt/anaconda3/bin/python")

x_mat <- as.matrix(seu_obj_subset@assays$RNA$data)
t_x_mat <- t(x_mat)
ad <-  AnnData(X=t_x_mat)
write_h5ad(anndata=ad, filename ="anndata_seuObj_Subset_EpithelialRV_allDonor.h5ad")

markers <- FindMarkers(seu_obj, ident.1="ciliatedinfected", ident.2 = "ciliatednon-infected", min.pct=0.1)
write.table(markers, "markers_ciliatedInfected_vs_CiliatedNonInfected_04292024.txt", quote=F)