source("Utils_scRNA.R")

data_dir <- 'filtered_feature_bc_matrix/'
list.files(data_dir)
expression_matrix <- Read10X(data.dir = data_dir)
seurat_obj <- CreateSeuratObject(counts = expression_matrix, project="Epithelial_rv", names.delim = "-", names.field = 2)

head(seurat_obj)

metadata_scrna <- read.table("MetadataSC.txt", sep = "\t", header = TRUE)
cellranger_agr <- read.table("aggr_file2.csv", sep = ",", header = TRUE)
cellranger_agr$orig.ident <- paste(1:10)

seurat_meta <- as.data.frame(seurat_obj@meta.data)

seurat_meta_cellranger <- data.frame(row.names = rownames(seurat_meta),
                                     "sample_id" = ifelse(seurat_meta$orig.ident == 1, "D1_12h", 
                                                          ifelse(seurat_meta$orig.ident == 2, "D1_16h",
                                                                 ifelse(seurat_meta$orig.ident == 3, "D1_20h", 
                                                                        ifelse(seurat_meta$orig.ident == 4, "D1_24h",
                                                                               ifelse(seurat_meta$orig.ident == 5, "D1_2h",
                                                                                      ifelse(seurat_meta$orig.ident == 6, "D1_42h", 
                                                                                             ifelse(seurat_meta$orig.ident == 7, "D1_8h",
                                                                                                    ifelse(seurat_meta$orig.ident == 8, "D2_24h", 
                                                                                                           ifelse(seurat_meta$orig.ident == 9, "D3_24h", "no-virus"))))))))))


#create metadata

raw_exprs <- expression_matrix
dim(raw_exprs)

cells <- colnames(raw_exprs)
nUMI <- Matrix::colSums(raw_exprs)
nGene <- Matrix::colSums(raw_exprs > 0)
mito_genes <- grep("^MT-", rownames(raw_exprs), value = TRUE, ignore.case = TRUE)
percent_mito <- Matrix::colSums(raw_exprs[mito_genes, ]) / Matrix::colSums(raw_exprs) * 100

meta_all <- data.frame("cell" = cells, "nUMI" = nUMI, "nGene" = nGene, "percent_mito" = percent_mito, stringsAsFactors = F) 
head(meta_all)

genes_ribo <- grep("^RPS\\d+|^RPL\\d+", row.names(raw_exprs), value = TRUE)
length(genes_ribo)
meta_all$percent.ribo <- Matrix::colSums(raw_exprs[genes_ribo, meta_all$cell]) / 
                                        Matrix::colSums(raw_exprs[, meta_all$cell]) * 100
head(meta_all$percent.ribo)

meta_all[,"sample_id"] <- seurat_meta_cellranger[rownames(meta_all),"sample_id"]

#filter 24h time point only

dim(meta_all)
meta_filt <- meta_all[which(meta_all$sample_id %in% c("no-virus", "D1_24h", "D2_24h", "D3_24h")),]
dim(meta_filt)

meta_all <- meta_filt
mRNA_exprs <- raw_exprs[,meta_all$cell]
dim(mRNA_exprs)

 
#filter data
#keep cells with nGene > 500 and percent mitochondrial lower than 20%

dim(meta_all)
meta_filt <- meta_all[-which(meta_all$nGene <= 500 | meta_all$percent_mito >= 20),]
dim(meta_filt)

meta_all <- meta_filt
mRNA_exprs <- raw_exprs[,meta_all$cell]
dim(mRNA_exprs)

exprs_norm <- mRNA_exprs %>% NormalizeDataSeurat()
dim(exprs_norm)

meta_all$sample <- rep("AS001", nrow(meta_all))
genes_exclude <- grep("^MT-|^RPL|^RPS|MALAT1|MIR-|^IG", rownames(exprs_norm), value=TRUE)
vargenes_df <- FindVariableGenesBatch(exprs_norm, meta_all, genes_exclude, 1000)
var_genes <- vargenes_df$gene

exprs_scaled_sigl <- exprs_norm[var_genes, ] %>% ScaleDataSeurat()
exprs_cosine_sigl <- cosineNorm(exprs_scaled_sigl,mode=c("matrix"))

 
#umap and clustering
pca_res <- irlba::prcomp_irlba(t(exprs_cosine_sigl), 20)
umap_res <- umap(pca_res$x, n_neighbors=30,metric="cosine", min_dist=.3,preserve.seed = TRUE)
meta_all$UMAP1 <- umap_res$layout[,1]
meta_all$UMAP2 <- umap_res$layout[,2]
dim(meta_all)
head(meta_all)
tempMeta <- cbind(meta_all, pca_res$x)
tempMeta[1:4,]
meta_all <- tempMeta
colnames(tempMeta)
colnames(meta_all)
#saveRDS(meta_all, "~/Documents/AADCRC/EpithelialRV/allDonor/results/D1-2-3-24h_qc_pca_umap.rds")

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)

ids_cos <- Reduce(cbind, lapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, modularity = 1,
                                   resolution = res_use, algorithm = 3, n.start = 10,
                                   n.iter = 10, random.seed = 0, print.output = FALSE,
                                   temp.file.location = NULL, edge.file.name = NULL)
}))

ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)
dim(ids_cos)
ids_cos[1:5, 1:5]

meta_all$res_0.20 <- ids_cos$res_0.20
meta_all$res_0.40 <- ids_cos$res_0.40
meta_all$res_0.60 <- ids_cos$res_0.60
meta_all$res_0.80 <- ids_cos$res_0.80
meta_all$res_1.00 <- ids_cos$res_1.00

length(unique(meta_all$res_0.20))
length(unique(meta_all$res_0.40))
length(unique(meta_all$res_0.60))
length(unique(meta_all$res_0.80))
length(unique(meta_all$res_1.00))

## #saveRDS(meta_all, "~/Documents/AADCRC/EpithelialRV/allDonor/results/D1-2-3-24h_qc_pca_umap_clust.rds")

harmony <- HarmonyMatrix(pca_res$x, meta_all, do_pca = FALSE,
                         c("donor","condition"),
                         epsilon.cluster = -Inf,
                         epsilon.harmony = -Inf,
                         max.iter.cluster = 50,
                         max.iter.harmony = 30,
                         plot_convergence = T)

colnames(harmony) <- paste0("harmonized_", colnames(harmony), sep="")
meta_all <- cbind(meta_all, harmony)
colnames(meta_all)
#saveRDS(meta_all, "~/Documents/AADCRC/EpithelialRV/allDonor/results/D1-2-3-24h_qc_harmony.rds")

#UMAP with 20 PCs
umap_res <- umap(harmony[, c(1:20)], n_neighbors = 30, metric = "cosine", min_dist = .1)
meta_all$harmonized_UMAP1 = umap_res$layout[,1]
meta_all$harmonized_UMAP2 = umap_res$layout[,2]
colnames(meta_all)

#Louvain calculation
snn_pcs <- BuildSNNSeurat(meta_all[, c(37:56)], nn.eps = .5)

resolution_list <- c(0.2, 0.4, 0.6, 0.8, 1.0)

ids_cos <- Reduce(cbind, mclapply(resolution_list, function(res_use) {
  Seurat:::RunModularityClustering(SNN = snn_pcs, modularity = 1,
                                   resolution = res_use, algorithm = 3, n.start = 10,
                                   n.iter = 10, random.seed = 0, print.output = FALSE,
                                   temp.file.location = NULL, edge.file.name = NULL)
}, mc.cores = min(16, length(resolution_list))))


ids_cos %<>% data.frame()
colnames(ids_cos) <- sprintf("res_%.2f", resolution_list)
dim(ids_cos)
ids_cos[1:5, 1:5]

meta_all$res_0.20 <- ids_cos$res_0.20
meta_all$res_0.40 <- ids_cos$res_0.40
meta_all$res_0.60 <- ids_cos$res_0.60
meta_all$res_0.80 <- ids_cos$res_0.80
meta_all$res_1.00 <- ids_cos$res_1.00

length(unique(meta_all$res_0.20))
length(unique(meta_all$res_0.40))
length(unique(meta_all$res_0.60))
length(unique(meta_all$res_0.80))
length(unique(meta_all$res_1.00))

#saveRDS(meta_all, "~/Documents/AADCRC/EpithelialRV/allDonor/results/D1-2-3-24h_qc_pca_umap_clust_harmony.rds")

geneListMarkers <- c("KRT5", "TP63", "KRT14", "KRT15", "S100A2", "BCAM", "DAPL1")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=3)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=3)

##Basal proliferating
 
geneListMarkers <- c("MKI67", "TOP2A", "CDK1", "HES1")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=3)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=3)

## Club markers

geneListMarkers <- c("SCGB1A1","KRT15","LYPD2")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=3)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=3)

 ## Deuterosomal markers

geneListMarkers <- c("DEUP1", "FOXJ1", "FOXN4", "CDC20B")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=3)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=3)

## Ciliated markers

geneListMarkers <- c("FOXJ1", "PIFO", "TPPP3", "SNTN", "FAM183A", "LRRIQ1", "DNAH12", "SNTN", "CAPS", "TUBB4B", "DNAH5", "TSPAN1")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=3)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=3)

## Secretory
 
geneListMarkers <- c("SCGB1A1", "MUC5AC", "MUC5B", "TFF3", "BPIFB1", "MSMB", "SLPI", "WFDC2")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=3)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=3)

## Neuroendocrine

geneListMarkers <- c("CHGA", "ASCL1","INSM1", "HOXB5")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=2)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=2)

## Tuft
geneListMarkers <- c("POU2F3", "AVIL", "GNAT3","TRPM5")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=3)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=3)

## Ionocyte

geneListMarkers <- c("FOXL1", "CFTR", "ASCL3", "INSM1", "HOXB5")

output1 <- plot_UMAP_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output1, ncol=3)

output2 <- plot_UMAP_harmonized_markers(geneListMarkers, exprs_norm, meta_all)

grid.arrange(grobs=output2, ncol=3)

#Cluster assignment
meta_all <- meta_all %>%
  mutate(cell_type = case_when(res_0.20 == 0 ~ "Transitional",
         res_0.20 == 1 ~ "Ciliated", 
         res_0.20 == 2 ~ "Secretory", 
         res_0.20 == 3 ~ "Basal", 
         res_0.20 == 4 ~ "Deuterosomal", 
         res_0.20 == 5 ~ "Deuterosomal", 
         TRUE ~ "Rare"))

colors <- c("#FC1CBF","#BDCDFF","#B5EFB5","#DEA0FD","#7ED7D1","#1C8356")

#saveRDS(meta_all, "~/Documents/AADCRC/EpithelialRV/allDonor/results/D1-2-3-24h_qc_pca_umap_clust_harmony_annotated.rds")

#create h5ad file for scDRS
seu_obj <- CreateSeuratObject(mRNA_exprs, project = "scRNA_24h", assay = "RNA")
seu_obj <- AddMetaData(seu_obj, meta_all)

SaveH5Seurat(seu_obj, filename = "~/Documents/AADCRC/EpithelialRV/allDonor/results/seu_obj.h5Seurat", overwrite = T,verbose = T)
Convert("~/Documents/AADCRC/EpithelialRV/allDonor/results/seu_obj.h5Seurat", assay ="RNA", dest = "h5ad",verbose = T,overwrite = T)

