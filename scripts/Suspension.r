# Scripts and functions:
  func.dir="/home/vincent/Documents/Rlib/functions/"
  source("/home/vincent/Documents/Rlib/functions/load_functions.r")
   
 
######################################################
### libraries
library(Seurat)
library(anndata)
library(matrixStats)




########################################################
## preprocessing - we are going to python in the first instance,
## to join the data today, and do embeddings for quick annotations


# Lazar2025
	x = readRDS("data/suspension/Lazar2025/HDCA_heart_sc_full_extended.rds")
	write.table(rownames(x), file="data/suspension/Lazar2025/features_Lazar2025.csv", sep=",", col.names=F, quote=F, row.names=F)
	y = anndata::AnnData(X = t(GetAssayData(x, layer="counts")))
	y$obs = x@meta.data
	anndata::write_h5ad(y, file="data/suspension/Lazar2025/Lazar2025.h5ad")


# Farah2024
	x = readRDS("data/suspension/Farah2024/overall_heart.rds")
	write.table(rownames(x), file="data/suspension/Farah2024/features_Farah2024.csv", sep=",", col.names=F, quote=F, row.names=F)
	y = anndata::AnnData(X = t(GetAssayData(x, layer="counts")))
	y$obs = x@meta.data
	anndata::write_h5ad(y, file="data/suspension/Farah2024/Farah2024.h5ad")


# DeBono2025 - already in 10X files, we can read in python
	
# Wang2020 - as mixed CSV matrices
	x = read.csv("data/suspension/Wang2020/GSE109816_normal_heart_umi_matrix.csv.gz", row.names = 1)
	x.cluster.info = read.csv("data/suspension/Wang2020/GSE109816_normal_heart_cell_cluster_info.txt", sep = "\t")
	x.qc.info = read.csv("data/suspension/Wang2020/GSE109816_normal_heart_cell_info.txt", sep = "\t")
	x.barcode.info = read.csv("data/suspension/Wang2020/GSE109816_metadata_barcodes_9994cells.txt", sep = "\t")
        rownames(x.cluster.info) = x.cluster.info[,1]
	rownames(x.qc.info) = x.qc.info[,1]
	rownames(x.barcode.info) = x.barcode.info[,1]

	# matrix conversion from data.frame
	x = as.matrix(x)

	z = read.csv("data/suspension/Wang2020/GSE121893_human_heart_sc_umi.csv.gz", row.names = 1)
	z.cluster.info = read.csv("data/suspension/Wang2020/GSE121893_all_heart_cell_cluster_info.txt.gz", sep = "\t")
	z.qc.info = read.csv("data/suspension/Wang2020/GSE121893_human_heart_sc_info.txt", sep = "\t")
	rownames(z.cluster.info) = z.cluster.info[,1]
	rownames(z.qc.info) = z.qc.info[,1]
	
	# matrix conversion from data.frame
	z = as.matrix(z)
	
	# METADATA
	# it looks like z.cluster.info is the correct metadata file - 
	# although, it IS missing 200 - 300 cells comparable to values given in publication
	# basically, the metadata is shit - missing cells, unknown barcodes, too large matrices, nothing is right. 
	# We'll cut using the z.cluster.info column, ID, then reapply QC filters as in published materials

	all.metadata = z.cluster.info
	x2 = x[,colnames(x) %in% all.metadata[,1]]
	z2 = z[,colnames(z) %in% all.metadata[,1]]
	
	disease_condition = gsub("^N.*", "healthy", all.metadata$condition)
	disease_condition[disease_condition != "healthy"] = "disease"
	
	all.metadata[,"disease_condition"] = disease_condition
	
	#### yeeesh now the variables are different... WHHHHYYY
	x2 = x2[rownames(z2),]
	
	# collect matrix columns
	y_matrix = cbind(x2, z2)

	# re-order metadata
	rownames(all.metadata) = all.metadata[,1]
	all.metadata = all.metadata[colnames(y_matrix),]
	
	# add QC metadata to metadata
	x_barcodes = colnames(x2)
	z_barcodes = colnames(z2)
	
	# mito
	all.metadata$mito = 0
	all.metadata[x_barcodes,"mito"] = x.qc.info[x_barcodes,"mito.perc"]
	all.metadata[z_barcodes,"mito"] = z.qc.info[z_barcodes,"mito.perc"]
	
	# Alignment
	all.metadata$alignment = 0
	all.metadata$total_alignment = 0
	all.metadata[x_barcodes,"alignment"] = x.qc.info[x_barcodes,"Aligned"]
	all.metadata[z_barcodes,"alignment"] = z.qc.info[z_barcodes,"Aligned"]
	
	all.metadata[x_barcodes,"total_alignment"] = x.qc.info[x_barcodes,"Total.Barcode.Alignments"]
	all.metadata[z_barcodes,"total_alignment"] = z.qc.info[z_barcodes,"Total.Barcode.Alignments"]
	
	
	# check QC parameters
	### Quality control! - I just reproduced what they did in their methods
	gene.filter = all.metadata$nGene > 500
	UMI.filter = all.metadata$nUMI > 10^(mean(log10(all.metadata$nUMI)) -
	2 * sd(log10(all.metadata$nUMI))) & all.metadata$nUMI < 10^(mean(log10(all.metadata$nUMI)) +
	2 * sd(log10(all.metadata$nUMI)))
	alignment.filter = (all.metadata$alignment/all.metadata$total_alignment) > 0.5
	mt.filter = mt.fraction < 0.72
	CM_digest.filter = all.metadata$nUMI > 10000
	CM_digest.filter[!all.metadata$group %in% c("CM")] = TRUE
	total.filter = gene.filter & UMI.filter & alignment.filter & mt.filter & CM_digest.filter
	y_matrix.filtered = y_matrix[,total.filter]

	x = y_matrix
	# EXPORT
	write.table(rownames(x), file="data/suspension/Wang2020/features_Wang2020.csv", sep=",", col.names=F, quote=F, row.names=F)
	y = anndata::AnnData(X = t(x))
	y$obs = all.metadata
	anndata::write_h5ad(y, file="data/suspension/Wang2020/features_Wang2020.h5ad")

	
# KnightSchrijver2022 - already as h5ad from czi


# Kuppe2022 - already as h5ad from czi



############### We need to process all datasets in python ################

  
  
