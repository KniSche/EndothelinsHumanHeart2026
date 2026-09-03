# Scripts and functions:
  func.dir="/home/vincent/Documents/Rlib/functions/"
  source("/home/vincent/Documents/Rlib/functions/load_functions.r")
   
 
######################################################
### libraries
library(Seurat)
library(anndata)
library(matrixStats)
library(DESeq2)
library(parallel)
library(BiocParallel)


  

# # Load all data
# adata = anndata::read_h5ad("data/suspension/anndata_combined.h5ad")
# cell_names = rownames(adata$X)

# Load samples from pseudobulk samples file:
bulk_samples = readRDS("data/suspension/bulk_samples_20.rds")
bulk_metadata = readRDS("data/suspension/bulk_metadata_20.rds")
bulk_size.factors = readRDS("bulk_size.factors_20.rds")
bulk_data = readRDS("bulk_data_20FULL.rds")


# we have 20 iterations... each has ~ 1400 psuedobulk samples
# We loop iteration, then sample


bulk_data = mclapply(
	1:20,
	function(i){
		bulk_data = sapply(
			1:length(bulk_samples[[i]]),
			function(condition){
				iteration_condition_samples = bulk_samples[[i]][[condition]]
				iteration_condition_bulkdata = colSums(adata[iteration_condition_samples,]$X)
				return(iteration_condition_bulkdata)
			}
		)
		return(bulk_data)
	}, mc.cores=5
)
for(i in 1:length(bulk_data)){
	colnames(bulk_data[[i]]) = rownames(bulk_metadata[[i]])

}
gc()
saveRDS(bulk_data, "bulk_data_20FULL.rds")


# now we have the bulk dataset... we can estimate the variables we need
# we can try to calculate this in parallel

	
	precalc_columns = c("baseMean",     "baseVar",      "allZero",      "dispGeneEst",  "dispGeneIter",
        			"dispFit",      "dispersion",   "dispIter",     "dispOutlier",  "dispMAP")
	
	celltypes2test = c("vCM", "aCM", "SMC", "Pericytes", "EC_Valves_Pericytes", "EC_Venous", "EC_Arterial", "EC_Capillary", "EC_Lymphatic",
				"EC_Mixed", "Fibroblasts", "Neurons_SchwannGlialCells", "Leukocytes")
				
	genes_EDN = c("EDN1", "EDN2", "EDN3", "EDNRA", "EDNRB", "ECE1", "ECE2")
	
	global_dispersions = list()
	global_results = list()
	
	for(i in 1:20){
			
			global_results_i = mclapply(1:length(celltypes2test),function(j){
				celltype = celltypes2test[j]
			
			# First create the DESeq object, adult and celltype only:

			# subset to adult only:
			adult_mask = bulk_metadata[[i]]$stage == "adult"

			# subset to celltype only:
			celltype_mask = bulk_metadata[[i]]$celltype == celltype

			# remove pseudobulk samples with fewer than 10 cells
			cell_number_mask = bulk_metadata[[i]]$cell_number >= 10

			# total mask
			DESeq_mask = adult_mask & celltype_mask & cell_number_mask

			DESeq_metadata = bulk_metadata[[i]][DESeq_mask,]
			DESeq_data = bulk_data[[i]][,DESeq_mask]
			DESeq_size.factors = bulk_size.factors[[i]][DESeq_mask]
				 DESeq_size.factors = DESeq_size.factors / exp(mean(log(DESeq_size.factors)))
				


			
			dds_full = DESeqDataSetFromMatrix(
				countData = round(DESeq_data),
				colData = DESeq_metadata,
				design = ~ cell_number_log + cell_or_nuclei + dataset + disease_subcategory
			)
			sizeFactors(dds_full) <- DESeq_size.factors
				
			# re level all to healthy
			dds_full$disease_subcategory <- relevel(dds_full$disease_subcategory, "healthy")


			dds_full = DESeq(dds_full, parallel=F)
			
			
			# moving to results
			dds_lrt <- nbinomLRT(dds_full,
			full = ~ cell_number_log + cell_or_nuclei + dataset + disease_subcategory,
			reduced = ~ cell_number_log + cell_or_nuclei + dataset
			)
	
			
		        # results from the LRT
		            res <- results(dds_lrt)
		            LRT_pvalue = res$padj 

                        # Get the unique conditions
                            conditions <- as.character(unique(dds_lrt$disease_subcategory))
                        
                        # Generate all unique combinations (pairs)
                        # This creates a matrix where each column is a pair
                            condition_pairs <- combn(conditions, 2)
                        
                        # Create a list to store the results
                            all_results <- list()
            
                            # Comparisons Loop through the pairs and run 'results()'
                                for(k in 1:ncol(condition_pairs)) {
                                    target <- condition_pairs[1, k]
                                    ref <- condition_pairs[2, k]
                                    
                                    res_name <- paste0(target, "_vs_", ref)
                                    
                                    # Extract results for this specific contrast
                                    all_results[[res_name]] <- results(
                                        dds_lrt,
                                        contrast = c("disease_subcategory", target, ref),
                                        test="Wald"
                                    )
                                }
                        
                        # add LRT to results
                            for(k in 1:length(all_results)){
                                all_results[[k]][,"LRT_padj"] = LRT_pvalue
                                all_results[[k]] = all_results[[k]][genes_EDN,]
                            }

			return(
				list(
					"dispersions" = mcols(dds_full)[,precalc_columns], 
					"results" = all_results
				)
			)
		}, mc.cores=13
		)
		
		
		
		global_results[[i]] = lapply(1:length(global_results_i), function(j){global_results_i[[j]][["results"]]})
		names(global_results[[i]]) = celltypes2test

		global_dispersions[[i]] = lapply(1:length(global_results_i), function(j){global_results_i[[j]][["dispersions"]]})
		names(global_dispersions[[i]]) = celltypes2test
		
		gc()
	}


	
	for(i in 1:length(global_results)){
		for(j in 1:length(global_results[[i]])){
			for(k in 1:length(global_results[[i]][[j]])){
				global_results[[i]][[j]][[k]] = global_results[[i]][[j]][[k]][genes_EDN,]
			}
		}
	}
	

	saveRDS(global_dispersions, "global_dispersions.rds")
	saveRDS(global_results, "global_results.rds")


# What about the Dataset Subsets?
# Kuppe2022:
	precalc_columns = c("baseMean",     "baseVar",      "allZero",      "dispGeneEst",  "dispGeneIter",
        			"dispFit",      "dispersion",   "dispIter",     "dispOutlier",  "dispMAP")
	
	celltypes2test = c("vCM", "aCM", "SMC", "Pericytes", "EC_Valves_Pericytes", "EC_Venous", "EC_Arterial", "EC_Capillary", "EC_Lymphatic",
				"EC_Mixed", "Fibroblasts", "Neurons_SchwannGlialCells", "Leukocytes")
				
	genes_EDN = c("EDN1", "EDN2", "EDN3", "EDNRA", "EDNRB", "ECE1", "ECE2")
	
	Kuppe2022_dispersions = list()
	Kuppe2022_results = list()
	

	for(i in 1:20){
			print(paste("--------------- ", i, " / 20: ", "", " -------------", sep=""))
			
			Kuppe2022_results_i = mclapply(1:length(celltypes2test),function(j){
				celltype = celltypes2test[j]
			
			# First create the DESeq object, adult and celltype only:

			# subset to adult only:
			adult_mask = bulk_metadata[[i]]$stage == "adult"

			# subset to celltype only:
			celltype_mask = bulk_metadata[[i]]$celltype == celltype

			# remove pseudobulk samples with fewer than 10 cells
			cell_number_mask = bulk_metadata[[i]]$cell_number >= 10

			# dataset only
			dataset_mask = bulk_metadata[[i]]$dataset == "Kuppe2022"
			#
			# total mask
			DESeq_mask = adult_mask & celltype_mask & cell_number_mask & dataset_mask 

			DESeq_metadata = bulk_metadata[[i]][DESeq_mask,]
			DESeq_data = bulk_data[[i]][,DESeq_mask]
			DESeq_size.factors = bulk_size.factors[[i]][DESeq_mask]
				 DESeq_size.factors = DESeq_size.factors / exp(mean(log(DESeq_size.factors)))
				


			
			dds_full = DESeqDataSetFromMatrix(
				countData = round(DESeq_data),
				colData = DESeq_metadata,
				design = ~ cell_number_log + disease_subcategory
			)
			sizeFactors(dds_full) <- DESeq_size.factors
				
			# re level all to healthy
			dds_full$disease_subcategory <- relevel(dds_full$disease_subcategory, "healthy")
			dds_full = DESeq(dds_full, parallel=F)
			
			# moving to results
			dds_lrt <- nbinomLRT(dds_full,
			full = ~ cell_number_log + disease_subcategory,
			reduced = ~ cell_number_log
			)

		        # results from the LRT
		            res <- results(dds_lrt)
		            LRT_pvalue = res$padj 

                        # Get the unique conditions
                            conditions <- as.character(unique(dds_lrt$disease_subcategory))
                        
                        # Generate all unique combinations (pairs)
                        # This creates a matrix where each column is a pair
                            condition_pairs <- combn(conditions, 2)
                        
                        # Create a list to store the results
                            all_results <- list()
            
                            # Comparisons Loop through the pairs and run 'results()'
                                for(k in 1:ncol(condition_pairs)) {
                                    target <- condition_pairs[1, k]
                                    ref <- condition_pairs[2, k]
                                    
                                    res_name <- paste0(target, "_vs_", ref)
                                    
                                    # Extract results for this specific contrast
                                    all_results[[res_name]] <- results(
                                        dds_lrt,
                                        contrast = c("disease_subcategory", target, ref),
                                        test="Wald"
                                    )
                                }
                        
                        # add LRT to results
                            for(k in 1:length(all_results)){
                                all_results[[k]][,"LRT_padj"] = LRT_pvalue
                                all_results[[k]] = all_results[[k]][genes_EDN,]
                            }

			
			return(
				list(
					"dispersions" = mcols(dds_full)[,precalc_columns], 
					"results" = all_results
				)
			)
		}, mc.cores=13
		)
		
		Kuppe2022_results[[i]] = lapply(1:length(Kuppe2022_results_i), function(j){Kuppe2022_results_i[[j]][["results"]]})
		names(Kuppe2022_results[[i]]) = celltypes2test
		
		Kuppe2022_dispersions[[i]] = lapply(1:length(Kuppe2022_results_i), function(j){Kuppe2022_results_i[[j]][["dispersions"]]})
		names(Kuppe2022_dispersions[[i]]) = celltypes2test
		
		gc()
	}


	saveRDS(Kuppe2022_dispersions, "Kuppe2022_dispersions.rds")
	saveRDS(Kuppe2022_results, "Kuppe2022_results.rds")





##### Wang2020
	precalc_columns = c("baseMean",     "baseVar",      "allZero",      "dispGeneEst",  "dispGeneIter",
        			"dispFit",      "dispersion",   "dispIter",     "dispOutlier",  "dispMAP")
	
	celltypes2test = c("vCM", "aCM", "SMC", "Pericytes", "EC_Valves_Pericytes", "EC_Venous", "EC_Arterial", "EC_Capillary", "EC_Lymphatic",
				"EC_Mixed", "Fibroblasts", "Neurons_SchwannGlialCells", "Leukocytes")
	
	wang_ignore_celltypes =c("Pericytes", "EC_Lymphatic", "Neurons_SchwannGlialCells")
				
	genes_EDN = c("EDN1", "EDN2", "EDN3", "EDNRA", "EDNRB", "ECE1", "ECE2")
	
	Wang2020_dispersions = list()
	Wang2020_results = list()
	
	
	for(i in 1:20){
			print(paste("--------------- ", i, " / 20: ", "", " -------------", sep=""))
			
			Wang2020_results_i = mclapply(1:length(celltypes2test),function(j){
				celltype = celltypes2test[j]
				if(!celltype %in% wang_ignore_celltypes){
				# First create the DESeq object, adult and celltype only:

				# subset to adult only:
				adult_mask = bulk_metadata[[i]]$stage == "adult"

				# subset to celltype only:
				celltype_mask = bulk_metadata[[i]]$celltype == celltype

				# remove pseudobulk samples with fewer than 10 cells
				cell_number_mask = bulk_metadata[[i]]$cell_number >= 10

				# dataset only
				dataset_mask = bulk_metadata[[i]]$dataset == "Wang2020"
				#
				# total mask
				DESeq_mask = adult_mask & celltype_mask & cell_number_mask & dataset_mask 

				DESeq_metadata = bulk_metadata[[i]][DESeq_mask,]
				DESeq_data = bulk_data[[i]][,DESeq_mask]
				DESeq_size.factors = bulk_size.factors[[i]][DESeq_mask]
					 DESeq_size.factors = DESeq_size.factors / exp(mean(log(DESeq_size.factors)))
					


				
				dds_full = DESeqDataSetFromMatrix(
					countData = round(DESeq_data),
					colData = DESeq_metadata,
					design = ~ cell_number_log + disease_subcategory
				)
				sizeFactors(dds_full) <- DESeq_size.factors
					
				# re level all to healthy
				dds_full$disease_subcategory <- relevel(dds_full$disease_subcategory, "healthy")
				dds_full = DESeq(dds_full, parallel=F)
				
				# moving to results
				dds_lrt <- nbinomLRT(dds_full,
				full = ~ cell_number_log + disease_subcategory,
				reduced = ~ cell_number_log
				)

				# results from the LRT
				    res <- results(dds_lrt)
				    LRT_pvalue = res$padj 

		                # Get the unique conditions
		                    conditions <- as.character(unique(dds_lrt$disease_subcategory))
		                
		                # Generate all unique combinations (pairs)
		                # This creates a matrix where each column is a pair
		                    condition_pairs <- combn(conditions, 2)
		                
		                # Create a list to store the results
		                    all_results <- list()
		    
		                    # Comparisons Loop through the pairs and run 'results()'
		                        for(k in 1:ncol(condition_pairs)) {
		                            target <- condition_pairs[1, k]
		                            ref <- condition_pairs[2, k]
		                            
		                            res_name <- paste0(target, "_vs_", ref)
		                            
		                            # Extract results for this specific contrast
		                            all_results[[res_name]] <- results(
		                                dds_lrt,
		                                contrast = c("disease_subcategory", target, ref),
		                                test="Wald"
		                            )
		                        }
		                
		                # add LRT to results
		                    for(k in 1:length(all_results)){
		                        all_results[[k]][,"LRT_padj"] = LRT_pvalue
		                        all_results[[k]] = all_results[[k]][genes_EDN,]
		                    }
		                
				
				return(
					list(
						"dispersions" = mcols(dds_full)[,precalc_columns], 
						"results" = all_results
					)
				)
			} else {
				return(list("results" = NULL, "dispersions" = NULL))
			}
		}, mc.cores=10
		)
		
		Wang2020_results[[i]] = lapply(1:length(Wang2020_results_i), function(j){Wang2020_results_i[[j]][["results"]]})
		names(Wang2020_results[[i]]) = celltypes2test
		
		Wang2020_dispersions[[i]] = lapply(1:length(Wang2020_results_i), function(j){Wang2020_results_i[[j]][["dispersions"]]})
		names(Wang2020_dispersions[[i]]) = celltypes2test
		
		gc()
	}


	saveRDS(Wang2020_dispersions, "Wang2020_dispersions.rds")
	saveRDS(Wang2020_results, "Wang2020_results.rds")


