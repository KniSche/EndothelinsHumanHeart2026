### Supplementary code for Endothelins in the human heart 2026

1. figures/
  - Figure panels as raw input to main figures, or supplementary figures produced by the scripts.
2. scripts/
  - scripts (exported `HTML`, `.pdf`, of R or Python script files associated with stages of the analysis).
3. source_data/ 
  - data generated in the scripts that tie with various figures in the article.
4. functions/
  - Some of the functions used throughout the analysis - not comprehensive as some private libraries were used for basic functions [https://github.com/KniSche/R_functions]
  
Files:
- `bulk_*_20.rds`: R data store objects for bulk dataset, metadata, sampleIDs, and size_factors
- `DESeq_results_*_20.rds`: Results from running DESeq2 on endothelin subset only over 20 iterations. 
- `*global_*`: global model (all datasets, including technical covariates)
- `*Kuppe2022_*`: Kuppe2022 dataset model (MI vs healthy)
- `*Wang2020_*`: Wang2020 dataset model (heart failure vs healthy)
- `*_results.rds`: DESeq2 results calculated from all genes, then subset by endothelin genes
- `*_dispersions.rds`: DESeq2 dispersions calculated from all genes, injected into DESeq2 when running on a low-performance computer.


Note: dispersions objects were too big to include on GitHub. However, `DESeq_results_*_20.rds` objects should contain some information and include the global results (see script). Please get in touch if the pre-calculated dispersions are desired.
