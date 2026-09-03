# Then we can feed this into a function...
  #  CellMatrix::
    plotspatialGene = function(
      image_h5ad,
      gene = NULL,
      rotation = 0,
      reflection = "none",
      opacity_filter = 0.9,
      sample_id,
      evaluate_image = T,
      feature_names = "default",
      nCounts_column = "n_counts",
      nGenes_column = "n_genes",
      QC_threshold = 1000
    ){
      # Image dataset (anndata object)  list.files(image_h5ads)

	if(class(image_h5ad)[1] == "character"){
	cat("\n 1. Loading Spatial Dataset: ", image_h5ad)
        	x <- read_h5ad(image_h5ad)
        } else {
        cat("\n 1. Loading Spatial Dataset: ", unique(image_h5ad$obs$sample))
		x = image_h5ad
        }
        sample_id_levels = levels(x$obs$sample)



cat("\n 2. Available sample IDs: ", sample_id_levels) 
  cat("\n ... processing sample: ", sample_id) 
         
    # Process sample, pull out the data
      image_adata <- x[x$obs$sample == sample_id,]

    # Point filtering (depth, standard QC for visium SD?)
    # This part could be changed depending on QC metrics and 
    # current standard practices
      spot_depths = image_adata$obs[, nCounts_column]
      spot_complexity = image_adata$obs[, nGenes_column]
      good_spots = which(spot_complexity > QC_threshold)
    
    if(evaluate_image == T){  
    cat("\n 3. \"evaluate_image\" set to TRUE:") 
    cat("\n ... Skipping spot matrix processing and gene expression writing")
    }
    
    if(evaluate_image != T){      
    cat("\n 3. \"evaluate_image\" set to FALSE:") 
    cat("\n ... Transforming spot matrix") 

    # filter to retain good spots only
      spot_matrix = image_adata[good_spots,]
      spot_matrix = spot_matrix$X
      gc()

    # Matrix  
      spot_matrix = as.matrix(spot_matrix)
      gc()
    # depth correction        
      spot_matrix = spot_matrix/rowSums(spot_matrix)*1e4
      gc()
    # log2 transformation, including a pseudocount of 1
      spot_matrix = log2(spot_matrix+1)
      gc()
    # round to 3 D.P. - this is a pragmatic space and time saving measure
      spot_matrix = round(spot_matrix,3)
      gc(full=T)
    
    # One of the challenges here is to keep the consistency between Matrix and Spatial 
    # We need to send the genes to their alphabetised positions, just like in the CellMatrix Code
      # Change column names, based on the adata SYMBOL, otherwise use supplied vector in argument
        if(feature_names == "default"){
          colnames(spot_matrix) = image_adata$var$feature_names
        } else {
          colnames(spot_matrix) = feature_names
        }
        
cat("\n ... Picking up the gene lookup table from: \"data/gene_lookup_table.csv\"" )
      # pull in the look up table as in CellMatrix suspension data
        gene_lookup_table = read.csv("data/gene_lookup_table.csv", head=F)
        gene_lookup_table = gene_lookup_table[gene_lookup_table[,1] %in% colnames(spot_matrix),]
      
      
      # Create gene lookup subset
      # (needed because gene_address will change between
      # the suspension, and spatial data based on features)
        gene_subset_spatial = gene_lookup_table
      # Check for existing file...
        if(!"gene_lookup_spatial.csv" %in% list.files("data/spatial/")){
          cat("\n ... Writing a spatial feature subset gene lookup table to: \"data/spatial/gene_lookup_spatial.csv\"")
    
          write.table(gene_subset_spatial, file="data/spatial/gene_lookup_spatial.csv", sep=",", quote=F, col.names=F, row.names=F)
        } else {
          cat("\n ... Picking up the spatial subset of gene lookup table from: \"data/spatial/gene_lookup_spatial.csv\"")
    
          gene_subset_spatial <- read.csv("data/spatial/gene_lookup_spatial.csv", head=F)
        }
          gene_groups = unique(gene_subset_spatial[,2])
      
cat("\n ... Writing the genes into their files: \n")
    
      # Write all gene groups to data/sample/  
        for(k in gene_groups){
          genes_k = gene_subset_spatial[gene_subset_spatial[,2] == k,1]
          b = spot_matrix[,genes_k]
          b = t(b)
          
          write.table(b,
            file=paste(
              "data/spatial/",
              sample_id,
              "/",
              "gene_matrices/",
              k,
              "_genes.csv", sep=""
            ),
            sep=",",
            quote=F,
            col.names=F,
            row.names=F
          )
cat(" ... Writing the genes into their files: ", k, "/", tail(gene_groups,1), sep="", fill=3)
    
        }
      } # END OF PROCESS GENE INFORMATION
    
    cat("\n 4. processing image(s).") 
    cat("\n ... pulling out image data from anndata file")
    
    # BEGIN process image information     
    # Parse adata to get the sample_id image
      img.idx <- grep(sample_id, names(image_adata$uns$spatial))

    # Pick the hi res version
      img_quality = "hires" # either hires, or lowres
      # Extract and process Image:
        image_data <- Image(
          image_adata$uns$spatial[[img.idx]][["images"]][[img_quality]]
        )
      # find the scale factor
        image_scale <- image_adata$uns$spatial[[img.idx]][["scalefactors"]][[paste("tissue_", img_quality, "_scalef", sep="")]]

      # Process Image:
        # Color mode
          colorMode(image_data) = "Color"
        
        # ROTATIONS ARE MANUAL!!!
          # Rotation...
          
          # native x
            x.native = image_adata$obsm[[grep("spatial" ,names(image_adata$obsm))]][good_spots,1] * image_scale
            y.native = image_adata$obsm[[grep("spatial" ,names(image_adata$obsm))]][good_spots,2] * image_scale
          
          # apply point rotation... for some reason the datasets have x and y mixed up?
          # I haven't looked into why, perhaps this is standard (also something to do with Fudicial alignment)
            x.coords <- y.native
            y.coords <- x.native
            
    cat("\n ... applying selected transformations")    
          # Orientation:
          # Rotation Clockwise
          # If rotation is not zero, put the image through the rotation functions below
            if(rotation != 0){
            
             # radians   
               radian_angle = (360-rotation) * (pi/180)
             
             # original origin (image centre)
               origin = dim(image_data)[1:2]/2
                
             # apply rotation using EBImage to the image
               image_data = rotate(image_data[,,1:3], rotation)
             
             # use Image centre as origin of data points
             # New origin (image centre)
               origin2 = dim(image_data)[1:2]/2
                      
              # rotation of datapoints around an origin:
              # difference to origin...
                x_r_a = x.coords - origin[1]
                y_r_a = y.coords - origin[2]
              
              # apply rotation after shifting to origin
                x_r =   x_r_a * cos(radian_angle) + y_r_a * sin(radian_angle)
                y_r = - x_r_a * sin(radian_angle) + y_r_a * cos(radian_angle)
              
              # add new origin 
                x_r = x_r + origin2[1] # origin[1]#*x.origin.scale_factor
                y_r = y_r + origin2[2] # origin[2]#*y.origin.scale_factor
              # update the coordinates      
                x.coords = x_r
                y.coords = y_r
            }
  
              
          # Reflection
            if(reflection == "vertical"){
              image_data = image_data[,ncol(image_data[,,1]):1,1:3]
              
              x.coords = x.coords
              y.coords = -y.coords+dim(image_data)[2]
            }
            if(reflection == "horizontal"){
              image_data = image_data[nrow(image_data[,,1]):1,,1:3]
              
              x.coords = -x.coords+dim(image_data)[1]
              y.coords = y.coords
            }
     
         
#         #### This commented section is just for checking the image orientation while constructing the function
##              pseudoGrey = (normalise((image_data[,,1]/3 + image_data[,,2]/3 + image_data[,,3]/3)))
##               #    display(channel(pseudoGrey*3, "grey"), method="raster", new=F)
##                  
##         
##            # Alter overlay image (tissue / opacity)
##              # construct a new image with opacity layer:
##                image_overlay = array(0,dim=c(dim(image_data)[1:2],4))
##                
##                # greyscale the image              
##                                  
##                image_overlay[,,1] = pseudoGrey*1
##                image_overlay[,,2] = pseudoGrey*1
##                image_overlay[,,3] = pseudoGrey*1
##                image_overlay[,,4]=1-(pseudoGrey*1)
##                
##                image_overlay = Image(image_overlay)
##                colorMode(image_overlay) = "Color"
##                
##                image_overlay[,,4][pseudoGrey > opacity_filter] = 1
##            
##          
##              display(channel(image_overlay, "grey")*light_or_dark, method="raster", new=F)
##      
##              # Construct the point underlay - this is purely for optimising the image during processing.
##              points(
##                  x.coords,
##                  y.coords,
##                  col="red",
##                  pch=16,
##                  cex=0.4#normalise(log(x+1)[plot.order])
##                ) # with colour
##           
##           
#           
    cat("\n 5. Generating grey scale images and mask")    
    ################ # plot 1
    ##################################################
            # We will work with the light palette first
            # Create a background Image (plotting under EBImage)
            #  image_background = Image(matrix(1, dim(image_data)[1], dim(image_data)[2]))
            #  display(image_background, method="raster", new=F)
            # Create the background image:
               pseudoGrey <- (normalise((image_data[,,1]/3 + image_data[,,2]/3 + image_data[,,3]/3)))
               #    display(channel(pseudoGrey*3, "grey"), method="raster", new=F)
                  
         
            # Alter overlay image (tissue / opacity)
              # construct a new image with opacity layer:
                image_overlay <- array(0,dim=c(dim(image_data)[1:2],4))
#                
                # greyscale the image                                             
                image_overlay[,,1] = pseudoGrey*1
                image_overlay[,,2] = pseudoGrey*1
                image_overlay[,,3] = pseudoGrey*1
                image_overlay[,,4]=1-(pseudoGrey*1)
                
                image_overlay = Image(image_overlay)
                colorMode(image_overlay) = "Color"
#                
                # had issues with this code during function calls, assignment of tmp objects               
#                image_overlay[,,4][pseudoGrey > opacity_filter] = 1
                # so I split it up:
                  tmp = matrix(as.vector(image_overlay[,,4]),dim(image_data)[1],dim(image_data)[2])
                  tmp[pseudoGrey > opacity_filter] = 1
                  image_overlay[,,4] = tmp
                  rm(tmp)
                  gc()
#            
          if(evaluate_image == T){
           cat("\n ... \"evaluate_image\" is TRUE, plotting images for evaluation")
        
           
            par(mfrow=c(1,2))
            # plot image
              display(channel(image_overlay, "grey")*1, method="raster", new=F)
    
              # Construct the point underlay - this is purely for optimising the image during processing.
              points(
                  x.coords,
                  y.coords,
                  col="red",
                  pch=16,
                  cex=0.4#normalise(log(x+1)[plot.order])
                ) # with colour
              
            # plot overlay
              display(image_overlay, method="raster", new=T)
        } 
    
    
        if(evaluate_image != T){   
          cat("\n 6. writing image 1 (light version): ", file.path("data/spatial", sample_id, "images", 1, "overlay.png"))
          
           
          # Overlay light
      
            writeImage(image_overlay, file.path("data/spatial", sample_id, "images", 1, "overlay.png"))
          # We can actually use the same image, but change the channel to grey scale
         }
       
       
       
    ################ # plot 2
    ############################################
        # DARK PLOT (black background)
          pseudoGrey <- 1-(normalise((image_data[,,1]/3 + image_data[,,2]/3 + image_data[,,3]/3)))
            
          
        # Alter overlay image (tissue / opacity)
          # construct a new image with opacity layer:
            image_overlay <- array(0,dim=c(dim(image_data)[1:2],4))
          
            # greyscale the image              
                              
              image_overlay[,,1] = pseudoGrey*1
              image_overlay[,,2] = pseudoGrey*1
              image_overlay[,,3] = pseudoGrey*1
              image_overlay[,,4]=1-(pseudoGrey*3)
              
              image_overlay = Image(image_overlay)
              colorMode(image_overlay) = "Color"
              
         # had issues with this code during function calls, assignment of tmp objects               
#                image_overlay[,,4][pseudoGrey > opacity_filter] = 1
                # so I split it up:
                  tmp = matrix(as.vector(image_overlay[,,4]),dim(image_data)[1],dim(image_data)[2])
                  tmp[pseudoGrey > opacity_filter] = 1
                  image_overlay[,,4] = tmp
                  rm(tmp)
                  gc()
                  
        # write images       
          if(evaluate_image == T){
          

            # plot image
              display(channel(image_overlay, "grey")*2, method="raster", new=F)
    
              # Construct the point underlay - this is purely for optimising the image during processing.
              points(
                  x.coords,
                  y.coords,
                  col="red",
                  pch=16,
                  cex=0.4#normalise(log(x+1)[plot.order])
                ) # with colour
              
            # plot overlay
              display(image_overlay, method="raster", new=T)
        }
        
        if(evaluate_image != T){  
        cat("\n ... writing image 2 (dark version): ", file.path("data/spatial", sample_id, "images", 2, "overlay.png") )
                
        
            writeImage(image_overlay, file.path("data/spatial", sample_id, "images", 2, "overlay.png"))
        
        cat("\n 7. writing datapoint coordinates to: ", file.path("data/spatial", sample_id, "coordinates.csv"))
        
         # We need to write point coordinates in the same way!
            write.table(cbind(x.coords, y.coords), file=file.path("data/spatial", sample_id, "coordinates.csv"), sep=",", quote=F, col.names=F, row.names=F)
      
         }
      
         
         gc(full=T)
         cat("\n")
         if(evaluate_image == T){
          cat("\n IMAGE EVALUATION: Please check the plots to ensure adequate overlap of images and spots, as well as any adjustments needed for image orientation. Once you are happy, re-run the function with the argument \"evaluate_image\" set to FALSE. This will then process the gene expression matrix, and write the data to the directory for the shiny app to handle.")
         } else {
          cat("\n IMAGE EVALUATION: no image evaluation was carried out. If transformations are needed, please re-run the function with the argument \"evaluate_image\" set to TRUE (default). This will skip the gene matrix processing and plot the images, re-run iteratively using the \"rotation\" and \"reflection\" arguments until the desired orientation is achieved.")
         }
          cat("\n")
          cat("\n")
    }   # END OF FUNCTION     
    
    
    
### process images?
#

#library(EBImage)
#library(anndata)

#func.dir="functions"
#source("functions/load_functions.r")

#  preprocessImage(
#    "../TEMP_VisiumStuff/processed_h5ads/Visium_4PCW_raw.h5ad",
#    rotation = 0,
#    reflection = "none",
#    opacity_filter = 0.9,
#    sample_id = "HCAHeartST13162332",
#    evaluate_image = T,
#    feature_names = "default"
#  )
  
#  preprocessImage(
#    image_h5ad = "../TEMP_VisiumStuff/processed_h5ads/Visium_4PCW_raw.h5ad"
#    rotation = 0
#    reflection = "none"
#    opacity_filter = 0.9
#    sample_id = "HCAHeartST13162332"
#    evaluate_image = F
#    feature_names = "default"
#  )
    
    
    
    

  # We can probably speed up the app by separating data collection from visualisation...
  # We will end this function, spitting out the variables of:
  # Output variables:
  
  # Plot 1: RGB
  VisiumData1 = function(
      genes_input,
      sample_id_list = sample_id_list
    ){
       
  
    # genes to collect
#      rgb.list = lapply(gene.list,
#        function(gene){
#          gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
#          gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
#          gene_address = paste("../CellMatrix/data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
#        
#          x = as.numeric(
#            scan(
#              gene_address, 
#              nlines=1, 
#              sep=",", 
#              what="character",
#              na.strings="NA",
#              skip = gene_position-1
#            )
#          )
#          return(x)
#        }
#      )
#      
      genes_found = genes_input
      genes_notFound = character()
      
      rgb.list_samples = lapply(sample_id_list,function(i){
        sample_id = i    
      # Reading in data as vectors
        gene.vector_all = 0
        if(length(genes_input) > 0){
            gene.vector_all = sapply(
              genes_input,
              function(gene){
                if(gene %in% gene_subset_spatial[,1]){
                  gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
                  gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
                  gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                  x = as.numeric(
                      scan(
                        gene_address, 
                        nlines=1, 
                        sep=",", 
                        what="character",
                        na.strings="NA",
                        skip = gene_position-1
                      )
                    )
                   } else {
                    genes_found <<- genes_found[!genes_found %in% gene]
                    genes_notFound <<- c(genes_notFound, gene)
                    
                    gene_lookup_group = gene_subset_spatial[1,2]
                    gene_position = 1
                    gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                    x = as.numeric(
                        scan(
                          gene_address, 
                          nlines=1, 
                          sep=",", 
                          what="character",
                          na.strings="NA",
                          skip = gene_position-1
                        )
                      )
                    x[x >= 0] = NA
                 }
                return(x)
              }  
            )
            gene.vector_all = rowMeans(gene.vector_all, na.rm=T) # average if multiple genes are present
          }
          gene.vector_all[is.na(gene.vector_all)] = 0
          
          gene.vector_G = 0
          gene.vector_B = 0
          
          # genes to collect
            rgb.list = list(
              gene.vector_all,
              gene.vector_G,
              gene.vector_B,
              genes_found,
              genes_notFound
            )
            
            names(rgb.list) = c(paste(genes_found, collapse=" "), paste("ignore", collapse=" "), paste("ignore", collapse=" "))

            return(rgb.list)
          }
        )
        names(rgb.list_samples) = sample_id_list
      return(rgb.list_samples)
    }
    
    
  # Plot 2: RGB
  VisiumData2 = function(
      gene.red,
      gene.green,
      gene.blue,
      sample_id_list = sample_id_list
    ){
       
  
    # genes to collect
#      rgb.list = lapply(gene.list,
#        function(gene){
#          gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
#          gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
#          gene_address = paste("../CellMatrix/data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
#        
#          x = as.numeric(
#            scan(
#              gene_address, 
#              nlines=1, 
#              sep=",", 
#              what="character",
#              na.strings="NA",
#              skip = gene_position-1
#            )
#          )
#          return(x)
#        }
#      )
#      
      rgb.list_samples = lapply(sample_id_list,function(i){
        sample_id = i    
      # Reading in data as vectors
        gene.vector_R = 0
        if(length(gene.red) > 0){
            gene.vector_R = sapply(
              gene.red,
              function(gene){
                if(gene %in% gene_subset_spatial[,1]){
                  gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
                  gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
                  gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                  x = as.numeric(
                      scan(
                        gene_address, 
                        nlines=1, 
                        sep=",", 
                        what="character",
                        na.strings="NA",
                        skip = gene_position-1
                      )
                    )
                   } else {
                    gene_lookup_group = gene_subset_spatial[1,2]
                    gene_position = 1
                    gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                    x = as.numeric(
                        scan(
                          gene_address, 
                          nlines=1, 
                          sep=",", 
                          what="character",
                          na.strings="NA",
                          skip = gene_position-1
                        )
                      )
                    x[x >= 0] = NA
                 }
                return(x)
              }  
            )
            gene.vector_R = rowMeans(gene.vector_R, na.rm=T) # average if multiple genes are present
          }
          gene.vector_R[is.na(gene.vector_R)] = 0
          
          gene.vector_G = 0
          if(length(gene.green) > 0){
            gene.vector_G = sapply(
              gene.green,
              function(gene){
              if(gene %in% gene_subset_spatial[,1]){
                  gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
                  gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
                  gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                  x = as.numeric(
                      scan(
                        gene_address, 
                        nlines=1, 
                        sep=",", 
                        what="character",
                        na.strings="NA",
                        skip = gene_position-1
                      )
                    )
                   } else {
                    gene_lookup_group = gene_subset_spatial[1,2]
                    gene_position = 1
                    gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                    x = as.numeric(
                        scan(
                          gene_address, 
                          nlines=1, 
                          sep=",", 
                          what="character",
                          na.strings="NA",
                          skip = gene_position-1
                        )
                      )
                    x[x >= 0] = NA
                 }
                return(x)
              }
            )
            gene.vector_G = rowMeans(gene.vector_G, na.rm=T) # average if multiple genes are present
          }
          gene.vector_G[is.na(gene.vector_G)] = 0
          
          gene.vector_B = 0
          if(length(gene.blue) > 0){
            gene.vector_B = sapply(
              gene.blue,
              function(gene){
               if(gene %in% gene_subset_spatial[,1]){
                  gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
                  gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
                  gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                  x = as.numeric(
                      scan(
                        gene_address, 
                        nlines=1, 
                        sep=",", 
                        what="character",
                        na.strings="NA",
                        skip = gene_position-1
                      )
                    )
                   } else {
                    gene_lookup_group = gene_subset_spatial[1,2]
                    gene_position = 1
                    gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                    x = as.numeric(
                        scan(
                          gene_address, 
                          nlines=1, 
                          sep=",", 
                          what="character",
                          na.strings="NA",
                          skip = gene_position-1
                        )
                      )
                    x[x >= 0] = NA
                 }
                return(x)
              }  
            )
            gene.vector_B = rowMeans(gene.vector_B, na.rm=T) # average if multiple genes are present
          }
          gene.vector_B[is.na(gene.vector_B)] = 0
           
          
          # genes to collect
            rgb.list = list(
              gene.vector_R,
              gene.vector_G,
              gene.vector_B
            )
            
            names(rgb.list) = c(paste(gene.red, collapse=" "), paste(gene.green, collapse=" "), paste(gene.blue, collapse=" "))

            return(rgb.list)
          }
        )
        names(rgb.list_samples) = sample_id_list
      return(rgb.list_samples)
    }

 
  
  # Plot.4 Ligand Receptors
 
  VisiumData4 = function(
      genes_L,
      genes_R,
      genes_selected,
      direction,
      sample_id_list = sample_id_list#, 
      #CPDB_list
    ){
       
   
    # identify direction of genes 
      if(direction == "send"){      
        gene.list_L = genes_L
        gene.list_R = genes_R[genes_R %in% genes_selected]
        gene.list = c(gene.list_L, gene.list_R)
      }
       
      if(direction == "receive"){      
        gene.list_R = genes_R
        gene.list_L = genes_L[genes_L %in% genes_selected]
        gene.list = c(gene.list_R, gene.list_L)
      }
    
          
    # genes to collect
#      rgb.list = lapply(gene.list,
#        function(gene){
#          gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
#          gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
#          gene_address = paste("../CellMatrix/data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
#        
#          x = as.numeric(
#            scan(
#              gene_address, 
#              nlines=1, 
#              sep=",", 
#              what="character",
#              na.strings="NA",
#              skip = gene_position-1
#            )
#          )
#          return(x)
#        }
#      )
#      
      rgb.list_samples = lapply(sample_id_list,function(i){
        sample_id = i    
      # Reading in data as vectors
            gene.vector_L = 0
        if(length(gene.list_L) > 0){
            gene.vector_L = sapply(
              gene.list_L,
              function(gene){
                if(gene %in% gene_subset_spatial[,1]){
                  gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
                  gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
                  gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                  x = as.numeric(
                      scan(
                        gene_address, 
                        nlines=1, 
                        sep=",", 
                        what="character",
                        na.strings="NA",
                        skip = gene_position-1
                      )
                    )
                   } else {
                    gene_lookup_group = gene_subset_spatial[1,2]
                    gene_position = 1
                    gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                    x = as.numeric(
                        scan(
                          gene_address, 
                          nlines=1, 
                          sep=",", 
                          what="character",
                          na.strings="NA",
                          skip = gene_position-1
                        )
                      )
                    x[x >= 0] = NA
                 }
                return(x)
              }  
            )
            gene.vector_L = rowMeans(gene.vector_L) # average if multiple genes are present
          }
          gene.vector_L[is.na(gene.vector_L)] = 0
          
            gene.vector_R = 0
          if(length(gene.list_R) > 0){
            gene.vector_R = sapply(
              gene.list_R,
              function(gene){
                if(gene %in% gene_subset_spatial[,1]){
                  gene_lookup_group = gene_subset_spatial[gene_subset_spatial[,1] == gene,2]
                  gene_position = which(gene_subset_spatial[gene_subset_spatial[,2] == gene_lookup_group,1] == gene) # reading
                  gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                  x = as.numeric(
                      scan(
                        gene_address, 
                        nlines=1, 
                        sep=",", 
                        what="character",
                        na.strings="NA",
                        skip = gene_position-1
                      )
                    )
                   } else {
                    gene_lookup_group = gene_subset_spatial[1,2]
                    gene_position = 1
                    gene_address = paste("data/spatial/", sample_id, "/", "gene_matrices/", gene_lookup_group, "_genes.csv", sep="") 
                    x = as.numeric(
                        scan(
                          gene_address, 
                          nlines=1, 
                          sep=",", 
                          what="character",
                          na.strings="NA",
                          skip = gene_position-1
                        )
                      )
                    x[x >= 0] = NA
                 }
                return(x)
              }  
            )
            gene.vector_R = rowMeans(gene.vector_R) # average if multiple genes are present
          }
          gene.vector_R[is.na(gene.vector_R)] = 0
          
          # genes to collect
            rgb.list = list(
              gene.vector_R,
              gene.vector_L
            )
            names(rgb.list) = c(paste(gene.list_R, collapse=" "), paste(gene.list_L, collapse=" "))

            return(rgb.list)
          }
        )
        names(rgb.list_samples) = sample_id_list
      return(rgb.list_samples)
    }






   # TEST FUNCTION
  plotVisiumEmpty = function(
      light_or_dark,
      mask_overlay,
      plot_size
      )
      {
        
        plotcex = 3.1*plot_size
        
        inverted = as.integer(light_or_dark)
        
        
        
        for(sample_id in sample_id_list){
          rgb.list = rgb.samples_list[[sample_id]]
          x.coords = image_coordinate_list[[sample_id]][,1]
          y.coords = image_coordinate_list[[sample_id]][,2]
          
        # Process opacity statements:          
        # colour processing:
         # if light, fade black to white
#           if(light_or_dark == 1){
#              # Transformation (black to white mapping)
#                plot.light()
#             }
#            
#            if(light_or_dark == 2){
#               plot.dark()
#             }
#                   
       
        # plot the image(s)
          par(mfrow=c(1, length(sample_id_list)))
          display(channel(image_overlay_list[[sample_id]][[inverted]], "grey")*inverted, method="raster", new=F)
          if(mask_overlay == T) {display(image_overlay_list[[sample_id]][[inverted]], method="raster", new=T)}
        
        
      image_dims = dim(image_overlay_list[[sample_id]][[inverted]])
      
      leg1.coords.y = image_dims[1]*0.85
      leg2.coords.y = image_dims[1]*0.9
      leg3.coords.y = image_dims[1]*0.95
      
  #         
      x.coords.1 = image_dims[2]*0.1
      x.coords.2 = x.coords.1 + image_dims[2]*0.1
      x.coords.mid = mean(c(x.coords.1, x.coords.2))
     

          text(image_dims[2]*0.03, image_dims[1]*0.03, paste(sample_id, "", sep="\n"), col=c("black", "white")[inverted], adj=0, cex=plotcex, xpd=NA)
          
#          legend("topleft",
#            legend=paste(sample_id, "16 pcw", sep=": "), 
#            col=NA, lty=0, pch=NA, bty="n", xpd=NA, text.col=c("black", "white")[inverted], ncol=1, cex=plotcex
#          )
      }
    }





   # TEST FUNCTION
  plotVisium1 = function(
      rgb.samples_list,
      light_or_dark,
      plot_size,
      mask_overlay,
      min.col,
      opacity_count_dependent,
      opacity,
      dot_or_not,
      sample_id_list = sample_id_list
      ){
        
        plotcex = 3.1*plot_size
        
        inverted = as.integer(light_or_dark)
        
        
        
        for(sample_id in sample_id_list){
          rgb.list = rgb.samples_list[[sample_id]]
          x.coords = image_coordinate_list[[sample_id]][,1]
          y.coords = image_coordinate_list[[sample_id]][,2]
          
          
                 
        # summation used for other calculations
        #  x = normalise(log(rgb.list[[1]]+0.1))
          x = rgb.list[[1]]
          
          x_colour = normalise(x)
          
          # NA control..
            x[is.na(x)] = 0
            x_colour[is.na(x_colour)]  = 0

        # Pick up initial colour palette:  
        # We aren't using RGB for colour, only for transparency here.
          plot.colours = plot.cols[[as.numeric(light_or_dark)]]
          plot_colours = plotcols(c(0,x), plot.colours[[1]][1:6])
          plot_colours = plot_colours[-1]
          
              
          leg_npoints = 20
          leg_cols = plotcols(seq(0, max(x), length=leg_npoints), plot.colours[[1]][1:6])
          leg_cols = col2rgb(leg_cols)
          leg_cols = leg_cols / 255
          
        # Process opacity statements:          
        # colour processing:
         # if light, fade black to white
           if(light_or_dark == 1){
              # Transformation (black to white mapping)
                plot.light()
                light_cols = col2rgb(plot_colours)
                light_cols = light_cols / 255
              if(opacity_count_dependent == 0){
                plot_colours = rgb(
                  light_cols[1,], 
                  light_cols[2,], 
                  light_cols[3,],
                  opacity
                  )
                leg_cols = rgb(
                  leg_cols[1,], 
                  leg_cols[2,], 
                  leg_cols[3,],
                  opacity
                  )
               } else {
                 plot_colours = rgb(
                   light_cols[1,], 
                   light_cols[2,], 
                   light_cols[3,],
                    x_colour
                 )
                 leg_cols = rgb(
                   leg_cols[1,], 
                   leg_cols[2,], 
                   leg_cols[3,],
                   seq(0, 1, length=leg_npoints)
                 )
               }
             }
            
            if(light_or_dark == 2){
               dark_cols = col2rgb(plot_colours)
               dark_cols = dark_cols / 255
               if(opacity_count_dependent == 0){   
                 plot_colours = rgb(
                   dark_cols[1,], 
                   dark_cols[2,], 
                   dark_cols[3,],
                   opacity
                 )
                 leg_cols = rgb(
                   leg_cols[1,], 
                   leg_cols[2,], 
                   leg_cols[3,],
                   opacity
                 )
               } else {
                 plot_colours = rgb(
                   dark_cols[1,], 
                   dark_cols[2,], 
                   dark_cols[3,],
                    x_colour
                 )
                 leg_cols = rgb(
                   leg_cols[1,], 
                   leg_cols[2,], 
                   leg_cols[3,],
                   seq(0, 1, length=leg_npoints)
                 )
               }
             }
                   
                
        # Cex processing:
          if(dot_or_not == 2){
            pointsize = plotcex * (x+0.1)/max((x+0.1), na.rm=T)
          } else {
            pointsize = plotcex * rep(1, length(x))
          }
        
        # other processing:
          plot.order = order(x)
        
        # plot the image(s)
          par(mfrow=c(1, length(sample_id_list)))
          display(channel(image_overlay_list[[sample_id]][[inverted]], "grey")*inverted, method="raster", new=F)
              points(
                x.coords[plot.order],
                y.coords[plot.order],
                col=plot_colours[plot.order],
                pch=16,
                cex=pointsize[plot.order]
              ) # with colour
          if(mask_overlay == T) {display(image_overlay_list[[sample_id]][[inverted]], method="raster", new=T)}
        
        
      image_dims = dim(image_overlay_list[[sample_id]][[inverted]])
      
      leg1.coords.y = image_dims[1]*0.85
      leg2.coords.y = image_dims[1]*0.9
      leg3.coords.y = image_dims[1]*0.95
      
  #         
      x.coords.1 = image_dims[2]*0.1
      x.coords.2 = x.coords.1 + image_dims[2]*0.1
      x.coords.mid = mean(c(x.coords.1, x.coords.2))
     

      
      
        if(sum(rgb.samples_list[[1]][[1]])>0){
          points(seq(x.coords.1,x.coords.2,length=leg_npoints), rep(leg1.coords.y, leg_npoints), cex=seq(min(pointsize),max(pointsize),length=leg_npoints), pch=16, col=
            leg_cols
          )
  #        
  #        if(length(genes_input) > 1){
  #                text(x.coords.mid, leg1.coords.y+2, "multiple genes", xpd=NA, cex=2)
  #                legend(x.coords.1-4, leg1.coords.y-3, legend=rgb.samples_list[[4]], bty="n", xpd=NA, ncol=4, pch=22, pt.bg=plotcols(gene.sum, c("black", "white"))[-1])
  #              } else {
  #                text(x.coords.mid, leg1.coords.y+2, gene, xpd=NA, cex=2)
  #              }
                
          text(x.coords.1, leg1.coords.y+image_dims[1]*0.015, signif(min(x),3), col=c("black", "white")[inverted], cex=plotcex)
          text(x.coords.2, leg1.coords.y+image_dims[1]*0.015, signif(max(x),3), col=c("black", "white")[inverted], cex=plotcex)
#          text(x.coords.2+image_dims[1]*0.015*plot_size, leg1.coords.y, paste(rgb.samples_list[[sample_id]][[4]], collapse=" "), adj=0, col=plot.colours[[1]][9])
          legend(x.coords.1-image_dims[1]*0.015*plot_size, leg2.coords.y, 
            legend=rgb.samples_list[[sample_id]][[4]], 
            col=NA, lty=0, pch=NA, bty="n", xpd=NA, text.col=plot.colours[[1]][9], ncol=3, cex=plotcex
          )
          
          
          text(image_dims[2]*0.03, image_dims[1]*0.03, paste(sample_id, "", sep="\n"), col=c("black", "white")[inverted], adj=0, cex=plotcex, xpd=NA)
          
#          legend("topleft",
#            legend=paste(sample_id, "16 pcw", sep=": "), 
#            col=NA, lty=0, pch=NA, bty="n", xpd=NA, text.col=c("black", "white")[inverted], ncol=1, cex=plotcex
#          )
        }
      }
    }

  
  
  
  # TEST FUNCTION
  plotVisium2 = function(
      rgb.samples_list,
      light_or_dark,
      plot_size,
      mask_overlay,
      min.col,
      red_on,
      green_on,
      blue_on,
      opacity_count_dependent,
      opacity,
      dot_or_not,
      sample_id_list = sample_id_list
      ){
        plotcex = 3.1*plot_size
        
        inverted = as.integer(light_or_dark)
        
        
        
        for(sample_id in sample_id_list){
          rgb.list = rgb.samples_list[[sample_id]]
          x.coords = image_coordinate_list[[sample_id]][,1]
          y.coords = image_coordinate_list[[sample_id]][,2]
          
       
        # summation used for other calculations
          #x = normalise(log(rgb.list[[1]] + rgb.list[[2]] + rgb.list[[3]]+0.1))
          #x = log(rgb.list[[1]] + rgb.list[[2]] + rgb.list[[3]]+0.1)
          #x = x/max(x)
          
          x.red = rgb.list[[1]] / max(rgb.list[[1]])
          x.green = rgb.list[[2]] / max(rgb.list[[2]])
          x.blue = rgb.list[[3]] / max(rgb.list[[3]])

         x.red = c(0, x.red)
         x.green = c(0, x.green)
         x.blue = c(0, x.blue)
            
         # x.red = normalise(rgb.list[[1]])
         # x.green =  normalise(rgb.list[[2]])
         # x.blue =  normalise(rgb.list[[3]])

          x.red[is.na(x.red)] = 0 
          x.green[is.na(x.green)] = 0           
          x.blue[is.na(x.blue)] = 0       
        
        
          x = (x.red+x.green+x.blue)/3
          x = x/max(x)

          x2 = (x.red+x.green+x.blue)
          x2[x2 > 1] = 1
          
          x[is.na(x)] = 0
          x2[is.na(x2)] = 0  

                    
        # colour processing:
#          if(opacity_count_dependent == 0){
#            plot_colours = rgb(
#              (min.col*(1-x) + red_on * (x.red)*(1-min.col*(1-x))),
#              (min.col*(1-x) + green_on * (x.green)*(1-min.col*(1-x))),
#              (min.col*(1-x) + blue_on * (x.blue)*(1-min.col*(1-x))),
#              opacity
#            )
#          } else {
#            plot_colours = rgb(
#              (min.col*(1-x) + red_on * (x.red)*(1-min.col*(1-x))),
#              (min.col*(1-x) + green_on * (x.green)*(1-min.col*(1-x))),
#              (min.col*(1-x) + blue_on * (x.blue)*(1-min.col*(1-x))),
#              x
#            )
#          }
        
          if(opacity_count_dependent == 0){
            plot_colours = rgb(
              (min.col*(1-x2) + red_on * (x.red)*(1-min.col*(1-x2))),
              (min.col*(1-x2) + green_on * (x.green)*(1-min.col*(1-x2))),
              (min.col*(1-x2) + blue_on * (x.blue)*(1-min.col*(1-x2))),
              opacity
            )
          } else {
            plot_colours = rgb(
              (min.col*(1-x2) + red_on * (x.red)*(1-min.col*(1-x2))),
              (min.col*(1-x2) + green_on * (x.green)*(1-min.col*(1-x2))),
              (min.col*(1-x2) + blue_on * (x.blue)*(1-min.col*(1-x2))),
              x2
            )
          }
        
        
         # if light, fade black to white
           if(light_or_dark == 1) {
                plot.light()
                light_cols = col2rgb(plot_colours)
                light_cols = light_cols / 255
                plot_colours = create_white_fade_colormap_alt(
                  light_cols[1,], 
                  light_cols[2,], 
                  light_cols[3,]
                  )
              if(opacity_count_dependent == 0){   
                light_cols = col2rgb(plot_colours)
                light_cols = light_cols / 255
                plot_colours = rgb(
                  light_cols[1,], 
                  light_cols[2,], 
                  light_cols[3,],
                  opacity
                  )
               } else {
                 light_cols = col2rgb(plot_colours)
                 light_cols = light_cols / 255
                 plot_colours = rgb(
                   light_cols[1,], 
                   light_cols[2,], 
                   light_cols[3,],
                   x2
                 )
               }
             }
                   
        plot_colours = plot_colours[-1]
            
        # Cex processing:
          if(dot_or_not == 2){
            pointsize = plotcex * x2[-1]
          } else {
            pointsize = plotcex * rep(1, length(x)-1)
          }
        
        # other processing:
          plot.order = order(x)
        
        # plot the image(s)
          par(mfrow=c(1, length(sample_id_list)))
          display(channel(image_overlay_list[[sample_id]][[inverted]], "grey")*inverted, method="raster", new=F)
              points(
                x.coords[plot.order],
                y.coords[plot.order],
                col=plot_colours[plot.order],
                pch=16,
                cex=pointsize[plot.order]
              ) # with colour
           
          if(mask_overlay == T) {display(image_overlay_list[[sample_id]][[inverted]], method="raster", new=T)}
        }
    
      image_dims = dim(image_overlay_list[[sample_id]][[inverted]])
      
      leg1.coords.y = image_dims[1]*0.85
      leg2.coords.y = image_dims[1]*0.9
      leg3.coords.y = image_dims[1]*0.95
      
  #         
      x.coords.1 = image_dims[2]*0.1
      x.coords.2 = x.coords.1 + image_dims[2]*0.1
      x.coords.mid = mean(c(x.coords.1, x.coords.2))
     

  
      
      leg_npoints = 20
      if(sum(rgb.samples_list[[1]][[1]])>0){
        points(seq(x.coords.1,x.coords.2,length=leg_npoints), rep(leg1.coords.y, leg_npoints), cex=seq(min(pointsize),max(pointsize),length=leg_npoints), pch=16, col=
          rgb(
              (min.col*(1-seq(0.0,1,length=leg_npoints)) + red_on * (seq(0.0,1,length=leg_npoints))*(1-min.col*(1-seq(0.0,1,length=leg_npoints)))),
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              opacity*seq(1-opacity_count_dependent,1,length=leg_npoints),
            )
        )
        # text(x.coords.1, leg1.coords.y+image_dims[1]*0.015, signif(min(rgb.list[[1]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.1, leg1.coords.y+image_dims[1]*0.015, 0, col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2, leg1.coords.y+image_dims[1]*0.015, signif(max(rgb.list[[1]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2+image_dims[1]*0.015*plot_size, leg1.coords.y, names(rgb.samples_list[[1]])[1], col="red", adj=0, cex=plotcex)
      }
      if(sum(rgb.samples_list[[1]][[2]])>0){
        points(seq(x.coords.1,x.coords.2,length=leg_npoints), rep(leg2.coords.y, leg_npoints), cex=seq(min(pointsize),max(pointsize),length=leg_npoints), pch=16, col=
          rgb(
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              (min.col*(1-seq(0.0,1,length=leg_npoints)) + green_on * (seq(0.0,1,length=leg_npoints))*(1-min.col*(1-seq(0.0,1,length=leg_npoints)))),
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              opacity*seq(1-opacity_count_dependent,1,length=leg_npoints),
            )
        )
        #text(x.coords.1, leg2.coords.y+image_dims[1]*0.015, signif(min(rgb.list[[2]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.1, leg2.coords.y+image_dims[1]*0.015, 0, col=c("black", "white")[inverted], cex=plotcex)  
        text(x.coords.2, leg2.coords.y+image_dims[1]*0.015, signif(max(rgb.list[[2]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2+image_dims[1]*0.015*plot_size, leg2.coords.y, names(rgb.samples_list[[1]])[2], col="green", adj=0, cex=plotcex)
      }
      if(sum(rgb.samples_list[[1]][[3]])>0){
        points(seq(x.coords.1,x.coords.2,length=leg_npoints), rep(leg3.coords.y, leg_npoints), cex=seq(min(pointsize),max(pointsize),length=leg_npoints), pch=16, col=
          rgb(
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              (min.col*(1-seq(0.0,1,length=leg_npoints)) + blue_on * (seq(0.0,1,length=leg_npoints))*(1-min.col*(1-seq(0.0,1,length=leg_npoints)))),
              opacity*seq(1-opacity_count_dependent,1,length=leg_npoints),
            )
        )
        #text(x.coords.1, leg3.coords.y+image_dims[1]*0.015, signif(min(rgb.list[[3]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.1, leg3.coords.y+image_dims[1]*0.015, 0, col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2, leg3.coords.y+image_dims[1]*0.015, signif(max(rgb.list[[3]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2+image_dims[1]*0.015*plot_size, leg3.coords.y, names(rgb.samples_list[[1]])[3], col="blue", adj=0, cex=plotcex)
      }
   
    text(image_dims[2]*0.03, image_dims[1]*0.03, paste(sample_id, "", sep="\n"), col=c("black", "white")[inverted], adj=0, cex=plotcex, xpd=NA)
          
#          legend("topleft",
#            legend=paste(sample_id, "16 pcw", sep=": "), 
#            col=NA, lty=0, pch=NA, bty="n", xpd=NA, text.col=c("black", "white")[inverted], ncol=1, cex=plotcex
#          )
    }





  # TEST FUNCTION
  plotVisium4 = function(
      rgb.samples_list,
      light_or_dark,
      plot_size,
      mask_overlay,
      min.col,
      red_on,
      blue_on,
      opacity_count_dependent,
      opacity,
      dot_or_not,
      sample_id_list = sample_id_list
      ){
        plotcex = 3.1*plot_size
        
        inverted = as.integer(light_or_dark)
        
        
        for(sample_id in sample_id_list){
          rgb.list = rgb.samples_list[[sample_id]]
          x.coords = image_coordinate_list[[sample_id]][,1]
          y.coords = image_coordinate_list[[sample_id]][,2]
          
          
                 
        # summation used for other calculations
          x = normalise(log(rgb.list[[1]] + rgb.list[[2]]+0.1))
          x.red = normalise(rgb.list[[1]])
          x.blue = normalise(rgb.list[[2]])
          
          # NA control..
            x[is.na(x)] = 0
            x.blue[is.na(x.blue)] = 0 
            x.red[is.na(x.red)] = 0 
          
              
        # colour processing:
          if(opacity_count_dependent == 0){
            plot_colours = rgb(
              (min.col*(1-x) + red_on * (x.red)*(1-min.col*(1-x))),
              min.col*(1-x),
              (min.col*(1-x) + blue_on * (x.blue)*(1-min.col*(1-x))),
              opacity
            )
          } else {
            plot_colours = rgb(
              (min.col*(1-x) + red_on * (x.red)*(1-min.col*(1-x))),
              min.col*(1-x),
              (min.col*(1-x) + blue_on * (x.blue)*(1-min.col*(1-x))),
              x
            )
          }
                
                
        # Cex processing:
          if(dot_or_not == 2){
            pointsize = plotcex * x
          } else {
            pointsize = plotcex * rep(1, length(x))
          }
        
        # other processing:
          plot.order = order(x)
        
        # plot the image(s)
          par(mfrow=c(1, length(sample_id_list)))
          display(channel(image_overlay_list[[sample_id]][[inverted]], "grey")*inverted, method="raster", new=F)
              points(
                x.coords[plot.order],
                y.coords[plot.order],
                col=plot_colours[plot.order],
                pch=16,
                cex=pointsize[plot.order]
              ) # with colour
          if(mask_overlay == T) {display(image_overlay_list[[sample_id]][[inverted]], method="raster", new=T)}
        }
    
      image_dims = dim(image_overlay_list[[sample_id]][[inverted]])
      
      leg1.coords.y = image_dims[1]*0.85
      leg2.coords.y = image_dims[1]*0.9
      leg3.coords.y = image_dims[1]*0.95
      
  #         
      x.coords.1 = image_dims[2]*0.1
      x.coords.2 = x.coords.1 + image_dims[2]*0.1
      x.coords.mid = mean(c(x.coords.1, x.coords.2))
     

  
      

          leg_npoints = 20
          
     
         
      if(sum(rgb.samples_list[[1]][[1]])>0){
        points(seq(x.coords.1,x.coords.2,length=leg_npoints), rep(leg1.coords.y, leg_npoints), cex=seq(min(pointsize),max(pointsize),length=leg_npoints), pch=16, col=
          rgb(
              (min.col*(1-seq(0.0,1,length=leg_npoints)) + red_on * (seq(0.0,1,length=leg_npoints))*(1-min.col*(1-seq(0.0,1,length=leg_npoints)))),
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              opacity*seq(1-opacity_count_dependent,1,length=leg_npoints),
            )
        )
        text(x.coords.1, leg1.coords.y+image_dims[1]*0.015, signif(min(rgb.list[[1]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2, leg1.coords.y+image_dims[1]*0.015, signif(max(rgb.list[[1]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2+image_dims[1]*0.015*plot_size, leg1.coords.y, names(rgb.list)[1], col="red", adj=0, cex=plotcex)
      }

      if(sum(rgb.samples_list[[1]][[2]])>0){
        points(seq(x.coords.1,x.coords.2,length=leg_npoints), rep(leg2.coords.y, leg_npoints), cex=seq(min(pointsize),max(pointsize),length=leg_npoints), pch=16, col=
          rgb(
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              (min.col*(1-seq(0.0,1,length=leg_npoints))),
              (min.col*(1-seq(0.0,1,length=leg_npoints)) + blue_on * (seq(0.0,1,length=leg_npoints))*(1-min.col*(1-seq(0.0,1,length=leg_npoints)))),
              opacity*seq(1-opacity_count_dependent,1,length=leg_npoints),
            )
        )
        text(x.coords.1, leg2.coords.y+image_dims[1]*0.015, signif(min(rgb.list[[2]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2, leg2.coords.y+image_dims[1]*0.015, signif(max(rgb.list[[2]]),3), col=c("black", "white")[inverted], cex=plotcex)
        text(x.coords.2+image_dims[1]*0.015*plot_size, leg2.coords.y, names(rgb.list)[2], col="blue", adj=0, cex=plotcex)
      }
      
   text(image_dims[2]*0.03, image_dims[1]*0.03, paste(sample_id, "", sep="\n"), col=c("black", "white")[inverted], adj=0, cex=plotcex, xpd=NA)
          
#          legend("topleft",
#            legend=paste(sample_id, "16 pcw", sep=": "), 
#            col=NA, lty=0, pch=NA, bty="n", xpd=NA, text.col=c("black", "white")[inverted], ncol=1, cex=plotcex
#          )
#       
#      points(seq(x.coords.1,x.coords.2,length=leg_npoints), rep(leg1.coords.y, leg_npoints), cex=seq(min(pointsize),max(pointsize),length=leg_npoints), pch=16, col=
#        rgb(
#            (min.col*(1-seq(0.0,1,length=leg_npoints)) + red_on * (seq(0.0,1,length=leg_npoints))*(1-min.col*(1-seq(0.0,1,length=leg_npoints)))),
#            (min.col*(1-seq(0.0,1,length=leg_npoints))),
#            (min.col*(1-seq(0.0,1,length=leg_npoints))),
#            opacity*seq(1-opacity_count_dependent,1,length=leg_npoints),
#          )
#      )
#      text(x.coords.1, leg1.coords.y+image_dims[1]*0.015, signif(min(rgb.list[[1]]),3), col="white", cex=3.3*plot_size)
#      text(x.coords.2, leg1.coords.y+image_dims[1]*0.015, signif(max(rgb.list[[1]]),3), col="white", cex=3.3*plot_size)
#      
#      points(seq(x.coords.1,x.coords.2,length=leg_npoints), rep(leg2.coords.y, leg_npoints), cex=seq(min(pointsize),max(pointsize),length=leg_npoints), pch=16, col=
#        rgb(
#            (min.col*(1-seq(0.0,1,length=leg_npoints))),
#            (min.col*(1-seq(0.0,1,length=leg_npoints))),
#            (min.col*(1-seq(0.0,1,length=leg_npoints)) + blue_on * (seq(0.0,1,length=leg_npoints))*(1-min.col*(1-seq(0.0,1,length=leg_npoints)))),
#            opacity*seq(1-opacity_count_dependent,1,length=leg_npoints),
#          )
#      )
#      text(x.coords.1, leg2.coords.y+image_dims[1]*0.015, signif(min(rgb.list[[2]]),3), col=c("black", "white")[light_or_dark], cex=3.3*plot_size)
#      text(x.coords.2, leg2.coords.y+image_dims[1]*0.015, signif(max(rgb.list[[2]]),3), col=c("white", "black")[light_or_dark], cex=3.3*plot_size)
#      
      

      }




 
# Lists for pre-selection
  
  gene_lookup_table <<- read.csv("data/gene_lookup_table.csv", head=F)
  gene_subset_spatial <<- read.csv("data/spatial/gene_lookup_spatial.csv", head=F)
  
  col.scale.1 = plotcols(1:9, c("black", verdant))
  col.scale.2 = plotcols(1:9, c("black", rev(Heat), tail(verdant,1)))
  col.scale.3 = plotcols(1:9, c("black",cambridge.sunrise[1:6], tail(verdant,1)))
  
# colour set ups
  plot.cols <<- list(
      "light" = list(
        rev(col.scale.1[-1]),
        rev(col.scale.2[-1]),
        rev(col.scale.3[-1]),
        "white"
      ),
      "dark" = list(
        col.scale.1,
        col.scale.2,
        col.scale.3,
        rgb(0.06,0.06,0.06)
      )
    )

   plot.dims <<- function(variable.list, annotations2keep, genes_input, plot_size, plot_type){
   
    if(length(variable.list) == 1){
        variable.1 = variable.list[1]
       # gene.list = unique(unlist(strsplit(genes_input, split=",|\t|\n| ")))
       # gene.list = gene.list[gene.list != ""]
       # gene.list = gsub(" ", "", gene.list)
       # gene.list = toupper(gene.list)
       # gene.list = gene.list[gene.list %in% gene_lookup_table[,1]]
      
        var1.names = read.csv(paste("data/gene_matrices/", variable.1, "/", variable.1, "/", "column_names.csv", sep=""), head=F)[,1]   
        
         if(sum(annotations2keep %in% var1.names) > 0){
          var1.names = var1.names[na.omit(match(annotations2keep, var1.names))]
         }
        
        matrix.dim1 <- length(genes_input)
        matrix.dim2 <- length(var1.names)
      }

      if(length(variable.list) == 2){
        variable.1 = variable.list[1]
        variable.2 = variable.list[2]

        var1.names = read.csv(paste("data/gene_matrices/", variable.1, "/", variable.1, "/", "column_names.csv", sep=""), head=F)[,1]
        var2.names = read.csv(paste("data/gene_matrices/", variable.1, "/", variable.2, "/", "row_names.csv", sep=""), head=F)[,1]



         if(sum(annotations2keep %in% var1.names) > 0){
          var1.names = var1.names[na.omit(match(annotations2keep, var1.names))]
         }
        
         if(sum(annotations2keep %in% var2.names) > 0){
          var2.names = var2.names[na.omit(match(annotations2keep, var2.names))]
         }
        
        matrix.dim1 <- length(var2.names)
        matrix.dim2 <- length(var1.names)
      }


      if(length(variable.list) == 3){
        variable.1 = variable.list[1]
        variable.2 = variable.list[2]
        variable.3 = variable.list[3]

        var1.names = read.csv(paste("data/gene_matrices/", variable.1, "/", variable.2, "/", "column_names.csv", sep=""), head=F)[,1]
        var2.names = read.csv(paste("data/gene_matrices/", variable.1, "/", variable.2, "/", "row_names.csv", sep=""), head=F)[,1]
        var3.names = read.csv(paste("data/gene_matrices/", variable.2, "/", variable.3, "/", "row_names.csv", sep=""), head=F)[,1]
        
        
         if(sum(annotations2keep %in% var1.names) > 0){
          var1.names = var1.names[na.omit(match(annotations2keep, var1.names))]
         }
        
         if(sum(annotations2keep %in% var2.names) > 0){
          var2.names = var2.names[na.omit(match(annotations2keep, var2.names))]
         }
         
         if(sum(annotations2keep %in% var3.names) > 0){
          var3.names = var3.names[na.omit(match(annotations2keep, var3.names))]
         }
         
        matrix.dim1 <- length(var2.names) + length(var3.names) + 2
        matrix.dim2 <- length(var1.names) + length(var3.names) + 2
      }

      
      if(plot_type == 2){
      
      
        variable.1 = variable.list[1]
        variable.2 = variable.list[2]

        var1.names = read.csv(paste("data/gene_matrices/", variable.1, "/", variable.1, "/", "column_names.csv", sep=""), head=F)[,1]
        var2.names = read.csv(paste("data/gene_matrices/", variable.1, "/", variable.2, "/", "row_names.csv", sep=""), head=F)[,1]



         if(sum(annotations2keep %in% var1.names) > 0){
          var1.names = var1.names[na.omit(match(annotations2keep, var1.names))]
         }
        
         if(sum(annotations2keep %in% var2.names) > 0){
          var2.names = var2.names[na.omit(match(annotations2keep, var2.names))]
         }
        
        matrix.dim1 <- length(var2.names)
        matrix.dim2 <- length(var1.names)
      }
      
      
      if(plot_type == 3){
        var1.names = read.csv(paste("data/gene_matrices/", "fine_grain", "/", "fine_grain", "/", "column_names.csv", sep=""), head=F)[,1]
        var2.names = read.csv(paste("data/gene_matrices/", "week", "/", "week", "/", "column_names.csv", sep=""), head=F)[,1]
#       
#         if(sum(annotations2keep %in% var1.names) > 0){
#          var1.names = var1.names[na.omit(match(annotations2keep, var1.names))]
#         }
#        
#         if(sum(annotations2keep %in% var2.names) > 0){
#          var2.names = var2.names[na.omit(match(annotations2keep, var2.names))]
#         }
#        
        matrix.dim1 <- length(var2.names) + 20 # add on the empty space
        matrix.dim2 <- length(var1.names) 
      }
        # Figure size...
        margin.size = 778
        square.size = 32
        between.size = 2
        axis_distance = 18 
      
      figdims = c(
        matrix.dim2*square.size + (matrix.dim2-1)*between.size + 2*margin.size + 2*axis_distance,
        matrix.dim1*square.size + (matrix.dim1-1)*between.size + 2*margin.size + 2*axis_distance
        )*plot_size
    return(figdims)
   
   }
   
   
   
   format_genes = function(genes_input, gene.red, gene.green, gene.blue, plot_type){
     if(plot_type == 1){
        gene.list = unique(unlist(strsplit(genes_input, split=",|\t|\n| ")))
        gene.list = gene.list[gene.list != ""]
        gene.list = gsub(" ", "", gene.list)
        gene.list = toupper(gene.list)
        gene.list = gene.list[gene.list %in% gene_lookup_table[,1]]
        
     }
     if(plot_type == 2){
             
       gene.red = toupper(gene.red)
       gene.red = gene.red[gene.red %in% gene_lookup_table[,1]]
       if(length(gene.red) == 0){gene.red=""}
      
       gene.green = toupper(gene.green)
       gene.green = gene.green[gene.green %in% gene_lookup_table[,1]]
       if(length(gene.green) == 0){gene.green=""}
      
       gene.blue = toupper(gene.blue)
       gene.blue = gene.blue[gene.blue %in% gene_lookup_table[,1]]
       if(length(gene.blue) == 0){gene.blue=""}
      
#      gene.red = "UPK3B"
#      gene.blue = "ITLN1"
#      gene.green = "CPB1"
      
       gene.list = c(gene.red, gene.green, gene.blue)
       gene.list = gene.list[gene.list != ""]
       
     }
     if(plot_type == 3){
        gene.list = unique(unlist(strsplit(genes_input, split=",|\t|\n| ")))
        gene.list = gene.list[gene.list != ""]
        gene.list = gsub(" ", "", gene.list)
        gene.list = toupper(gene.list)
        gene.list = gene.list[gene.list %in% gene_lookup_table[,1]]
        
     }
     return(gene.list)
   }
   
   collate_genes = function(area, red, green, blue){
    return(c(area, red, green, blue))
   }
   
   


# Colour map?
   create_white_fade_colormap_alt <<- function(r_values, g_values, b_values) {
  #    """
  #    Alternative implementation: Each channel represents intensity from white
  #    Higher values move that channel away from white (1) towards 0
  #    """
      
      # Validate inputs (same as main function)
      if (!is.numeric(r_values) || !is.numeric(g_values) || !is.numeric(b_values)) {
        stop("All inputs must be numeric vectors")
      }
      if (length(r_values) != length(g_values) || length(r_values) != length(b_values)) {
        stop("All input vectors must have the same length")
      }
      if (any(r_values < 0 | r_values > 1) || 
          any(g_values < 0 | g_values > 1) || 
          any(b_values < 0 | b_values > 1)) {
        stop("All values must be between 0 and 1")
      }
      
      # Simple approach: intensity directly controls how much to reduce from white
      final_r <- 1 - (1 - r_values) * (g_values + b_values > 0)
      final_g <- 1 - (1 - g_values) * (r_values + b_values > 0)  
      final_b <- 1 - (1 - b_values) * (r_values + g_values > 0)
      
      # Alternative: more intuitive approach
      # Each channel intensity reduces the OTHER channels
      final_r <- pmax(r_values, 1 - g_values - b_values)
      final_g <- pmax(g_values, 1 - r_values - b_values)
      final_b <- pmax(b_values, 1 - r_values - g_values)
      
      # Ensure values stay within 0-1 range
      final_r <- pmax(0, pmin(1, final_r))
      final_g <- pmax(0, pmin(1, final_g))
      final_b <- pmax(0, pmin(1, final_b))
      
      colors <- rgb(final_r, final_g, final_b)
      return(colors)
    }   
     
   
############ Hierarchy list, for plotting and annotation mapping   ###################   
    hierarchy_list = list(
        "Cardiomyocytes" = list(
           "AtrialCardiomyocytes" = list(
              "AtrialCardiomyocytesLeft",
              "AtrialCardiomyocytesRight",
              "AtrialCardiomyocytesCycling"
            ),
            "VentricularCardiomyocytes" = list(
              "VentricularCardiomyocytesLeftCompact",
              "VentricularCardiomyocytesRightCompact",
              "VentricularCardiomyocytesLeftTrabeculated",
              "VentricularCardiomyocytesRightTrabeculated",
              "VentricularCardiomyocytesCycling"
            ),
            "CardiacConductionSystem" = list(
              "SinoatrialNodePacemakerCells",
              "AtrioventricularNodePacemakerCells",
              "VentricularConductionSystemProximal",
              "VentricularConductionSystemDistal"
            )
        ),
        "Mesenchymal" = list(
          "Fibroblasts" = list(
            "GreatVesselAdventitialFibroblasts",  
            "CoronaryVesselAdventitialFibroblasts",
            "MyocardialInterstitialFibroblasts",
            "SubEpicardialFibroblasts",   
            "Myofibroblasts",                            
            "LymphNodeFibroblasticReticularCells",       
            "ValveInterstitialCells"
          ),                     
          "MuralCells" = list(
            "GreatVesselSmoothMuscleCells",              
            "CoronarySmoothMuscleCells",                 
            "DuctusArteriosusSmoothMuscleCells",         
            "CoronaryPericytes" 
          ),
          "PericardialCells" = list(
            "PericardialCellsIntermediate",              
            "PericardialCellsFibrous",                   
            "PericardialCellsParietal" 
          )
        ),
        "Endothelium" = list(
          "BloodVesselEndothelialCells" = list(
            "GreatVesselArterialEndothelialCells",       
            "GreatVesselVenousEndothelialCells",         
            "CoronaryArterialEndothelialCells",          
            "CoronaryVenousEndothelialCells",            
            "CoronaryCapillaryEndothelialCells" 
          ),
          "EndocardialCells" = list(
            "EndocardialCells",                  
            "EndocardialCushionCells",                   
            "ValveEndothelialCells"
          ),
          "LymphaticEndothelialCells" = list(
            "LymphaticEndothelialCells"
          )
        ),
        "Epicardium" = list(
          "EpicardialCells" = list(
            "MesothelialEpicardialCells",                
            "EpicardiumDerivedCells"
          )
        ),
        "Neural" = list(
          "Neurons" = list(
            "NeuronPrecursors",                          
            "ChromaffinCells",                           
            "SympatheticNeurons",                        
            "ParasympatheticNeurons"
          ),
          "Glia" = list(
            "SchwannCellPrecursors",                     
            "SchwannCells"
          )
        ),
        "Leukocytes" = list(
          "MyeloidCells" = list(
            "MonocytesMPOpos",                           
            "Monocytes",                                 
            "MonocyteDerivedCells",                      
            "MacrophagesCX3CR1pos",                      
            "MacrophagesTIMD4pos",                       
            "MacrophagesLYVE1pos",                       
            "MacrophagesATF3pos",                        
            "DendriticCellsType1",                       
            "DendriticCellsMature",                      
            "PlasmacytoidDendriticCells",                
            "MastCells",                                 
            "Megakaryocytes"
          ),
          "LymphoidCells" = list(
            "TCellsCD4pos",                              
            "TCellsCD8pos",                              
            "TregsCD4pos",                               
            "ProBCells",                                 
            "BCells",                                    
            "BCellsMS4A1pos",                            
            "NaturalKillerCells",                        
            "InnateLymphoidCells"
          )
        )
      )
      
   
    hierarchy_list2 = list(
        "Cardiomyocytes" = list(
           "aCMs" = list(
              "aCMs - Left",
              "aCMs - Right",
              "aCMs - Cycling"
            ),
            "vCMs" = list(
              "vCMs - LeftCompact",
              "vCMs - RightCompact",
              "vCMs - LeftTrabeculated",
              "vCMs - RightTrabeculated",
              "vCMs - Cycling"
            ),
            "CardiacConductionSystem" = list(
              "SAN - PacemakerCells",
              "AVN - PacemakerCells",
              "vConductionSystemProximal",
              "vConductionSystemDistal"
            )
        ),
        "Mesenchymal" = list(
          "Fibroblasts" = list(
            "GreatVesselAdventitialFibroblasts",  
            "CoronaryVesselAdventitialFibroblasts",
            "MyocardialInterstitialFibroblasts",
            "SubEpicardialFibroblasts",   
            "Myofibroblasts",                            
            "LymphNodeFibroblasticReticularCells",       
            "ValveInterstitialCells"
          ),                     
          "MuralCells" = list(
            "GreatVesselSmoothMuscleCells",              
            "CoronarySmoothMuscleCells",                 
            "DuctusArteriosusSmoothMuscleCells",         
            "CoronaryPericytes" 
          ),
          "PericardialCells" = list(
            "PericardialCellsIntermediate",              
            "PericardialCellsFibrous",                   
            "PericardialCellsParietal" 
          )
        ),
        "Endothelium" = list(
          "BloodVesselEndothelialCells" = list(
            "GreatVesselArterialEndothelialCells",       
            "GreatVesselVenousEndothelialCells",         
            "CoronaryArterialEndothelialCells",          
            "CoronaryVenousEndothelialCells",            
            "CoronaryCapillaryEndothelialCells" 
          ),
          "EndocardialCells" = list(
            "EndocardialCells",                  
            "EndocardialCushionCells",                   
            "ValveEndothelialCells"
          ),
          "LymphaticEndothelialCells" = list(
            "LymphaticEndothelialCells"
          )
        ),
        "Epicardium" = list(
          "EpicardialCells" = list(
            "MesothelialEpicardialCells",                
            "EpicardiumDerivedCells"
          )
        ),
        "Neural" = list(
          "Neurons" = list(
            "NeuronPrecursors",                          
            "ChromaffinCells",                           
            "SympatheticNeurons",                        
            "ParasympatheticNeurons"
          ),
          "Glia" = list(
            "SchwannCellPrecursors",                     
            "SchwannCells"
          )
        ),
        "Leukocytes" = list(
          "MyeloidCells" = list(
            "MonocytesMPOpos",                           
            "Monocytes",                                 
            "MonocyteDerivedCells",                      
            "MacrophagesCX3CR1pos",                      
            "MacrophagesTIMD4pos",                       
            "MacrophagesLYVE1pos",                       
            "MacrophagesATF3pos",                        
            "DendriticCellsType1",                       
            "DendriticCellsMature",                      
            "PlasmacytoidDendriticCells",                
            "MastCells",                                 
            "Megakaryocytes"
          ),
          "LymphoidCells" = list(
            "TCellsCD4pos",                              
            "TCellsCD8pos",                              
            "TregsCD4pos",                               
            "ProBCells",                                 
            "BCells",                                    
            "BCellsMS4A1pos",                            
            "NaturalKillerCells",                        
            "InnateLymphoidCells"
          )
        )
      )
     
   # Load in all images.... (light and dark plot variations too!)
  
  rgb.samples_list = list()
  
  image_overlay_list = lapply(sample_id_list, function(sample_id){
      image_overlays = list(
        readImage(file.path("data/spatial", sample_id, "images", 1, "overlay.png")), # light
        readImage(file.path("data/spatial", sample_id, "images", 2, "overlay.png"))  # dark
      )
      return(image_overlays)
    }
  )
  names(image_overlay_list) = sample_id_list
  
  image_coordinate_list = lapply(sample_id_list, function(sample_id){
      sample_coordinates = read.csv(file.path("data/spatial", sample_id, "coordinates.csv"), head=F)
      return(sample_coordinates)
    }
  )
  names(image_coordinate_list) = sample_id_list
   
    
