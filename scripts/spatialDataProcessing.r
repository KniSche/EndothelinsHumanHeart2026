# Then we can feed this into a function...
  #  CellMatrix::
    preprocessImage = function(
      image_h5ad,
      rotation = 0,
      reflection = "none",
      opacity_filter = 0.9,
      sample_id,
      evaluate_image = T,
      feature_names = "default",
      nCounts_column = "n_counts",
      nGenes_column = "n_genes",
      QC_threshold = 600,
      contrast_mod = 1
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
      good_spots = which(spot_d > QC_threshold)
    
    if(evaluate_image == T){  
    cat("\n 3. \"evaluate_image\" set to TRUE:") 
    cat("\n ... Skipping spot matrix processing and gene expression writing")
    }
    
    if(evaluate_image != T){
    mkdir("data/spatial/",
              sample_id,
              "/",
              "gene_matrices/")
                    
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
        if(feature_names[1] == "default"){
          colnames(spot_matrix) = image_adata$var$SYMBOL
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
               if(contrast_mod != 1){
               	pseudoGrey = pseudoGrey-0.5
               	pseudoGrey = pseudoGrey * contrast_mod
               	pseudoGrey = pseudoGrey+0.5
               }
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
      	mkdir(file.path("data/spatial", sample_id, "images", 1))
            writeImage(image_overlay, file.path("data/spatial", sample_id, "images", 1, "overlay.png"))
          # We can actually use the same image, but change the channel to grey scale
         }
       
       
       
    ################ # plot 2
    ############################################
        # DARK PLOT (black background)
          pseudoGrey <- 1-(normalise((image_data[,,1]/3 + image_data[,,2]/3 + image_data[,,3]/3)))
             if(contrast_mod != 1){
               	pseudoGrey = pseudoGrey-0.5
               	pseudoGrey = pseudoGrey * contrast_mod
               	pseudoGrey = pseudoGrey+0.5
               }
          
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
                
            mkdir(file.path("data/spatial", sample_id, "images", 2))
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
    
    
    
    
    
    
 
