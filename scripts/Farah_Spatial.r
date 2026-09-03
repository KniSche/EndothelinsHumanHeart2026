# Scripts and functions:
  func.dir="/home/vincent/Documents/Rlib/functions/"
  source("/home/vincent/Documents/Rlib/functions/load_functions.r")
   
 
######################################################



######################################################
######################################################
######################################################

# Dataset

# FARAH2024
datasets = c("x13PCW_1", "x13PCW_2", "x13PCW_3", "x15PCW")
n_squares = 100
figure_dir = "figures/spatial/Farah2024/"
data_dir = "data/spatial/Farah2024/"	
unique.genes = gsub("\\.csv", "", list.files(paste(data_dir, "x13PCW_1/", sep="")))
  	
	
	

### FIGURES
n_squares = 100

for(dataset in datasets){
	# gene selection
	  gene.r = "MYH7"
	  gene.g = "EDNRA"
	  gene.b = "MYH6"
	 
	# gene prep
	  if(!gene.r %in% unique.genes){
	    gene.r = "empty_gene"
	  }
	  if(!gene.g %in% unique.genes){
	    gene.g = "empty_gene"
	  }
	  if(!gene.b %in% unique.genes){
	    gene.b = "empty_gene"
	  }   
	
	threshold = 1
	dataset_meta = read.csv(paste(data_dir, dataset,"/plotdata.csv", sep=""))
	
          # get values
	  gene.r.cells = read.csv(paste(data_dir, dataset,"/", gene.r, ".csv", sep=""), head=F)[,1]
	  gene.g.cells = read.csv(paste(data_dir, dataset,"/", gene.g, ".csv", sep=""), head=F)[,1]
	  gene.b.cells = read.csv(paste(data_dir, dataset,"/", gene.b, ".csv", sep=""), head=F)[,1]
	  # scale  
	  gene.r.cells = gene.r.cells / max(gene.r.cells)
	  gene.g.cells = gene.g.cells / max(gene.g.cells)
	  gene.b.cells = gene.b.cells / max(gene.b.cells)
	  # remove NA
	  gene.r.cells[is.na(gene.r.cells)] = 0
	  gene.g.cells[is.na(gene.g.cells)] = 0
	  gene.b.cells[is.na(gene.b.cells)] = 0
	  # create RGB space
	  col2plot = rgb(
	    gene.r.cells,
	    gene.g.cells,
	    gene.b.cells
	  )
	  
	png(paste(figure_dir, "red", gene.r, "_green", gene.g, "_blue", gene.b, "_", dataset, ".png", sep=""), width=500, height=500)
	  	plot.dark()
	  	par("bg" = rgb(0.08, 0.08, 0.08))
		  plot(
		    -dataset_meta[,"global_y"],
		    dataset_meta[,"global_x"],
		    pch=".",
		    cex=dataset_meta[,"volume"]/20000,
		    col=col2plot,
		    asp=T,
		    xlab=dataset,
		    ylab="",
		    bty="n",
		    xaxt="n",
		    yaxt="n"
		  )
	mtext(gsub("empty_gene", "", gene.r), adj=0, col="red")
	mtext(gsub("empty_gene", "", gene.g), adj=0.5, col="green")
	mtext(gsub("empty_gene", "", gene.b), adj=1, col="blue")
	dev.off()

	
	bin.dist = c(
		diff(range(dataset_meta[,"global_y"]))/n_squares,
		diff(range(dataset_meta[,"global_y"]))/n_squares
	)
	
	x.r = collect.points(
	 -dataset_meta[,"global_y"],
	 dataset_meta[,"global_x"],
	 colours=alt.cols,
	 weights = gene.r.cells,
	 method="derp",
	 bin.dist=bin.dist
	 )
	x.g = collect.points(
	 -dataset_meta[,"global_y"],
	 dataset_meta[,"global_x"],
	 colours=alt.cols,
	 weights = gene.g.cells,
	 method="derp",
	 bin.dist=bin.dist
	 )
	x.b = collect.points(
	 -dataset_meta[,"global_y"],
	 dataset_meta[,"global_x"],
	 colours=alt.cols,
	 weights = gene.b.cells,
	 method="derp",
	 bin.dist=bin.dist
	 )
	
	 col2plot = rgb(
	    x.r[,3],
	    x.g[,3],
	    x.b[,3]
	  )
	  
 	pdf(paste(figure_dir, "red", gene.r, "_green", gene.g, "_blue", gene.b, "_", dataset, "_CollectedDensity.pdf", sep=""), height=bin.dist[1]/5.5, width=bin.dist[1]/5.5)
		plot.dark()
		plot(x.r[,1], x.r[,2], col=col2plot, pch=".", cex=bin.dist[1]/10, asp=T,xaxt="n", yaxt="n", ylab="", xlab="")
		mtext(gsub("empty_gene", "", gene.r), adj=0.1, col="red")
		mtext(gsub("empty_gene", "", gene.g), adj=0.5, col="green")
		mtext(gsub("empty_gene", "", gene.b), adj=0.9, col="blue")
		mtext(dataset, 1, adj=0)
		#gradient.legend(x.r, colours=c("black", "red"), scale=0.75)
	dev.off()
}



# checking for overlap with other genes:
for(dataset in datasets){
	gene = "EDNRA"
	threshold = 1
	dataset_meta = read.csv(paste(data_dir, dataset,"/plotdata.csv", sep=""))
	
        gene.values = read.csv(paste(data_dir, dataset,"/", gene, ".csv", sep=""), head=F)[,1]

#	Atrial_mean.exp = mean(gene.values[AtrialClusters==2])
#	Atrial_fraction_pos = round(sum(gene.values[AtrialClusters==2]>threshold)/sum(AtrialClusters==2)*100,1)

#	Ventricular_mean.exp = mean(gene.values[VentricularClusters==2])
#	Ventricular_fraction_pos = round(sum(gene.values[VentricularClusters==2]>threshold)/sum(VentricularClusters==2)*100,1)

#	 gene.r.cells = read.csv(paste(data_dir, "x13PCW_2/", gene.r, ".csv", sep=""), head=F)[,1]
#	  gene.g.cells = read.csv(paste(data_dir, "x13PCW_2/", gene.g, ".csv", sep=""), head=F)[,1]
#	  gene.b.cells = read.csv(paste(data_dir, "x13PCW_2/", gene.b, ".csv", sep=""), head=F)[,1]
#	  x13PCW_2_meta = read.csv(data_dir, "x13PCW_2/plotdata.csv")
	  
#	  gene.r.cells = gene.r.cells / max(gene.r.cells)
#	  gene.g.cells = gene.g.cells / max(gene.g.cells)
#	  gene.b.cells = gene.b.cells / max(gene.b.cells)
	   
#	  gene.r.cells[is.na(gene.r.cells)] = 0
#	  gene.g.cells[is.na(gene.g.cells)] = 0
#	  gene.b.cells[is.na(gene.b.cells)] = 0
	  
#	  col2plot = rgb(
#	    gene.r.cells,
#	    gene.g.cells,
#	    gene.b.cells
#	  )
  
  
	
	png(paste(figure_dir, gene, "_", dataset, ".png",sep=""), height=1200, width=1200)
	plot.dark()
	par("bg" = rgb(0.08, 0.08, 0.08))
	  plot(
	    -dataset_meta[,"global_y"],
	    dataset_meta[,"global_x"],
	    pch=16,
	    cex=0.5*dataset_meta[,"volume"]/20000,
	    col=plotcols(gene.values),
	    asp=T,
	    xlab=dataset,
	    ylab="",
	    bty="n",
	    xaxt="n",
	    yaxt="n"
	  )
	  mtext(dataset, adj=0, cex=4)
	# mtext(paste("expression: ", signif(Atrial_mean.exp,3), sep=""), adj=0.1, col="white")
	# mtext(paste("+ve fraction: ", signif(Atrial_fraction_pos,3), sep=""), adj=0.9, col="white")
	dev.off() 
	  
	bin.dist = c(
		diff(range(dataset_meta[,"global_y"]))/n_squares,
		diff(range(dataset_meta[,"global_y"]))/n_squares
	)
	
	x = collect.points(
	 -dataset_meta[,"global_y"],
	 dataset_meta[,"global_x"],
	 colours=alt.cols,
	 weights = gene.values,
	 method="derp",
	 bin.dist=bin.dist
	 )
	
	pdf(paste(figure_dir, gene, "_", dataset, ".pdf",sep=""), height=bin.dist[1]/5.5, width=bin.dist[1]/5.5)
		plot.dark()
        	plot(x[,1], x[,2], col=x[,4], pch=".", cex=bin.dist/10, asp=T, xaxt="n", yaxt="n", ylab="", xlab="")
        	mtext(gene, adj=1)
        	mtext(dataset, adj=0)
        	gradient.legend(gene.values, colours=alt.cols, scale=0.75)
	dev.off()
}



# checking for overlap with other genes:
for(dataset in datasets){
	gene = "ECE1"
	threshold = 1
	dataset_meta = read.csv(paste(data_dir, dataset,"/plotdata.csv", sep=""))
	
        gene.values = read.csv(paste(data_dir, dataset,"/", gene, ".csv", sep=""), head=F)[,1]

#	Atrial_mean.exp = mean(gene.values[AtrialClusters==2])
#	Atrial_fraction_pos = round(sum(gene.values[AtrialClusters==2]>threshold)/sum(AtrialClusters==2)*100,1)

#	Ventricular_mean.exp = mean(gene.values[VentricularClusters==2])
#	Ventricular_fraction_pos = round(sum(gene.values[VentricularClusters==2]>threshold)/sum(VentricularClusters==2)*100,1)

#	 gene.r.cells = read.csv(paste(data_dir, "x13PCW_2/", gene.r, ".csv", sep=""), head=F)[,1]
#	  gene.g.cells = read.csv(paste(data_dir, "x13PCW_2/", gene.g, ".csv", sep=""), head=F)[,1]
#	  gene.b.cells = read.csv(paste(data_dir, "x13PCW_2/", gene.b, ".csv", sep=""), head=F)[,1]
#	  x13PCW_2_meta = read.csv(data_dir, "x13PCW_2/plotdata.csv")
	  
#	  gene.r.cells = gene.r.cells / max(gene.r.cells)
#	  gene.g.cells = gene.g.cells / max(gene.g.cells)
#	  gene.b.cells = gene.b.cells / max(gene.b.cells)
	   
#	  gene.r.cells[is.na(gene.r.cells)] = 0
#	  gene.g.cells[is.na(gene.g.cells)] = 0
#	  gene.b.cells[is.na(gene.b.cells)] = 0
	  
#	  col2plot = rgb(
#	    gene.r.cells,
#	    gene.g.cells,
#	    gene.b.cells
#	  )
  
  
	
	png(paste(figure_dir, gene, "_", dataset, ".png",sep=""), height=1200, width=1200)
	plot.dark()
	par("bg" = rgb(0.08, 0.08, 0.08))
	  plot(
	    -dataset_meta[,"global_y"],
	    dataset_meta[,"global_x"],
	    pch=16,
	    cex=0.5*dataset_meta[,"volume"]/20000,
	    col=plotcols(gene.values),
	    asp=T,
	    xlab=dataset,
	    ylab="",
	    bty="n",
	    xaxt="n",
	    yaxt="n"
	  )
	  mtext(dataset, adj=0, cex=4)
	# mtext(paste("expression: ", signif(Atrial_mean.exp,3), sep=""), adj=0.1, col="white")
	# mtext(paste("+ve fraction: ", signif(Atrial_fraction_pos,3), sep=""), adj=0.9, col="white")
	dev.off() 
	  
	bin.dist = c(
		diff(range(dataset_meta[,"global_y"]))/n_squares,
		diff(range(dataset_meta[,"global_y"]))/n_squares
	)
	
	x = collect.points(
	 -dataset_meta[,"global_y"],
	 dataset_meta[,"global_x"],
	 colours=alt.cols,
	 weights = gene.values,
	 method="derp",
	 bin.dist=bin.dist
	 )
	
	pdf(paste(figure_dir, gene, "_", dataset, ".pdf",sep=""), height=bin.dist[1]/5.5, width=bin.dist[1]/5.5)
		plot.dark()
        	plot(x[,1], x[,2], col=x[,4], pch=".", cex=bin.dist/10, asp=T, xaxt="n", yaxt="n", ylab="", xlab="")
        	mtext(gene, adj=1)
        	mtext(dataset, adj=0)
        	gradient.legend(gene.values, colours=alt.cols, scale=0.75)
	dev.off()
}




# Other genes??
for(dataset in datasets){
# gene selection
	  gene.r = "HEY2"
	  gene.g = "PRRX1"
	  gene.b = "CNN1"
	 
	# gene prep
	  if(!gene.r %in% unique.genes){
	    gene.r = "empty_gene"
	  }
	  if(!gene.g %in% unique.genes){
	    gene.g = "empty_gene"
	  }
	  if(!gene.b %in% unique.genes){
	    gene.b = "empty_gene"
	  }   
	
	threshold = 1
	dataset_meta = read.csv(paste(data_dir, dataset,"/plotdata.csv", sep=""))
	
          # get values
	  gene.r.cells = read.csv(paste(data_dir, dataset,"/", gene.r, ".csv", sep=""), head=F)[,1]
	  gene.g.cells = read.csv(paste(data_dir, dataset,"/", gene.g, ".csv", sep=""), head=F)[,1]
	  gene.b.cells = read.csv(paste(data_dir, dataset,"/", gene.b, ".csv", sep=""), head=F)[,1]
	  # scale  
	  gene.r.cells = gene.r.cells / max(gene.r.cells)
	  gene.g.cells = gene.g.cells / max(gene.g.cells)
	  gene.b.cells = gene.b.cells / max(gene.b.cells)
	  # remove NA
	  gene.r.cells[is.na(gene.r.cells)] = 0
	  gene.g.cells[is.na(gene.g.cells)] = 0
	  gene.b.cells[is.na(gene.b.cells)] = 0
	  # create RGB space
	  col2plot = rgb(
	    gene.r.cells,
	    gene.g.cells,
	    gene.b.cells
	  )
	  
	png(paste(figure_dir, "red", gene.r, "_green", gene.g, "_blue", gene.b, "_", dataset, ".png", sep=""), width=500, height=500)
	  	plot.dark()
	  	par("bg" = rgb(0.08, 0.08, 0.08))
		  plot(
		    -dataset_meta[,"global_y"],
		    dataset_meta[,"global_x"],
		    pch=".",
		    cex=dataset_meta[,"volume"]/20000,
		    col=col2plot,
		    asp=T,
		    xlab=dataset,
		    ylab="",
		    bty="n",
		    xaxt="n",
		    yaxt="n"
		  )
	mtext(gsub("empty_gene", "", gene.r), adj=0, col="red")
	mtext(gsub("empty_gene", "", gene.g), adj=0.5, col="green")
	mtext(gsub("empty_gene", "", gene.b), adj=1, col="blue")
	dev.off()

	
	bin.dist = c(
		diff(range(dataset_meta[,"global_y"]))/n_squares,
		diff(range(dataset_meta[,"global_y"]))/n_squares
	)
	
	x.r = collect.points(
	 -dataset_meta[,"global_y"],
	 dataset_meta[,"global_x"],
	 colours=alt.cols,
	 weights = gene.r.cells,
	 method="derp",
	 bin.dist=bin.dist
	 )
	x.g = collect.points(
	 -dataset_meta[,"global_y"],
	 dataset_meta[,"global_x"],
	 colours=alt.cols,
	 weights = gene.g.cells,
	 method="derp",
	 bin.dist=bin.dist
	 )
	x.b = collect.points(
	 -dataset_meta[,"global_y"],
	 dataset_meta[,"global_x"],
	 colours=alt.cols,
	 weights = gene.b.cells,
	 method="derp",
	 bin.dist=bin.dist
	 )
	
	 col2plot = rgb(
	    x.r[,3],
	    x.g[,3],
	    x.b[,3]
	  )
	  
 	pdf(paste(figure_dir, "red", gene.r, "_green", gene.g, "_blue", gene.b, "_", dataset, "_CollectedDensity.pdf", sep=""), height=bin.dist[1]/5.5, width=bin.dist[1]/5.5)
		plot.dark()
		plot(x.r[,1], x.r[,2], col=col2plot, pch=".", cex=bin.dist[1]/10, asp=T,xaxt="n", yaxt="n", ylab="", xlab="")
		mtext(gsub("empty_gene", "", gene.r), adj=0.1, col="red")
		mtext(gsub("empty_gene", "", gene.g), adj=0.5, col="green")
		mtext(gsub("empty_gene", "", gene.b), adj=0.9, col="blue")
		mtext(dataset, 1, adj=0)
		#gradient.legend(x.r, colours=c("black", "red"), scale=0.75)
	dev.off()
}	














	
	
	
	
	
	 
  
