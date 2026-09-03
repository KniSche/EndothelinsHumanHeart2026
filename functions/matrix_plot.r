# matrixplot function, basic form without breaks

black.divergent = c(
  "#ffec69",
  "#ffcd00",
  "#ff9100",
  "#000000",
  "#005fff",
  "#0090ff",
  "#76faff"
  )

# First we need this to auto determine the rounding space for the matrix.

  decimalplaces_vec <- function(x) {

    vector <- c()
    for (i in 1:length(x)){

      if(!is.na(as.numeric(x[i]))){

        if ((as.numeric(x[i]) %% 1) != 0) {
          vector <- c(vector, nchar(strsplit(sub('0+$', '', as.character(x[i])), ".", fixed=TRUE)[[1]][[2]]))


        }else{
          vector <- c(vector, 0)
        }
      }else{
        vector <- c(vector, NA)
      }
    }
    return(max(vector))
  }
    
    
    
    
   
 # matrix2plot = t(  t(  as.matrix(foetal.merged.cm$table)  )  /  (colSums(as.matrix(foetal.merged.cm$table))+1))
 
 
 
matrix.plot <- function(matrix2plot, col.scale=Heat, cex.multiplier=3, y.axis=T, x.axis=T, pch=15, cex.axis=1,  round.n = "default", generate.scale=F, aspect=T, bg.col = par("bg"), ...){ 
  # Identify colors for the matrix - very basic
  # Range of colours (equal on either side of zero)
  # First, find how many decimal places to round at
  if(round.n == "default"){
   round.n = round(mean(sapply(sample(1:length(matrix2plot), round(0.01*length(matrix2plot))), function(i){decimalplaces_vec(matrix2plot[i])})))
    if(round.n > 4){round.n = 4}
  }
  # if both tails exist
    matrix.values = matrix2plot[!is.na(matrix2plot)]
    range.col = 0
    if(1 %in% sign(range(matrix.values)) & -1 %in% sign(range(matrix.values))){
      range.col=c(
        seq(-max(abs(matrix.values)),0,by=1*10^-round.n),
        seq(0, max(abs(matrix.values)),by=1*10^-round.n)
      )
    }
   
  # if only a positive tail exists...
    if(1 %in% sign(range(matrix.values)) & !-1 %in% sign(range(matrix.values))){
      range.col=seq(0, max(matrix.values),by=1*10^-round.n)
    }
    
  # if only a negative tail exists...
    if(!1 %in% sign(range(matrix.values)) & -1 %in% sign(range(matrix.values))){
      range.col=seq(-max(matrix.values),0,by=1*10^-round.n)
    }
# Colour mapping matrix     
  # Dataframe for colours and values
    col.df= data.frame(
      colorRampPalette(col.scale)(length(range.col)),
      range.col
      )
      
     if(generate.scale==T){
      return(col.df)
     }
  # Match the rounded value of the matrix2plot to the colour gradient 
  #  col.df[match(round(matrix2plot,round.n), col.df[,2]),1]




#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  
# Colour matrix
    if(class(col.scale)[1] == "matrix"){matrix2plot.col = col.scale} else {
    # Create new matrix of the colours for our function 
      matrix2plot.col = matrix(
        col.df[match(as.character(round(matrix2plot,round.n-1)), as.character(round(col.df[,2],round.n-1))),1],
    #    col.df[match(round(matrix2plot), col.df[,2]),1], 
        dim(matrix2plot)[1],
        dim(matrix2plot)[2]
      )
   }


  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  # 
  # plotting
  # Identify cex values (pt size)                    
    matrix2plot.cex = matrix(1, dim(matrix2plot)[1],dim(matrix2plot)[2])

  # make initial plot
    plot(c(1, dim(matrix2plot)[2])+2, c(1, dim(matrix2plot)[1]), pch=NA, ... , xaxt="n", yaxt="n", xaxs="i", yaxs="i", ylim=c(0,dim(matrix2plot)[1]+1), xlim=c(0,dim(matrix2plot)[2]+1), bty="n", asp=aspect)
    rect(-200,-200, 100+dim(matrix2plot)[2], 100+dim(matrix2plot)[1], col=bg.col, lty=0, xpd=NA)
    
  # plot points from matrix.
     points(rep(1:dim(matrix2plot)[2], each=dim(matrix2plot)[1]),  rep(dim(matrix2plot)[1]:1,dim(matrix2plot)[2]),
      cex = matrix2plot.cex*cex.multiplier,
      col = matrix2plot.col,
      pch=pch
     )

#  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  #  
# plot axes and labels
if(y.axis==T){
  axis(2,
      at=dim(matrix2plot)[1]:1,
      labels=rownames(matrix2plot),
      las=2,
      cex.axis=cex.axis
    ) 
}
if(x.axis==T){    
  axis(3,
      at=1:dim(matrix2plot)[2],
      labels=colnames(matrix2plot),
      las=2,
      cex.axis=cex.axis
    ) 
}
# end?
}



## Run function (custom matrix plot function (see library))
#  matrix.plot(matrix2plot, Heat, ylab="", xlab="")

## Calculate per-class mean recall and precision (ignoring NAs)  
#  mean.recall = sapply(1:12, function(i){mean(recall.df[i,6:12][!is.na(recall.df[i,6:12])])})
#  mean.precision = sapply(1:12, function(i){mean(precision.df[i,6:12][!is.na(precision.df[i,6:12])])})
#   
## Plot the mean recall and precision under column and rows respectively
#  points(1:12, rep(-1, 12), col=colorRampPalette(alt.cols)(101)[round(100*mean.recall)+1], pch=15, cex=3, xpd=NA)
#    axis(2, at=-1, label="Recall", xpd=NA, las=2)
#  points(rep(14, 12), 1:12, col=colorRampPalette(alt.cols)(101)[round(100*mean.precision)+1], pch=15, cex=3, xpd=NA)
#    axis(3, at=14, label="Precision", xpd=NA, las=2)

## legend
#  mtext("Predicted", 2, padj=-10)
#  mtext("Reference", 3, padj=-10)

#  gradient.legend(seq(0,1,0.1), colours=alt.cols, distance=2, scale=0.8)

## Two matrices
## one of numerical value or test results 
## 
#  mem = mem
#  mem.p = mem.col
## Three matrices - one main split into colour and cex values  
#  matrix2plot = mem
#  
#  
## Then, devise splits for matrix - if no splits just plot matrix - TO DO

#  plot.matrix = function(matrix2plot, col = c(alt.cols)){
#    matrix2plot = mem
#    
#  
#  
#  }
#  
#    
#  matrix2plot.col = mem.col
#  matrix2plot.cex = mem.cex
#  
#  matrix2plot.sub = matrix2plot[9:13,,drop=F]
#  matrix2plot.col.sub = matrix2plot.col[9:13,,drop=F]
#  matrix2plot.cex.sub = matrix2plot.cex[9:13,,drop=F]
#  
#  points(rep(1:dim(matrix2plot.sub)[2], each=dim(matrix2plot.sub)[1]),  rep((meta.ylim-10):(meta.ylim-14),dim(matrix2plot.sub)[2]),
#    cex = matrix2plot.cex.sub,
#    col = matrix2plot.col.sub,
#    pch=15
#   )
#   
#  axis(2,
#    at=c(meta.ylim, y2p)[9:13],
#    labels=rownames(matrix2plot)[9:13],
#    las=2,
#    cex.axis=0.4
#  ) 

#  
#  plot.type="dot"
#  
#  
#  
#  
#  
#  
#  
#  if(plot.type=="dot"){
#  
#  
#  par(oma=c(10,7,2,6))
#    plot(c(1, dim(mem)[2]), c(1, dim(mem)[1]), pch=NA, xlab="", ylab="cluster", xaxt="n", yaxt="n", xaxs="i", yaxs="i", ylim=c(0,dim(mem)[1]+1), xlim=c(0,dim(mem)[2]+1), bty="n")
#     points(rep(1:dim(mem)[2], each=dim(mem)[1]),  rep(dim(mem)[1]:1,dim(mem)[2]),
#      cex = mem.cex,
#      col = mem.col,
#      pch=16
#     )
#    
#    axis(2, at=dim(mem)[1]:1, label = rownames(mem), las = 2, cex.axis=0.6)
#    axis(4, at=dim(mem)[1]:1, label = rownames(mem), las = 2, cex.axis=0.6)
#    axis(1, at=1:dim(mem)[2], label = colnames(mem), las = 2, cex.axis=0.6)
#    axis(3, at=1:dim(mem)[2], label = colnames(mem), las = 2, cex.axis=0.6)
#    cluster.labels = clus.man.lab$man.lab
#    for(k in levels(Idents(all.subsample.integrated))[cid.cluster$order]){
#      axis(1, at=range(which(featureset.df[,2] == k)), tcl=0.2, labels=c("", ""), line=3.7)
#      axis(1, at=mean(range(which(featureset.df[,2] == k))), tcl=-0.5, labels=paste(clus.man.lab$man.lab[as.numeric(k)+1], "-", k), line=3.7, las=2)
#    }
#    mtext("Unique top 5 cluster gene markers", 3, padj=-8)
#   
#  gradient.legend(c(0,1), colours=alt.cols, scale=0.7, cex=0.7, label = "Gene-scaled mean log2 expression", cex.axis=0.6, distance=0.25)
#  legend(dim(mem)[2]*1.05,dim(mem)[2]*0.28, 
#    legend=c(0.1, 0.33, 0.75, 1), pch=16, pt.cex=c(0.1, 0.33, 0.75, 1)+0.2, xpd=NA, bty="n",
#    ncol=1, cex=0.6, title = "Fraction in cluster"
#    )
