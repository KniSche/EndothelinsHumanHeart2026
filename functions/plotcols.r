# Simple numeric scale to colour grade...


plotcols = function(data, colours = alt.cols, zero=F){
    x = as.numeric(cut(data,101))
    if(sum(x==51) == length(x) & zero==T){
      x[x==51] = 1
    }
    colorRampPalette(color=colours)(101)[x]
  }
  
