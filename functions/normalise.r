
    normalise <- function(data, data2, full){
    #data is data to be normalised
    #data2 is data to be normalised to (example, you have points in [,1] and error bars in [,2], data = [,2] and data2 [,1] to normalise the error bars to the data points.
    #full is a logical saying to transpose to (0 - > 1) or just (x / max(x))
    if(missing(full)){full = T}
    if(missing(data2)){d2.check = F} else {d2.check=T}
     if(d2.check==F){

      if(full == T) {output = (data - min(as.numeric(na.omit(data))))/(max(as.numeric(na.omit(data)))-min(as.numeric(na.omit(data))))} else {output = data / max(as.numeric(na.omit(data)))}
    }
    if(d2.check==T){
     if(full == T) {output =  (data - min(as.numeric(na.omit(data2))))/(max(as.numeric(na.omit(data2)))-min(as.numeric(na.omit(data2))))} else {output =  data / max(as.numeric(na.omit(data2)))}
    }
      output
    }

  
