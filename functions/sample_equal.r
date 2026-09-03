sample.equal <- function(vector, classes, n, override = F){
  n2sub <-  min(table(classes))
  if(n < n2sub){a = n} else {a = n2sub}
  
  equal.sample.vector <- unlist(lapply(levels(classes), function(k) {
      if(override==T){sample(vector[classes==k], min(c(length(vector[classes==k]), n)), replace=F)} else {
        sample(vector[classes==k], a, replace=F)
      }
    }))
  
  return(equal.sample.vector)  
}

# & length(vector[classes==k]) >= n

#sample.equal <- function(vector, classes, n, override = F){
#  n2sub <-  min(table(classes))
#  if(n < n2sub){a = n} else {a = n2sub}
#  if(override==T){a = min(n)}
#  sub.vector <- unlist(lapply(levels(classes), function(k) {
#    sample(which(classes==k), a, replace=F)
#    }))
#    
#  equal.sample.vector <- vector[sub.vector]
#  
#  return(equal.sample.vector)  
#}


