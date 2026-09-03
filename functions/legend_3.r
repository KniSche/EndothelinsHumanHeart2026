# plot legend on the y axis right hand side
     legend.3 <- function(...,xjust=0.5){
        legend(mean(par("usr")[1:2]), 1.4*par("usr")[4], xjust=xjust, xpd=NA, bty="n", ...)
      }   
      
