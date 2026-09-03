# plot legend on the y axis right hand side
     legend.4 <- function(...,yjust=0.5){
        legend(par("usr")[2], mean(par("usr")[3:4]), yjust=yjust, xpd=NA, ...)
      }   
