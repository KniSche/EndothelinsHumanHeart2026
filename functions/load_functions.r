# No example plots! This is a generic command that I personally pass to my own plot functions 
#       - otherwise one may be inundated with multiple plots spawning upon function loading
        no.example.plx = T

# Vector of files to exclude from loading (such as this one)
        file.exceptions <- c("load_functions.r", "volcano_plot_manual")
        file.exceptions <- c(file.exceptions, list.files(func.dir)[grep("\\.r$", list.files(func.dir), invert=T)])
# default linux directory , at least for me...
#        func.dir <- "~/Documents/Rlib/functions"
# Begin functions
        load.functions <- function(directory){

        # locate names
                dir.ls <- list.files(directory)

        # Exclude functions and non-R. files
                functions.list <- dir.ls[!(dir.ls %in% file.exceptions) & 1:(length(dir.ls)) %in% grep("\\.r$", dir.ls)]

        # load files
                derpaderp <- list()
                for(i in functions.list){
                derpaderp[[i]] <- tryCatch(source(file.path(func.dir,i)), warning = function(w) {paste(i)}, error = function(e) {paste(i, "was not loaded properly, please check the code")})
                }
                 return(cat("The following functions contained errors and were not loaded properly, please check the code:\n\n", paste(names(which(sapply(functions.list, function(i){class(derpaderp[[i]])}) == "character")), collapse = "  "),"\n\n", sep=""))
        }

 load.functions(func.dir)
