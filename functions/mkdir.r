# This is a really cool function that I made, it creates directories recursively in R when you give it a filepath
# Of note: it will NOT overwrite a directory that already exists as checked in line 8. If you're concerned do a few tests.
# The overwrite=F is not included yet.

	mkdir <- function(filepath, ..., overwrite=F){
			directories <- unlist(strsplit(c(filepath, ...),"/"))
				for(i in 1:length(directories)){
					filepath.temp <- paste(directories[1:i],collapse="/")
		            ifelse(filepath.temp != "" & !dir.exists(file.path(filepath.temp)), dir.create(file.path(filepath.temp)),FALSE)
		    	}
		}



# usage:

# simple
#	mkdir("path/to/my/folder")
#	mkdir(file.path("path","to","my","folder"))
#
# embed string object 
#	my.temporary.path = file.path("path","to","my","folder")
# include multiple objects using the ... argument.
#	mkdir(my.temporary.path,"additional/path")
