rm(list = ls(all = TRUE)) 
gc()
#install.packages("devtools")
library("devtools")
#install_github("klutometis/roxygen")
library("roxygen2")
#setwd("/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg")
#create("trap17")
working_dir <- "/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg/trap17-pkg"
setwd("/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg/trap17-pkg/trap17")
document()
setwd("..")
install("trap17")
library("trap17")
dirs <- set_dirs(working_dir = working_dir)
#data
trapdata <- process_data(filename = "TRAP17.csv",
                         dirs = dirs,
                         return_data = TRUE)
setwd("/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg/trap17-pkg/trap17")
#use_data(trapdata) 
document()
setwd("..")
install("trap17")
