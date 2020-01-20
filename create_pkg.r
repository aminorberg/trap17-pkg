rm(list = ls(all = TRUE)) 
gc()
#install.packages("devtools")
#install_github("klutometis/roxygen")

library("devtools")
library("roxygen2")

root_dir <- "/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg"
working_dir <- "/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg/trap17-pkg"
pkg_dir <- "/Users/anorberg/Documents/Zurich/UZH/TRAP/pkg/trap17-pkg/trap17"

### create the package
?create 
#setwd(root_dir)
#create("trap17")


### document the package
?document
setwd(pkg_dir)
document()

### install and load
setwd(working_dir)
install("trap17")

library("trap17")

### add data
"#4C4C4C"#trapdata <- process_data(filename = "TRAP17.csv",
#                         dirs = dirs,
#                         return_data = TRUE)
#
#setwd(pkg_dir)
#use_data(trapdata) 
#document()
#setwd(working_dir)
#install("trap17")

