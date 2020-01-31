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
#setwd(root_dir)
#create("trap17")


### document the package
setwd(pkg_dir)
document()

### install and load
setwd(working_dir)
install("trap17")

library("trap17")

### add data
#dirs <- set_dirs(working_dir = working_dir, raw_data = TRUE)
#trapdata <- process_data(filename = "TRAP17.csv",
#                  dirs = dirs,
#                  save_data = TRUE,
#                  return_data = TRUE)
#setwd(pkg_dir)
#use_data(trapdata, overwrite = TRUE) 
#document()
#setwd(working_dir)
#install("trap17")

### unit tests
#setwd(pkg_dir)
#usethis::use_testthat()
#usethis::use_test("basic-tests")
