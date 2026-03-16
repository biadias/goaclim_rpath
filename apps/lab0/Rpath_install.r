# 563 LAB Spring 2020 first setup for Rpath (Ecopath and Ecosim)
# Kerim.Aydin@noaa.gov

# STEP 1: Install/update R (Version 3.6.3 preferred) and needed R libraries
#
# I tend to use Rstudio (GUI for R), you can use Rstudio or plain R
#
# Packages to Install in R (start R or Rstudio and enter the following)
# (you should only have to do this once in the class)
#
install.packages("devtools")     #devtools is a big package, may take a while
install.packages("Rcpp")
install.packages("data.table")
install.packages("shiny")


# STEP 2 Install or update Rpath library from github using branch for the lab
# (will do this before each lab session to install new examples)
#
library(devtools)
install_github('NOAA-EDAB/Rpath', ref="lab_563_spring2020") 


# STEP 3 Running the Rpath code test
# (load these libraries before every Rpath script (or shiny script).
  library(Rpath)
  library(shiny)

# Load an example Ecosystem 
# Find example ecosystem files for Gulf of Alaska (GOA)
  GOA_base <- system.file("extdata", "GOA_tiny_base.csv", package="Rpath", mustWork=T)
  GOA_diet <- system.file("extdata", "GOA_tiny_diet.csv", package="Rpath", mustWork=T)
  GOA_ped  <- system.file("extdata", "GOA_tiny_ped.csv",  package="Rpath", mustWork=T)

# Read in an unbalanced Ecopath model from csv file inputs
  path.unbal <- read.rpath.params(GOA_base, GOA_diet, GOA_ped)
  
# Make sure everything is installed and the files are loaded correctly 
# with this command
  check.rpath.params(path.unbal)
# You should get the following output:  'Rpath parameter file is functional'

# Now check if shiny app works and example is installed
runRpathShiny("lab0")

# This should open a browser window using your default browser 
# (e.g. Chrome, Firefox, etc.) and the R terminal should read: 
# 'Listening on http://127.0.0.1:'
#
# The browser window should display a drop-down list of functional groups,
# and display two bar graphs (mortality and diet) that change when the dropdown
# is changed.

# To exit browser mode, close the browser, or click back to the R terminal
# and hit ESC key









