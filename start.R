# The file to run on my local computer.
# Goes throught the different created food webs in \Input and runs to create stability vs. synchrony calculations
files <- list.files("/Volumes/Macintosh\ HD\ 2/Projects/ecological/ED/Input")
#files <- c("S5C05n_V1_P2.rda", "S5C05n_V1_P5.rda", "S5C05n_V11_P2.rda","S5C05n_V11_P5.rda")
for (the.names in files){
  commandArgs <- function() the.names
  source("runIt.R")
}
