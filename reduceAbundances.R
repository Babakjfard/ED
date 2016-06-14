#================== 6/7/2016
# to average different variations of one condition into one file
theFiles <- list.files(path = "./Input")

# put them into collections, so must end up having 5x5x4 = 100 different cases
# 5 no.species <- {5, 10, 20, 50, 100}, 5 connectances <-{0.05, 0.1, 0.18, 0.25, 0.3}, 4 patch.numbers <- {2,  5, 15, 50}
S <-  as.character(c(5))  # No. of species
C <- c("05", "1", "18", "25", "3")   # connectance values
P <- as.character(c(2, 5))   # number of patches
combs <- expand.grid(S, C, P)

#puts the same variations in the same group, so in the next step each group will be averaged into one file
groups <- lapply(1:dim(combs)[1], function(i) grep(paste0("S",combs[i,1],"C",combs[i,2],"n_V[[:digit:]]+_P",combs[i,3],".rda"), theFiles, value = TRUE))

# each collections contains arrays with equal sizes (not more than 30 arrays in each collection)
# the array of mean values for each collection is calculated!!# This file averages the abundances of different variations of the same scenario in 1 file
# so there will be (5species x 5 connectance x 4 patches)= 100 files

