# theResults <- list of processed files(each file is the average of different variations of the same food web) in Output folder
# size[3] of all abundances should be the same, which is different dispersal values (for the first test it is 10, then 100 in final runs)
# open the first file
load(theResults[1])
time.steps <- dim(abundances)[1]
dispersal.steps <- dim(abundances)[3]
analytics <- array(dim = c(time.steps, dispersal.steps, 4), dimnames = list(NULL, NULL, c("Stability.Local", "Synchrony.Local", "Stability.Regional", "Synchrony.Regional", )))

for (i in theResults){
   load(i)
   # add the name to the first column of all 4 tables, so S, C, P will be clear
   for (d in 1:dispersal.steps){
      # calculate regional stability
      # add it to Stability.Regional for corresponding column for d
      # calculate local stability
      # add it to Stability.Local
      # calculate Regional Synchrony
      # add it to Synchrony.Regional
      # calculate local synchrony
      # add it to Synchrony.Local
   }	
}
