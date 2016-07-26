# 5/31/2016 the first instance to be staged into Git/Github
# Start the clock!
#require("RevoUtilsMath")
#setMKLthreads(4)
library(deSolve)

source("analysis.R")

args <- commandArgs()
#args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0){
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1){
  args[2] = "params.csv"
}


load(paste0("Input/",args[1]))
par <- read.csv(args[2], header=TRUE)
load("Init_B.rda")
print(args[1])

no.webs <- length(webs)
no.species <- unlist(lapply(webs, FUN = nrow))
max.species <- max(no.species)
init <- B[1:max.species, 1:no.webs]

#parameters of Yidsiz model
a<- matrix(nrow = max.species, ncol=max.species, par$a)  #attack rate
h <- matrix(nrow = max.species, ncol=max.species, par$h)    #handling time
d <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 0.001, .01, .1, .4, .5)    #dispersal rates
m <- rep(.3,max.species) # m <- runif(max.species)*.4            #mortality. same 

### The required parameters for mechanistic competition
basals <- matrix(nrow=max.species, ncol=no.webs,0)
ff <- matrix(nrow=max.species, ncol = no.webs, 0)
for (i in 1:no.webs){
    adj <- webs[[i]]
    bas <- which(colSums(adj)==0)
    basals[bas,i] <- 1  #vector of species with 1 for basals and 0 for others
    
    ff.this <- 1/colSums(adj)
    ff[1:length(ff.this),i] <- ff.this
}
ff[is.infinite(ff)] <- 0


# Creating neighbors maps. for now just using periodic boundaries.
# Can use the other two boundary conditions also, to compare the results
neighs <- find.neighs(no.webs)
dispers.maps <- dispersal.maps(neighs$neighs.per, no.species)

# Scenario 1: the same array in each food web represents the same species in all webs.
# using this way only one matrix for each parameter will be enough for covering all
# required data for the webs
parms <- list(webs=webs, a=a, h=h, basals=basals, no.species=no.species, ff=ff, lambda=par$lambda, alpha=par$alpha, beta=par$beta,
               neighs=neighs$neighs.per, dispers.maps=dispers.maps, R=par$R)

#------------------------------------------------------------------------------------------------------------

max.time =3000 
nreps <- 3
# putting all inputs into vector
initials <- c(init)


abundances <- array(0, dim=c(max.time+1, (max.species)*no.webs+1, length(d)))

for (i in 1:length(d)){
    
    parms$d <- d[i]
    temp.abund <- list()
    for (j in 1:nreps){
        temp.abund[[j]] <- ode(y=initials, func=solve.model, times=0:max.time, parms=parms)
         
    
    }
    abundances[,,i] <- Reduce('+', temp.abund)/length(temp.abund)
    
}
fileName <- args[1]
theFile <- paste0("Output/out_",fileName)
save(abundances, file = theFile) 
