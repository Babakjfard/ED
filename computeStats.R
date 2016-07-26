# this file contains different statistical calculation on the output data resulting from the analyses


library(vegan)


# Returns a matrix to be used for calculating the regional(global) 
# properties of the species. showing whcih columns represent same species in different patches.
to.regional <- function(no.species, no.webs){
    b <- matrix(seq(from=1, to=no.species*no.webs), ncol=no.webs)
    return(b)
}

# calculates the global stability values for each species.
# reg.Map is calculated by toRegional function. matrix statistics value
# in regional scale vs. time
reg.stability <- function(abundances, reg.map){
    a <- abundances
    b <- reg.map
    no.species <- nrow(reg.map)
    
    
    reg.Mean <- sapply(1:no.species, function(i) rowMeans(a[,b[i,]]))
    #reg.Sd <- sapply(1:no.species, function(i) apply(a[,b[i,]],1,sd))
    temp.sd <- apply(reg.Mean, MARGIN = 2, sd)
    
    stability <- colMeans(reg.Mean)/temp.sd
    return(stability)
}

#Calculates the local stability values for each species
loc.stability <- function(abundances, reg.map){
    a <- abundances
    b <- reg.map
    loc.mean <- colMeans(a)
    loc.sd <- apply(a, MARGIN = 2, sd)
    indv.stab <- loc.mean/loc.sd
    no.species <- nrow(reg.map)
    
    loc.stab <- sapply(1:no.species, function(i) mean(indv.stab[b[i,]]))
    return(loc.stab)
}


# Calculates Kendall's W for synchrony, uses formulation for no ties
# from Book Biostatistical Analysis, Data are entered in a matrix, where
# each column represents a parameter and the number of rows are the numbers
# of data samples for each parameter (number of ranks)
Kendall.w.noties <- function (theData){
    M <- ncol(theData) #number of variables
    n <- nrow(theData) #number of rankers
    R <- rowSums(apply(birds, MARGIN = 2, rank))
    
    w <- (sum(R^2)-(sum(R))^2/n)/(M^2*(n^3-n)/12)
    return(w)
}

# calculates the intraspecific synchrony for a community of patches
# Input is a matrix, for that each row represents the abundances throughout
# the whole community (in one time step).
# Then each row is set into a matrix for which each column represents the whole 
# abundance in one site. so there will be no.web columns
synchrony.Intra <- function(temp.data, no.species, no.webs){
    d <- array(c(t(abund)), c(no.species,no.webs,nrow(abund)))
    w <- apply(d, MARGIN = 3, FUN = kendall.global)
    return(w)
}

# calculates the interspecific synchrony. uses the same logic as in regional stability- 12/15/2015
# array of size (time, no.web, no.species)
inter.sync <- function(abundances, reg.map){
    a <- abundances
    b <- reg.map
    no.species <- nrow(reg.map)
    
    w <- sapply(1:no.species, function(i) kendall.global(a[,b[i,]]))
    return()
}

#calculates the intraspecific synchrony
# array of size (time, no.species, no.web)
intra.sync <- function(abundances, reg.map){
    a <- abundances
    b <- reg.map
    no.webs <- ncol(reg.map)
    
    w <- sapply(1:no.webs, function(i) kendall.global(a[,b[,i]]))
}

## abundaned code------------------------------------------
#### Compute the stats
#compute.stats <- function(timeseries, d, no.webs, max.species){
#    d <- d
#    R <- as.matrix(timeseries[,(max.species*no.webs+2):ncol(timeseries)])
#    N <- as.matrix(timeseries[,2:(max.species*no.webs+1)])
#    
#    R.stab <- mean(rowMeans(R))/sd(rowMeans(R))
#    global.stab <- mean(rowMeans(N))/sd(rowMeans(N))
#    
#    loc.stab <- matrix(nrow=nrow(N),ncol=no.webs)    
#    for (count in 0:(no.webs-1)){
#        loc.stab[,count+1] <- rowSums(N[,(count*max.species+1):((count+1)*max.species)])        
#    }
#    local.stab <- mean(apply(loc.stab, 2, FUN = function(x) mean(x)/sd(x)))
#    
#    
#    return(c(d, local.stab, global.stab, R.stab))
#}
#
#local.stability <- function(abundance.Matrix, start){
#    abund <- abundance.Matrix[start:nrow(abundance.Matrix),] #take initial values out of statistical computations
#    #result <- colMeans(abund)/apply(abund,MARGIN = 2,FUN = sd)
#    theMeans <- colMeans(abund)
#    sds <- apply(abund, MARGIN = 2, FUN = sd)
#    return(list(means=theMeans, Sds=sds))
#}