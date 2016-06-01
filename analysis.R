#----------------------- food web models
library(deSolve)
# ---- niche-model ----
# Niche food web model
niche <- function (S, C) {
  
  adj <- matrix(nrow=S, ncol=S, 0)
  niches <- runif(S) # Random niches n
  radius <- rbeta(S, 1, 1/(2*C)-1) * niches # Random radii r
  # Feed between [r/2, n)
  center <- runif(S, min = radius/2, max = niches)  
  for (i in 1:S) {
    consumed <- which(niches > (center[i] - radius[i]/2) & 
                        niches < (center[i] + radius[i]/2))
    adj[consumed, i] <- 1
  }
  return (adj)
}

# ---- random-model ----
# Random food web model
random <- function (S, C) {
    adj <- matrix(nrow=S, ncol=S, 0)
    rand <- matrix(nrow=S, ncol=S, runif(S*S))
    adj[rand < C] <- 1
    return (adj)
}

# ---- cascade-model ----
# Cascade food web model
cascade <- function (S, C = NULL) {
    adj <- matrix(nrow=S, ncol=S, 0)
    niches <- runif(S) # Random niches n
    if (is.null(C)) {
        prob <- 1
    }
    else {
        prob <- 2*C*S/(S-1)
    }
    for (i in 1:S) {
        consumed <- which(niches < niches[i] & runif(1) < prob)
        adj[consumed, i] <- 1
    }
    return (adj)
}


# no.speciesx no.patches
#This file contains the functions which are used for analysing the model(s)

extract.Patches <- function(abundances, max.species, no.webs){
    patches <- array(0, dim=c(nrow(abundances), max.species, no.webs))
    
    for (i in 0:(no.webs-1)){
        patches[,,i+1] <- abundances[, ((i*max.species)+1):((i+1)*max.species)]
        
    }
    return(patches)    
}

#Code taken from lab 4 Solutions
remove.disconnected <- function(adj) {
    # Q: has a bug. if a maatrix which disconnected removed before, it returns empty matrix!!
    adj <- as.matrix(adj)
    # Remove species that neither consume nor are consumed
    empty <- which(rowSums(adj) == 0 & colSums(adj) == 0)
    if (length(empty)>0){
        adj <- adj[-empty,-empty]
    }
     
    return(adj)
}
#----------------------------------------------------------

# Code from lab-5 Solutions
find.neighs <- function(npops) {
    neighs.abs <- matrix(nrow = npops, ncol = 3)
    neighs.abs[, 1] <- 1:npops
    colnames(neighs.abs) <- c("pop.number", "neighb.left", "neighb.right")
    neighs.abs[2:(npops - 1), 2:3] <- c((2:(npops - 1)) - 1, (2:(npops -
                                                                     1)) + 1)
    neighs.per <- neighs.abs
    neighs.ref <- neighs.abs
    # Periodic boundaries
    neighs.per[1, 2:3] <- c(npops, 2)
    neighs.per[npops, 2:3] <- c(npops - 1, 1)
    # Reflecting boundaries
    neighs.ref[1, 2:3] <- c(1, 2)
    neighs.ref[npops, 2:3] <- c(npops-1, npops)
    # Absorbing boundaries
    neighs.abs[1, 2:3] <- c(NA, 2)
    neighs.abs[npops, 2:3] <- c(npops-1, NA)
    return(list(neighs.abs = neighs.abs, neighs.per = neighs.per, neighs.ref = neighs.ref))
}
#----------------------------------------------------------------

# calculates the three required maps for immigration and emigration 
# among arbitrary number of food webs with different number of species.
# works for the scenarios of each food web having one left and one right in
# connection (with three boundary conditions)
# neighbors.map: the map of connections of food web. It has nrow=no.wes,ncol=3
# no.species: a vector representing number of species in each web
dispersal.maps <- function(neighbors.map, no.species){
    left.immig <- matrix(nrow = max(no.species), ncol=nrow(neighbors.map), 0)
    right.immig <- matrix(nrow = max(no.species), ncol=nrow(neighbors.map), 0)
    emigration <- matrix(nrow = max(no.species), ncol=nrow(neighbors.map), 0)
    
    for (i in 1:length(no.species)){
        left.immig[1:no.species[neighbors.map[i,2]],i] <- 0.5
        right.immig[1:no.species[neighbors.map[i,3]],i] <- 0.5
        s <- rep(1, max(no.species))
        step <- ifelse(no.species[neighbors.map[i,2]]==no.species[neighbors.map[i,3]], 0, 1)
        s[(min(no.species[neighbors.map[i,2]], no.species[neighbors.map[i,3]])+step):length(s)]<-2
        if (no.species[neighbors.map[i,2]]==no.species[neighbors.map[i,3]]){
            s[no.species[neighbors.map[i,2]]]=1
        }
        left.immig[,i] <- left.immig[,i]*s
        right.immig[,i] <- right.immig[,i] *s
        #emigration[1:max(no.species[neighbors.map[i,2]], no.species[neighbors.map[i,3]]),i] <-1
        emigration <- left.immig+right.immig
    }
    return(list(left.immig=left.immig, right.immig=right.immig, emigration=emigration))
}
#--------------------------------------------------------------------------------

# the ode function to be solved over the webs. 
# Later? develop also a model for asynchronic models which can create a 
# whole new effects and research issue.
solve.model <- function(t, y, parms){
    
    y[y<0] <-0
    with(parms,{
        # return from vector form into matrix form for calculations
        #(R <- as.matrix(y[(max(no.species)*length(no.species)+1):length(y)]))
        (N <- matrix(y[1:(max(no.species)*length(no.species))], ncol=length(no.species)))
        
        ddB <- matrix(nrow=max(no.species), ncol=length(no.species))
        #dy2 <- matrix(nrow=length(no.species), ncol=1)
        
        
        no.webs <- length(no.species)
        for (i in 1:no.webs){
            #B <- N[1:no.species[i],i]
            B <- y[((i-1)*max(no.species)+1):(i*max(no.species) )]
            adj <- webs[[i]]
            
            a.this <- a[1:no.species[i],i]  #rate of attack for a predator is the same for all preys
            ff.this <- ff[1:no.species[i],i]
            h.this <- h[1:no.species[i],i]
            basals.this <- basals[1:no.species[i],i]
            
            ff.this <- ff[1:no.species[i],i]
            h.this <- h[1:no.species[i], i]
            
            # Consumers
            r <- B%*%adj
            ressource <- lambda*a.this*ff.this*B*r/(1+a.this*h.this*ff.this*r)
            
            ff.pred <- adj%*%ff.this
            p <- adj%*%(ff.this*B)
            
            preddation <- a.this*p*B/(1+a.this*h.this*ff.pred*B)
            no.basals <- rep(1,no.species[i])-basals.this
            
            ddB[,i] <- (lambda*a.this*R/(1+a.this*h.this*R))*basals.this*B+t(ressource)-preddation-alpha*B*no.basals-beta*B^2
            
        }
        
        #Calculating dispersals .they can be easily replaced
        # by adjacency maps of connections between food webs arbitrarily!
        disp.left <- N*d*dispers.maps$left.immig
        disp.left <- disp.left[,neighs[,2]]
        
        disp.right <- N*d*dispers.maps$right.immig
        disp.right <- disp.right[,neighs[,3]]
        
        emig <- N*d*dispers.maps$emigration
        mortality <- m*N
        
        
        dy1 <- ddB+disp.left+disp.right-emig
        #dy1 <- dy1+disp.left+disp.right-emig-mortality
        return(list(c(dy1)))
    })
}


