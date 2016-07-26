#form matrix of predictors
x <- cbind(punt$R_Strength, punt$L_Strength, punt$R_Flexibility, punt$L_Flexibility, punt$O_Strength
           )
cor(x)

#Compute multicolinearity
R2.table <- matrix(NA, 5, 2) # to hold the results
dimnames(R2.table) <- list(c("R_Strength", "L_Strength", "R_Flexibility", "L_Flexibility","O_Strength"), c("R2", "VIF"))

n <- nrow(x)

for (j in 1:5){
    tmp <- lsfit(x[,-j], x[,j])
    SSTot <- (n-1)*var(x[,j])
    SSErr <- sum(tmp$res^2)
    R2 <- 1 - SSErr / SSTot
    VIF <- 1 / (1-R2)
    R2.table[j, 1] <- R2
    R2.table[j,2] <- VIF
}

#====================================================
# 11/24/2015
# calculating kendall W for intraspecific synchrony:
# from page 444 on the book Biostatistical analysis
birds <- c(10.4, 10.8, 11.1, 10.2, 10.3, 10.2, 10.7, 10.5, 10.8, 11.2, 10.6, 11.4, 7.4, 7.6, 7.9, 7.2, 7.4, 7.1, 7.4, 7.2, 7.8, 7.7, 7.8, 8.3, 17, 17, 20, 14.5, 15.5, 13, 19.5, 16, 21, 20, 18, 22)
birdsmat <- matrix(birds, ncol=3)

#sample of creating an array
(a <- matrix(8, 2, 3))
(b <- matrix(9, 2, 3))
array(data = c(a,b), dim = c(2,3,4))

#calculate the global abundance of species, to be used for interpecific synchrony
samps <- sapply(1:5, function(i) rowSums(samp[,reg.map[i,]]))
#=================================================
#3/21/2016
webs <- vector("list", 1)
webs[[1]] <- adj
max.species <- dim(adj)[1]
N <- y
a <- matrix(nrow = max.species, ncol=no.webs, parms$a)

h <- matrix(rep(parms$h,max.species^2), ncol=max.species)
basals <- parms$basals


parms <- list(webs=webs, a=a, b=b, h=h, m=m, basals=basals, mu=mu, Y=Y, K=K, no.species=no.species, ff=ff, lambda=.65, alpha=.3, beta=.5,
              flow=flow,S=S, neighs=neighs$neighs.per, dispers.maps=dispers.maps, R= 500)

#=========================
abunds <- list()
for (i in 1:10){
    load(paste0("P5S20_results_difAdj",i,".rda"))
    abunds[[i]] <- abundances
}
#===========================
p <- 2
S <- 5
theFile <- "07_11_abundances.pdf"
pdf(file = theFile)

for (k in 1:dim(abunds)[3]){
 par(mfrow=c(2,2), oma=c(0,0,2,0))
    toPlot <- abunds[,,k]
   
    for (j in 1:p){
    
    matplot(toPlot[,1], toPlot[,((j-1)*S+2):(j*S+1)], t="l", lty=1, xlab="time", ylab = "Biomass")
    
    }
        title(main=paste0("model",k, "\n h=",hs[k]), outer = TRUE)
    
}
dev.off()
#=================================4/6/2016
for (i in 1:10){
    load(paste0("P5S20_results_difAdj",k,".rda"))
    NS[,,i]<-N
}
par(mfrow=c(2,5))
for (i in 1:10){
    barplot(1:20 , NS[,1:5,i], t="l", lty=1)
}
#============================
max.time <- 1000
d <- c(0)
h <- seq(.1, 1, length.out = 100)
for (k in 1:10){
    
    initials <- c(N)
    abundances <- array(0, dim=c(max.time, length(N)+1, length(d)))
    for (i in 1:length(h)){
        
        temp.abund <- list()
        for (j in 1:nreps){
            parms$h <- h[i]
            temp.abund[[j]]<- ode(y = initials, parms = parms, times = seq(from = 0, to = 500, len = 1000),
                                  method = "ode45", func = solve.model)
            
        }
        
        abundances[,,i] <- (temp.abund[[1]]+temp.abund[[2]]+temp.abund[[3]])/3
        
    }
    fileName <- paste0("P",ncol(N),"S", nrow(N),"_variable_h",k,".rda")
    save(d,abundances, file = fileName )
}
#=============================================================source("analyses3_30.R")
max.time = 1000
start.time = 150
nreps <- 3
patches <- 5
# putting all inputs into vector

#r <- seq(from = 1e-7,to=.1, length.out = 100)
d<- c(0,1e-7,.1,.5)
# required statistics

for (k in 1:10){
    
     initials <- c(NS[,,k])
    abundances <- array(0, dim=c(max.time, length(N)+1, length(d)))
    for (i in 1:length(d)){
        
        temp.abund <- list()
        for (j in 1:nreps){
            parms$d <- d[i]
            temp.abund[[j]]<- ode(y = initials, parms = parms, times = seq(from = 0, to = 500, len = 1000),
                                  method = "ode45", func = solve.model)
            
        }
        
        abundances[,,i] <- (temp.abund[[1]]+temp.abund[[2]]+temp.abund[[3]])/3
        
    }
    fileName <- paste0("P",ncol(N),"S", nrow(N),"_results_difBsameAdj",k,".rda")
    save(d,abundances, file = fileName )
}
#===================================================================
# computing connectedness
for (i in 1:dim(adj)[3]){
    adjs <- adj[,,i]
    l <- dim(adj)[1]
    (con[i] <- sum(adjs)/(l*(l-1)))
}
#===================================================================
# computing connectedness
bas <- vector()
for (i in 1:dim(adj)[3]){
    b <- length(which(colSums(adj[,,i])==0))
    
    (bas[i] <- b/20)
}
#===================================================================
# to keep adj and B-init constant and change h
# needs initials and parms be in the memory to get used.
# hs <- c(...)
patches <- 5
species <- 20
max.time <- 1000
nreps <- 3
hs <- c(.1,.15,.17,.2,.25,.27,.3,.35,.37,.4,.45,.47,.5,.55,.57,.6,.65,.67,.7,.8,1,10,20,30)

initials <- N


for (adjs in dim(adj)[3]){
    abunds <- array(dim = c(max.time, patches*species+1, length(hs)))
    for(webs in 1:5){parms$webs[[webs]] <- adj[,,adjs]}
    
for (i in 1:length(hs)){
    parms$h <- hs[i]
    
    temp.abund <- list()
    for (j in 1:nreps){
        parms$d <- 0
        temp.abund[[j]]<- ode(y = initials, parms = parms, times = seq(from = 0, to = 500, len = 1000),
                              method = "ode45", func = solve.model)
        
    }
    abunds[,,i] <- (temp.abund[[1]]+temp.abund[[2]]+temp.abund[[3]])/3
    
}

fileName <- paste0("dif_H-P",ncol(N),"S", nrow(N),"_",adjs,".rda")
save(N,hs, parms, abunds, file = fileName )
}
#====================
source("analyses3_30.R")
set.seed(sample(1:100,1))
repeat{
  
  # repeat{
  (samp <- niche(species, C[i]))
  samp[diag(samp)] <-0     #no cannibalism
  (samp <- remove.disconnected(samp))
  if (nrow(as.matrix(samp))==species) break
}
samp
#===================

#calculating sum of each of the adjacency matrices. each column represents a connectedness.
(links <- sapply(1:length(adj.variations) , function(i) apply(X = adj.variations[[i]], MARGIN = 3, FUN = sum)))
(t <- !(apply(links, MARGIN = 2 , FUN = duplicated)))

load("adjacencies3.rda")
colSums(unique.adjs)

a <- matrix(-999, nrow = 30, ncol=30)
for (i in 1:30){
  for (j in 1:i) {
    a[j,i] <- sum(abs(test[,,j]-test[,,i]))
    
  }
}
#================== 5/10/2016
# to put each (unique) adjacency matrix in one file!!
setwd("/Volumes/Macintosh HD 2/ecological/codes/data")
thefiles <- c("niche_adjs10.rda", "niche_adjs20.rda", "niche_adjs50.rda")
for (theAdjs in thefiles) {
  load(theAdjs)
  thenames <- names(adj.variations)
  test1 <- substr(thenames[5],3,4)
  test2 <- substr(thenames[5],8,11)
  print(test1)
  print(test2)
}

#================== 6/15/2016
#for statistics of the model
for (i in 1:length(final.abundances)){
  abundance <- final.abundances[[i]]
  thename <- names(final.abundances)[i]
  no.species <- strtoi(substr(thename, start = regexpr("S", thename)[1]+1, stop = regexpr("C", thename)[1]-1))
  connectance <- substr(thename, start = regexpr("C", thename)[1]+1, stop = regexpr("n", thename)[1]-1)
  no.webs <- strtoi(substr(thename, start = regexpr("P", thename)[1]+1, stop = regexpr("r", thename)[1]-2))

  reg.map <- to.regional(no.species, no.webs)
  reg.stab <- sapply(1:dim(abundance)[3], function(i) reg.stability(abundance[,,i], reg.map))
  loc.stab <- sapply(1:dim(abundance)[3], function(i) loc.stability(abundance[,,i], reg.map))

  #######!!!!!!!!!!!!!!!!!!The following two lines need improvements
  intra.synch <-  sapply(1:dim(abundance)[3], function(i) intra.sync(abundance[,,i], reg.map)) 
  inter.synch <-  sapply(1:dim(abundance)[3], function(i) inter.sync(abundance[,,1], reg.map))
}

#=======================7/12
f <- list.files("Output")
theFile <- "07_19_Reg_stabiilities.pdf"
pdf(file = theFile)
for (fl in f){
  load(paste0("Output/",fl))
  specs <- model.spec.by.name(fl)
  reg.map <- to.regional(no.species = strtoi(specs[1]), no.webs = strtoi(specs[2]))
  reg.stab <- sapply(1:dim(abundances)[3], function(i) reg.stability(abundances = abundances[2500:3001,2:dim(abundances)[2],i],reg.map = reg.map))
  reg.stab <- t(reg.stab)
  d <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 0.001, .01, .1, .4, .5)    #dispersal rates
  matplot(d, reg.stab, t="l", lty=1, cex = d, xlab = "d", ylab="Regional Stability")
  title(main = paste0("Species=",specs[1], ", Patches=", specs[2], " Connectance=",specs[3]))
}
dev.off()

