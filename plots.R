# 20/7/2016
source("computeStats.R")

model.spec.by.name <- function(thename){
  no.species <- strtoi(substr(thename, start = regexpr("S", thename)[1]+1, stop = regexpr("C", thename)[1]-1))
  connectance <- substr(thename, start = regexpr("C", thename)[1]+1, stop = regexpr("n", thename)[1]-1)
  no.webs <- strtoi(substr(thename, start = regexpr("P", thename)[1]+1, stop = regexpr("r", thename)[1]-2))
  v <- substr(thename, start = regexpr("V", thename)[1]+1, stop=regexpr("V", thename)[1]+2)
  return(c(no.species = no.species, no.webs = no.webs, connectance=connectance, v=v))
}


#To plot regional stability vs. d diagrams
plot.reg.stab <- function(folder, output.file="reg_Stab_plots2.pdf"){
  f <- list.files(folder)
  theFile <- output.file
  pdf(file = theFile)
  for (fl in f){
    print(fl)
    load(paste0(folder,"/",fl))
    specs <- model.spec.by.name(fl)
    reg.map <- to.regional(no.species = strtoi(specs[1]), no.webs = strtoi(specs[2]))
    reg.stab <- sapply(1:dim(abundances)[3], function(i) reg.stability(abundances = abundances[2500:3001,2:dim(abundances)[2],i],reg.map = reg.map))
    reg.stab <- t(reg.stab)
    d <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 0.001, .01, .02, .03, .04)    #dispersal rates
    matplot(d, reg.stab, t="l", lty=1, cex = d, xlab = "d", ylab="Regional Stability")
    title(main = paste0("Species=",specs[1], ", Patches=", specs[2], ", Connectance=",specs[3], ", V",specs[4]))
  }
  dev.off()
  
}

#To plot local stability vs. d diagrams (exactly the same code as for reg.stab. to be merged into one function!!!!!!!!!!!!!!!!!!!!!!!!!!!)
plot.loc.stab <- function(folder, output.file="loc_Stab_plots2.pdf"){
  f <- list.files(folder)
  theFile <- output.file
  pdf(file = theFile)
  for (i in 1:length(f)){
    print(paste0(i,": ",f[i]))
    fl <- f[i]
    load(paste0(folder,"/",fl))
    specs <- model.spec.by.name(fl)
    reg.map <- to.regional(no.species = strtoi(specs[1]), no.webs = strtoi(specs[2]))
    loc.stab <- sapply(1:dim(abundances)[3], function(i) loc.stability(abundances = abundances[2500:3001,2:dim(abundances)[2],i],reg.map = reg.map))
    loc.stab <- t(loc.stab)
    d <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 0.001, .01, .02, .03, .04)    #dispersal rates
    matplot(d, loc.stab, t="l", lty=1, cex = d, xlab = "d", ylab="Regional Stability")
    title(main = paste0("Species=",specs[1], ", Patches=", specs[2], ", Connectance=",specs[3], ", V",specs[4]))
  }
  dev.off()
  
}

plot.sync.inter <- function(folder, output.file="sync_Inter.pdf"){
  f <- list.files(folder)
  theFile <- output.file
  pdf(file = theFile)
  for (fl in f){
    print(fl)
    load(paste0(folder,"/",fl))
    specs <- model.spec.by.name(fl)
    reg.map <- to.regional(no.species = strtoi(specs[1]), no.webs = strtoi(specs[2]))
    result <- sapply(1:dim(abundances)[3], function(i) inter.sync(abundances = abundances[2500:3001,2:dim(abundances)[2],i],reg.map = reg.map))
    result <- t(result)
    d <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 0.001, .01, .02, .03, .04)    #dispersal rates
    matplot(d, result, t="l", lty=1, cex = d, xlab = "d", ylab="Interspecies sunchrony")
    title(main = paste0("Species=",specs[1], ", Patches=", specs[2], ", Connectance=",specs[3], ", V",specs[4]))
  }
  dev.off()
}
  

#To plot local stability vs. d diagrams (exactly the same code as for reg.stab. to be merged into one function!!!!!!!!!!!!!!!!!!!!!!!!!!!)
plot.sync.intra <- function(folder, output.file="sync.intra.pdf"){
  f <- list.files(folder)
  theFile <- output.file
  pdf(file = theFile)
  for (i in 1:length(f)){
    print(paste0(i,": ",f[i]))
    fl <- f[i]
    load(paste0(folder,"/",fl))
    specs <- model.spec.by.name(fl)
    reg.map <- to.regional(no.species = strtoi(specs[1]), no.webs = strtoi(specs[2]))
    result <- sapply(1:dim(abundances)[3], function(i) intra.sync(abundances = abundances[2500:3001,2:dim(abundances)[2],i],reg.map = reg.map))
    result <- t(result)
    d <- c(1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 0.001, .01, .02, .03, .04)    #dispersal rates
    matplot(d, sapply(1:len), t="l", lty=1, cex = d, xlab = "d", ylab="Intraspecies synchrony")
    title(main = paste0("Species=",specs[1], ", Patches=", specs[2], ", Connectance=",specs[3], ", V",specs[4]))
  }
  dev.off()
  
}



