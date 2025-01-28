# DIC for Geometric
# fix theta at intervals between 0.1 and 0.9

dic.geom <- function(nn=nn, xtot=xtot, aa=aa, bb=bb){
  # approximated value of DIC for Geometric
  # just need to pass on sum(X_i) and sample size n
  
  alltot <- nn+aa + xtot+bb
  atot <- nn + aa
  btot <- xtot+bb
  
  t1 <- -2*nn*(log(atot) - log(alltot))
  
  t2.log <- log(btot) -log(atot) -log(alltot)
  t2 <- 2*nn*exp(t2.log)
  
  t3.log <- log(btot) - log(alltot)
  t3 <- -2*xtot*t3.log
  
  t4.log <- log(atot) + log(xtot) - log(btot) -log(alltot)
  t4 <- 2*exp(t4.log)
  
  dic <- t1+t2+t3+t4
  return(dic/nn)  # DIC/n
}

dic.limit <- function(theta){ 
  # theoretical limit of DIC.n
  aa <- exp(log(1-theta) - log(theta))*log(1-theta) + log(theta)
  return(-2*aa)
}


# settings, initialization
# can use range of priors for which theta has high to low probabilities

thetavec <- c(rep(c(1,10,100),3),rep(1,3),rep(10,3), rep(100,3)) #arbitrary
thetapriormat <- matrix(thetavec,length(thetavec)/2,2)
nreps <- 10 #no of replicate data sets under each setting
nsamlist <- 10^seq(2,7, by=1)  #memory issue after 10^8! also prob not needed

#  simulate to have a uniform range of thetas,
# (may have to change priors so in each case have those which give high prior mass to values?)

simthetavec <- seq(0.1,0.9, by=0.1)
nsettings.theta <- length(simthetavec)
Xtots <- vector("list", nsettings.theta)

# testing model should use different priors, but group results by single correct theta

dic.limit(simthetavec)
origtheta <- as.data.frame(cbind(1:nsettings.theta,thetapriormat,simthetavec,dic.limit(simthetavec))) # hyperparameters, theta, dic.limit
names(origtheta) <- c("index","a","b","theta","dic.lim.geom")
origtheta
save(origtheta,file="origtheta_fixed_geom.obj")


# dic has to be for nsam varying and different settings
nlist <- length(nsamlist)
dic.n <-  vector("list",nlist)
for(j in 1:nlist){
  dic.n[[j]] <- vector("list",nsettings.theta)
  for(i in 1:nsettings.theta){
    dic.n[[j]][[i]] <- vector("list",nsettings.theta)
  }
}


# main evaluation loop
for(i in 1:nsettings.theta){
  cat("theta setting ", i,"\n")
  
for(j in 1:nlist){ # for each n

    nsam <- nsamlist[j] #no of samples per set
     # create a "nsam"-sized dataset for each value of theta, and n replicates each
    
    Xlist <- replicate(n=nreps, expr=rgeom(n=nsam, prob=simthetavec[i]), simplify=F ) #dataset for theta i
    Xtots <- lapply(Xlist,sum)
  
    # calculate value of DIC in each case (1 for each replicate)
    # do for each combination of hyperparameters
    
    for(k in 1:nsettings.theta){ # over all priors
    
    dic.n[[j]][[i]][[k]] <- lapply(Xtots,dic.geom, nn=nsam, aa=thetapriormat[k,1], bb=thetapriormat[k,2])
    }
    # Do for increasing values of nsam 
  }
}

# replicates 1 to nreps each
# columns are n, alpha, beta, replicate no, DIC

longtable <- NULL

for(j in 1:nlist){
  nn <- nsamlist[j]

    init.table <- NULL #prior and data settings
    for(k in 1:nsettings.theta){
      aalpha <- thetapriormat[k,1]; bbeta <- thetapriormat[k,2]
      tmptable <- matrix(rep(c(nn,aalpha, bbeta),nreps), nreps, 3, byrow=T)
      init.table <- rbind(init.table, tmptable)
    }
  
    for(i in 1:nsettings.theta){
     ivec <- rep(i,nreps*nsettings.theta)
     dicvec <- unlist(dic.n[[j]][[i]]) 
     tmptable <- cbind(ivec, init.table,rep(1:nreps,nsettings.theta), dicvec)
     longtable <- rbind(longtable, tmptable)
  }
}

dicgeomtable <- as.data.frame(longtable)
names(dicgeomtable) <- c("index","n","alpha", "beta", "repl", "dic.geom")
# compare with dic limits, indexed by i

#Add column with true theta values, needed for legend
dicgeomtable$theta <- simthetavec[dicgeomtable$index]

write.csv(dicgeomtable, "dic-geom-fixThetanew.csv", row.names=F) 

#plots
# some settings if not running the whole file

# dicgeomtable <- read.csv("dic-geom-fixThetanew.csv")
# load("origtheta_fixed_geom.obj")
simthetavec <- seq(0.1,0.9, by=0.1)
nsettings.theta <- length(simthetavec)

library(ggplot2)
library(scales) # to access break formatting functions
library(RColorBrewer)

mycols <- brewer.pal(n = nsettings.theta, name = 'Paired')

#unlogged DIC_n: plot by index

my.labs <- vector("list",nsettings.theta)
for(i in 1:nsettings.theta){
  value <- simthetavec[i]
  my.labs[[i]] <- bquote(theta==.(value))
}
my.labs

pdf("DIC_n-geom_fixedTheta.pdf",width=10)
p <- ggplot(dicgeomtable, aes(n,dic.geom))
p + geom_point(aes(colour=factor(theta))) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_brewer(palette = "Paired",labels=my.labs) +
  geom_hline(yintercept = origtheta$dic.lim.geom, col=mycols[1:9]) +
  #scale_color_manual(labels = theta.true, values = mycols[1:9])+
  labs(col="", y="DIC.n")+
  theme_bw() +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16),
        legend.text=element_text(size=rel(1.5)))
dev.off()



