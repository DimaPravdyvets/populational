#Mutation dynamics
#recurrent
par(mfrow=c(3,4),bg='powderblue')

mut_dyn <- function(p0,mu,tgen){
  plot(0,0,xlim = c(0,tgen),ylim = c(0,1), main = paste("Recurrent mutation – mu = ",mu), type = "l",xlab=paste("Number of generations: ",tgen),ylab="q")
  #q = 1 - p0*(1- mu)^tgen
  freqvector <- c()
  freqvector[1]<- 1-p0
  for (i in 1:tgen){
    freqvector[i+1]<-freqvector[i] + mu*(1-freqvector[i])
  if(freqvector[i+1]==1){
    print(paste('q reached 1 at generation: ',i))
  }
  }
  #print(freqvector)
  lines(freqvector, type = "l")
}
mut_dyn(0.1,0.001,10000)
mut_dyn(0.4,0.001,10000)
mut_dyn(0.8,0.001,10000)
mut_dyn(1,0.001,10000)
mut_dyn(0.1,0.0001,10000)
mut_dyn(0.4,0.0001,10000)
mut_dyn(0.8,0.0001,10000)
mut_dyn(1,0.0001,10000)
mut_dyn(0.1,0.001,100000)
mut_dyn(0.4,0.001,100000)
mut_dyn(0.8,0.001,100000)
mut_dyn(1,0.001,100000)



#recurrent retromutation
mut_dyn_retro <- function(p0,mu,nu,tgen)
{plot(0,0,xlim = c(0,tgen),ylim = c(0,1),type = "l", main = paste("Mutation and retromutaion – mu = ",mu,"nu =",nu),
      xlab="Generation",ylab="q")
  freqvector <- c()
  freqvector[1]<- 1-p0
  for (i in 1:tgen){freqvector[i+1]<-freqvector[i] + mu*(1-freqvector[i])- nu*freqvector[i]}
  lines(freqvector, type = "l")
}
mut_dyn_retro(1,0.000001,0.001,10000)
