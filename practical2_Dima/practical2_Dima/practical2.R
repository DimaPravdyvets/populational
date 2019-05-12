############################## Script 1 ############################## 
#Genetic Drift
#Simulation of random drift trajectory of allele frequency for one locus with two alleles 
# Ne = effective population size, tgen= number of generation simulated, p_0 = initial frequency allele A)
#A single trajectory (a linear random walk)
#using sample function
genetic_drift<-function(p_0,Ne,tgen) {
  alelos<-c(1,0) #1 A and 0 al (numerical for estimating frequency directly)
  pvector<-c();
  pvector[1]<-p_0;
  p_ = p_0;
  xaxis<-c(0:tgen)
  
  for (i in 2:(tgen+1)){
    s_i =  sample(alelos,2*Ne, c(p_,1-p_), replace = T)
    p_i=sum(s_i)/length(s_i)
    pvector[i]<-p_i
    #print(pvector[i])
    p_=p_i
  }
  
plot(xaxis,pvector,xlim= c(0,tgen),ylim= c(0,1),type = "l", main = paste("Genetic Drift -Population size = ",Ne), xlab="Generation",ylab="Frequency")
} 

genetic_drift(0.3,5,100)
genetic_drift(0.3,25,100)
genetic_drift(0.3,50,100)
genetic_drift(0.3,100,100)
genetic_drift(0.3,250,100)
genetic_drift(0.3,1000,100)
genetic_drift(0.3,5000,100)
genetic_drift(0.3,50000,100)


############################## Script 2 ##############################
genetic_drift2<-function(p_0,Ne,tgen) {
  alelos<-c(1,0) #1 A and 0 al (numerical for estimating frequency directly)
  pvector<-c();
  pvector[1]<-p_0;
  p_ = p_0;
  xaxis<-c(0:tgen)
  mean_freq=c()
  mean_freq[1]<-p_0
  for (i in 2:(tgen+1)){
    s_i =  sample(alelos,2*Ne, c(p_,1-p_), replace = T)
    p_i=sum(s_i)/length(s_i)
    pvector[i]<-p_i
    p_=p_i
    mean_freq[i]<-mean(pvector)
  }
  

#print(mean_freq)
plot(mean_freq,xlab = "Generation",ylab = "frequency")
print(paste('mean: ', mean(mean_freq)))
print(paste('var: ', var(mean_freq)))
  
#plot()
}
genetic_drift2(0.3,5,1000)
genetic_drift2(0.3,25,1000)
genetic_drift2(0.3,50,1000)
genetic_drift2(0.3,100,1000)
genetic_drift2(0.3,250,1000)
genetic_drift2(0.3,1000,1000)
genetic_drift2(0.3,5000,1000)
genetic_drift2(0.3,50000,1000)

############################## Script 3 ##############################
genetic_drift_bottle<-function(p_0,Ne,tgen) {
  alelos<-c(1,0) #1 A and 0 al (numerical for estimating frequency directly)
  pvector<-c();
  pvector[1]<-p_0;
  p_ = p_0;
  xaxis<-c(0:tgen)
  bn<-c()
  for (i in 2:(tgen+1)){
    if (i%%9==0){
      #print(i)
      Ne=floor(Ne%/%10)
      bn<-c(bn,i)
    print(paste('bottleneck',Ne))
    s_i =  sample(alelos,2*Ne, c(p_,1-p_), replace = T)
    p_i=sum(s_i)/length(s_i)
    pvector[i]<-p_i
    p_=p_i
    #print(pvector[i])
      }
    else{
      Ne=floor(Ne*1.5)
    print(Ne)
      s_i =  sample(alelos,2*Ne, c(p_,1-p_), replace = T)
      p_i=sum(s_i)/length(s_i)
      pvector[i]<-p_i
      #print(pvector[i])
      p_=p_i
    }
    
  
  } 
  print(bn)
  
  plot(xaxis,pvector,xlim= c(0,tgen),ylim= c(0,1),type = "l", main = paste("Genetic Drift -Population size = ",Ne), xlab="Generation",ylab="Frequency")
  for(i in bn){
    
    abline(v=i,col='forestgreen')
  }
  legend('topleft',c('bottleneck'),col=c('forestgreen'),pch=15,cex=.8,box.lty=0)
  }
genetic_drift_bottle(0.3,5,30)
genetic_drift_bottle(0.5,10,100)
genetic_drift_bottle(0.3,5,100)
genetic_drift_bottle(0.9,50,100)
genetic_drift_bottle(0.3,20,100)


############################## Script 4 ##############################

install.packages('dae')
library(dae)
genetic_drift_varNe<-function(p_0,Ne_range,tgen) {
  alelos<-c(1,0) #1 A and 0 al (numerical for estimating frequency directly)
  pvector<-c();
  pvector[1]<-p_0;
  p_ = p_0;
  xaxis<-c(0:tgen)
  NeS<-c()
  
  
  for (i in 1:(tgen+1)){
    Ne=runif(1, Ne_range[1],Ne_range[2])
    NeS[i]<-Ne
    s_i =  sample(alelos,2*Ne, c(p_,1-p_), replace = T)
    p_i=sum(s_i)/length(s_i)
    pvector[i]<-p_i
    #print(pvector[i])
    p_=p_i
  }
  #print(NeS)
  plot(xaxis,pvector,xlim= c(0,tgen),ylim= c(0,1),type = "l", main = paste("Genetic Drift -Population size = ",Ne), xlab="Generation",ylab="Frequency")
  print(paste('Harmonic mean: ', harmonic.mean(NeS)))
}

genetic_drift_varNe(0.3,c(5,300),50)
genetic_drift(0.3,56.0048396951463,50)

genetic_drift_varNe(0.3,c(5,300),100)
genetic_drift(0.3,55.1927358828573,100)

genetic_drift_varNe(0.3,c(2,3000),50)

genetic_drift_varNe(0.3,c(1,300000),50)

genetic_drift_varNe(0.3,c(2,5),50)

############################## Script 5 ##############################

genetic_drift_binomial<-function(p_0,Ne,tgen) {
  alelos<-c(1,0) #1 A and 0 al (numerical for estimating frequency directly)
  pvector<-c();
  pvector[1]<-p_0;
  p_ = p_0;
  xaxis<-c(0:tgen)

  par(bg='aquamarine1')
  
  for (i in 2:(tgen+1)){
    p_=rbinom(1,2*Ne,p_)/(2*Ne)
    pvector=c(pvector,p_)
    }
  
  plot(xaxis,pvector,xlim= c(0,tgen),ylim= c(0,1),type = "l", main = paste("Genetic Drift -Population size = ",Ne), xlab="Generation",ylab="Frequency")
} 

genetic_drift_binomial(0.3,50,100)

############################## Script 6 ##############################


#Genetic Drift#Simulation of random drift trajectory of allele frequency for one locus with two alleles 
# Ne = effective population size, tgen= number of generation simulated, p_0 = initial frequency allele A)
#Multiple trajectories (multiple replicas)
library(RColorBrewer)

genetic_drift_big<-function(nrepl,p_0,Ne,tgen) {
  alelos<-c(1,0) #1 A and 0 al (numerical for estimating frequency directly)
  xaxis<-c(0:tgen)
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  par(bg='darkorchid4')
  
  plot(0,0,xlim = c(0,tgen),ylim= c(0,1),main = paste("Genetic Drift -Population size = ",Ne), type = "l",xlab="Generation",ylab="Frequency")
  mean_freq=c()
  mean_freq[1]<-p_0
  for (n in 1:nrepl){
    
    p_=p_0
    pvector<-c();
    pvector[1]<-p_0;
    
    for (i in 2:(tgen+1)){
      s_i =  sample(alelos,2*Ne, c(p_,1-p_), replace = T)
      p_i=sum(s_i)/length(s_i)
      pvector[i]<-p_i
      p_=p_i
      mean_freq[i]<-mean(pvector)
    }
    color=sample(col_vector,1)
    lines(pvector,col=color)
  }
  print(mean_freq)
  print(paste('mean: ', mean(mean_freq)))
  print(paste('var: ', var(mean_freq)))
}
genetic_drift_big(10,0.2,500,1000)


