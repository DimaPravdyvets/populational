###################################EX 1 ###################################


#ESTIMATION OF GENOTYPE AND ALLELE FREQUENCIES FROM A GENOTYPING SURVEY

gene_freq<-function(N11,N12,N22){
  N=N11+N12+N22;#N total sample
  p=(N11+N12/2)/N; #p and q allele frequency X1 and X2 respectively
  q=1-p;
  f11=N11/N;
  f12=N11/N;
  f22=N22/N;
  print(paste(p,"=p",q,'=q',f11,'=f11',f12,'=f12',f22,'=f22'));
}
gene_freq(1787,3037,1305)

#optionsl function for N numbers
install.packages('rlist')
library(rlist)


#input should be in this format: (list(N11,N12,N13,N1N,N21,N22,N23,N2N,ETC),#of alleles)
#remembre that 12 and 21 is the same, so in case you need to, just divide those combinations by 2
# Constructing Quadratic Formula

Ngen_greq <- function(L,n){
  N=0
  for (i in L){
    N=N+i
  }
  #print(floor(result(1,1,-2*length(L))))
  falafel = vector("list")
  step = 1
  while (step<=length(L)){
    p = (L[[step]]/N)
    falafel=c(falafel,sqrt(p))
    step=step+n+1
  }
  freq = vector("list")
  for (j in L){
    w = (j/N)
    freq=c(freq,w)
  }
  #print(do.call(sum,freq))
  #print(do.call(sum,fallel))
  print("frequencies genotypes:")
  print(paste(freq))
  print(paste("Sum of genotype frequencies: ",Reduce('+',freq)))
  
  print("frequencies alleles:")
  print(paste(falafel))
  print(paste("sum of alleles frequencies:",Reduce('+',falafel)))
}
Ngen_greq(list(1787,3037/2,3037/2,1305),2)
Ngen_greq(list(1787,3037/2,1305/2,1234/2,2343,3325/2,1567,2235/2,3321),3)
Ngen_greq(list(1787,3037/2,1305/2,1234/2,2343/2,3325,1567/2,2235/2,3321/2,3254,1458/2,2541/2,4444/2,4565/2,1211/2,1111),4)
Ngen_greq(list(11,12/2,13/2,14/2,15/2,21/2,22,23/2,24/2,25/2,31/2,32/2,33,34/2,35/2,41/2,42/2,43/2,44,45/2,51/2,52/2,53/2,54/2,55),5)


###################################EX 2 ###################################



#Series of scripts for the random mating exercises

#MENDELIAN POPULATION
#RANDOM MATING -HWE
#One-gene two-alleles (X1 and X2). Genotypes: X11, X12 and X22
#One-gene two-alleles (X1 and X2). Genotypes: X11, X12 and X22
#f11, f12 and f22 genotype probabilities (frequencies) initial generation genotypes X1X1, X1X2 and X2X2, respectively
random_mating<-function(f11,f12,f22){
  p=f11+f12/2;
  q=1-p #p and q allele frequency X1 and X2 respectively
  #Mating pairs (frequencies) probabilities
  f11_11=f11*f11
  f11_12=2*f11*f12
  f11_22=2*f11*f22
  f12_12=f12*f12
  f12_22=2*f12*f22
  f22_22=f22*f22
  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  f11_=f11_11+f11_12/2+f12_12/4
  f12_=f11_12/2+f11_22+f12_12/2+f12_22/2
  f22_=f12_12/4+f12_22/2+f22_22
  #allele probability of progeny p_
  p_=f11_+f12_/2;
  q_=1-p_
  #output
  print(paste("Allele prob. initial generation (p and q) ",p,q))
  print(paste("Genotype prob. initial generation (f11, f12, f22) ",f11,f12,f22))
  print(paste("Allele prob. next generation (p' and q') ",p_,q_))
  print(paste("Genotype prob. next generation (f'11, f'12, f'22) ",f11_,f12_,f22_))
  print(paste("Genotype prob. next generation according HWE ",p_*p_,2*p_*q_,q_*q_))
}
random_mating(0.25,0.5,0.25)


#Complete assortative positive mating among genotypes

random_mating1<-function(f11,f12,f22){
  p=f11+f12/2;
  q=1-p #p and q allele frequency X1 and X2 respectively
  #Mating pairs (frequencies) probabilities
  ProbAllMatings = f11*f11 + f12*f12 + f22*f22
  f11_11 = f11*f11/ ProbAllMatings
  f11_12 = 0
  f11_22 = 0
  f12_12 = f12*f12/ ProbAllMatings
  f12_22 = 0
  f22_22 = f22*f22/ ProbAllMatings
  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  f11_=f11_11+f11_12/2+f12_12/4
  f12_=f11_12/2+f11_22+f12_12/2+f12_22/2
  f22_=f12_12/4+f12_22/2+f22_22
  #allele probability of progeny p_
  p_=f11_+f12_/2;
  q_=1-p_
  # output
  print(paste("Allele prob. initial generation (p and q) ",p,q))
  print(paste("Genotype prob. initial generation (f11, f12, f22) ",f11,f12,f22))
  print(paste("Allele prob. next generation (p' and q') ",p_,q_))
  print(paste("Genotype prob. next generation (f'11, f'12, f'22) ",f11_,f12_,f22_))
  print(paste("Genotype prob. next generation according HWE ",p_*p_,2*p_*q_,q_*q_))
}
random_mating1(0.25,0.5,0.25)
random_mating1(0.333333333333333,0.333333333333333,0.333333333333333)


#Complete assortative positive mating among phenotypes
random_mating2<-function(f11,f12,f22){
  p=f11+f12/2;
  q=1-p #p and q allele frequency X1 and X2 respectively
  #Mating pairs (frequencies) probabilities
  ProbAllMatings = f11*f11 + f12*f12 + f22*f22
  f11_11 = f11*f11/ProbAllMatings
  f11_12 = f11*f12/ProbAllMatings
  f11_22 = 0
  f12_12 = f12*f12/ProbAllMatings
  f12_22 = 0
  f22_22 = f22*f22/ProbAllMatings
  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  f11_=f11_11+f11_12/2+f12_12/4
  f12_=f11_12/2+f11_22+f12_12/2+f12_22/2
  f22_=f12_12/4+f12_22/2+f22_22
  #allele probability of progeny p_
  p_=f11_+f12_/2;
  q_=1-p_
  # output
  print(paste("Allele prob. initial generation (p and q) ",p,q))
  print(paste("Genotype prob. initial generation (f11, f12, f22) ",f11,f12,f22))
  print(paste("Allele prob. next generation (p' and q') ",p_,q_))
  print(paste("Genotype prob. next generation (f'11, f'12, f'22) ",f11_,f12_,f22_))
  print(paste("Genotype prob. next generation according HWE ",p_*p_,2*p_*q_,q_*q_))
}
random_mating2(0.25,0.5,0.25)

#Complete assortative negative mating among genotypes (a
#given genotype never mate with its same genotype and mate at
#random with other genotypes)

random_mating3<-function(f11,f12,f22){
  p=f11+f12/2;
  q=1-p #p and q allele frequency X1 and X2 respectively
  #Mating pairs (frequencies) probabilities
  ProbAllMatings = f11*f11 + f12*f12 + f22*f22
  f11_11 = 0
  f11_12 = f11*f12/ ProbAllMatings
  f11_22 = f11*f22/ ProbAllMatings
  f12_12 = 0
  f12_22 = f12*f22/ ProbAllMatings
  f22_22 = 0
  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  f11_=f11_11+f11_12/2+f12_12/4
  f12_=f11_12/2+f11_22+f12_12/2+f12_22/2
  f22_=f12_12/4+f12_22/2+f22_22
  #allele probability of progeny p_
  p_=f11_+f12_/2;
  q_=1-p_
  # output
  print(paste("Allele prob. initial generation (p and q) ",p,q))
  print(paste("Genotype prob. initial generation (f11, f12, f22) ",f11,f12,f22))
  print(paste("Allele prob. next generation (p' and q') ",p_,q_))
  print(paste("Genotype prob. next generation (f'11, f'12, f'22) ",f11_,f12_,f22_))
  print(paste("Genotype prob. next generation according HWE ",p_*p_,2*p_*q_,q_*q_))
}
random_mating3(0.25,0.5,0.25)


#Complete assortative negative mating among phenotypes (a
#given phenotypes never mate with its same phenotype)
random_mating4<-function(f11,f12,f22){
  p=f11+f12/2;
  q=1-p #p and q allele frequency X1 and X2 respectively
  #Mating pairs (frequencies) probabilities
  ProbAllMatings = f11*f11 + f12*f12 + f22*f22
  f11_11 = 0
  f11_12 = 0
  f11_22 = f11*f22/ ProbAllMatings
  f12_12 = 0
  f12_22 = f12*f22/ ProbAllMatings
  f22_22 = 0
  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  f11_=f11_11+f11_12/2+f12_12/4
  f12_=f11_12/2+f11_22+f12_12/2+f12_22/2
  f22_=f12_12/4+f12_22/2+f22_22
  #allele probability of progeny p_
  p_=f11_+f12_/2;
  q_=1-p_
  # output
  print(paste("Allele prob. initial generation (p and q) ",p,q))
  print(paste("Genotype prob. initial generation (f11, f12, f22) ",f11,f12,f22))
  print(paste("Allele prob. next generation (p' and q') ",p_,q_))
  print(paste("Genotype prob. next generation (f'11, f'12, f'22) ",f11_,f12_,f22_))
  print(paste("Genotype prob. next generation according HWE ",p_*p_,2*p_*q_,q_*q_))
}
random_mating4(0.25,0.5,0.25)


#segregation#did not understand the exercise


#x-linked# not working, just gave it a try

mating_xlinked<-function(xx11,xx12,xx22,x1,x2){
  p = (xx11+xx12/2+x1)/(sum(xx11,xx12,xx22,x1,x2))
  q=1-p
  #Mating pairs (frequencies) probabilities
  xx11_11=xx11*xx11
  xx11_12=2*xx11*xx12
  xx11_22=2*xx11*xx22
  xx12_12=xx12*xx12
  xx12_22=2*xx12*xx22
  xx22_22=xx22*xx22
  xx11_1=xx11*x1
  xx11_2=1.5*xx11*x2
  xx12_1=1.5*xx12*x1
  xx12_2=1.5*xx12*x2
  xx22_1=1.5*xx22*x1
  xx22_2=xx22*x2

  #Genotype probabilities of progeny (offspring) f11_, f12_, f22_
  xx11_=xx11_11+xx11_12/2+xx12_12/4
  xx12_=xx11_12/2+xx11_22+xx12_12/2+xx12_22/2
  xx22_=xx12_12/4+xx12_22/2+xx22_22
  #allele probability of progeny p_
  p_=xx11_+xx12_/2;
  q_=1-p_
  #output
  print(paste("Allele prob. initial generation (p and q) ",p,q))
  print(paste("Genotype prob. initial generation (f11, f12, f22) ",xx11,xx12,xx22))
  print(paste("Allele prob. next generation (p' and q') ",p_,q_))
  print(paste("Genotype prob. next generation (f'11, f'12, f'22) ",xx11_,xx12_,xx22_))
  print(paste("Genotype prob. next generation according HWE ",p_*p_,2*p_*q_,q_*q_))
  
  
  print(paste(p,q))
}

mating_xlinked(12,23,22,1,4)



###################################EX 3 ###################################



dima_square <- function(O,E){
  d2=0
  for (i in 1:length(O)){
    d2=d2+(((O[[i]]-E[[i]])**2)/E[[i]])
    
  }
  return(d2)
}

#HWE CHI-SQUARE TEST ON COUNTS FROM A GENOTYPING SURVEY
#One-gene two-alleles (X1 and X2). Genotypes: X11, X12 and X22. Sample size genotype ij -> Nij
HWE_test <- function(N11,N12,N22){
  N=N11+N12+N22; #N total sample
  #Allele frequencies
  p=(N11+N12/2)/N; #p and q allele frequency X1 and X2 respectively
  q=1-p;
  #Genotype frequencies
  f11=N11/N #f11, f12 and f22 genotype frequency genotypes X1X1, X1X2 and X2X2, respectively
  f12=N12/N
  f22=N22/N
  #Expected HWE
  E11=N*p^2
  E12=N*2*p*q
  E22=N-E11-E12
  #Chi-square value and probability
  chicuadrado <- dima_square(list(N11,N12,N22),list(E11,E12,E22))
  #print(chicuadrado)
  pval_chi <- pchisq(chicuadrado,df=2,lower=T)
  #output
  print(paste("p-value of chi-square: ",pval_chi))
  if(pval_chi<=0.05){
    print('Data is in hwe')
  }
  else{
    print('Data is not in hwe')
  } 
  print(paste('p: ',p, 'q: ', q))
  print(paste('Genotype frequencies: ','f11: ',f11, 'f12: ', f12, 'f22: ', f22))
}

HWE_test(1787,3037,1305)




###################################EX 4 ###################################



#ESTIMATION OF ALLELE AND GENOTYPE FREQUENCIES FROM COUNTS OF DOMINANT AND
#RECESSIVE PHENOTYPES ASSUMING HWE
#one-locus two-alleles with dominance A>a (i.e. Rh group)
Freq_dominance <- function(NA_,Naa){
  N=NA_+Naa; #N total sample
  
  #Phenotypic frequencies
  FNA_ = NA_/N
  FNaa = 1-FNA_
  
  #Genotype frequencies according HWE
  q=sqrt(FNaa)
  p=1-q
  #output
  print(paste('Phenotypic frequencies: ',FNA_, FNaa))
  print(paste('p and q: ',p,q))
  print(paste('AA~Aa~aa: ',p**2, 2*p*q, q**2))
  print(paste('Expected AA: ',NA_*p**2, 'Expected Aa: ', NA_*2*p*q, 'Expected aa: ', FNaa*q**2))
  }
Freq_dominance(115,19)

