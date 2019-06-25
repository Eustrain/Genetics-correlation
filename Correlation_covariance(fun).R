library(lme4)
data <- read.csv("data.csv",head=T)

names(data)
Traits <- c("ears" ,"len"  ,   "weight"  ,"yield" )
co_gp(Traits,data$Parents,data$rep,data)

co_gp<- function(Traits,Entry,Rep,data){

  traits <- Traits
  geno<-as.factor(Entry)
  rep <-as.factor(Rep)
###########################################  
 nrep <- length(levels(rep))
leng_traits<- length(traits) # number of traits
  
###########################################
  
  # Creating the matrix  G and P
  
 
  G <- matrix(nrow = leng_traits, ncol = leng_traits) 
  P <- matrix(nrow = leng_traits, ncol = leng_traits) 
  # Estimation of  variance each traits.

    for (i in 1: leng_traits) {
      y <- data[,traits[i]]
      fm <- lmer(y ~ (1|geno)+(1|rep))
      vc <- VarCorr(fm,comp="Variance")
      G[i, i] <- vc$geno[1]
      P[i, i] <- vc$geno[1] + attr(vc, "sc")^2/nrep 
    }

  ####Sum de each variable  and  estimation of  variance components
  ###Example X1+X2=X1X2, X1+X3=X1X3 etc
  
      for (i in 1:( leng_traits - 1)) {
        for (j in (i + 1): leng_traits) {
       y<-data[,traits[i]] + data[,traits[j]]
          fm <-lmer(y ~ (1|  geno) + (1| rep))
        varcor <-VarCorr(fm) 
       G[i, j] <- G[j, i] <- (varcor$geno[1] - G[i, i] - G[j, j]) / 2
        P[i, j] <- P[j, i] <- (varcor$geno[1] + attr(varcor, "sc")^2 / nrep - P[i, i] - P[j, j]) / 2
       
        }
      }
  
  ####################Estimation of correlation and covariance
    diag_G <- diag(diag(G)^{-0.5},  leng_traits,  leng_traits)
    diag_P<- diag(diag(P)^{-0.5},  leng_traits,  leng_traits)
    GC <- diag_G %*% G %*% diag_G# Genotypic correlation matrix
    PC <- diag_P %*% P %*% diag_P # Phenotypic correlation matrix
####################Names of  matrix    
    row.names(G) <- Traits
    colnames(G) <- Traits
    row.names(P) <- Traits
    colnames(P) <- Traits
    row.names(GC) <- Traits
    colnames(GC) <- Traits
    row.names(PC) <- Traits
    colnames(PC) <- Traits
 
    G <- round(G,2)
    P <- round(P,2)
    GC <- round(GC,2)
    PC <- round(PC,2)
 # results
  
  results<- list(Genetic_Cov = G, Pheno_Cov = P, Genetic_Cor = GC, Pheno_Cor = PC)
  print(results)

}

##Singh, R.K. and Chaudhary, B.D. (1987) Biometrical Methods in Quantitative Genetic Analysis. Kalyani Publishers, New Delhi, Ludhiana, India, 318.
