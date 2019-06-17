names(data)
traits <- c("ears" ,"len"  ,   "weight"  ,"yield" )
co_ge(Traits,data$Parents,data$rep,data)
co_gp<- function(Traits,Entrada,Rep,data){
  names(traits)
  traits <- Traits
  geno<-Entrada
  rep <- rep
  names(traits)
###########################################
  
  geno <- as.character(data$Parents)
  rep <- as.character(data$rep)
  
  # Creating the matrix  G and P
  
  leng_traits<- length(traits) # number of traits
  G <- matrix(nrow = leng_traits, ncol = leng_traits) 
  P <- matrix(nrow = leng_traits, ncol = leng_traits) 
  l_g <-length(levels(as.factor(data$Parents))) 
  nrep <- length(levels(as.factor(data$rep))) 
  
  # Estimation of components of variance using the REML method.


    for (i in 1:nt) {
      y <- data[,traits[i]]
      fm <- lme4::lmer(y ~ (1|geno) + (1|rep))
      vc <- lme4::VarCorr(fm)
      G[i, i] <- vc$g[1]
      P[i, i] <- vc$g[1] + attr(vc, "sc")^2 /nrep
    }

  ####Sum de each variable  and  estimation of  variance components
  ###Example X1+X2=X1X2, X1+X3=X1X3 etc
  
      for (i in 1:(nt - 1)) {
        for (j in (i + 1):nt) {
       y<-data[,traits[i]] + data[,traits[j]]
          fm <- lme4::lmer(y ~ (1|g) + (1|r))
        vcz <- lme4::VarCorr(fm) 
       G[i, j] <- G[j, i] <- (vcz$g[1] - G[i, i] - G[j, j]) / 2
        P[i, j] <- P[j, i] <- (vcz$g[1] + attr(vcz, "sc")^2 / nrep - P[i, i] - P[j, j]) / 2
       
        }
      }
  
    d1 <- diag(diag(G)^{-0.5}, nt, nt)
    d2 <- diag(diag(P)^{-0.5}, nt, nt)
    GC <- d1 %*% G %*% d1 # Genotypic correlation matrix
    PC <- d2 %*% P %*% d2 # Phenotypic correlation matrix

  # results
  
  results<- list(Genetic_Cov = G, Pheno_Cov = P, Genetic_Cor = GC, Pheno_Cor = PC)
  print(results)

}