#########correlation and  genetics covariance
#########
library(tidyverse)
library(dplyr)
library(lme4)
#########################################################
###############################
######Load data, we need see how many variables we have
data <- read.csv("data.csv",head=T)
head(data)
str(data)
names(data)[1] <- "Parents"
data$Parents <- as.factor(data$Parents)
data$rep <- as.factor(data$rep)
names(data) 
#######################################################################
## First steps we do the sum of each variblae;for example; x1+x2=x1x2, x1+x3,. etc..
#####################################################################
########################

co1 <- names(data[4:6])
co2 <- names(data[5:6])
co3 <- names(data[6])
re <- NULL;
re2 <- NULL;
re3 <- NULL;
for (i in co1) {
  for (j in co2) {
    for (k in co3 ) {
      
      t <- data[3]+data[i]
      tt <- data[4]+data[j]
      ttt <- data[5]+data[k]
      re[[i]] <- t
      re2[[j]] <- tt 
      re3[[k]] <- ttt
    }
  }
}

########################################################################
#######################################################################
## Now we extract the value for  each varible
cor1 <- as.data.frame(re)
names(cor1) <- c("x1x2","x1x3","x1x4")
cor2 <- as.data.frame(re2 )
names(cor2) <- c("x2x3","x2x4")
cor3 <- as.data.frame(re3)
names(cor3) <- c("x3x4")
#########################################################################
#####put togheter the varibles, sums, repetition and variables
#### with this data frame we are going to work
sum_prod<- cbind(data,cor1,cor2,cor3)
head(sum_prod)

################################################
##############################################

valor_genotipo <-c()## <- save the genotype variance and covariance
valor_error<-c()####### <- save the error
y <-NULL;####################save  genotyoe, error , covariance and the phenotypic value
varlist <- names(sum_prod)[3:12] 
##########################################################

for (i in varlist){
  form <- reformulate(c("(1|Parents)","(1|rep)"),response=i)
  lmer_results <- lmer(form, data=sum_prod)
  lmer_summary <- summary(lmer_results)
  varc <- VarCorr(lmer_results,comp="Variance")
   geno <- as.vector(varc$Parents)
   err <- attr(varc,'sc')^2
   valor_genotipo[[i]] <-c(geno)
   valor_error[[i]] <- c(err)
   y <- rbind(y,c(geno,err,(geno+err)))
}

resul <- as.data.frame(y)
names(resul) <- c("Genotipo","Error","Phenotype")
rownames(resul) <-varlist
resul#########################the first  3 value  belong to  components of variance to our variables
################################the next are  variance componets from sums of products

#######Estimation de genetics correlations

###############################################
gen_var <- valor_genotipo[1:4]
cov_var <- valor_genotipo[5:10]
#################################################
##############
cov_x1 <- NULL
corr_x1 <- NULL
x1 <- NULL
for (i in 1:3) {
  v1 <- gen_var[2:4]
  v2 <- gen_var[3:4]
  cova_genict<- (cov_var[i]-gen_var[1]-v1[i])/2
  cov_x1[[i]]<-cova_genict
  corr_genetics<-cov_x1[[i]]/sqrt(gen_var[1]*v1[i])
 corr_x1[i]<-corr_genetics
  x1[[i]]<- cbind(cova_genict,corr_genetics)
  }
########################################################
cov_x2 <-NULL
corr_x2 <- NULL
x2 <-  NULL


for (i in 1:2) {
  vc <- cov_var[4:5]
  v2 <- gen_var[3:4]
  cova_genict<- (vc[i]-gen_var[2]-v2[i])/2
  cov_x2[[i]]<-cova_genict
  cor_gen <- cov_x2[[i]]/sqrt(gen_var[2]*v2[i])
  corr_x2[[i]] <-cor_gen 
  x2[[i]]<- cbind(cova_genict,cor_gen)
  }
##########################################################
cov_x3 <- (cov_var[6]-gen_var[3]-gen_var[4])/2
corr_x3 <- cov_x2/sqrt(gen_var[3]*gen_var[4])
#######################################################
######################################################
#the correlations  and covarince gentics among ears/len  ears/weight  ears/yield :
corr_x1 
cov_x1 
### the correlations  among le/weight  le /yield
corr_x2 
cov_x2
##### and the correlations  weight/ yield is:
corr_x3
cov_x3 
###
##Singh, R.K. and Chaudhary, B.D. (1987) Biometrical Methods in Quantitative Genetic Analysis. Kalyani Publishers, New Delhi, Ludhiana, India, 318.
