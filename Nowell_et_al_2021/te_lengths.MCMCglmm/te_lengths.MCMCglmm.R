setwd("./")

## libraries
library(MCMCglmm)
library(phangorn)

## READ DATA
df<-read.table("ectopic.mean_superfam.csv", header=T, sep=",", colClasses=c('factor','factor','factor','integer','numeric','numeric','factor','factor'))
str(df)

## READ IN TREE
tree<-read.tree("rotifers_tree.newick")
## root tree
tree<-midpoint(tree)
## plot tree
plot(tree, underscore=T)

## make ultrametric and read back in to correct class
tree<-chronos(tree,lambda=1)
tree<-read.tree(text=write.tree(tree))
tree$node.label<-NULL
## convert tree into covariance matrix for mcmcglmm
inv.phylo<-inverseA(tree,nodes="TIPS",scale=TRUE)
## make prior
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))

## TEST FOR LENGTH DIFFERENCES BETWEEN BDELLOIDS & MONOGONONTS
## H0: there are no significant differences in TEs between bdelloids and monogononts

## fit DNA transposons model
mod1.DNA<-MCMCglmm(length.mean.log ~ is.bdelloid, random=~ID, data=subset(df, df$type=="DNA"),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.DNA, file="mod1.DNA.RData")

## check for convergence in posterior output
plot(mod1.DNA)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.DNA$VCV)[2, , ])
diag(autocorr(mod1.DNA$Sol)[2, , ])

## summary
summary(mod1.DNA)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 7287.062 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.1307  0.03064   0.2743     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units    0.7457   0.7066   0.7847     2000
# 
# Location effects: length.mean.log ~ is.bdelloid 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     4.6326   3.9337   5.3507     2000 <5e-04 ***
# is.bdelloid1    0.3323  -0.5673   1.3097     2000   0.45    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## fit LTR transposons model
mod1.LTR<-MCMCglmm(length.mean.log ~ is.bdelloid, random=~ID, data=subset(df, df$type=="LTR"),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.LTR, file="mod1.LTR.RData")

## check for convergence in posterior output
plot(mod1.LTR)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.LTR$VCV)[2, , ])
diag(autocorr(mod1.LTR$Sol)[2, , ])

## summary
summary(mod1.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 1219.206 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     0.119 0.002684      0.4     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     1.084   0.9321    1.234     2000
# 
# Location effects: length.mean.log ~ is.bdelloid 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     5.1908   4.4519   5.8715     2132 <5e-04 ***
# is.bdelloid1    0.1586  -0.8049   1.0529     2000  0.663    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## fit LINE transposons model
mod1.LINE<-MCMCglmm(length.mean.log ~ is.bdelloid, random=~ID, data=subset(df, df$type=="LIN"),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.LINE, file="mod1.LINE.RData")

## check for convergence in posterior output
plot(mod1.LINE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.LINE$VCV)[2, , ])
diag(autocorr(mod1.LINE$Sol)[2, , ])

## summary
summary(mod1.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 3112.116 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID   0.05541 0.002244   0.1762     1836
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units      1.12    1.025    1.216     2000
# 
# Location effects: length.mean.log ~ is.bdelloid 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     4.8464   4.3570   5.3285     2000 <5e-04 ***
# is.bdelloid1    0.2569  -0.3829   0.9120     1855  0.343    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



## TEST FOR LENGTH DIFFERENCES BETWEEN DESICCATING & NONDESICCATING BDELLOID SPECIES
## H0: no significant differences in TE lengths between desiccating and nondesiccating bdelloids

## bdelloid tree  
tree.bdelloid<-drop.tip(tree,tree$tip.label[grep("B",tree$tip.label)])
plot(tree.bdelloid, underscore = T)

## convert tree into covariance matrix for mcmcglmm
inv.phylo<-inverseA(tree.bdelloid,nodes="TIPS",scale=TRUE)

## fit DNA transposon model
mod2.DNA<-MCMCglmm(length.mean.log ~ is.desiccating, random=~ID, data=subset(df, df$is.bdelloid==1 & df$type=="DNA"),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.DNA, file="mod2.DNA.RData")

## check for convergence in posterior output
plot(mod2.DNA)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.DNA$VCV)[2, , ])
diag(autocorr(mod2.DNA$Sol)[2, , ])
## summary
summary(mod2.DNA)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 6165.208 
# 
# G-structure:  ~ID
# 
#   post.mean l-95% CI u-95% CI eff.samp
# ID   0.06501  0.01235   0.1532     1578
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units    0.7948   0.7492   0.8383     2000
# 
# Location effects: length.mean.log ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)       4.92406  4.53246  5.34050     2000 <5e-04 ***
# is.desiccating1   0.06235 -0.19216  0.34265     2000  0.605    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit PLE transposon model
mod2.PLE<-MCMCglmm(length.mean.log ~ is.desiccating, random=~ID, data=subset(df, df$is.bdelloid==1 & df$type=="PLE"),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.PLE, file="mod2.PLE.RData")

## check for convergence in posterior output
plot(mod2.PLE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.PLE$VCV)[2, , ])
diag(autocorr(mod2.PLE$Sol)[2, , ])
## summary
summary(mod2.PLE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 19.45105 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID   0.07614 0.005479    0.208     1865
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units   0.07122  0.04645  0.09777     2389
# 
# Location effects: length.mean.log ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)       5.79505  5.33269  6.21470     2176 <5e-04 ***
# is.desiccating1   0.05808 -0.29186  0.33688     2165  0.662    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit LTR transposon model
mod2.LTR<-MCMCglmm(length.mean.log ~ is.desiccating, random=~ID, data=subset(df, df$is.bdelloid==1 & df$type=="LTR"),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.LTR, file="mod2.LTR.RData")

## check for convergence in posterior output
plot(mod2.LTR)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.LTR$VCV)[2, , ])
diag(autocorr(mod2.LTR$Sol)[2, , ])
## summary
summary(mod2.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 1053.097 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID   0.05638 0.002153    0.188     1866
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     1.167   0.9875    1.341     2219
# 
# Location effects: length.mean.log ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)        5.2051   4.7681   5.6551     2000 <5e-04 ***
# is.desiccating1    0.2054  -0.1357   0.5375     2000  0.199    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit LINE transposon model
mod2.LINE<-MCMCglmm(length.mean.log ~ is.desiccating, random=~ID, data=subset(df, df$is.bdelloid==1 & df$type=="LIN"),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.LINE, file="mod2.LINE.RData")

## check for convergence in posterior output
plot(mod2.LINE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.LINE$VCV)[2, , ])
diag(autocorr(mod2.LINE$Sol)[2, , ])
## summary
summary(mod2.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 2647.869 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID   0.03873 0.002497   0.1173     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     1.232     1.11    1.345     2000
# 
# Location effects: length.mean.log ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)       5.13957  4.78399  5.48862     2000 <5e-04 ***
# is.desiccating1  -0.02667 -0.29145  0.23584     2000  0.837    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
