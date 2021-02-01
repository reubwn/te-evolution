setwd("./")

##LIBS
library(phangorn)
library(MCMCglmm)

##READ IN DATA
df.te_load<-read.table("te_load.txt", head=T, sep="\t")
##drop megabubbles samples
df.te_load <- df.te_load[!grepl("10x_m",df.te_load$ID), ]
##convert ID, Desiccation to factor
df.te_load$ID <- as.factor(df.te_load$ID)
df.te_load$Desiccating <- as.factor(df.te_load$Desiccating)

##check what we have
str(df.te_load)
head(df.te_load)
levels(df.te_load$ID) ##IDs

##READ IN ROTIFER TREE (FIG 2A)
tree <- read.tree("rotifers_tree.newick")
##SUBSET BDELLOID TIPS
tree <- drop.tip(tree,tree$tip.label[grep("B",tree$tip.label)]) ## removes monogononts

##root tree
tree<-midpoint(tree)
##plot tree
plot(tree, underscore=T)

##make ultrametric and read back in to correct class
tree<-chronos(tree, lambda=1)
tree<-read.tree(text=write.tree(tree))
tree$node.label<-NULL
##convert tree into covariance matrix for mcmcglmm
inv.phylo<-inverseA(tree, nodes="TIPS", scale=TRUE)
##default prior
prior<-list(G=list(G1=list(V=1, nu=0.02)), R=list(V=1, nu=0.02))

##TESTS FOR EFFECT OF DESICCATION
##H0 = no effect of desiccation on repeat load
##hypothesis may apply to any interspersed repeated sequence in genome, not just TEs

##fit model for DNA transposons
mod.DNA<-MCMCglmm(log(DNA.transposons) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                  family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                  nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.DNA, file="mod.DNA.RData")

##check for convergence in posterior output
plot(mod.DNA)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.DNA$VCV)[2, , ])
diag(autocorr(mod.DNA$Sol)[2, , ])
##summary
summary(mod.DNA)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -81.68841 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.3319   0.1248   0.6296     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.004212 0.001874 0.007361     2000
# 
# Location effects: log(DNA.transposons) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)     0.5929  -0.2156   1.5098     2000 0.154
# Desiccating1    0.2702  -0.2685   0.8309     2000 0.325

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.DNA <- mod.DNA$VCV[,'ID']/(mod.DNA$VCV[,'ID']+mod.DNA$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.DNA)
posterior.mode(lambda.DNA)
HPDinterval(lambda.DNA)

##fit model for RC
mod.RC<-MCMCglmm(log(Rolling.circles) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                 family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                 nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.RC, file="mod.RC.RData")

##check for convergence in posterior output
plot(mod.RC)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.RC$VCV)[2, , ])
diag(autocorr(mod.RC$Sol)[2, , ])
##summary
summary(mod.RC)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 13.32131 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.1044 0.002066   0.3285     2051
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units   0.07577  0.03832   0.1215     2163
# 
# Location effects: log(Rolling.circles) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)   -1.55757 -2.13438 -1.01609     2000 <5e-04 ***
# Desiccating1   0.31031 -0.05095  0.75115     2000    0.1    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.RC <- mod.RC$VCV[,'ID']/(mod.RC$VCV[,'ID']+mod.RC$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.RC)
posterior.mode(lambda.RC)
HPDinterval(lambda.RC)

##fit model for PLEs
mod.PLE<-MCMCglmm(log(Penelope) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                  family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                  nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.PLE, file="mod.PLE.RData")

##check for convergence in posterior output
plot(mod.PLE)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.PLE$VCV)[2, , ])
diag(autocorr(mod.PLE$Sol)[2, , ])
##summary
summary(mod.PLE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -64.27166 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.4219    0.108   0.8595     1818
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.006917 0.002294   0.0122     1880
# 
# Location effects: log(Penelope) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)    -0.7065  -1.7797   0.1957     2000 0.148
# Desiccating1    0.2534  -0.3722   0.8787     2000 0.403

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.PLE <- mod.PLE$VCV[,'ID']/(mod.PLE$VCV[,'ID']+mod.PLE$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.PLE)
posterior.mode(lambda.PLE)
HPDinterval(lambda.PLE)

##fit model for LTRs
mod.LTR<-MCMCglmm(log(LTRs) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                  family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                  nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.LTR, file="mod.LTR.RData")

##check for convergence in posterior output
plot(mod.LTR)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.LTR$VCV)[2, , ])
diag(autocorr(mod.LTR$Sol)[2, , ])
##summary
summary(mod.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -14.12179 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     1.438 0.007867    3.206     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units    0.0381 0.002064  0.09064     2338
# 
# Location effects: log(LTRs) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)    -0.1102  -2.4048   2.0681     2000 0.920
# Desiccating1   -0.2529  -0.8634   0.4939     1797 0.427

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.LTR <- mod.LTR$VCV[,'ID']/(mod.LTR$VCV[,'ID']+mod.LTR$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.LTR)
posterior.mode(lambda.LTR)
HPDinterval(lambda.LTR)

##fit model for LINEs
mod.LINE<-MCMCglmm(log(LINEs) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.LINE, file="mod.LINE.RData")

##check for convergence in posterior output
plot(mod.LINE)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.LINE$VCV)[2, , ])
diag(autocorr(mod.LINE$Sol)[2, , ])
##summary
summary(mod.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -81.5403 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.6947   0.2469    1.212     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.004203 0.001456 0.008369     2000
# 
# Location effects: log(LINEs) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)     0.2278  -0.9042   1.6127     2000 0.693
# Desiccating1   -0.2107  -1.0420   0.5694     2302 0.579

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.LINE <- mod.LINE$VCV[,'ID']/(mod.LINE$VCV[,'ID']+mod.LINE$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.LINE)
posterior.mode(lambda.LINE)
HPDinterval(lambda.LINE)[2]

##fit model for Satellite
mod.Satellite<-MCMCglmm(log(Satellite) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                   family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                   nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.Satellite, file="mod.Satellite.RData")

##check for convergence in posterior output
plot(mod.Satellite)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.Satellite$VCV)[2, , ])
diag(autocorr(mod.Satellite$Sol)[2, , ])
##summary
summary(mod.Satellite)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 88.86726 
# 
# G-structure:  ~ID
# 
# post.mean l-95% CI u-95% CI eff.samp
# ID     0.977 0.002835    2.972     2000
# 
# R-structure:  ~units
# 
# post.mean l-95% CI u-95% CI eff.samp
# units    0.6853   0.3323    1.116     2000
# 
# Location effects: log(Satellite) ~ Desiccating 
# 
# post.mean l-95% CI u-95% CI eff.samp pMCMC   
# (Intercept)    -3.5558  -5.1470  -1.7948     2140 0.004 **
# Desiccating1   -0.5309  -1.7510   0.7145     2000 0.336
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.Satellite <- mod.Satellite$VCV[,'ID']/(mod.Satellite$VCV[,'ID']+mod.Satellite$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.Satellite)
posterior.mode(lambda.Satellite)
HPDinterval(lambda.Satellite)[2]

##fit model for Simple
mod.Simple<-MCMCglmm(log(Simple) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.Simple, file="mod.Simple.RData")

##check for convergence in posterior output
plot(mod.Simple)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.Simple$VCV)[2, , ])
diag(autocorr(mod.Simple$Sol)[2, , ])
##summary
summary(mod.Simple)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -53.35795 
# 
# G-structure:  ~ID
# 
# post.mean l-95% CI u-95% CI eff.samp
# ID    0.3737   0.1037   0.7902     2000
# 
# R-structure:  ~units
# 
# post.mean l-95% CI u-95% CI eff.samp
# units   0.00962 0.003464  0.01682     2000
# 
# Location effects: log(Simple) ~ Desiccating 
# 
# post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)    -0.2836  -1.2347   0.6222     2000 0.527
# Desiccating1    0.4180  -0.1951   1.0049     2000 0.145

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.Simple <- mod.Simple$VCV[,'ID']/(mod.Simple$VCV[,'ID']+mod.Simple$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.Simple)
posterior.mode(lambda.Simple)
HPDinterval(lambda.Simple)[2]

##fit model for Low.complexity
mod.Low.complexity<-MCMCglmm(log(Low.complexity) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                     family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                     nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.Low.complexity, file="mod.Low.complexity.RData")

##check for convergence in posterior output
plot(mod.Low.complexity)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.Low.complexity$VCV)[2, , ])
diag(autocorr(mod.Low.complexity$Sol)[2, , ])
##summary
summary(mod.Low.complexity)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -95.368 
# 
# G-structure:  ~ID
# 
# post.mean l-95% CI u-95% CI eff.samp
# ID    0.2105   0.0831   0.3974     1853
# 
# R-structure:  ~units
# 
# post.mean l-95% CI u-95% CI eff.samp
# units  0.003049 0.001368 0.005374     1827
# 
# Location effects: log(Low.complexity) ~ Desiccating 
# 
#              post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
# (Intercept)  -1.843510 -2.509065 -1.158930     2000 <5e-04 ***
# Desiccating1  0.461565 -0.003357  0.893357     2009  0.053 .  

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.Low.complexity <- mod.Low.complexity$VCV[,'ID']/(mod.Low.complexity$VCV[,'ID']+mod.Low.complexity$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.Low.complexity)
posterior.mode(lambda.Low.complexity)
HPDinterval(lambda.Low.complexity)[2]

##fit model for Unclassified
mod.Unclassified<-MCMCglmm(log(Unclassified) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                     family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                     nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.Unclassified, file="mod.Unclassified.RData")

##check for convergence in posterior output
plot(mod.Unclassified)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.Unclassified$VCV)[2, , ])
diag(autocorr(mod.Unclassified$Sol)[2, , ])
##summary
summary(mod.Unclassified)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -47.02112 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.5676   0.1962    1.086     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units   0.01089  0.00476  0.01811     2000
# 
# Location effects: log(Unclassified) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     2.9438   1.8333   4.1192     2221 <5e-04 ***
# Desiccating1    0.2585  -0.5686   0.9662     2000  0.464    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.Unclassified <- mod.Unclassified$VCV[,'ID']/(mod.Unclassified$VCV[,'ID']+mod.Unclassified$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.Unclassified)
posterior.mode(lambda.Unclassified)
HPDinterval(lambda.Unclassified)[2]

##fit model for Class I
mod.I<-MCMCglmm(log(SUM.known.I) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                nitt=420000, burnin=20000, thin=200, verbose=T)
##save model
save(mod.I, file="mod.I.RData")

##check for convergence in posterior output
plot(mod.I)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.I$VCV)[2, , ])
diag(autocorr(mod.I$Sol)[2, , ])
##summary
summary(mod.I)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -55.01899 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.4857   0.0822    1.063     1598
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units   0.00929 0.002753  0.01691     1653
# 
# Location effects: log(SUM.known.I) ~ Desiccating 
# 
#              post.mean  l-95% CI  u-95% CI eff.samp pMCMC  
# (Intercept)   0.987520  0.002477  2.100567     2000 0.063 .
# Desiccating1 -0.084265 -0.786097  0.605240     2000 0.773  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.I <- mod.I$VCV[,'ID']/(mod.I$VCV[,'ID']+mod.I$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.I)
posterior.mode(lambda.I)
HPDinterval(lambda.I)

##fit model for Class II
mod.II<-MCMCglmm(log(SUM.known.II) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                 family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                 nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.II, file="mod.II.RData")

##check for convergence in posterior output
plot(mod.II)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.II$VCV)[2, , ])
diag(autocorr(mod.II$Sol)[2, , ])
##summary
summary(mod.II)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -72.47144 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.2925   0.1038   0.5503     2253
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
  # units  0.005529 0.002521 0.009655     2000
# 
# Location effects: log(SUM.known.II) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)     0.6839  -0.2228   1.4213     2238 0.112
# Desiccating1    0.2962  -0.2018   0.8468     2215 0.267

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.II <- mod.II$VCV[,'ID']/(mod.II$VCV[,'ID']+mod.II$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.II)
posterior.mode(lambda.II)
HPDinterval(lambda.II)

##fit model for Class I and II
mod.known.all<-MCMCglmm(log(SUM.known.all) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                  family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                  nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.known.all, file="mod.known.all.RData")

##check for convergence in posterior output
plot(mod.known.all)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.known.all$VCV)[2, , ])
diag(autocorr(mod.known.all$Sol)[2, , ])
##summary
summary(mod.known.all)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -75.54045 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.2895   0.1052   0.5596     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.005051 0.002167 0.008471     2000
# 
# Location effects: log(SUM.known.all) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC    
# (Intercept)     1.5013   0.7321   2.3472     2651 0.001 ***
# Desiccating1    0.1561  -0.3541   0.6736     2162 0.545    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.known.all <- mod.known.all$VCV[,'ID']/(mod.known.all$VCV[,'ID']+mod.known.all$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.known.all)
posterior.mode(lambda.known.all)
HPDinterval(lambda.known.all)

##fit model for all repeats combined
mod.all<-MCMCglmm(log(SUM.TOTAL) ~ Desiccating, random=~ID, data=subset(df.te_load, df.te_load$Bdelloid==1),
                  family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                  nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.all, file="mod.all.RData")

##check for convergence in posterior output
plot(mod.all)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.all$VCV)[2, , ])
diag(autocorr(mod.all$Sol)[2, , ])
##summary
summary(mod.all)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -60.81562 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.3897   0.1184   0.7436     1766
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.007529 0.003365  0.01281     2000
# 
# Location effects: log(SUM.TOTAL) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     3.1887   2.3122   4.1591     2000 <5e-04 ***
# Desiccating1    0.2678  -0.3141   0.8600     2000  0.345    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##estimate lambda, the posterior probability of the phylogenetic signal
lambda.all <- mod.all$VCV[,'ID']/(mod.all$VCV[,'ID']+mod.all$VCV[,'units'])
##mean and HPD of lambda 
mean(lambda.all)
posterior.mode(lambda.all)
HPDinterval(lambda.all)

##PLOT LAMBDA
type <- c("DNA transposons","Rolling circles","PLEs","LTRs","LINEs","Satellite","Simple","Low complex","Unclassified","Class I TEs","Class II TEs","Class I and II TEs","All repeats")
mean <- c(mean(lambda.DNA), mean(lambda.RC), mean(lambda.PLE), mean(lambda.LTR), mean(lambda.LINE), mean(lambda.Satellite), mean(lambda.Simple), mean(lambda.Low.complexity), mean(lambda.Unclassified), mean(lambda.I), mean(lambda.II), mean(lambda.known.all),mean(lambda.all))
lower <- c(HPDinterval(lambda.DNA)[1],HPDinterval(lambda.RC)[1],HPDinterval(lambda.PLE)[1],HPDinterval(lambda.LTR)[1],HPDinterval(lambda.LINE)[1],HPDinterval(lambda.Satellite)[1],HPDinterval(lambda.Simple)[1],HPDinterval(lambda.Low.complexity)[1],HPDinterval(lambda.Unclassified)[1],HPDinterval(lambda.I)[1],HPDinterval(lambda.II)[1],HPDinterval(lambda.known.all)[1],HPDinterval(lambda.all)[1])
upper <- c(HPDinterval(lambda.DNA)[2],HPDinterval(lambda.RC)[2],HPDinterval(lambda.PLE)[2],HPDinterval(lambda.LTR)[2],HPDinterval(lambda.LINE)[2],HPDinterval(lambda.Satellite)[2],HPDinterval(lambda.Simple)[2],HPDinterval(lambda.Low.complexity)[2],HPDinterval(lambda.Unclassified)[2],HPDinterval(lambda.I)[2],HPDinterval(lambda.II)[2],HPDinterval(lambda.known.all)[2],HPDinterval(lambda.all)[2])
df.lambda <- data.frame(type, mean, lower, upper)
str(df.lambda)

## graphics
par(mar=c(6,3,2,1), oma=c(2,1,0,0), las=3)
par(tcl=-0.25)
par(mgp=c(2, 0.6, 0))
## bdelloid
plot(df.lambda$mean, ylim=c(0,1), xaxt="n", xlab="", ylab="Lambda (mean ± HPDi)", type="n")
rect(9.5,-0.5,13.5,1.5, col="grey95", border=NA); abline(v=9.5, lty=2)
arrows(x0=c(1:length(df.lambda$mean)),y0=df.lambda$lower, x1=c(1:length(df.lambda$upper)),y1=df.lambda$upper, code=3, angle=90, length=0.1)
points(df.lambda$mean, pch=21, bg="blue")
axis(1, at=c(1:length(df.lambda$mean)), cex.axis=1, labels=df.lambda$type);box()
