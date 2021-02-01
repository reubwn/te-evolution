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
##for just bdelloid df
str(df.te_load.bdelloids)
head(df.te_load.bdelloids)
levels(df.te_load.bdelloids$ID) ##IDs

##READ IN PROTOSTOME TREE
tree <- read.tree("protostome_tree.newick")
##root tree
tree<-midpoint(tree)
##plot tree
plot(tree, underscore=T, cex=0.8)

##make ultrametric and read back in to correct class
tree<-chronos(tree, lambda=1)
tree<-read.tree(text=write.tree(tree))
tree$node.label<-NULL
##convert tree into covariance matrix for mcmcglmm
inv.phylo<-inverseA(tree, nodes="TIPS", scale=TRUE)
##default prior
prior<-list(G=list(G1=list(V=1, nu=0.02)), R=list(V=1, nu=0.02))

##TESTS FOR EFFECT OF DESICCATION ACROSS PROTOSTOME TREE

##fit model for DNA transposons
mod.DNA.proto<-MCMCglmm(log(DNA.transposons) ~ Desiccating, random=~ID, data=df.te_load,
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.DNA.proto, file="mod.DNA.proto.RData")

##check for convergence in posterior output
plot(mod.DNA.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.DNA.proto$VCV)[2, , ])
diag(autocorr(mod.DNA.proto$Sol)[2, , ])
##summary
summary(mod.DNA.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -99.68614 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     1.276   0.7945    1.803     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.006516 0.001536   0.0138     2000
# 
# Location effects: log(DNA.transposons) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC  
# (Intercept)     0.9525   0.1777   1.7024     2000 0.011 *
# Desiccating1   -0.2264  -0.7037   0.2411     2000 0.351  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for RC
mod.RC.proto<-MCMCglmm(log(Rolling.circles) ~ Desiccating, random=~ID, data=df.te_load,
                       family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                       nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.RC.proto, file="mod.RC.proto.RData")

##check for convergence in posterior output
plot(mod.RC.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.RC.proto$VCV)[2, , ])
diag(autocorr(mod.RC.proto$Sol)[2, , ])
##summary
summary(mod.RC.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -13.76622 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     2.676    1.519    4.069     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units   0.02884 0.004958  0.06007     2000
# 
# Location effects: log(Rolling.circles) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)    -2.3779  -3.5119  -1.3019     1917 <5e-04 ***
# Desiccating1    0.3869  -0.3020   1.0932     2000  0.268    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for PLEs
mod.PLE.proto<-MCMCglmm(log(Penelope) ~ Desiccating, random=~ID, data=df.te_load,
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.PLE.proto, file="mod.PLE.proto.RData")

##check for convergence in posterior output
plot(mod.PLE.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.PLE.proto$VCV)[2, , ])
diag(autocorr(mod.PLE.proto$Sol)[2, , ])
##summary
summary(mod.PLE.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -45.18097 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     4.724    2.993     6.78     2171
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units    0.0174 0.002176   0.0394     2000
# 
# Location effects: log(Penelope) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC  
# (Intercept)    -1.9711  -3.4703  -0.5299     2000 0.011 *
# Desiccating1   -0.6759  -1.5825   0.2246     2000 0.152  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for LTRs
mod.LTR.proto<-MCMCglmm(log(LTRs) ~ Desiccating, random=~ID, data=df.te_load,
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=T)
##save model
save(mod.LTR.proto, file="mod.LTR.proto.RData")

##check for convergence in posterior output
plot(mod.LTR.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.LTR.proto$VCV)[2, , ])
diag(autocorr(mod.LTR.proto$Sol)[2, , ])
##summary
summary(mod.LTR.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 20.5099 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     1.531   0.5223    2.581     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units   0.05689 0.009506   0.1201     1861
# 
# Location effects: log(LTRs) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC  
# (Intercept)     0.5207  -0.3211   1.3216     2000 0.219  
# Desiccating1   -0.6330  -1.2351  -0.1082     2000 0.024 * [NS after correction for multiple testing]
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for LINEs
mod.LINE.proto<-MCMCglmm(log(LINEs) ~ Desiccating, random=~ID, data=df.te_load,
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=T)
##save model
save(mod.LINE.proto, file="mod.LINE.proto.RData")

##check for convergence in posterior output
plot(mod.LINE.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.LINE.proto$VCV)[2, , ])
diag(autocorr(mod.LINE.proto$Sol)[2, , ])
##summary
summary(mod.LINE.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -101.7482 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     1.858    1.174    2.624     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.005966 0.001527  0.01281     2000
# 
# Location effects: log(LINEs) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)     0.6933  -0.1600   1.6573     2000 0.137
# Desiccating1   -0.4223  -0.9636   0.1849     2000 0.164

##fit model for SINEs
mod.SINE.proto<-MCMCglmm(log(SINEs) ~ Desiccating, random=~ID, data=df.te_load,
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.SINE.proto, file="mod.SINE.proto.RData")

##check for convergence in posterior output
plot(mod.SINE.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.SINE.proto$VCV)[2, , ])
diag(autocorr(mod.SINE.proto$Sol)[2, , ])
##summary
summary(mod.SINE.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 146.4534 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     7.577    3.645    12.49     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units    0.4317   0.1543   0.7362     2000
# 
# Location effects: log(SINEs) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)   -4.75022 -6.82227 -3.01611     2006 <5e-04 ***
# Desiccating1   0.03396 -1.36058  1.30050     2000  0.954    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for Satellite
mod.Satellite.proto<-MCMCglmm(log(Satellite) ~ Desiccating, random=~ID, data=df.te_load,
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.Satellite.proto, file="mod.Satellite.proto.RData")

##check for convergence in posterior output
plot(mod.Satellite.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.Satellite.proto$VCV)[2, , ])
diag(autocorr(mod.Satellite.proto$Sol)[2, , ])
##summary
summary(mod.Satellite.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 147.8226 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     3.117   0.9857     5.74     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units    0.4933   0.2117   0.8238     2200
# 
# Location effects: log(Satellite) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)   -3.81071 -4.97003 -2.60720     2000 <5e-04 ***
# Desiccating1  -0.89365 -1.89301 -0.01265     2000  0.056 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for Simple
mod.Simple.proto<-MCMCglmm(log(Simple) ~ Desiccating, random=~ID, data=df.te_load,
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.Simple.proto, file="mod.Simple.proto.RData")

##check for convergence in posterior output
plot(mod.Simple.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.Simple.proto$VCV)[2, , ])
diag(autocorr(mod.Simple.proto$Sol)[2, , ])
##summary
summary(mod.Simple.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -90.21578 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     1.103   0.6644    1.586     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.007415 0.002183  0.01524     1748
# 
# Location effects: log(Simple) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)     0.5462  -0.0972   1.2259     2000 0.104
# Desiccating1   -0.2345  -0.6680   0.2401     2000 0.309

##fit model for Low.complexitys
mod.Low.complexity.proto<-MCMCglmm(log(Low.complexity) ~ Desiccating, random=~ID, data=df.te_load,
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.Low.complexity.proto, file="mod.Low.complexity.proto.RData")

##check for convergence in posterior output
plot(mod.Low.complexity.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.Low.complexity.proto$VCV)[2, , ])
diag(autocorr(mod.Low.complexity.proto$Sol)[2, , ])
##summary
summary(mod.Low.complexity.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -118.6479 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.6651   0.4089   0.9321     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.004477 0.001539   0.0085     2000
# 
# Location effects: log(Low.complexity) ~ Desiccating 
# 
#              post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
# (Intercept)  -1.226297 -1.772522 -0.693046     1854 <5e-04 ***
# Desiccating1  0.001419 -0.383256  0.335347     2000  0.989    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for Unclassifieds
mod.Unclassified.proto<-MCMCglmm(log(Unclassified) ~ Desiccating, random=~ID, data=df.te_load,
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.Unclassified.proto, file="mod.Unclassified.proto.RData")

##check for convergence in posterior output
plot(mod.Unclassified.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.Unclassified.proto$VCV)[2, , ])
diag(autocorr(mod.Unclassified.proto$Sol)[2, , ])
##summary
summary(mod.Unclassified.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -59.17599 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     1.223   0.7652    1.789     2277
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units   0.01265 0.002874  0.02494     2000
# 
# Location effects: log(Unclassified) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     2.6277   1.9036   3.3467     2000 <5e-04 ***
# Desiccating1    0.3283  -0.1333   0.8272     2000  0.175    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for Class I
mod.I.proto<-MCMCglmm(log(SUM.known.I) ~ Desiccating, random=~ID, data=df.te_load,
                      family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                      nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.I.proto, file="mod.I.proto.RData")

##check for convergence in posterior output
plot(mod.I.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.I.proto$VCV)[2, , ])
diag(autocorr(mod.I.proto$Sol)[2, , ])
##summary
summary(mod.I.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -76.94774 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     1.237   0.7799    1.804     2139
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.009452 0.002229   0.0192     2000
# 
# Location effects: log(SUM.known.I) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)    1.70537  0.94624  2.41743     2000 <5e-04 ***
# Desiccating1  -0.46331 -0.92704  0.03153     2308  0.064 .  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for Class II
mod.II.proto<-MCMCglmm(log(SUM.known.II) ~ Desiccating, random=~ID, data=df.te_load,
                       family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                       nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.II.proto, file="mod.II.proto.RData")

##check for convergence in posterior output
plot(mod.II.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.II.proto$VCV)[2, , ])
diag(autocorr(mod.II.proto$Sol)[2, , ])
##summary
summary(mod.II.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -98.1002 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     1.057    0.666    1.485     1842
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.006648 0.001752  0.01345     2317
# 
# Location effects: log(SUM.known.II) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC   
# (Intercept)    1.02542  0.33036  1.70694     1775 0.008 **
# Desiccating1  -0.06847 -0.56145  0.32964     2000 0.776   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for Class I and II
mod.all.known.proto<-MCMCglmm(log(SUM.known.all) ~ Desiccating, random=~ID, data=df.te_load,
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.all.known.proto, file="mod.all.known.proto.RData")

##check for convergence in posterior output
plot(mod.all.known.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.all.known.proto$VCV)[2, , ])
diag(autocorr(mod.all.known.proto$Sol)[2, , ])
##summary
summary(mod.all.known.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -97.95909 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.9152   0.5529    1.299     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.006599 0.001949  0.01406     2000
# 
# Location effects: log(SUM.known.all) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     2.3978   1.7737   3.0399     2000 <5e-04 ***
# Desiccating1   -0.2733  -0.6470   0.1477     1807  0.183    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##fit model for all repeats
mod.all.proto<-MCMCglmm(log(SUM.TOTAL) ~ Desiccating, random=~ID, data=df.te_load,
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=F)
##save model
save(mod.all.proto, file="mod.all.proto.RData")

##check for convergence in posterior output
plot(mod.all.proto)
##plot and check convergence (should be << 0.2)
diag(autocorr(mod.all.proto$VCV)[2, , ])
diag(autocorr(mod.all.proto$Sol)[2, , ])
##summary
summary(mod.all.proto)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: -80.48099 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID    0.6002   0.3474   0.8864     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units  0.008866 0.002948  0.01729     2006
# 
# Location effects: log(SUM.TOTAL) ~ Desiccating 
# 
#              post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)    3.48892  2.94007  4.00164     2000 <5e-04 ***
# Desiccating1   0.07506 -0.25569  0.42251     2000  0.683    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
