setwd("./")

## libraries
library(phangorn)
library(MCMCglmm)

## READ IN DATA
## 'span' represents the number of bp assigned to 'feature' in a 50 kb window drawn around BUSCO/TE of type 'type'
df<-read.table("neighbourhood.25kb.txt", col.names=c("ID","species","chrom","start","end","window","name","type","feature","span"), stringsAsFactors=T)

## add asexual and desiccation factors
df$is.bdelloid <- as.factor(ifelse(df$ID %in% c('Bc_PSC1','Bp_HYR1'), 0, 1))
df$is.desiccating <- as.factor(ifelse(df$ID %in% c('Bc_PSC1','Bp_HYR1','Rg_MAG1','Rg_MAG2','Rg_MAG3','Rg_RM9','Rg_RM15','Rs_AK11','Rs_AK15','Rs_AK16','Rs_AK27','Rs_RS1','Rc_Rc2018'), 0, 1))
## add BUSCO as factor
df$is.BUSCO <- as.factor(ifelse(df$type=="BUSCO", 1, 0))
## calculate density
df$density <- df$span/df$window
## have a look
str(df)

## READ IN TREE
tree<-read.tree("rotifers_tree.newick")
## drop tips not in the data
tree<-drop.tip(tree,tree$tip.label[c(grep(paste(levels(df$ID),collapse="|"),tree$tip.label,invert=T))])
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

## TEST FOR 'NEIGHBOURHOOD' DIFFERENCES BETWEEN BUSCO GENEs and TEs

## fit model for feature==genes
mod3.genes<-MCMCglmm(log(span+1) ~ is.BUSCO, random=~ID, data=subset(df, df$feature=="gene"),
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod3.genes, file="mod1.genes.PLE.RData")
## check for convergence in posterior output
plot(mod3.genes)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod3.genes$VCV)[2, , ])
diag(autocorr(mod3.genes$Sol)[2, , ])
## summary
summary(mod3.genes)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 191474 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     184.3     94.7    289.1     1723
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units      4.51     4.45     4.57     2000
# 
# Location effects: log(span + 1) ~ is.BUSCO 
# 
#             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     3.505  -13.800   21.270     2000  0.699    
# is.BUSCO1       5.808    5.760    5.853     2000 <5e-04 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## fit model for feature==genes
mod3.TEs<-MCMCglmm(log(span+1) ~ is.BUSCO, random=~ID, data=subset(df, df$feature=="TE"),
                     family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                     nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod3.TEs, file="mod1.TEs.PLE.RData")
## check for convergence in posterior output
plot(mod3.TEs)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod3.TEs$VCV)[2, , ])
diag(autocorr(mod3.TEs$Sol)[2, , ])
## summary
summary(mod3.TEs)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 196388.7 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     370.8      190    586.6     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     5.042    4.975    5.105     2000
# 
# Location effects: log(span + 1) ~ is.BUSCO 
# 
#             post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     8.190  -16.301   33.314     2250  0.498    
# is.BUSCO1      -3.110   -3.157   -3.061     2000 <5e-04 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## fit model for feature==genes
mod3.telo<-MCMCglmm(log(span+1) ~ is.BUSCO, random=~ID, data=subset(df, df$feature=="telo"),
                     family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                     nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod3.telo, file="mod1.telo.PLE.RData")
## check for convergence in posterior output
plot(mod3.telo)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod3.telo$VCV)[2, , ])
diag(autocorr(mod3.telo$Sol)[2, , ])
## summary
summary(mod3.telo)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 158715.8 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     148.7    83.49    236.9     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     2.144    2.116    2.172     2075
# 
# Location effects: log(span + 1) ~ is.BUSCO 
# 
# post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)    1.3461 -14.2730  15.6969     2000  0.857    
# is.BUSCO1      0.4495   0.4180   0.4843     2000 <5e-04 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1