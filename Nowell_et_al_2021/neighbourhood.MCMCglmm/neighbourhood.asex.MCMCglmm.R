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
## have a look
str(df)

## READ IN TREE
tree<-read.tree("rotifers_tree.newick")
## drop data not in tips
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

## TEST FOR 'NEIGHBOURHOOD' DIFFERENCES BETWEEN BDELLOIDS & MONOGONONTS
## H0: there are no differences in the density of genomic features (genes, TEs or telomeric repeats) surrounding TEs of type PLE, LTR or LINE
## for monogononts versus bdelloids

## fit model for feature==genes & type==LTR
mod1.genes.LTR<-MCMCglmm(log(span+1) ~ is.bdelloid, random=~ID, data=subset(df, df$feature=="gene" & df$type=="LTR"),
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.genes.LTR, file="mod1.genes.LTR.RData")
## check for convergence in posterior output
plot(mod1.genes.LTR)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.genes.LTR$VCV)[2, , ])
diag(autocorr(mod1.genes.LTR$Sol)[2, , ])
## summary
summary(mod1.genes.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 16900.73 
# 
# G-structure:  ~ID
# 
#    spost.mean l-95% CI u-95% CI eff.samp
# ID     645.6    270.4     1146     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     11.21     10.7    11.78     2000
# 
# Location effects: log(span + 1) ~ is.bdelloid 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)      5.628  -42.127   59.418     2124 0.851
# is.bdelloid1    -2.117  -69.036   67.469     2401 0.945

## fit model for feature==genes & type==LINE
mod1.genes.LINE<-MCMCglmm(log(span+1) ~ is.bdelloid, random=~ID, data=subset(df, df$feature=="gene" & df$type=="LINE"),
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.genes.LINE, file="mod1.genes.LINE.RData")
## check for convergence in posterior output
plot(mod1.genes.LINE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.genes.LINE$VCV)[2, , ])
diag(autocorr(mod1.genes.LINE$Sol)[2, , ])
## summary
summary(mod1.genes.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 18537.49 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID       448    184.2    770.6     2453
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     13.94    13.27     14.6     2000
# 
# Location effects: log(span + 1) ~ is.bdelloid 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)      5.921  -32.337   46.343     2000 0.771
# is.bdelloid1    -1.519  -57.411   49.961     1754 0.969

##
##
## fit model for feature==TEs & type==LTR
mod1.TEs.LTR<-MCMCglmm(log(span+1) ~ is.bdelloid, random=~ID, data=subset(df, df$feature=="TE" & df$type=="LTR"),
                      family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                      nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.TEs.LTR, file="mod1.TEs.LTR.RData")
## check for convergence in posterior output
plot(mod1.TEs.LTR)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.TEs.LTR$VCV)[2, , ])
diag(autocorr(mod1.TEs.LTR$Sol)[2, , ])
## summary
summary(mod1.TEs.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 10727.62 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     388.7    205.6    636.4     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units      1.64    1.556    1.718     2000
# 
# Location effects: log(span + 1) ~ is.bdelloid 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)     6.9254 -31.7001  44.8291     1521 0.731
# is.bdelloid1    0.2385 -50.6834  53.0974     2000 0.973

## fit model for feature==TEs & type==LINE
mod1.TEs.LINE<-MCMCglmm(log(span+1) ~ is.bdelloid, random=~ID, data=subset(df, df$feature=="TE" & df$type=="LINE"),
                       family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                       nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.TEs.LINE, file="mod1.TEs.LINE.RData")
## check for convergence in posterior output
plot(mod1.TEs.LINE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.TEs.LINE$VCV)[2, , ])
diag(autocorr(mod1.TEs.LINE$Sol)[2, , ])
## summary
summary(mod1.TEs.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 11308.04 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     273.8    139.6    438.1     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     1.643    1.564    1.722     2000
# 
# Location effects: log(span + 1) ~ is.bdelloid 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)     7.1645 -21.2234  36.6993     2000  0.64
# is.bdelloid1    0.2393 -41.4227  38.7547     2000  0.97

##
##
## fit model for feature==telos & type==LTR
mod1.telos.LTR<-MCMCglmm(log(span+1) ~ is.bdelloid, random=~ID, data=subset(df, df$feature=="telo" & df$type=="LTR"),
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.telos.LTR, file="mod1.telos.LTR.RData")
## check for convergence in posterior output
plot(mod1.telos.LTR)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.telos.LTR$VCV)[2, , ])
diag(autocorr(mod1.telos.LTR$Sol)[2, , ])
## summary
summary(mod1.telos.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 12087.04 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     62.11    17.29    118.1     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     2.507    2.392    2.633     2147
# 
# Location effects: log(span + 1) ~ is.bdelloid 
# 
#              post.mean  l-95% CI  u-95% CI eff.samp pMCMC
# (Intercept)    1.26986 -15.02951  16.09221     2000 0.871
# is.bdelloid1  -0.06356 -19.24534  21.61893     2000 0.992

## fit model for feature==telos & type==LINE
mod1.telos.LINE<-MCMCglmm(log(span+1) ~ is.bdelloid, random=~ID, data=subset(df, df$feature=="telo" & df$type=="LINE"),
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod1.telos.LINE, file="mod1.telos.LINE.RData")
## check for convergence in posterior output
plot(mod1.telos.LINE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod1.telos.LINE$VCV)[2, , ])
diag(autocorr(mod1.telos.LINE$Sol)[2, , ])
## summary
summary(mod1.telos.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 12731.77 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     60.78    16.57    113.7     1753
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     2.504    2.388    2.625     1973
# 
# Location effects: log(span + 1) ~ is.bdelloid 
# 
#              post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)     0.9758 -13.9522  15.1782     2000 0.878
# is.bdelloid1    0.8695 -18.9849  20.7808     2092 0.919