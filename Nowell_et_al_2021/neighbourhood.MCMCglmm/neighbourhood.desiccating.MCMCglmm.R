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
## drop tips not in the data
tree<-drop.tip(tree,tree$tip.label[c(grep(paste(levels(df$ID),collapse="|"),tree$tip.label,invert=T))])
## drop monogononts  
tree<-drop.tip(tree,tree$tip.label[grep("B",tree$tip.label)])
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

## TEST FOR 'NEIGHBOURHOOD' DIFFERENCES BETWEEN DESICCATING AND NONDESICCATING BDELLOIDS
## H0: there are no differences in the density of genomic features (genes, TEs or telomeric repeats) surrounding TEs of type PLE, LTR or LINE
## for desiccating versus nondesiccating bdelloids

## fit model for feature==genes & type==PLE
mod2.genes.PLE<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="gene" & df$type=="PLE" & is.bdelloid==1),
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.genes.PLE, file="mod2.genes.PLE.RData")
## check for convergence in posterior output
plot(mod2.genes.PLE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.genes.PLE$VCV)[2, , ])
diag(autocorr(mod2.genes.PLE$Sol)[2, , ])
## summary
summary(mod2.genes.PLE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 23111.57 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     109.5    45.33    193.3     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     11.35    10.89    11.85     2000
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)         4.678  -11.017   20.584     1843 0.550
# is.desiccating1    -1.258  -11.058    8.563     2000 0.789

## fit model for feature==genes & type==LTR
mod2.genes.LTR<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="gene" & df$type=="LTR" & is.bdelloid==1),
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.genes.LTR, file="mod2.genes.LTR.RData")
## check for convergence in posterior output
plot(mod2.genes.LTR)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.genes.LTR$VCV)[2, , ])
diag(autocorr(mod2.genes.LTR$Sol)[2, , ])
## summary
summary(mod2.genes.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 16851.31 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     234.9    92.45    407.7     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     11.21    10.66    11.74     2000
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                 post.mean  l-95% CI  u-95% CI eff.samp pMCMC
# (Intercept)       4.11233 -17.62408  27.71392     1747 0.707
# is.desiccating1  -0.04882 -15.88090  13.30788     1855 0.996

## fit model for feature==genes & type==LINE
mod2.genes.LINE<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="gene" & df$type=="LINE" & is.bdelloid==1),
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.genes.LINE, file="mod2.genes.LINE.RData")
## check for convergence in posterior output
plot(mod2.genes.LINE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.genes.LINE$VCV)[2, , ])
diag(autocorr(mod2.genes.LINE$Sol)[2, , ])
## summary
summary(mod2.genes.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 15804.58 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     173.4     69.4    311.1     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     13.99    13.26    14.67     2000
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)         5.497  -13.457   25.962     2000 0.565
# is.desiccating1    -1.444  -13.754   10.497     2152 0.831

##
##
## fit model for feature==TEs & type==PLE
mod2.TEs.PLE<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="TE" & df$type=="PLE" & is.bdelloid==1),
                      family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                      nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.TEs.PLE, file="mod2.TEs.PLE.RData")
## check for convergence in posterior output
plot(mod2.TEs.PLE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.TEs.PLE$VCV)[2, , ])
diag(autocorr(mod2.TEs.PLE$Sol)[2, , ])
## summary
summary(mod2.TEs.PLE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 13591.98 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     97.24    52.26    159.7     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     1.292    1.237    1.347     2000
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)        6.8151  -7.4972  20.5427     2000 0.346
# is.desiccating1    0.2125  -8.7502   9.7142     1775 0.961

## fit model for feature==TEs & type==LTR
mod2.TEs.LTR<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="TE" & df$type=="LTR" & is.bdelloid==1),
                      family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                      nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.TEs.LTR, file="mod2.TEs.LTR.RData")
## check for convergence in posterior output
plot(mod2.TEs.LTR)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.TEs.LTR$VCV)[2, , ])
diag(autocorr(mod2.TEs.LTR$Sol)[2, , ])
## summary
summary(mod2.TEs.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 10700.7 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     146.8    76.06    239.6     2004
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units      1.64    1.562    1.721     2064
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                  post.mean   l-95% CI   u-95% CI eff.samp pMCMC
# (Intercept)       6.448664 -12.810510  23.300584     2000 0.471
# is.desiccating1  -0.007022 -11.613784  10.991086     2244 0.989

## fit model for feature==TEs & type==LINE
mod2.TEs.LINE<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="TE" & df$type=="LINE" & is.bdelloid==1),
                       family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                       nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.TEs.LINE, file="mod2.TEs.LINE.RData")
## check for convergence in posterior output
plot(mod2.TEs.LINE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.TEs.LINE$VCV)[2, , ])
diag(autocorr(mod2.TEs.LINE$Sol)[2, , ])
## summary
summary(mod2.TEs.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 8923.535 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     111.1    60.11    180.3     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     1.284     1.22     1.35     2021
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)        6.7736  -9.0250  22.0109     2000 0.412
# is.desiccating1   -0.2032 -10.1552  10.3523     2000 0.961

##
##
## fit model for feature==telos & type==PLE
mod2.telos.PLE<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="telo" & df$type=="PLE" & is.bdelloid==1),
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.telos.PLE, file="mod2.telos.PLE.RData")
## check for convergence in posterior output
plot(mod2.telos.PLE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.telos.PLE$VCV)[2, , ])
diag(autocorr(mod2.telos.PLE$Sol)[2, , ])
## summary
summary(mod2.telos.PLE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 17832.85 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     14.21    4.491    27.12     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     3.407    3.268    3.548     1715
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)        0.6535  -5.5091   6.1954     1698 0.813
# is.desiccating1    1.2265  -2.3931   5.0568     2000 0.482

## fit model for feature==telos & type==LTR
mod2.telos.LTR<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="telo" & df$type=="LTR" & is.bdelloid==1),
                        family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                        nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.telos.LTR, file="mod2.telos.LTR.RData")
## check for convergence in posterior output
plot(mod2.telos.LTR)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.telos.LTR$VCV)[2, , ])
diag(autocorr(mod2.telos.LTR$Sol)[2, , ])
## summary
summary(mod2.telos.LTR)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 12056.78 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     20.44     5.64    38.24     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     2.511    2.398    2.642     1807
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)        0.2393  -6.6311   6.6129     2134 0.955
# is.desiccating1    1.3398  -3.3628   5.3009     2000 0.505

## fit model for feature==telos & type==LINE
mod2.telos.LINE<-MCMCglmm(log(span+1) ~ is.desiccating, random=~ID, data=subset(df, df$feature=="telo" & df$type=="LINE" & is.bdelloid==1),
                         family="gaussian", ginverse=list(ID=inv.phylo$Ainv), prior=prior,
                         nitt=420000, burnin=20000, thin=200, verbose=T)
## save model
save(mod2.telos.LINE, file="mod2.telos.LINE.RData")
## check for convergence in posterior output
plot(mod2.telos.LINE)
## plot and check convergence (should be << 0.2)
diag(autocorr(mod2.telos.LINE$VCV)[2, , ])
diag(autocorr(mod2.telos.LINE$Sol)[2, , ])
## summary
summary(mod2.telos.LINE)
# Iterations = 20001:419801
# Thinning interval  = 200
# Sample size  = 2000 
# 
# DIC: 11014.42 
# 
# G-structure:  ~ID
# 
#    post.mean l-95% CI u-95% CI eff.samp
# ID     21.91    7.702       43     2000
# 
# R-structure:  ~units
# 
#       post.mean l-95% CI u-95% CI eff.samp
# units     2.656    2.525    2.796     2000
# 
# Location effects: log(span + 1) ~ is.desiccating 
# 
#                 post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)        0.3691  -6.9754   7.1610     2000 0.910
# is.desiccating1    1.6973  -2.7362   6.1875     2000 0.441