setwd("./")

## libraries
library(phytools)
library(RColorBrewer)
library(cowplot)

## get data
rotifers.simple<-read.table("repeats_tab.simple.txt", head=T)
rotifers.relative<-read.table("repeats_tab.relative.txt", head=T)
tree<-read.tree("SUPERMATRIX.aln.mod_Nov2020.fasta.treefile.rooted.newick")
plot(tree, font=1, underscore=T, cex=0.8)
add.scale.bar(x=0, y=25, cex=0.8)
## reorder rows to correspond to tip labels
rotifers.simple <- rotifers.simple[order(factor(rotifers.simple$ID, levels=tree$tip.label)),]
rotifers.relative <- rotifers.relative[order(factor(rotifers.relative$ID, levels=tree$tip.label)),]

## protostome data
protostome.simple<-read.table("protostomes/protostome_tab.simple.txt", head=T)
protostome.relative<-read.table("protostomes/protostome_tab.relative.txt", head=T)
## reverse order rows for 'beside=T'
protostome.simple <- protostome.simple[order(factor(rev(protostome.simple$Species), levels=protostome.simple$Species)),]
protostome.relative <- protostome.relative[order(factor(rev(protostome.relative$Species), levels=protostome.relative$Species)),]

## set up colours for simple
cols.simple<-c("#F34C37","#FF2987","#FFCE6EFF","#3EA85A","#4594C7","#AE217E","grey50","grey60","grey70","grey80")

## set up colours for rotifers relative
dna.cols=colorRampPalette(rev(brewer.pal(n=9,name="Reds")))(length(grep("DNA", colnames(rotifers.relative))))
ltr.cols=colorRampPalette(rev(brewer.pal(n=9,name="Greens")))(length(grep("LTR", colnames(rotifers.relative))))
lines.cols=colorRampPalette(rev(brewer.pal(n=9,name="Blues")))(length(grep("LINE", colnames(rotifers.relative))))
## put into vector
cols.relative<-c(dna.cols,"#FF2987","#FFCE6EFF",ltr.cols,lines.cols,"#AE217E")

## set up colours for protostomes relative
protostome.relative.dna.cols=colorRampPalette(rev(brewer.pal(n=9,name="Reds")))(length(grep("DNA", colnames(protostome.relative))))
protostome.relative.ltr.cols=colorRampPalette(rev(brewer.pal(n=9,name="Greens")))(length(grep("LTR", colnames(protostome.relative))))
protostome.relative.lines.cols=colorRampPalette(rev(brewer.pal(n=9,name="Blues")))(length(grep("LINE", colnames(protostome.relative))))
## put into vector
cols.protostome.relative<-c(protostome.relative.dna.cols,"#FF2987","#FFCE6EFF",protostome.relative.ltr.cols,protostome.relative.lines.cols,"#AE217E")

## desiccation tip cols
des<-"#FFB380"
non<-"#AACCFF"

#######################
## multiplot with tree
#######################

p1<-~{
  par(mfrow=c(1,3))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(mar=c(0,0,4,1),las=1)
  
  ## plot tree
  plot(tree, font=1, cex=1.2, underscore=T, show.tip.label=F, x.lim=c(0,1.1))
  add.scale.bar(x=0, y=25, cex=1)
  ## plot legend
  legend("bottomleft", y.intersp=0.9,
         c("Bdelloid (desiccating)","Bdelloid (nondesiccating)","Monogonont (nondesiccating)"),
         pch=c(21,21,22), pt.bg=c(des,non,non), bg=NA, box.col=NA, cex=1, pt.cex=1.8
  )
  ## tip symbols
  tiplabels(tip=grep("B",tree$tip.label), pch=22, cex=1.8, bg=non, offset=0.02)
  tiplabels(tip=grep("Ar",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("As",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Av",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Dc",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Rg",tree$tip.label), pch=21, cex=1.8, bg=non, offset=0.02)
  tiplabels(tip=grep("Rs",tree$tip.label), pch=21, cex=1.8, bg=non, offset=0.02)
  tiplabels(tip=grep("Rc",tree$tip.label), pch=21, cex=1.8, bg=non, offset=0.02)
  tiplabels(tip=grep("Rd",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Rw",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Rp",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  ## plot simple
  barplot(t(as.matrix(rotifers.simple[,-1])), xlim=c(0,50),
    horiz=T, col=cols.simple, border=cols.simple, ylim=c(0.7,50), xaxt="n", 
    names.arg=tree$tip.label, cex.names=1.2)
  ## top axis
  axis(3, cex.axis=1.2, cex.lab=1.4, line=-0.5)
  mtext("Proportion of genome (%)", side=3, line=1.5, cex=1)
  ## plot legend
  legend(30,37, title="Repeat type", y.intersp=0.9,
    c("DNA","RC","PLE","LTR","LINE","SINE","Satellite","Simple","Low complex","Unclassified"),
    fill=cols.simple,border=cols.simple, bg="grey95", box.col="grey95", cex=1
  )
  ## plot relative
  barplot(t(as.matrix(rotifers.relative[,-1])),
          horiz=T, col=cols.relative, border=cols.relative, ylim=c(0.7,50), xaxt="n", yaxt="n")
  ## top axis
  axis(3, cex.axis=1.2, cex.lab=1.4, line=-0.5)
  mtext("Relative proportion of total TEs (%)", side=3, line=1.5, cex=1)
}

## plot with same ylim limits
# p2<-~{
#   par(mfrow=c(1,3))
#   par(tcl=-0.25)
#   par(mgp=c(2, 0.6, 0))
#   par(mar=c(4,0,1,1),las=1)
#   
#   ## blank plot
#   plot(1, type="n", axes=F, xlab="", ylab="", main="")
#   ## plot simple
#   barplot(t(as.matrix(protostome.simple[,-1])), xlim=c(0,50),
#           horiz=T, col=cols.simple, border=cols.simple, ylim=c(0.7,50), xlab="Proportion of genome (%)", cex.axis=1.2, cex.lab=1.4, 
#           names.arg=rev(protostome.simple$Species), cex.names=1.2)
#   ## plot relative
#   barplot(t(as.matrix(protostome.relative[,-1])),
#           horiz=T, col=cols.protostome.relative, border=cols.protostome.relative, ylim=c(0.7,50), xlab="Relative proportion of total TEs (%)", cex.axis=1.2, cex.lab=1.4)
# }

p2<-~{
  par(mfrow=c(1,3))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(mar=c(4,0,1,1),las=1)
  
  ## blank plot
  plot(1, type="n", axes=F, xlab="", ylab="", main="")
  ## plot simple
  b<-barplot(t(as.matrix(protostome.simple[,-1])), xlim=c(0,70),
          horiz=T, col=cols.simple, border=cols.simple, xaxt="n", 
          names.arg=protostome.simple$Species, cex.names=1.2, font=3)
  ## bottom axis
  axis(1, cex.axis=1.2, cex.lab=1.4)
  mtext("Proportion of genome (%)", side=1, line=2, cex=1)
  ## mark desiccating species
  text(colSums(t(as.matrix(protostome.simple[,-1])))+1, b, c(rep("",11),"D","ND","","D","ND",""), cex=1, font=2, adj=0, col=c(rep(NA,11),des,non,NA,des,non,NA))
  ## plot relative
  barplot(t(as.matrix(protostome.relative[,-1])),
          horiz=T, col=cols.protostome.relative, border=cols.protostome.relative, xaxt="n", yaxt="n")
  ## bottom axis
  axis(1, cex.axis=1.2, cex.lab=1.4)
  mtext("Relative proportion of total TEs (%)", side=1, line=2, cex=1)
}

## make plot
plot_grid(p1,p2, nrow=2, rel_heights=c(2.2,1), labels="AUTO")


#############
## KNOWN ONLY
#############
##

p1<-~{
  par(mfrow=c(1,3))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(mar=c(0,0,4,1),las=1)
  
  ## plot tree
  plot(tree, font=1, cex=1.2, underscore=T, show.tip.label=F, x.lim=c(0,1.16))
  ## plot legend
  legend("bottomleft", y.intersp=0.9,
         c("Bdelloid (desiccating)","Bdelloid (nondesiccating)","Monogonont (nondesiccating)"),
         pch=c(21,21,22), pt.bg=c(des,non,non), bg=NA, box.col=NA, cex=1, pt.cex=1.8
  )
  ## tip symbols
  tiplabels(tip=grep("B",tree$tip.label), pch=22, cex=1.8, bg=non, offset=0.02)
  tiplabels(tip=grep("Ar",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("As",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Av",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Dc",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Rg",tree$tip.label), pch=21, cex=1.8, bg=non, offset=0.02)
  tiplabels(tip=grep("Rs",tree$tip.label), pch=21, cex=1.8, bg=non, offset=0.02)
  tiplabels(tip=grep("Rc",tree$tip.label), pch=21, cex=1.8, bg=non, offset=0.02)
  tiplabels(tip=grep("Rd",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Rw",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  tiplabels(tip=grep("Rp",tree$tip.label), pch=21, cex=1.8, bg=des, offset=0.02)
  ## plot simple
  barplot(t(as.matrix(rotifers.simple[,-c(1,8:11)])), xlim=c(0,12),
          horiz=T, col=cols.simple, border=cols.simple, ylim=c(0.7,50), xaxt="n", 
          names.arg=tree$tip.label, cex.names=1.2)
  ## top axis
  axis(3, cex.axis=1.2, cex.lab=1.4, line=-0.5)
  mtext("Proportion of genome (%)", side=3, line=1.5, cex=1)
  ## plot legend
  legend(30,37, title="Repeat type", y.intersp=0.9,
         c("DNA","RC","PLE","LTR","LINE","SINE","Satellite","Simple","Low complex","Unclassified"),
         fill=cols.simple,border=cols.simple, bg="grey95", box.col="grey95", cex=1
  )
  ## plot relative
  barplot(t(as.matrix(rotifers.relative[,-1])),
          horiz=T, col=cols.relative, border=cols.relative, ylim=c(0.7,50), xaxt="n", yaxt="n")
  ## top axis
  axis(3, cex.axis=1.2, cex.lab=1.4, line=-0.5)
  mtext("Relative proportion of total TEs (%)", side=3, line=1.5, cex=1)
}

p2<-~{
  par(mfrow=c(1,3))
  par(tcl=-0.25)
  par(mgp=c(2,0.6,0))
  par(mar=c(4,0,1,1),las=1)
  
  ## blank plot
  plot(1, type="n", axes=F, xlab="", ylab="", main="")
  ## plot simple
  b<-barplot(t(as.matrix(protostome.simple[,-c(1,8:11)])), xlim=c(0,55),
          horiz=T, col=cols.simple, border=cols.simple, xaxt="n", 
          names.arg=protostome.simple$Species, cex.names=1.2)
  ## bottom axis
  axis(1, cex.axis=1.2, cex.lab=1.4)
  mtext("Proportion of genome (%)", side=1, line=2, cex=1)
  ## mark desiccating species
  text(colSums(t(as.matrix(protostome.simple[,-c(1,8:11)])))+1, b, c(rep("",11),"D","ND","","D","ND",""), cex=1, font=2, adj=0, col=c(rep(NA,11),des,non,NA,des,non,NA))
  ## plot relative
  barplot(t(as.matrix(protostome.relative[,-1])),
          horiz=T, col=cols.protostome.relative, border=cols.protostome.relative, xaxt="n", yaxt="n")
  ## bottom axis
  axis(1, cex.axis=1.2, cex.lab=1.4)
  mtext("Relative proportion of total TEs (%)", side=1, line=2, cex=1)
}

## make plot
plot_grid(p1,p2, nrow=2, rel_heights=c(2.2,1), labels="auto")
