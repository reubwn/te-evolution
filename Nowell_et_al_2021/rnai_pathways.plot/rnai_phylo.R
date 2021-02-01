setwd("./")

## libraries
library(phytools)
library(RColorBrewer)
library(cowplot)
## cols
As.col<-"#66C2A3"
Rs.col<-"#FA8C61"
Rg.col<-"#8C9ECA"
Rd.col<-"#E68AC2"
Rw.col<-"#A6D751"
Rp.col<-"#FFD92E"
fake.col<-rgb(0,0,0, max=255, alpha=0)
par(family="Ubuntu Light", ps=10, las=1)

## read trees
tree.ago<-read.tree("argonaute.newick")
tree.rdrp<-read.tree("rdrp.newick")
tree.dicer<-read.tree("dicer.newick")
tab<-read.table("species_avg.table", head=T)
## cols
ago.col<-brewer.pal(3,"Set1")[1]
dicer.col<-brewer.pal(3,"Set1")[3]
rdrp.col<-brewer.pal(3,"Set1")[2]

plot.genes<-~{
  par(mar=c(4,6,0,2), oma=c(0,0,1,1))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  barplot(t(as.matrix(tab[,-1])), beside=T, space=c(0,0),
          xlab="Number of genes",
          col=c(ago.col,rdrp.col,dicer.col), border=NA, names.arg=tab[,1], horiz=T, las=1); grid(nx=NULL, ny=NA)
  legend("topright", title="Gene name",
         c("Argonaute","RdRP","Dicer"),
         fill=c(ago.col,rdrp.col,dicer.col), border=NA, xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.7, pt.cex=1); box()
}

plot.ago<-~{
  ## plot tree
  plot(tree.ago, type="u", show.tip.label=F, no.margin=T, rotate.tree=90)
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("CAEEL",tree.ago$tip.label), pch=15, col="white") ## to provide a bit of white padding
  tiplabels(tip=grep("DROME",tree.ago$tip.label), pch=15, col="white")
  tiplabels(tip=grep("HUMAN",tree.ago$tip.label), pch=15, col="white")
  tiplabels(tip=grep("CAEEL",tree.ago$tip.label), text="C", bg=NA, frame="none", cex=0.6)
  tiplabels(tip=grep("DROME",tree.ago$tip.label), text="D", bg=NA, frame="none", cex=0.6)
  tiplabels(tip=grep("HUMAN",tree.ago$tip.label), text="H", bg=NA, frame="none", cex=0.6)
  ## plot tip cols for rotifers
  tiplabels(tip=grep("PSC1",tree.ago$tip.label), pch=23, bg="blue")
  tiplabels(tip=grep("BRAPC",tree.ago$tip.label), pch=23, bg="blue")
  tiplabels(tip=grep("As",tree.ago$tip.label), pch=21, bg=As.col)
  tiplabels(tip=grep("Rd",tree.ago$tip.label), pch=21, bg=Rd.col)
  tiplabels(tip=grep("Rw",tree.ago$tip.label), pch=21, bg=Rw.col)
  ## clade labels
  text(5.75,0.2, labels="Piwi", cex=8/10)
  text(5.5,4.5, labels="WAGO", cex=8/10)
  text(5.25,3.5, labels="BDAGO I", cex=8/10, font=2)
  text(0.5,4.5, labels="BDAGO II", cex=8/10, font=2)
  ## legend
  legend("topright", title="Species ID",
         c("As","Rd","Rp","Monog."), pch=c(21,21,21,23),
         pt.bg=c(As.col,Rd.col,Rw.col,"blue"), border=NA, xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.7, pt.cex=1)
}

plot.rdrp<-~{
  plot(tree.rdrp, type="u", show.tip.label=F, no.margin=T, rotate.tree=0)
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("CAEEL",tree.rdrp$tip.label), pch=15, col="white")
  tiplabels(tip=grep("ARATH",tree.rdrp$tip.label), pch=15, col="white")
  tiplabels(tip=grep("NEUCR",tree.rdrp$tip.label), pch=15, col="white")
  tiplabels(tip=grep("CAEEL",tree.rdrp$tip.label), text="C", bg=NA, frame="none", cex=0.6)
  tiplabels(tip=grep("ARATH",tree.rdrp$tip.label), text="A", bg=NA, frame="none", cex=0.6)
  tiplabels(tip=grep("NEUCR",tree.rdrp$tip.label), text="N", bg=NA, frame="none", cex=0.6)
  ## plot tip cols for rotifers
  tiplabels(tip=grep("PSC1",tree.rdrp$tip.label), pch=23, bg="blue")
  tiplabels(tip=grep("BRAPC",tree.rdrp$tip.label), pch=23, bg="blue")
  tiplabels(tip=grep("As",tree.rdrp$tip.label), pch=21, bg=As.col)
  tiplabels(tip=grep("Rd",tree.rdrp$tip.label), pch=21, bg=Rd.col)
  tiplabels(tip=grep("Rw",tree.rdrp$tip.label), pch=21, bg=Rw.col)
  ## clade labels
  text(1.25,1.25, labels=expression("RDR"*alpha), cex=8/10)
  text(4.25,3, labels=expression("RDR"*beta), cex=8/10)
  text(0.75,2.75, labels=expression("RDR"*gamma), cex=8/10)
  text(4.25,0.25, labels="RDR I", cex=8/10, font=2)
  text(1.8,3.4, labels="RDR II", cex=8/10, font=2)
}

plot.dicer<-~{
  plot(tree.dicer, type="u", show.tip.label=F, no.margin=T, rotate.tree=0)
  add.scale.bar(font=1, cex=6/10)
  ## plot tip labels for ref sequences
  tiplabels(tip=grep("CAEEL",tree.dicer$tip.label), pch=15, col="white")
  tiplabels(tip=grep("DROME",tree.dicer$tip.label), pch=15, col="white")
  tiplabels(tip=grep("HUMAN",tree.dicer$tip.label), pch=15, col="white")
  tiplabels(tip=grep("ARATH",tree.dicer$tip.label), pch=15, col="white")
  tiplabels(tip=grep("NEUCR",tree.dicer$tip.label), pch=15, col="white")
  tiplabels(tip=grep("CAEEL",tree.dicer$tip.label), text="C", bg=NA, frame="none", cex=0.6)
  tiplabels(tip=grep("DROME",tree.dicer$tip.label), text="D", bg=NA, frame="none", cex=0.6)
  tiplabels(tip=grep("HUMAN",tree.dicer$tip.label), text="H", bg=NA, frame="none", cex=0.6)
  tiplabels(tip=grep("ARATH",tree.dicer$tip.label), text="A", bg=NA, frame="none", cex=0.6)
  tiplabels(tip=grep("NEUCR",tree.dicer$tip.label), text="N", bg=NA, frame="none", cex=0.6)
  ## plot tip cols for rotifers
  tiplabels(tip=grep("PSC1",tree.dicer$tip.label), pch=23, bg="blue")
  tiplabels(tip=grep("BRAPC",tree.dicer$tip.label), pch=23, bg="blue")
  tiplabels(tip=grep("As",tree.dicer$tip.label), pch=21, bg=As.col)
  tiplabels(tip=grep("Rd",tree.dicer$tip.label), pch=21, bg=Rd.col)
  tiplabels(tip=grep("Rw",tree.dicer$tip.label), pch=21, bg=Rw.col)
  ## clade labels
  text(3,2, labels="Dicer1", cex=8/10)
  text(4,5, labels="Dicer2", cex=8/10)
  text(0.3,3, labels="Plant Dicers", cex=8/10)
}

## plot top
top<-plot_grid(plot.genes,plot.ago,plot.rdrp,plot.dicer, nrow=2,ncol=2, labels = c("A","B","C","D"))