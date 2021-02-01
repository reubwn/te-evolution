setwd("./")

## libraries
library(cowplot)
library(RColorBrewer)
## cols
As.col<-"#66C2A3"
Rs.col<-"#FA8C61"
Rg.col<-"#8C9ECA"
Rd.col<-"#E68AC2"
Rw.col<-"#A6D751"
Rp.col<-"#FFD92E"
## cols
ago.col<-brewer.pal(3,"Set1")[1]
dicer.col<-brewer.pal(3,"Set1")[3]
rdrp.col<-brewer.pal(3,"Set1")[2]

## Bp_HYR1 data 
BRACP<-read.table("pfam_Bp_HYR1.table", head=T)
BRACP$log.UniProt<-log(BRACP$UniProt)
BRACP$log.BRACP<-log(BRACP$BRACP)

plot.BRACP<-~{
  par(mar=c(3,3,2,1), oma=c(0,0,0,0))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  plot(BRACP$log.BRACP~BRACP$log.UniProt, col="grey", xlim=c(-2,8),ylim=c(-2,8),
       xlab="Ln Count eukaryote", ylab="Ln Count monogonont")
  abline(a=0,b=1,col="black",lty=2)
  points(BRACP$log.BRACP[grep("PF02170", BRACP$PFAM)]~BRACP$log.UniProt[grep("PF02170", BRACP$PFAM)], pch=2, col=ago.col)
  points(BRACP$log.BRACP[grep("PF02171", BRACP$PFAM)]~BRACP$log.UniProt[grep("PF02171", BRACP$PFAM)], pch=6, col=ago.col)
  points(BRACP$log.BRACP[grep("PF03368", BRACP$PFAM)]~BRACP$log.UniProt[grep("PF03368", BRACP$PFAM)], pch=1, col=dicer.col)
  points(BRACP$log.BRACP[grep("PF05183", BRACP$PFAM)]~BRACP$log.UniProt[grep("PF05183", BRACP$PFAM)], pch=0, col=rdrp.col)
}

## 10x pseudohap data
pseudohap<-read.table("pfam_10x_p.table", head=T)
pseudohap$log.UniProt<-log(pseudohap$UniProt)
pseudohap$log.pseudohap<-log(pseudohap$pseudohap)

plot.p<-~{
  par(mar=c(3,3,2,1), oma=c(0,0,0,0))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  plot(pseudohap$log.pseudohap~pseudohap$log.UniProt, col="grey", xlim=c(-2,8),ylim=c(-2,8),
       xlab="Ln Count eukaryote", ylab="Ln Count bdelloid")
  abline(a=0,b=1,col="black",lty=2)
  points(pseudohap$log.pseudohap[grep("PF02170", pseudohap$PFAM)]~pseudohap$log.UniProt[grep("PF02170", pseudohap$PFAM)], pch=24, bg=ago.col)
  points(pseudohap$log.pseudohap[grep("PF02171", pseudohap$PFAM)]~pseudohap$log.UniProt[grep("PF02171", pseudohap$PFAM)], pch=25, bg=ago.col)
  points(pseudohap$log.pseudohap[grep("PF03368", pseudohap$PFAM)]~pseudohap$log.UniProt[grep("PF03368", pseudohap$PFAM)], pch=21, bg=dicer.col)
  points(pseudohap$log.pseudohap[grep("PF05183", pseudohap$PFAM)]~pseudohap$log.UniProt[grep("PF05183", pseudohap$PFAM)], pch=22, bg=rdrp.col)
}

## take ratio
pseudohap$ratio <- pseudohap$pseudohap/pseudohap$UniProt
pseudohap$ratio.corr <- (pseudohap$pseudohap/2)/pseudohap$UniProt
pseudohap$log.ratio <- log(pseudohap$pseudohap/pseudohap$UniProt)
pseudohap$log.ratio.corr <- log((pseudohap$pseudohap/2)/pseudohap$UniProt)
head(pseudohap)

plot.histo<-~{
  par(mar=c(3,3,2,1), oma=c(0,0,0,0))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  ## make histogram
  hist(pseudohap$log.ratio.corr, breaks=50, 
       xlim=c(-4,4), ylim=c(0,800),
       col=NA, border=NA, main="", xlab="Ln Abundance score"); grid(nx=NA, ny=NULL)
  ## get quantiles
  q<-quantile(pseudohap$log.ratio.corr, probs=c(0.05,0.95)); q
  ## plot 95% CI rect and bars
  rect(q[1],0,q[2],800, col="grey95", border=NA)
  hist(pseudohap$log.ratio.corr, breaks=50, add=T, col="black", border="white", main="", xlab="", ylab="")
  ## plot lines
  abline(v=pseudohap$log.ratio.corr[grep("PF03368", pseudohap$PFAM)], col=dicer.col)
  abline(v=pseudohap$log.ratio.corr[grep("PF02170", pseudohap$PFAM)], col=ago.col)
  abline(v=pseudohap$log.ratio.corr[grep("PF02171", pseudohap$PFAM)], col=ago.col)
  abline(v=pseudohap$log.ratio.corr[grep("PF05183", pseudohap$PFAM)], col=rdrp.col)
  points(pseudohap$log.ratio.corr[grep("PF03368",pseudohap$PFAM)],800, pch=21, cex=1, xpd=T, bg=dicer.col)
  points(pseudohap$log.ratio.corr[grep("PF02170",pseudohap$PFAM)],800, pch=24, cex=1, xpd=T, bg=ago.col)
  points(pseudohap$log.ratio.corr[grep("PF02171",pseudohap$PFAM)],800, pch=25, cex=1, xpd=T, bg=ago.col)
  points(pseudohap$log.ratio.corr[grep("PF05183",pseudohap$PFAM)],800, pch=22, cex=1, xpd=T, bg=rdrp.col)
  legend("topleft", title="Domain",
         c("PAZ","PIWI","RdRP","Dicer"),
         pch=c(24,25,22,21), pt.bg=c(ago.col,ago.col,rdrp.col,dicer.col),
         xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.7, pt.cex=1);box()
}

## check position of quantiles
q<-quantile(pseudohap$log.ratio.corr, probs=c(0.05,0.95)); q
pseudohap$log.ratio.corr[grep("PF02170", pseudohap$PFAM)]
pseudohap$log.ratio.corr[grep("PF02171", pseudohap$PFAM)]
pseudohap$log.ratio.corr[grep("PF03368", pseudohap$PFAM)]
pseudohap$log.ratio.corr[grep("PF05183", pseudohap$PFAM)]

## plot
bottom<-plot_grid(plot.BRACP,plot.p,plot.histo, nrow=1, rel_widths=c(.66,.66,1), labels=c("B","C","D"))
