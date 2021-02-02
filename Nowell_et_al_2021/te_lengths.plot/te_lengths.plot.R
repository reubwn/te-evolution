setwd("~/Google_Drive/Bdelloid/results/pop_genomics/te_evolution/RepeatMasker/RM_full/outfiles/one_code/ectopic.plot/")

## libraries
library(cowplot)

## define cols
As.col<-"#66C2A3"
Rs.col<-"#FA8C61"
Rg.col<-"#8C9ECA"
Rd.col<-"#E68AC2"
Rw.col<-"#A6D751"
Rp.col<-"#FFD92E"

## get data
df1<-read.table("te_lengths.txt", head=T, colClasses=c("factor","character","integer","integer","integer","character","factor"), comment.char="")
str(df1)

## plot length distribution
p1<-~{
  par(mar=c(3,3,2,1), oma=c(1,1,0,0))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  ## plot the grid
  boxplot(log10(length/1000) ~ ID*type, data=df1, 
          at=(1:47)[-c(12,24,36)],
          xaxt="n", xlab="", ylab="Log10 TE length (kb)", ylim=c(-2,3),
          col=NA, outline=F, border=NA); grid(nx=NA, ny=NULL)## plot the boxes
  boxplot(log10(length/1000) ~ ID*type, data=df1, 
          at=(1:47)[-c(12,24,36)],
          xaxt="n", yaxt="n",
          col=c("blue","black",As.col,"#F9A65A","#CD7058",Rd.col,Rp.col,Rw.col,"#599AD3",Rg.col,Rs.col),
          boxcol=c("blue","black",As.col,"#F9A65A","#CD7058",Rd.col,Rp.col,Rw.col,"#599AD3",Rg.col,Rs.col),
          outline=T, outcex=0.5, outpch=16, outcol="grey",
          medcol="white", whisklty=1, staplecol=NA, add=T)
  axis(1, at=c(0,12,24,36,48), lab=F)
  ## xlab
  mtext(c("DNA","LINE","LTR","PLE"), at=c(6.5,18.5,30.5,42.5), side=1, line=0.5, las=1)
  ## legend
  legend("topright", title="Sample ID", ncol=3,
         c("Bc","Dc","As","Av","Ar","Rd","Rp","Rw","Rc","Rg","Rs"),
         pch=16, col=c("blue","black",As.col,"#F9A65A","#CD7058",Rd.col,Rp.col,Rw.col,"#599AD3",Rg.col,Rs.col),
         xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.7, pt.cex=1);box()
}

## get data
df2 <- read.table("te_counts.txt", head=T, colClasses=c("factor","character","integer","numeric"), comment.char="")
str(df2)

## PLOT
p2<-~{
  par(mar=c(3,3,2,1), oma=c(1,1,0,0))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  
  plot(subset(df2$mean/1000,df3$ID=="Bc_PSC1"), subset(df2$count,df3$ID=="Bc_PSC1"), type="n", xlim=c(0,3), ylim=c(0,1000), xlab="Mean TE length (kb)", ylab="TE count");grid()
  rect(0,-500,0.6,1500, col="grey95", border=NA); abline(v=0.6, lty=2) ## plot rectangle at 'threshold'
  points(subset(df2$mean/1000,df2$ID=="Bc_PSC1"), subset(df2$count,df2$ID=="Bc_PSC1"), pch=23, bg="blue")
  points(subset(df2$mean/1000,df2$ID=="Dc_DCAR706"), subset(df2$count,df2$ID=="Dc_DCAR706"), pch=16, col="black")
  points(subset(df2$mean/1000,df2$ID=="As_ASTE804"), subset(df2$count,df2$ID=="As_ASTE804"), pch=16, col=As.col)
  points(subset(df2$mean/1000,df2$ID=="Av_Av2013"), subset(df2$count,df2$ID=="Av_Av2013"), pch=3, col="#F9A65A")
  points(subset(df2$mean/1000,df2$ID=="Ar_Ar2018"), subset(df2$count,df2$ID=="Ar_Ar2018"), pch=3, col="#CD7058")
  points(subset(df2$mean/1000,df2$ID=="Rd_RSOR408"), subset(df2$count,df2$ID=="Rd_RSOR408"), pch=16, col=Rd.col)
  points(subset(df2$mean/1000,df2$ID=="Rp_RPSE411"), subset(df2$count,df2$ID=="Rp_RPSE411"), pch=16, col=Rp.col)
  points(subset(df2$mean/1000,df2$ID=="Rw_RSIL801"), subset(df2$count,df2$ID=="Rw_RSIL801"), pch=16, col=Rw.col)
  points(subset(df2$mean/1000,df2$ID=="Rc_Rc2018"), subset(df2$count,df2$ID=="Rc_Rc2018"), pch=3, col="#599AD3")
  points(subset(df2$mean/1000,df2$ID=="Rg_MAG1"), subset(df2$count,df2$ID=="Rg_MAG1"), pch=16, col=Rg.col)
  points(subset(df2$mean/1000,df2$ID=="Rs_AK11"), subset(df2$count,df2$ID=="Rs_AK11"), pch=16, col=Rs.col)
  legend("topright", title="Sample ID", ncol=2,
         c("Bc","Av","Ar","Rc","Dc","As","Rd","Rp","Rw","Rg","Rs"),
         pch=c(23,3,3,3,16,16,16,16,16,16,16),
         col=c("black","#F9A65A","#CD7058","#599AD3","black",As.col,Rd.col,Rp.col,Rw.col,Rg.col,Rs.col),
         pt.bg=c("blue",NA,NA,NA,NA,NA,NA,NA,NA,NA,NA),
         xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.7, pt.cex=1);box()
}

## plot grid
plot_grid(p1,p2, labels="AUTO")
