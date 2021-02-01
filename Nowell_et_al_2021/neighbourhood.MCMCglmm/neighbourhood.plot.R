setwd("~/Google_Drive/Bdelloid/results/pop_genomics/te_evolution/type_I_RT_survey/RT_2020_Av_TEs_Athena_P_full/hmm1e-5/neighbourhood_analysis/neighbourhood.MCMCglmm/")

library(cowplot)

##READ IN DATA
df<-read.table("neighbourhood.25kb.txt", col.names=c("ID","species","chrom","start","end","window","name","type","feature","span"), stringsAsFactors=T)
## explanation: 'span' represents the number of bp assigned to 'feature' in a 50 kb window drawn around BUSCO/TE of type 'type'
## add asexual and desiccation factors
df$is.bdelloid <- as.factor(ifelse(df$ID %in% c('Bc_PSC1','Bp_HYR1'), 0, 1))
df$is.desiccating <- as.factor(ifelse(df$ID %in% c('Bc_PSC1','Bp_HYR1','Rg_MAG1','Rg_MAG2','Rg_MAG3','Rg_RM9','Rg_RM15','Rs_AK11','Rs_AK15','Rs_AK16','Rs_AK27','Rs_RS1','Rc_Rc2018'), 0, 1))
## add BUSCO as factor
df$is.BUSCO <- as.factor(ifelse(df$type=="BUSCO", 1, 0))
## have a look
str(df)

## pull out subset for plotting (one per species)
df.subset<-df[df$ID %in% c('Bc_PSC1','Dc_DCAR706','Av_Av2013','Ar_ARIC003','As_ASTE805','Rg_MAG3','Rs_AK11','Rc_Rc2018','Rd_RSOR408','Rw_RSIL806','Rp_RPSE503'),]
df.subset<-droplevels(df.subset)
##reorder
df.subset$ID<-factor(df.subset$ID, levels=c('Bc_PSC1','Dc_DCAR706','Av_Av2013','Ar_ARIC003','As_ASTE805','Rg_MAG3','Rs_AK11','Rc_Rc2018','Rd_RSOR408','Rw_RSIL806','Rp_RPSE503'))
str(df.subset)

##PLOT RAW DATA

## define cols
As.col<-"#66C2A3"
Av.col<-"#F9A65A"
Ar.col<-"#CD7058"
Rs.col<-"#FA8C61"
Rg.col<-"#8C9ECA"
Rd.col<-"#E68AC2"
Rw.col<-"#A6D751"
Rp.col<-"#FFD92E"
Rc.col<-"#599AD3"

## plot
p1<-~{
  par(mar=c(2,3,2,1), oma=c(1,1,0,0))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  boxplot(span/1000 ~ ID*type, data=subset(df.subset,feature=="gene"),
          at=(1:47)[-c(12,24,36)], ylim=c(0,50),
          xaxt="n", xlab="", ylab="Other genes (kb)", cex.axis=0.8,
          col=NA, outline=T, border=NA); grid(nx=NA, ny=NULL)## plot the boxes
  boxplot(span/1000 ~ ID*type, data=subset(df.subset,feature=="gene"),
          at=(1:47)[-c(12,24,36)],
          xaxt="n", yaxt="n",
          col=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
          boxcol=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
          outline=T, outcex=0.5, outpch=16, outcol="grey",
          medcol="white", whisklty=1, staplecol=NA,
          add=T)
  axis(1, at=c(0,12,24,36,48), lab=F)
  ## xlab
  mtext(c("BUSCO","LINE","LTR","PLE"), at=c(6.5,18.5,30.5,42.5), side=1, line=0.5, las=1, cex=0.9)
  ## legend
  legend("topright", title="Sample ID", ncol=3,
         c("Bc","Dc","Av","Ar","As","Rg","Rs","Rc","Rd","Rw","Rp"),
         pch=16, col=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
         xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.7, pt.cex=1);box()
  
}

p2<-~{
  par(mar=c(2,3,2,1), oma=c(1,1,0,0))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  boxplot(span/1000 ~ ID*type, data=subset(df.subset,feature=="TE"),
          at=(1:47)[-c(12,24,36)], ylim=c(0,50),
          xaxt="n", xlab="", ylab="Other TEs (kb)", cex.axis=0.8,
          col=NA, outline=T, border=NA); grid(nx=NA, ny=NULL)## plot the boxes
  boxplot(span/1000 ~ ID*type, data=subset(df.subset,feature=="TE"),
          at=(1:47)[-c(12,24,36)],
          xaxt="n", yaxt="n",
          col=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
          boxcol=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
          outline=T, outcex=0.5, outpch=16, outcol="grey",
          medcol="white", whisklty=1, staplecol=NA,
          add=T)
  axis(1, at=c(0,12,24,36,48), lab=F)
  ## xlab
  mtext(c("BUSCO","LINE","LTR","PLE"), at=c(6.5,18.5,30.5,42.5), side=1, line=0.5, las=1, cex=0.9)
  ## legend
  legend("topright", title="Sample ID", ncol=3,
         c("Bc","Dc","Av","Ar","As","Rg","Rs","Rc","Rd","Rw","Rp"),
         pch=16, col=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
         xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.7, pt.cex=1);box()
  
}
p3<-~{
  par(mar=c(2,3,2,1), oma=c(1,1,0,0))
  par(tcl=-0.25)
  par(mgp=c(2, 0.6, 0))
  boxplot(span/1000 ~ ID*type, data=subset(df.subset,feature=="telo"),
          at=(1:47)[-c(12,24,36)], ylim=c(0,0.5),
          xaxt="n", xlab="", ylab="Telomeric repeat (kb)", cex.axis=0.8,
          col=NA, outline=T, border=NA); grid(nx=NA, ny=NULL)## plot the boxes
  boxplot(span/1000 ~ ID*type, data=subset(df.subset,feature=="telo"),
          at=(1:47)[-c(12,24,36)],
          xaxt="n", yaxt="n",
          col=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
          boxcol=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
          outline=T, outcex=0.5, outpch=16, outcol="grey",
          medcol="white", whisklty=1, staplecol=NA,
          add=T)
  axis(1, at=c(0,12,24,36,48), lab=F)
  ## xlab
  mtext(c("BUSCO","LINE","LTR","PLE"), at=c(6.5,18.5,30.5,42.5), side=1, line=0.5, las=1, cex=0.9)
  ## legend
  legend("topright", title="Sample ID", ncol=3,
         c("Bc","Dc","Av","Ar","As","Rg","Rs","Rc","Rd","Rw","Rp"),
         pch=16, col=c("blue","black",Av.col,Ar.col,As.col,Rg.col,Rs.col,Rc.col,Rd.col,Rw.col,Rp.col),
         xjust=1, y.intersp=0.9, bg="grey95", box.col=NA, cex=0.7, pt.cex=1);box()
  
}
## normal scale
plot_grid(p1,p2,p3, ncol=3, labels=c("B","C","D"))
