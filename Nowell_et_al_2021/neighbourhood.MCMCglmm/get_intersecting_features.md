## Pipeline for neighbourhood analysis
- Code is the same for LINEs, LTRs and PLEs; BUSCOs slightly different (see below).
- Not all species could be included in for loops if filenames didn't match, but the logic is the same

### Putative RTEs (LINEs, LTRs, PLEs)

Working from a dir containing lists of putative RT genes (one per species, e.g. `Ar_ARIC003.LINEs.list` etc):
1. Grep gene coordinates from annotation GFF and convert to bed format:
```bash
for p in {Ar_ARIC003,Av_Av2013,As_ASTE804,As_ASTE805,As_ASTE806,Bc_PSC1,Dc_DCAR505,Dc_DCAR706,Rc_Rc2018,Rd_RSOR408,Rd_RSOR410,Rd_RSOR504,Rg_MAG1,Rg_MAG2,Rg_MAG3,Rg_RM9,Rg_RM15,Rp_RPSE411,Rp_RPSE503,Rp_RPSE809,Rp_RPSE812,Rs_AK11,Rs_AK15,Rs_AK16,Rs_AK27,Rs_RS1,Rw_RSIL801,Rw_RSIL802,Rw_RSIL804,Rw_RSIL806}; do 
 export P=$p; echo $p; 
 grep -Ff ${p}.LINEs.list ./annotation/augustus_${p}.fixed.gff3 \
 | grep "transcript\|mRNA" \
 | perl -lane 'if(/ID=(\w+\|\w+\.\w\w)\;/){print join("\t",$F[0],$F[3],$F[4],"ID=$1;type=LINE")}' \
 > ${p}.LINEs.list.bed; 
done
```

2. Make 50 kb window around each gene:
```bash
for p in {Ar_ARIC003,As_ASTE804,As_ASTE805,As_ASTE806,Bc_PSC1,Dc_DCAR505,Dc_DCAR706,Rd_RSOR408,Rd_RSOR410,Rd_RSOR504,Rg_MAG1,Rg_MAG2,Rg_MAG3,Rg_RM9,Rg_RM15,Rp_RPSE411,Rp_RPSE503,Rp_RPSE809,Rp_RPSE812,Rs_AK11,Rs_AK15,Rs_AK16,Rs_AK27,Rs_RS1,Rw_RSIL801,Rw_RSIL802,Rw_RSIL804,Rw_RSIL806}; do 
 echo $p; 
 bedtools slop \
  -i ${p}.LINEs.list.bed \
  -g ./genome_files/scaffolds_${p}_maxhap_genomic.fna.genome.txt \
  -b 25000 \
 > ${p}.LINEs.flanking.25kb.bed; 
done
```

3. Find intersecting features:
```bash
for p in {Ar_ARIC003,Av_Av2013,As_ASTE804,As_ASTE805,As_ASTE806,Bc_PSC1,Bp_HYR1,Dc_DCAR505,Dc_DCAR706,Rc_Rc2018,Rd_RSOR408,Rd_RSOR410,Rd_RSOR504,Rg_MAG1,Rg_MAG2,Rg_MAG3,Rg_RM9,Rg_RM15,Rp_RPSE411,Rp_RPSE503,Rp_RPSE809,Rp_RPSE812,Rs_AK11,Rs_AK15,Rs_AK16,Rs_AK27,Rs_RS1,Rw_RSIL801,Rw_RSIL802,Rw_RSIL804,Rw_RSIL806}; do 
 export P=$p; echo $p; 
 bedtools intersect \
  -a ${p}.LINEs.flanking.25kb.bed \
  -b <(bedtools intersect -a ./annotation/${p}/augustus_${p}.gff3.genes.bed -b ./one_code/${p}.out.one_code.bed -v) ./one_code/${p}.out.one_code.bed ./bedfiles/${p}.telo.filtered.bed \
  -wo -names gene TE telo \
 | perl -lane 'print join("\t",$F[0],$F[1],$F[2],$F[3],$F[4],$F[-1])' \
 | sort -k1,1n -k4,4n -k5,5n \
 | bedtools groupby -g 1-5 -c 6 -o sum \
 | perl -lane '$h1{(join("\t",$F[0],$F[1],$F[2],$F[3]))}{$F[4]}=$F[5];END{foreach $k1 (keys %h1){%h2=%{$h1{$k1}}; if($h2{gene}){print join("\t",$k1,"gene",$h2{gene})}else{print join("\t",$k1,"gene",0)}; if($h2{TE}){print join("\t",$k1,"TE",$h2{TE})}else{print join("\t",$k1,"TE",0)}; if($h2{telo}){print join("\t",$k1,"telo",$h2{telo})}else{print join("\t",$k1,"telo",0)}}}' \
 | sort -k1,1n -k4,4n -k5,5n \
 | perl -lane 'print join("\t",$ENV{P},substr($ENV{P},0,2),$F[0],$F[1],$F[2],($F[2]-$F[1]),$F[3],+(split /=/,$F[3])[2],$F[4],$F[5])' \
 > ${p}.LINEs.flanking.25kb.bed.span; 
done
```

### BUSCOs
For BUSCO genes code is slightly different.

1. Make 50 kb window around each BUSCO:
```bash
for p in {Ar_ARIC003,As_ASTE804,As_ASTE805,As_ASTE806,Bc_PSC1,Dc_DCAR505,Dc_DCAR706,Rd_RSOR408,Rd_RSOR410,Rd_RSOR504,Rg_MAG1,Rg_MAG2,Rg_MAG3,Rg_RM9,Rg_RM15,Rp_RPSE411,Rp_RPSE503,Rp_RPSE809,Rp_RPSE812,Rs_AK11,Rs_AK15,Rs_AK16,Rs_AK27,Rs_RS1,Rw_RSIL801,Rw_RSIL802,Rw_RSIL804,Rw_RSIL806}; do 
 echo $p; 
 bedtools slop \
  -i full_table_metazoa_odb9_scaffolds_${p}_maxhap_genomic.fna.tsv.bed \
  -g ./genome_files/scaffolds_${p}_maxhap_genomic.fna.genome.txt \
  -b 25000 \
 > ${p}.BUSCO.flanking.25kb.bed; 
done
```

2. Find intersecting features:
```bash
for p in {Ar_ARIC003,Av_Av2013,As_ASTE804,As_ASTE805,As_ASTE806,Bc_PSC1,Bp_HYR1,Dc_DCAR505,Dc_DCAR706,Rc_Rc2018,Rd_RSOR408,Rd_RSOR410,Rd_RSOR504,Rg_MAG1,Rg_MAG2,Rg_MAG3,Rg_RM9,Rg_RM15,Rp_RPSE411,Rp_RPSE503,Rp_RPSE809,Rp_RPSE812,Rs_AK11,Rs_AK15,Rs_AK16,Rs_AK27,Rs_RS1,Rw_RSIL801,Rw_RSIL802,Rw_RSIL804,Rw_RSIL806}; do 
 export P=$p; echo $p; 
 bedtools intersect \
  -a ${p}.BUSCO.flanking.25kb.bed \
  -b <(bedtools intersect -a ./annotation/${p}/augustus_${p}.gff3.genes.bed -b ./one_code/${p}.out.one_code.bed -v) ./one_code/${p}.out.one_code.bed ./telomere_out/bedfiles/${p}.telo.filtered.bed \
  -wo -names gene TE telo \
 | perl -lane 'print join("\t",$F[0],$F[1],$F[2],$F[3],$F[5],$F[-1])' \
 | sort -k1,1n -k4,4n -k5,5n \
 | bedtools groupby -g 1-5 -c 6 -o sum \
 | perl -lane '$h1{(join("\t",$F[0],$F[1],$F[2],$F[3]))}{$F[4]}=$F[5];END{foreach $k1 (keys %h1){%h2=%{$h1{$k1}}; if($h2{gene}){print join("\t",$k1,"gene",$h2{gene})}else{print join("\t",$k1,"gene",0)}; if($h2{TE}){print join("\t",$k1,"TE",$h2{TE})}else{print join("\t",$k1,"TE",0)}; if($h2{telo}){print join("\t",$k1,"telo",$h2{telo})}else{print join("\t",$k1,"telo",0)}}}' \
 | sort -k1,1n -k4,4n -k5,5n \
 | perl -lane 'print join("\t",$ENV{P},substr($ENV{P},0,2),$F[0],$F[1],$F[2],($F[2]-$F[1]),$F[3],"BUSCO",$F[4],$F[5])' \
 > ${p}.BUSCO.flanking.25kb.bed.span; 
done
```

### Concatenate
Cat the whole lot together:
```bash
cat BUSCOs/*flanking.25kb.bed.span LINEs/*flanking.25kb.bed.span LTRs/*flanking.25kb.bed.span PLEs/*flanking.25kb.bed.span > neighbourhood.25kb.txt
```
