#!/bin/bash

Group=GROUPGOESHERE
subgroup=SUBgroupGOESHERE
Histone=HISTONEGOESHERE

cd /projects/lu_lab/Tom/LSD_Op/ChIP_seq/$Group/$Histone/$subgroup

cd Raw_Data/
gunzip *.gz
trim_galore *.fastq

cd ..

FILES=$PWD/Raw_Data/*.fq
FQ=.fq
SAM=.sam
LOG=.log
BAM=.bam
SORT=_sort.bam
UNI=_unique.bam
PREBED=_pre.bed
BED=.bed
EXT=_extend.bed
EXPWIN=_extend_promoter_win.bed
EXG4000=_extend_geno_4000.bed
EXPCOR=_extend_promoter_nor.bed
EXGCOR=_extend_geno_4000_nor.bed
G100=_geno_win_100.bed
SORTG100=_geno_win_100_sort.bed
BedG=.bedGraph
NEBw=.bw

mkdir Aligned_BAM
mkdir Aligned_SAM
mkdir SORT_BAM
mkdir UNIQUE_BAM
mkdir BED
mkdir EXTEND_BED
mkdir EXT_PROMOTER_WIN
mkdir EXT_GENO_WIN
mkdir Correlation
mkdir MACS
mkdir GenoWin100Sort
mkdir BedGraph
mkdir Nor_Ext_BW

input_length=$(wc -l < /ChIP_seq/Input/$Group/$Histone/$subgroup/input.bed )

for fn in $FILES
do
echo `basename "$fn"`
f=`basename "${fn%.*}"`
echo $f

bowtie2 -p 16 -x /ChIP_seq/Data/mm10/Sequence/Bowtie2Index/mm10 -U $PWD/Raw_Data/$f$FQ -S $PWD/Aligned_SAM/$f$SAM 2>$PWD/Aligned_SAM/$f$LOG

samtools view -bq 10 $PWD/Aligned_SAM/$f$SAM > $PWD/Aligned_BAM/$f$BAM

samtools sort $PWD/Aligned_BAM/$f$BAM -o $PWD/SORT_BAM/$f$SORT

samtools index $PWD/SORT_BAM/$f$SORT

samtools rmdup -s $PWD/SORT_BAM/$f$SORT $PWD/UNIQUE_BAM/$f$UNI

bedtools bamtobed -i $PWD/UNIQUE_BAM/$f$UNI > $PWD/BED/$f$PREBED

bedtools subtract -a $PWD/BED/$f$PREBED -b /ChIP_seq/Data/mm10/mm10.blacklist.bed > $PWD/BED/$f$BED

bedToBam -i $PWD/BED/$f$BED -g /ChIP_seq/Data/bed_winbed/mm10/mm10genome.bed > $PWD/BED/$f$BAM

bedtools slop -i $PWD/BED/$f$BED -g /ChIP_seq/Data/bed_winbed/mm10/mm10genome.bed -b 100 > $PWD/EXTEND_BED/$f$EXT

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /ChIP_seq/Data/Promoter/mm10/mm10Promoter.bed > $PWD/EXT_PROMOTER_WIN/$f$EXPWIN

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /ChIP_seq/Data/bed_winbed/mm10/mm10genome_win_4000.bed > $PWD/EXT_GENO_WIN/$f$EXG4000

bedtools coverage -counts -b $PWD/EXTEND_BED/$f$EXT -a /ChIP_seq/Data/bed_winbed/mm10/mm10genome_win_100.bed > $PWD/EXT_GENO_WIN/$f$G100

sort -k1,1 -k2,2g -u -o $PWD/GenoWin100Sort/$f$SORTG100 $PWD/EXT_GENO_WIN/$f$G100

ChIP_length=$(wc -l < $PWD/BED/$f$BED)

paste $PWD/EXT_PROMOTER_WIN/$f$EXPWIN /ChIP_seq/Input/$Group/$Histone/$subgroup/input_extend_promoter_win.bed | awk -v OFS="\t" '{print $4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Correlation/$f$EXPCOR

paste $PWD/EXT_GENO_WIN/$f$EXG4000 /ChIP_seq/Input/$Group/$Histone/$subgroup/input_extend_genome_4000win.bed | awk -v OFS="\t" '{print $4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Correlation/$f$EXGCOR

macs2 callpeak -t $PWD/BED/$f$BED -c /ChIP_seq/Input/$Group/$Histone/$subgroup/input.bed -f BED -g mm -n $f -q 0.05 --outdir $PWD/MACS

paste $PWD/EXT_GENO_WIN/$f$G100 /LSD_Op/ChIP_seq/Input/$Group/$Histone/$subgroup/input_extend_genome_100win.bed | awk -v OFS="\t" '{print $1,$2,$3,$4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/BedGraph/$f$BedG

bedSort $PWD/BedGraph/$f$BedG $PWD/BedGraph/$f$BedG

bedGraphToBigWig $PWD/BedGraph/$f$BedG /ChIP_seq/Data/bed_winbed/mm10/mm10genome.bed $PWD/Nor_Ext_BW/$f$NEBw

done

/ChIP_seq/Summary.sh