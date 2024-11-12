#!/bin/bash

cd /ChIP_seq/Input/V_V/FEMALE/NAC/

gunzip *.gz
trim_galore *.fastq

bowtie2 -x /ChIP_seq/Data/mm10/Sequence/Bowtie2Index/mm10 -U /ChIP_seq/Input/V_V/FEMALE/NAC/VEH-VEH-Input-NAC-Female_trimmed.fq -S input.sam 2> alignSumm.log

samtools view -bq 1 ./input.sam > ./input.bam

samtools sort input.bam -o input_sort.bam

samtools index input_sort.bam

samtools rmdup -s input_sort.bam input_unique.bam

bedtools bamtobed -i input_unique.bam > input_pre.bed

bedtools subtract -a ./input_pre.bed -b /ChIP_seq/Data/mm10/mm10.blacklist.bed > input.bed

bedtools slop -i input.bed -g /ChIP_seq/Data/bed_winbed/mm10/mm10genome.bed -b 100 > input_extend.bed

bedtools coverage -counts -b input_extend.bed -a /ChIP_seq/Data/Promoter/mm10/mm10Promoter.bed > input_extend_promoter_win.bed

bedtools coverage -counts -b input_extend.bed -a /ChIP_seq/Data/bed_winbed/mm10/mm10genome_win_100.bed > input_extend_genome_100win.bed

bedtools coverage -counts -b input_extend.bed -a /ChIP_seq/Data/bed_winbed/mm10/mm10genome_win_4000.bed > input_extend_genome_4000win.bed
