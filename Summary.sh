

FOLDER=$1

wc -l $FOLDER/MACS/*_summits.bed > Peaks.txt
wc -l $FOLDER/Raw_Data/*.fastq > RawReads.txt
wc -l $FOLDER/Raw_Data/*.fq > Trimmedreads.txt
wc -l $FOLDER/BED/*.bed > Unique.txt
paste RawReads.txt Unique.txt Peaks.txt Trimmedreads.txt | awk -v OFS="\t" '{print $2,$1/4,$7/4,$3,($3/($1/4)*100),$5}' > Summary.txt


