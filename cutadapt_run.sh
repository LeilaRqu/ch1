#!/bin/bash

for FWDpath in /home/ler4794/ch1/ch1raw/demux/*_R1_001.fastq; do
	FWD=$(basename $FWDpath) #pull file name for fwd read
	REV=$(echo $FWD | sed 's/_R1_001.fastq/_R2_001.fastq/') #pull file name for rev read
	cutadapt -g AACTTTYRRCAAYGGATCWCT -G AGCCTCCGCTTATTGATATGCTTAART --discard-untrimmed --overlap 5 -e 0.15 --max-n=0 -o /home/ler4794/ch1/trimmed/trimmed_$FWD -p /home/ler4794/ch1/trimmed/trimmed_$REV /home/ler4794/ch1/ch1raw/demux/$FWD /home/ler4794/ch1/ch1raw/demux/$REV
done
