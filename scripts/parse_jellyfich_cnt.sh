#!/bin/sh

scriptfolder="/Users/kxs624/Documents/workspace/alignment_evaluation/scripts"
datafolder="~/Dropbox/Work/projects/strobealign/journal2022/REVISION_GenomeBiology/data/v0.7/seed_analysis"

for k in 20 30 40 50 60 70 80 90 100
do
	python scriptfolder/jellyfish_kmer_cnt_parse.py $datafolder/jellyfish_kmers_hg38/histo.$k.tsv $k hg38 
done

for k in 20 30 40 50 60 70 80 90 100
do
	python scriptfolder/jellyfish_kmer_cnt_parse.py $datafolder/jellyfish_kmers_chm13/histo.$k.tsv $k CHM13 
done