#!/bin/bash


# RUN scripts e.g. as: ./mini_analysis.sh /Users/kxs624/Documents/data/genomes/human/ /Users/kxs624/tmp/STROBEALIGN/ 

if [ $# -lt 2 ]; then
    # TODO: print usage
    echo "./mini_analysis.sh <ref_path> <outroot>"
    exit 1
fi

refs=$1
outroot=$2
eval_script_dir="/Users/kxs624/Documents/workspace/alignment_evaluation/evaluation_PE"
echo $outroot

mkdir -p $outroot


 # python get_accuracy.py --truth /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads.sam --predicted_sam ~/tmp/STROBEALIGN/strobealign_chr18_100k_PE_reads.sam 

# /usr/bin/time -l minimap2 -t 1 --eqx -ax sr /Users/kxs624/Documents/data/genomes/human/hg38_chr18.fa  /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads_L.fq /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads_R.fq > ~/tmp/STROBEALIGN/minimap2_chr18_100k_PE_reads_f0.00001.sam

# /usr/bin/time -l ./strobealign -t 1  -o ~/tmp/STROBEALIGN/strobealign_chr18_100k_PE_reads.sam /Users/kxs624/Documents/data/genomes/human/hg38_chr18.fa /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads_L.fq /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads_R.fq




echo -n  "tool","ref","%-aligned","accuracy"$'\n'



for read_lengh in 100 150 200 
do 
    for chr_id in hg38_chr21 #hg38_chrX hg38_chr18 hg38_chr15 hg38_chr1
    do
        # mkdir -p $outroot/$chr_id/
        # mason_simulator -ir $refs/$chr_id.fa -n 100000 --illumina-read-length $read_lengh --fragment-mean-size 300 -o $outroot/$chr_id/$read_lengh.L.fq -or $outroot/$chr_id/$read_lengh.R.fq -oa $outroot/$chr_id/$read_lengh.sam

        /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_lengh.L.fq $outroot/$chr_id/$read_lengh.R.fq 1> $outroot/$chr_id/$read_lengh.minimap2.sam 2>  $outroot/$chr_id/$read_lengh.minimap2.stderr
        echo -n $chr_id,minimap2,$read_lengh,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_lengh.sam --predicted_sam $outroot/$chr_id/$read_lengh.minimap2.sam

        /usr/bin/time -l strobealign -t 1  -o $outroot/$chr_id/$read_lengh.srobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_lengh.L.fq $outroot/$chr_id/$read_lengh.R.fq &>  $outroot/$chr_id/$read_lengh.strobealign.stderr
        echo -n $chr_id,strobealign,$read_lengh,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_lengh.sam --predicted_sam $outroot/$chr_id/$read_lengh.strobealign.sam
    done
done



