#!/bin/bash


# RUN scripts e.g. as: ./mini_analysis.sh /Users/kxs624/Documents/data/genomes/human/ /Users/kxs624/tmp/STROBEALIGN/ 

if [ $# -lt 2 ]; then
    # TODO: print usage
    echo "./mini_analysis.sh <ref_path> <outroot>"
    exit 1
fi

refs=$1
outroot=$2
eval_script_dir="/Users/kxs624/Documents/workspace/alignment_evaluation/scripts"
strobealign_dev_dir="/Users/kxs624/Documents/workspace/StrobeAlign"
echo $outroot

mkdir -p $outroot


 # python get_accuracy.py --truth /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads.sam --predicted_sam ~/tmp/STROBEALIGN/strobealign_chr18_100k_PE_reads.sam 

# /usr/bin/time -l minimap2 -t 1 --eqx -ax sr /Users/kxs624/Documents/data/genomes/human/hg38_chr18.fa  /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads_L.fq /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads_R.fq > ~/tmp/STROBEALIGN/minimap2_chr18_100k_PE_reads_f0.00001.sam

# /usr/bin/time -l ./strobealign -t 1  -o ~/tmp/STROBEALIGN/strobealign_chr18_100k_PE_reads.sam /Users/kxs624/Documents/data/genomes/human/hg38_chr18.fa /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads_L.fq /Users/kxs624/Documents/workspace/StrobeAlign/data/hg38_chr18_100k_PE_reads_R.fq


# USE BITCOUNT n_BITS as parameter to software! bitcount 8 for < 100nt, 16 for 150-200, 32 for 200-250, 64 for >=300

# Test read length 50-70: 
# DECIDED k 20,l -3, u 0,c 8

# read length 100:
# DECIDED k 20, l -3, u 3, c 8

# read length 150: 
# DECIDED  k=20, l 0, u 7, c 8 

# read length 200-250:
# DECIDED k=22, l 2, u 10, c 8

# read length 300:
# DECIDED  k=23, l 2, u 10, c 8

# It is possible we are overfitting on smaller chromosomes, and that choosing 
#slightly larger k would be better on whole hg38! In that case, test on larger 
# joint ref with several chromosomes! Set k one or two larger

# Create one large dataset of 1M reads and chrX, chrY, chr18,chr15, chr6
# python ~/Documents/workspace/uLTRA/scripts/filter_fasta.py --keep_only chr6 chr15 chr18 chrX chrY --outfile ../hg38_chr6_15_18_X_Y.fa GRCh38.p13.genome.fa
#

# for k in 17 18 19 20 
# for l in -2 -1 0 1
# for u in 3 5 7 9

echo -n  "tool","ref","%-aligned","accuracy,time(sec),Mem(MB)"$'\n'

# for read_lengh in 100 #100 #150 #200 250 300 
# do 
#     for chr_id in hg38_chr21 hg38_chr1 hg38_chr15 hg38_chr18 hg38_chrX  #hg38_chr1_2 hg38_chr6_15_18_X_Y # hg38_chr1 hg38_chr15 hg38_chr18 hg38_chr21 hg38_chrX #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
#     do
#         # mkdir -p $outroot/$chr_id/
#         # mason_simulator -ir $refs/$chr_id.fa -n 100000 --illumina-read-length $read_lengh --fragment-mean-size 300 -o $outroot/$chr_id/$read_lengh.L.fq -or $outroot/$chr_id/$read_lengh.R.fq -oa $outroot/$chr_id/$read_lengh.sam

#         # minimap2 stats
#         /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_lengh.L.fq $outroot/$chr_id/$read_lengh.R.fq 1> $outroot/$chr_id/$read_lengh.minimap2.sam 2>  $outroot/$chr_id/$read_lengh.minimap2.stderr
#         echo -n $chr_id,$read_lengh,minimap2,align,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_lengh.sam --predicted_sam $outroot/$chr_id/$read_lengh.minimap2.sam --time_mem $outroot/$chr_id/$read_lengh.minimap2.stderr

#         /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -k 22 -l 0 -u 10 -o $outroot/$chr_id/$read_lengh.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_lengh.L.fq $outroot/$chr_id/$read_lengh.R.fq &>  $outroot/$chr_id/$read_lengh.strobealign.stderr
#         echo -n $chr_id,$read_lengh,strobealign,align,8,22,0,10,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_lengh.sam --predicted_sam $outroot/$chr_id/$read_lengh.strobealign.sam --time_mem $outroot/$chr_id/$read_lengh.strobealign.stderr

#         for bc in 8 #16 # 32
#         do
#             for k in 19 20 #18 19 20 #18 19 20 21
#             do
#                 for l in -3 -2
#                 do
#                     for u in 3 #3 5 #9 #5 7 9
#                     do

#                         /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -c $bc -k $k -l $l -u $u -o $outroot/$chr_id/$read_lengh.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_lengh.L.fq $outroot/$chr_id/$read_lengh.R.fq &>  $outroot/$chr_id/$read_lengh.strobealign.stderr
#                         echo -n $chr_id,$read_lengh,strobealign,align,$bc,$k,$l,$u,
#                         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_lengh.sam --predicted_sam $outroot/$chr_id/$read_lengh.strobealign.sam --time_mem $outroot/$chr_id/$read_lengh.strobealign.stderr

#                         # Mapping stats
#                         # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -k $k -x -o $outroot/$chr_id/$read_lengh.strobealign.paf $refs/$chr_id.fa $outroot/$chr_id/$read_lengh.L.fq $outroot/$chr_id/$read_lengh.R.fq &>  $outroot/$chr_id/$read_lengh.strobealign.stderr
#                         # echo -n $chr_id,$read_lengh,strobealign_k20,map,
#                         # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_lengh.sam --predicted_paf $outroot/$chr_id/$read_lengh.strobealign.paf --time_mem $outroot/$chr_id/$read_lengh.strobealign.stderr
#                     done
#                 done
#             done
#             echo
#         done
#         echo "NEW CHR"
#         echo
#     done
# done


for read_lengh in 100 150 200 250 300 
do 
    for chr_id in hg38_chr1_2 hg38_chr6_15_18_X_Y # hg38_chr1 hg38_chr15 hg38_chr18 hg38_chr21 hg38_chrX #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
    do

        # minimap2 stats
        /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_lengh.L.fq $outroot/$chr_id/$read_lengh.R.fq 1> $outroot/$chr_id/$read_lengh.minimap2.sam 2>  $outroot/$chr_id/$read_lengh.minimap2.stderr
        echo -n $chr_id,$read_lengh,minimap2,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_lengh.sam --predicted_sam $outroot/$chr_id/$read_lengh.minimap2.sam --time_mem $outroot/$chr_id/$read_lengh.minimap2.stderr

        /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -r $read_lengh -o $outroot/$chr_id/$read_lengh.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_lengh.L.fq $outroot/$chr_id/$read_lengh.R.fq &>  $outroot/$chr_id/$read_lengh.strobealign.stderr
        echo -n $chr_id,$read_lengh,strobealign,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_lengh.sam --predicted_sam $outroot/$chr_id/$read_lengh.strobealign.sam --time_mem $outroot/$chr_id/$read_lengh.strobealign.stderr

    done
done


