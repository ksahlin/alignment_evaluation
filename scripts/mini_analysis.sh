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
snv_rate="0.005"
small_indel_rate="0.005"
indel_rate="0.00001"

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

# for read_length in 100 150 200 250 300 #100 #150 #200 250 300 
# do 
#     for chr_id in sim_20contigs #hg38_chrY # hg38_chr21 hg38_chr1 hg38_chr15 hg38_chr18 hg38_chrX  #hg38_chr1_2 hg38_chr6_15_18_X_Y # hg38_chr1 hg38_chr15 hg38_chr18 hg38_chr21 hg38_chrX #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
#     do
#         if [[ ! -f $outroot/$chr_id/$read_length.sam ]]
#         then
#             # echo "SIMULATING VARIANTS"
#             mkdir -p $outroot/$chr_id/
#             # mason_variator --sv-indel-rate 0.000005 --snp-rate 0.001 --small-indel-rate 0.0001 --max-small-indel-size 50   -ir $refs/$chr_id.fa -ov $refs/$chr_id.vcf &> /dev/null
#             echo "SIMULATING READS"
#             # mason_simulator -ir $refs/$chr_id.fa -iv $refs/$chr_id.vcf -n 50000 --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
            
#             # FOR CONTIGS
#           if  ((read_length >= 250));
#             then  
#             echo mason_simulator -ir $refs/$chr_id.fa -n 100000 --illumina-read-length $read_length --fragment-mean-size 700 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
#             mason_simulator -ir $refs/$chr_id.fa -n 100000 --illumina-read-length $read_length --fragment-mean-size 700 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
#             else
#             echo mason_simulator -ir $refs/$chr_id.fa -n 100000 --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
#             mason_simulator -ir $refs/$chr_id.fa -n 100000 --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
#             fi
#         fi

#         # # minimap2 stats
#         # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2.sam 2>  $outroot/$chr_id/$read_length.minimap2.stderr
#         # # /usr/bin/time -l minimap2 -t 1 --eqx -ax map-hifi $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2.sam 2>  $outroot/$chr_id/$read_length.minimap2.stderr
#         # echo -n $chr_id,$read_length,minimap2,align,
#         # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2.sam --time_mem $outroot/$chr_id/$read_length.minimap2.stderr

#         # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         # echo -n $chr_id,$read_length,strobealign,align,
#         # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr


#         # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -k 22 -l 0 -u 10 -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         # echo -n $chr_id,$read_length,strobealign,align,8,22,0,10,
#         # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

#         # for bc in 8 #16 # 32
#         # do
#         #     for k in 15 17 20 #19 20 #18 19 20 #18 19 20 21
#         #     do
#         #         for l in -2
#         #         do
#         #             for u in 3 7 #3 5 #9 #5 7 9
#         #             do

#         #                 /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -c $bc -k $k -l $l -u $u -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         #                 echo -n $chr_id,$read_length,strobealign,align,$bc,$k,$l,$u,
#         #                 python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

#         #                 # Mapping stats
#         #                 # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -k $k -x -o $outroot/$chr_id/$read_length.strobealign.paf $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         #                 # echo -n $chr_id,$read_length,strobealign_k20,map,
#         #                 # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_paf $outroot/$chr_id/$read_length.strobealign.paf --time_mem $outroot/$chr_id/$read_length.strobealign.stderr
#         #             done
#         #         done
#         #     done
#         #     echo
#         # done
#         # echo "NEW CHR"
#         # echo
#     done
# done



for read_length in 100 150 200 250 300 # 
do 
    for chr_id in  sim_repeat_genome500 hg38_chr21 hg38_chr18 hg38_chrX hg38_chr1 hg38_chr15 # sim_50contigs # hg38_chr1_2 hg38_chr6_15_18_X_Y  #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
    do
        if [[ ! -f $outroot/$chr_id.vcf ]]
        then
            echo "SIMULATING VARIANTS"
            mkdir -p $outroot/$chr_id/
            # SIM4 rates
            echo   mason_variator --sv-indel-rate $indel_rate --snp-rate $snv_rate --small-indel-rate $small_indel_rate --max-small-indel-size 100   -ir $refs/$chr_id.fa -ov $outroot/$chr_id.vcf
            mason_variator --sv-indel-rate $indel_rate --snp-rate $snv_rate --small-indel-rate $small_indel_rate  --max-small-indel-size 100   -ir $refs/$chr_id.fa -ov $outroot/$chr_id.vcf &> /dev/null
        fi

        if [[ ! -f $outroot/$chr_id/$read_length.sam ]]
        then
            echo "SIMULATING READS"
          if  ((read_length >= 250));
            then  
            echo mason_simulator -ir $refs/$chr_id.fa -n 100000 -iv $outroot/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 700 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
            mason_simulator -ir $refs/$chr_id.fa -n 100000 -iv $outroot/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 700 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
            else
            echo mason_simulator -ir $refs/$chr_id.fa -n 100000 -iv $outroot/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
            mason_simulator -ir $refs/$chr_id.fa -n 100000 -iv $outroot/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
            fi
        fi

        # HIGH ERROR

        # # minimap2 stats
        # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/high_error.$read_length.L.fq $outroot/$chr_id/high_error.$read_length.R.fq 1> $outroot/$chr_id/high_error.$read_length.minimap2.sam 2>  $outroot/$chr_id/high_error.$read_length.minimap2.stderr
        # echo -n $chr_id,$read_length,minimap2,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/high_error.$read_length.sam --predicted_sam $outroot/$chr_id/high_error.$read_length.minimap2.sam --time_mem $outroot/$chr_id/high_error.$read_length.minimap2.stderr

        # /usr/bin/time   -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length -o $outroot/$chr_id/high_error.$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/high_error.$read_length.L.fq $outroot/$chr_id/high_error.$read_length.R.fq &>  $outroot/$chr_id/high_error.$read_length.strobealign.stderr
        # echo -n $chr_id,$read_length,strobealign,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/high_error.$read_length.sam --predicted_sam $outroot/$chr_id/high_error.$read_length.strobealign.sam --time_mem $outroot/$chr_id/high_error.$read_length.strobealign.stderr
        

        # NORMAL ERROR
        # # minimap2 stats
        # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/sim_500contigs.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2.sam 2>  $outroot/$chr_id/$read_length.minimap2.stderr

        # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2.sam 2>  $outroot/$chr_id/$read_length.minimap2.stderr
        # echo -n $chr_id,$read_length,minimap2,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2.sam --time_mem $outroot/$chr_id/$read_length.minimap2.stderr

        # /usr/bin/time   -l $strobealign_dev_dir/./strobealign3 -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
        # echo -n $chr_id,$read_length,strobealign-0.3,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr
        

        # /usr/bin/time   -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/sim_500contigs.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr

        ################################
        ######### PAIRED END ###########
        # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2_PE.sam 2>  $outroot/$chr_id/$read_length.minimap2_PE.stderr
        # echo -n $chr_id,$read_length,minimap2_PE,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2_PE.sam --time_mem $outroot/$chr_id/$read_length.minimap2_PE.stderr

        # /usr/bin/time  -l strobealign-v0.4 -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign_0.4.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign_0.4.stderr
        # echo -n $chr_id,$read_length,strobealign-0.4,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.4.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.4.stderr

        /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_0.6_PM.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.6_PM.stderr
        echo -n $chr_id,$read_length,strobealign-0.6.2,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.6_PM.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.6_PM.stderr

        /usr/bin/time  -l $strobealign_dev_dir/./strobealign-v0.6.1 -t 1 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_0.6-full_ssw.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.6-full_ssw.stderr
        echo -n $chr_id,$read_length,strobealign-0.6.1,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.6-full_ssw.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.6-full_ssw.stderr

        /usr/bin/time  -l strobealign-v0.6 -t 1 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_0.6-full_ssw.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.6-full_ssw.stderr
        echo -n $chr_id,$read_length,strobealign-0.6,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.6-full_ssw.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.6-full_ssw.stderr

        # /usr/bin/time  -l strobealign-v0.5 -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign_0.5.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign_0.5.stderr
        # echo -n $chr_id,$read_length,strobealign-0.5,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.5.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.5.stderr


        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -N 5 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_0.6_MM.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.6_MM.stderr
        # echo -n $chr_id,$read_length,strobealign-0.6_MM,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.6_MM.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.6_MM.stderr

        ################################
        ################################

        echo 

        # # #####################################
        # # # ########### SINGLE END ###########

        # # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.minimap2_SE.sam 2>  $outroot/$chr_id/$read_length.minimap2_SE.stderr
        # # echo -n $chr_id,$read_length,minimap2_SE,align,
        # # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2_SE.sam --time_mem $outroot/$chr_id/$read_length.minimap2_SE.stderr

        # /usr/bin/time  -l strobealign-v0.4 -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign_0.4.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq &>  $outroot/$chr_id/$read_length.strobealign_0.4.stderr
        # echo -n $chr_id,$read_length,strobealign-0.4-SE,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.4.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.4.stderr

        # /usr/bin/time  -l strobealign-v0.5 -t 1 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.strobealign_0.5.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.5.stderr
        # echo -n $chr_id,$read_length,strobealign-0.5-SE,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.5.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.5.stderr

        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.strobealign_0.6_SE_PM.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.6_SE_PM.stderr
        # echo -n $chr_id,$read_length,strobealign-0.6-SE_PM,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.6_SE_PM.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.6_SE_PM.stderr

        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -N 5 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.strobealign_0.6_SE_MM.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.6_SE_MM.stderr
        # echo -n $chr_id,$read_length,strobealign-0.6-SE_MM,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.6_SE_MM.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.6_SE_MM.stderr

        # # #################################
        # # #################################
        
        echo

        # /usr/bin/time   -l $strobealign_dev_dir/./strobealign -t 1 -x -r $read_length -o $outroot/$chr_id/$read_length.strobealign.paf $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
        # echo -n $chr_id,$read_length,strobealign-0.4,map,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_paf $outroot/$chr_id/$read_length.strobealign.paf --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

                

        # python $eval_script_dir/get_erroneous.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam_method1 $outroot/$chr_id/$read_length.minimap2.sam \
        #                                          --predicted_sam_method2 $outroot/$chr_id/$read_length.strobealign.sam --fq1 $outroot/$chr_id/$read_length.L.fq --fq2 $outroot/$chr_id/$read_length.R.fq --om $outroot/$chr_id/$read_length.misaligned \
        #                                          --ou $outroot/$chr_id/$read_length.misaligned

        # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -R 0 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
        # echo -n $chr_id,$read_length,strobealign,NO_RESC,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

        # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -x -r $read_length -o $outroot/$chr_id/$read_length.strobealign.paf $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
        # echo -n $chr_id,$read_length,strobealign,map,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_paf $outroot/$chr_id/$read_length.strobealign.paf --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

        # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -x -R 0 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.paf $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
        # echo -n $chr_id,$read_length,strobealign,map_NO_RESC,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_paf $outroot/$chr_id/$read_length.strobealign.paf --time_mem $outroot/$chr_id/$read_length.strobealign.stderr
    done
    echo
    echo
done


###### FOR BUGFIXING #########

# for read_length in 100 #150 200 250 300 
# do 
#     for chr_id in hg38_chr18 # hg38_chr1 hg38_chr15 hg38_chr18 hg38_chr21 hg38_chrX # hg38_chr1_2 hg38_chr6_15_18_X_Y  #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
#     do

#         # minimap2 stats
#         /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.misaligned_1.fq $outroot/$chr_id/$read_length.misaligned_2.fq 1> $outroot/$chr_id/$read_length.minimap2.sam 2>  $outroot/$chr_id/$read_length.minimap2.stderr
#         echo -n $chr_id,$read_length,minimap2,align,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2.sam --time_mem $outroot/$chr_id/$read_length.minimap2.stderr

#         /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.misaligned_1.fq $outroot/$chr_id/$read_length.misaligned_2.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         echo -n $chr_id,$read_length,strobealign,align,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

#     done
#     echo
# done


# #### FOR LONG READS ############

# for read_length in 1000 
# do 
#     for chr_id in hg38_chr18 #hg38_chrX hg38_chr21 hg38_chr1 hg38_chr15 # hg38_chr1_2 hg38_chr6_15_18_X_Y  #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
#     do

#         if [[ ! -f $outroot/$chr_id/$read_length.sam ]]
#         then
#             echo "SIMULATING VARIANTS"
#             mkdir -p $outroot/$chr_id/
#             mason_variator --sv-indel-rate 0.000005 --snp-rate 0.001 --small-indel-rate 0.0002 --max-small-indel-size 50   -ir $refs/$chr_id.fa -ov $refs/$chr_id.vcf &> /dev/null
#             echo "SIMULATING READS"
#               mason_simulator -ir $refs/$chr_id.fa -n 10000 --seq-technology 454 \
#                 --454-read-length-mean $read_length --454-read-length-stddev 100  \
#                 --fragment-mean-size 2000 --fragment-size-std-dev 200 \
#                 -o $outroot/$chr_id/$read_length.fq -oa $outroot/$chr_id/$read_length.sam  
#             # echo mason_simulator --seq-technology sanger -ir $refs/$chr_id.fa -iv $refs/$chr_id.vcf --sanger-prob-mismatch-scale 2.0 -n 100000 --sanger-read-length-min 500 --sanger-read-length-max 3000 --sanger-read-length-mean $read_length -o $outroot/$chr_id/$read_length.fq -oa $outroot/$chr_id/$read_length.sam  
#             # mason_simulator --seq-technology sanger -ir $refs/$chr_id.fa -iv $refs/$chr_id.vcf --sanger-prob-mismatch-scale 2.0 -n 100000 --sanger-read-length-min 500 --sanger-read-length-max 3000 --sanger-read-length-mean $read_length -o $outroot/$chr_id/$read_length.fq -oa $outroot/$chr_id/$read_length.sam  
#         fi

#         # HIGH ERROR

#         # minimap2 stats
#         /usr/bin/time -l minimap2 -t 1 -ax map-hifi --eqx $refs/$chr_id.fa $outroot/$chr_id/$read_length.fq 1> $outroot/$chr_id/$read_length.minimap2.sam 2>  $outroot/$chr_id/$read_length.minimap2.stderr
#         echo -n $chr_id,$read_length,minimap2,align,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2.sam --time_mem $outroot/$chr_id/$read_length.minimap2.stderr

#         /usr/bin/time   -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         echo -n $chr_id,$read_length,strobealign,align,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr      
#     done
# done
