#!/bin/sh
#SBATCH -A snic2020-5-651
#SBATCH -N 1
#SBATCH -C mem128GB
#SBATCH -c 20
#SBATCH --time=02-00:00:00
# Memory per node specification is in MB. It is optional. 
# The default limit is 3000MB per core.
#SBATCH --job-name="sa_300"
#SBATCH --mail-user=ksahlin@kth.se
#SBATCH --mail-type=ALL

#conda init bash
set -o errexit

########################
### EDIT THESE LINES ###

# 100
read_length="100"
bc_sizes=`seq 8 8 16`
k_sizes=$(seq 18 21)
offsets=$(seq 1 2)
spans=$(seq 4 3 7)

# # 150
# read_length="150"
# bc_sizes=`seq 8 8 16`
# k_sizes=$(seq 18 21)
# offsets=`seq 3 2 5`
# spans=`seq 6 3 9`

# #200
# read_length="200"
# bc_sizes=`seq 8 8 16`
# k_sizes=$(seq 21 24)
# offsets=`seq 7 2 9`
# spans=`seq 6 3 9`

# #250
# read_length="200"
# bc_sizes=`seq 8 8 16`
# k_sizes=$(seq 21 24)
# offsets=`seq 7 2 9`
# spans=`seq 6 3 9`

# #300
# read_length="300"
# bc_sizes=`seq 8 8 16`
# k_sizes=$(seq 20 23)
# offsets=`seq 5 2 7`
# spans=`seq 6 3 9`
########################
########################

outroot="/proj/snic2020-16-138/strobemap_eval/tmp/"
eval_script_dir="/home/kris/source/alignment_evaluation/scripts"
echo $outroot

mkdir -p $outroot

hg38="/proj/snic2020-16-138/ultra_eval/genomes/Homo_sapiens.GRCh38.dna.primary_assembly_modified_headers.fa"
reads_dir="/proj/snic2020-16-138/strobemap_eval/reads_PE"


# /usr/bin/time -v strobealign -t 1 -r 250 -o /proj/snic2020-16-138/strobemap_eval/tmp/SIM3_250.sam /proj/snic2020-16-138/ultra_eval/genomes/Homo_sapiens.GRCh38.dna.primary_assembly_modified_headers.fa 
#                                       /proj/snic2020-16-138/strobemap_eval/reads_PE/SIM3/250_L.fq /proj/snic2020-16-138/strobemap_eval/reads_PE/SIM3/250_R.fq 2>&1 | tee /proj/snic2020-16-138/strobemap_eval/tmp/SIM3_250.out

echo -n  "tool","ref","%-aligned","accuracy,time(sec),Mem(MB)"$'\n'

for dataset in SIM1 SIM2 SIM3 
do 
    echo
    echo $dataset
    echo
    mkdir -p $outroot/$dataset/$read_length
    strobealign_pred=$outroot/$dataset/$read_length/v0.2.strobealign.sam
    truth=$reads_dir/$dataset/$read_length.sam

    /usr/bin/time -v strobealign -t 8 -r $read_length -o $strobealign_pred $hg38 $reads_dir/$dataset/${read_length}_L.fq $reads_dir/$dataset/${read_length}_R.fq &>  $outroot/$dataset/$read_length/v0.2.strobealign.stderr
    echo -n default,$read_length,strobealign,align,
    python $eval_script_dir/get_stats_linux.py --truth $truth --predicted_sam $strobealign_pred --time_mem $outroot/$dataset/$read_length/v0.2.strobealign.stderr
    rm $strobealign_pred
    echo


    for k in $k_sizes
    do
        for offset in $offsets
        do
            l=$(($offset - $k/5))
            for span in $spans
            do
                u=$(($l + $span))
                for bc in $bc_sizes
                do
                    # echo $k,$l,$u,$bc 
                    strobealign_pred=$outroot/$dataset/$read_length/$k.$l.$u.$bc.strobealign.sam
                    /usr/bin/time -v strobealign -t 8 -k $k -l $l -u $u -c $bc -o $strobealign_pred $hg38 $reads_dir/$dataset/${read_length}_L.fq $reads_dir/$dataset/${read_length}_R.fq &>  $outroot/$dataset/$read_length/$k.$l.$u.$bc.strobealign.stderr
                    echo -n $k,$l,$u,$bc,$read_length,strobealign,align,
                    python $eval_script_dir/get_stats_linux.py --truth $truth --predicted_sam $strobealign_pred --time_mem $outroot/$dataset/$read_length/$k.$l.$u.$bc.strobealign.stderr 2> $outroot/$dataset/$read_length/$k.$l.$u.$bc.analysis.stderr
                    rm $strobealign_pred
                done
            done
        done
        echo
    done

done


# for read_length in 100 150 #200 250 300 
# do 
#     for chr_id in hg38_chr18 hg38_chrX hg38_chr21 hg38_chr1 hg38_chr15 # hg38_chr1_2 hg38_chr6_15_18_X_Y  #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
#     do

#         # minimap2 stats
#         # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2.sam 2>  $outroot/$chr_id/$read_length.minimap2.stderr
#         # echo -n $read_length,minimap2,align,
#         # python $eval_script_dir/get_stats_linux.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2.sam --time_mem $outroot/$chr_id/$read_length.minimap2.stderr

#         /usr/bin/time   -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         echo -n $chr_id,$read_length,strobealign,align,
#         python $eval_script_dir/get_stats_linux.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr
#         # python $eval_script_dir/get_erroneous.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam_method1 $outroot/$chr_id/$read_length.minimap2.sam \
#         #                                          --predicted_sam_method2 $outroot/$chr_id/$read_length.strobealign.sam --fq1 $outroot/$chr_id/$read_length.L.fq --fq2 $outroot/$chr_id/$read_length.R.fq --om $outroot/$chr_id/$read_length.misaligned \
#         #                                          --ou $outroot/$chr_id/$read_length.misaligned

#         # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -R 0 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         # echo -n $chr_id,$read_length,strobealign,NO_RESC,
#         # python $eval_script_dir/get_stats_linux.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

#         # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -x -r $read_length -o $outroot/$chr_id/$read_length.strobealign.paf $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         # echo -n $chr_id,$read_length,strobealign,map,
#         # python $eval_script_dir/get_stats_linux.py --truth $outroot/$chr_id/$read_length.sam --predicted_paf $outroot/$chr_id/$read_length.strobealign.paf --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

#         # /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -x -R 0 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.paf $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         # echo -n $chr_id,$read_length,strobealign,map_NO_RESC,
#         # python $eval_script_dir/get_stats_linux.py --truth $outroot/$chr_id/$read_length.sam --predicted_paf $outroot/$chr_id/$read_length.strobealign.paf --time_mem $outroot/$chr_id/$read_length.strobealign.stderr
#     done
#     echo
# done


###### FOR BUGFIXING #########

# for read_length in 100 #150 200 250 300 
# do 
#     for chr_id in hg38_chr18 # hg38_chr1 hg38_chr15 hg38_chr18 hg38_chr21 hg38_chrX # hg38_chr1_2 hg38_chr6_15_18_X_Y  #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
#     do

#         # minimap2 stats
#         /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.misaligned_1.fq $outroot/$chr_id/$read_length.misaligned_2.fq 1> $outroot/$chr_id/$read_length.minimap2.sam 2>  $outroot/$chr_id/$read_length.minimap2.stderr
#         echo -n $chr_id,$read_length,minimap2,align,
#         python $eval_script_dir/get_stats_linux.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2.sam --time_mem $outroot/$chr_id/$read_length.minimap2.stderr

#         /usr/bin/time -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.misaligned_1.fq $outroot/$chr_id/$read_length.misaligned_2.fq &>  $outroot/$chr_id/$read_length.strobealign.stderr
#         echo -n $chr_id,$read_length,strobealign,align,
#         python $eval_script_dir/get_stats_linux.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign.sam --time_mem $outroot/$chr_id/$read_length.strobealign.stderr

#     done
#     echo
# done


