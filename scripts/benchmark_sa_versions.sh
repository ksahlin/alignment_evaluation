#!/bin/bash

# RUN scripts e.g. as: ./benchmark_sa_versions.sh /Users/ksahlin/prefix/data/genomes/subset/ /Users/ksahlin/tmp/STROBEALIGN/

if [ $# -lt 2 ]; then
    # TODO: print usage
    echo "./benchmark_sa_versions.sh <ref_path> <outroot>"
    exit 1
fi

refs=$1
outroot=$2
eval_script_dir="/Users/ksahlin/prefix/source/alignment_evaluation/scripts"
echo $outroot
mkdir -p $outroot


snv_rate="0.001"
small_indel_rate="0.0005"
indel_rate="0.00001"


echo -n  "tool","ref","%-aligned","accuracy,time(sec),Mem(MB)"$'\n'


# TBD: Several E-coli strains? 

# Human PE-reads

for read_length in 75 100 125 #150 250 
do 
    for chr_id in  hg38_chr1_2_3 hg38_chr15 # #hg38_chr21 hg38_chrX # hg38_chr18  hg38_chr15 hg38_chr1
    do
        if [[ ! -f $outroot/$chr_id.vcf ]]
        then
            echo "SIMULATING VARIANTS"
            mkdir -p $outroot/$chr_id/
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
   

        ################################
        ######### PAIRED END ###########
        
        if [[ ! -f $refs/$chr_id.fa.bwt ]]
        then
            echo "BWA MEM INDEXING"
            bwa index $refs/$chr_id.fa
        fi

        /usr/bin/time  -l strobealign-master-4735d37 -t 1 $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_master-4735d37.sam 2>  $outroot/$chr_id/$read_length.strobealign_master-4735d37.stderr
        /usr/bin/time  -l strobealign-v0.7.1 -t 1 $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_0.7.1.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.7.1.stderr
        /usr/bin/time  -l strobealign-master-4735d37 -r $read_length -i $outroot/$chr_id/$read_length.strobealign_master-4735d37.sti -t 1 $refs/$chr_id.fa 2>  $outroot/$chr_id/$read_length.strobealign_master-4735d37-indexing.stderr
        /usr/bin/time  -l strobealign-master-4735d37  -t 1 $outroot/$chr_id/$read_length.strobealign_master-4735d37.sti $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_master-4735d37-preindexed.sam 2>  $outroot/$chr_id/$read_length.strobealign_master-4735d37-preindexed.stderr


        # /usr/bin/time  -l strobealign-master-4735d37 -i $outroot/$chr_id/$read_length.strobealign_master-4735d37.sti -t 1 $refs/$chr_id.fa 2>  $outroot/$chr_id/$read_length.strobealign_master-4735d37-indexing.stderr
        # /usr/bin/time  -l strobealign-master-4735d37 -t 1 $outroot/$chr_id/$read_length.strobealign_master-4735d37.sti $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_master-4735d37-preindexed.sam 2>  $outroot/$chr_id/$read_length.strobealign_master-4735d37-preindexed.stderr
        # /usr/bin/time -l  bwa mem -t 1 -o $outroot/$chr_id/$read_length.bwa_mem_PE.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 2>  $outroot/$chr_id/$read_length.bwa_mem_PE.stderr
        # /usr/bin/time  -l strobealign-master-4735d37 -t 1 $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_master-4735d37.sam 2>  $outroot/$chr_id/$read_length.strobealign_master-4735d37.stderr
        # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2_PE.sam 2>  $outroot/$chr_id/$read_length.minimap2_PE.stderr
        # /usr/bin/time  -l strobealign-v0.7.1 -t 1 $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_0.7.1.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.7.1.stderr

        # echo -n $chr_id,$read_length,bwa_mem,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.bwa_mem_PE.sam --time_mem $outroot/$chr_id/$read_length.bwa_mem_PE.stderr
        # echo -n $chr_id,$read_length,minimap2_PE,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2_PE.sam --time_mem $outroot/$chr_id/$read_length.minimap2_PE.stderr
        
        echo -n $chr_id,$read_length,strobealign-v0.7.1,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.7.1.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.7.1.stderr
        echo -n $chr_id,$read_length,strobealign-master,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_master-4735d37.sam --time_mem $outroot/$chr_id/$read_length.strobealign_master-4735d37.stderr
        echo -n $chr_id,$read_length,strobealign-master-preindexed,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_master-4735d37-preindexed.sam --time_mem $outroot/$chr_id/$read_length.strobealign_master-4735d37-preindexed.stderr

        ################################
        ################################

        echo 

    done
    echo
    echo
done


# # Human SE-reads


# for read_length in 75 100 125 150 250 # 
# do 
#     for chr_id in hg38_chr1_2_3 hg38_chr21 hg38_chr18 hg38_chrX hg38_chr15 hg38_chr1
#     do
#         if [[ ! -f $outroot/$chr_id.vcf ]]
#         then
#             echo "SIMULATING VARIANTS"
#             mkdir -p $outroot/$chr_id/
#             mason_variator --sv-indel-rate $indel_rate --snp-rate $snv_rate --small-indel-rate $small_indel_rate  --max-small-indel-size 100   -ir $refs/$chr_id.fa -ov $outroot/$chr_id.vcf &> /dev/null
#         fi

#         if [[ ! -f $outroot/$chr_id/$read_length.sam ]]
#         then
#             echo "SIMULATING READS"
#           if  ((read_length >= 250));
#             then  
#             echo mason_simulator -ir $refs/$chr_id.fa -n 100000 -iv $outroot/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 700 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
#             mason_simulator -ir $refs/$chr_id.fa -n 100000 -iv $outroot/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 700 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
#             else
#             echo mason_simulator -ir $refs/$chr_id.fa -n 100000 -iv $outroot/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
#             mason_simulator -ir $refs/$chr_id.fa -n 100000 -iv $outroot/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
#             fi
#         fi
 

#         #####################################
#         # ########### SINGLE END ###########

#         /usr/bin/time  -l strobealign-v0.7.1 -t 1 $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.strobealign_0.7.1_SE.sam 2>  $outroot/$chr_id/$read_length.strobealign_0.7.1_SE.stderr
#         echo -n $chr_id,$read_length,strobealign-v0.7.1_SE,align,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_0.7.1_SE.sam --time_mem $outroot/$chr_id/$read_length.strobealign_0.7.1_SE.stderr

#         /usr/bin/time  -l strobealign-master-4735d37 -t 1 $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.strobealign_master-4735d37_SE.sam 2>  $outroot/$chr_id/$read_length.strobealign_master-4735d37_SE.stderr
#         echo -n $chr_id,$read_length,strobealign-master_SE,align,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_master-4735d37_SE.sam --time_mem $outroot/$chr_id/$read_length.strobealign_master-4735d37_SE.stderr

#         /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.minimap2_SE.sam 2>  $outroot/$chr_id/$read_length.minimap2_SE.stderr
#         echo -n $chr_id,$read_length,minimap2_SE,align,
#         python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2_SE.sam --time_mem $outroot/$chr_id/$read_length.minimap2_SE.stderr


#         #################################
#         #################################
        
#         echo


#     done
#     echo
#     echo
# done


