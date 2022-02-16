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


echo -n  "tool","ref","%-aligned","accuracy,time(sec),Mem(MB)"$'\n'

chr_id="sim_repeat_genome500"
if [[ ! -f $refs/$chr_id.vcf ]]; then
    echo "SIMULATING VARIANTS"
    mkdir -p $outroot/$chr_id/
    echo   mason_variator --sv-indel-rate 0.00005 --snp-rate 0.005 --small-indel-rate 0.005 --max-small-indel-size 50   -ir $refs/$chr_id.fa -ov $refs/$chr_id.vcf
    mason_variator --sv-indel-rate 0.00005 --snp-rate 0.005 --small-indel-rate 0.005 --max-small-indel-size 50   -ir $refs/$chr_id.fa -ov $refs/$chr_id.vcf &> /dev/null
fi

bcftools sort -Oz $refs/$chr_id.vcf -o $refs/$chr_id.sorted.vcf.gz
bcftools index $refs/$chr_id.sorted.vcf.gz

true_indel=$(grep "sim_small_indel" $refs/$chr_id.vcf | wc -l)
true_snp=$(grep "sim_snp" $refs/$chr_id.vcf | wc -l)
echo "TRUE INDELS" "${true_indel}"
echo "TRUE SNPS" "${true_snp}"

for read_length in 100 150 200 250 300 # 
do 
    for chr_id in sim_repeat_genome500 # sim_repeat_genome5 #hg38_chr21 hg38_chr18 hg38_chrX hg38_chr1 hg38_chr15 # sim_repeat_genome500 # sim_50contigs # hg38_chr1_2 hg38_chr6_15_18_X_Y  #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
    do

        if [[ ! -f $outroot/$chr_id/$read_length.sam ]]
        then    
            
          echo "SIMULATING READS"
          if  ((read_length >= 250));
            then  
                echo mason_simulator -ir $refs/$chr_id.fa -n 2000000 -iv $refs/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 700 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
                mason_simulator -ir $refs/$chr_id.fa -n 2000000 -iv $refs/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 700 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
            else
                echo mason_simulator -ir $refs/$chr_id.fa -n 2000000 -iv $refs/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
                mason_simulator -ir $refs/$chr_id.fa -n 2000000 -iv $refs/$chr_id.vcf --illumina-read-length $read_length --fragment-mean-size 300 -o $outroot/$chr_id/$read_length.L.fq -or $outroot/$chr_id/$read_length.R.fq -oa $outroot/$chr_id/$read_length.sam
          fi
        fi


        ################################
        ######### PAIRED END ###########

        # MINIMAP2
        /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2_PE.sam 2>  $outroot/$chr_id/$read_length.minimap2_PE.stderr
        echo -n $chr_id,$read_length,minimap2_PE,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2_PE.sam --time_mem $outroot/$chr_id/$read_length.minimap2_PE.stderr

        # sam to sorted indexed bam
        samtools view -u $outroot/$chr_id/$read_length.minimap2_PE.sam | samtools sort -o $outroot/$chr_id/$read_length.minimap2_PE.bam &> /dev/null
        rm $outroot/$chr_id/$read_length.minimap2_PE.sam
        samtools index $outroot/$chr_id/$read_length.minimap2_PE.bam &> /dev/null

        # call variants
        # echo bcftools mpileup -O v --fasta-ref $refs/$chr_id.fa $outroot/$chr_id/$read_length.minimap2_PE.bam > $outroot/$chr_id/$read_length.minimap2_PE.vcf
        bcftools mpileup -O z --fasta-ref $refs/$chr_id.fa $outroot/$chr_id/$read_length.minimap2_PE.bam > $outroot/$chr_id/$read_length.minimap2_PE.vcf.gz 2> /dev/null
        bcftools call -v -c -O v $outroot/$chr_id/$read_length.minimap2_PE.vcf.gz > $outroot/$chr_id/$read_length.minimap2_PE_variants.vcf 2> /dev/null

        bcftools sort -Oz $outroot/$chr_id/$read_length.minimap2_PE_variants.vcf -o $outroot/$chr_id/$read_length.minimap2_PE_variants.sorted.vcf.gz &> /dev/null
        bcftools index $outroot/$chr_id/$read_length.minimap2_PE_variants.sorted.vcf.gz &> /dev/null
        bcftools isec -O u $refs/$chr_id.sorted.vcf.gz  $outroot/$chr_id/$read_length.minimap2_PE_variants.sorted.vcf.gz -p $outroot/$chr_id/$read_length.minimap2_PE_variants_ovl &> /dev/null

        pred_snp="$(grep -E "\t[ACGT]\t[ACGT]\t" $outroot/$chr_id/$read_length.minimap2_PE_variants.vcf | wc -l)"
        pred_indel="$(grep "INDEL" $outroot/$chr_id/$read_length.minimap2_PE_variants.vcf  | wc -l)"
        ovl_total="$(cat $outroot/$chr_id/$read_length.minimap2_PE_variants_ovl/sites.txt | wc -l)"
        ovl_snp="$(grep -E "\t[ACGT]\t[ACGT]\t" $outroot/$chr_id/$read_length.minimap2_PE_variants_ovl/sites.txt| wc -l)"
        ovl_indel=$(( ovl_total - ovl_snp))
        echo $chr_id,$read_length,minimap2_PE,align,$pred_snp,$pred_indel,$ovl_snp,$ovl_indel



        # # STROBEALIGN v0.4
        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign-0.4 -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign_PE_0.4.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq &>  $outroot/$chr_id/$read_length.strobealign_PE_0.4.stderr
        # echo -n $chr_id,$read_length,strobealign-0.4,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_PE_0.4.sam --time_mem $outroot/$chr_id/$read_length.strobealign_PE_0.4.stderr

        # # sam to sorted indexed bam
        # samtools view -u $outroot/$chr_id/$read_length.strobealign_PE_0.4.sam | samtools sort -o $outroot/$chr_id/$read_length.strobealign_PE_0.4.bam
        # samtools index $outroot/$chr_id/$read_length.strobealign_PE_0.4.bam

        # # call variants
        # bcftools mpileup -O v --fasta-ref $refs/$chr_id.fa $outroot/$chr_id/$read_length.strobealign_PE_0.4.bam > $outroot/$chr_id/$read_length.strobealign_PE_0.4.vcf
        # bcftools call -v -c -O v $outroot/$chr_id/$read_length.strobealign_PE_0.4.vcf > $outroot/$chr_id/$read_length.strobealign_PE_0.4_variants.vcf



        # STROBEALIGN v0.5
        /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_PE_0.5.sam 2>  $outroot/$chr_id/$read_length.strobealign_PE_0.5.stderr
        echo -n $chr_id,$read_length,strobealign-0.5_PE,align,
        python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_PE_0.5.sam --time_mem $outroot/$chr_id/$read_length.strobealign_PE_0.5.stderr

        # sam to sorted indexed bam
        samtools view -u $outroot/$chr_id/$read_length.strobealign_PE_0.5.sam | samtools sort -o $outroot/$chr_id/$read_length.strobealign_PE_0.5.bam &> /dev/null
        rm $outroot/$chr_id/$read_length.strobealign_PE_0.5.sam
        samtools index $outroot/$chr_id/$read_length.strobealign_PE_0.5.bam &> /dev/null

        # call variants
        bcftools mpileup -O z --fasta-ref $refs/$chr_id.fa $outroot/$chr_id/$read_length.strobealign_PE_0.5.bam > $outroot/$chr_id/$read_length.strobealign_PE_0.5.vcf.gz 2> /dev/null
        bcftools call -v -c -O v $outroot/$chr_id/$read_length.strobealign_PE_0.5.vcf.gz > $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants.vcf 2> /dev/null

        bcftools sort -Oz $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants.vcf -o $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants.sorted.vcf.gz &> /dev/null
        bcftools index $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants.sorted.vcf.gz &> /dev/null
        bcftools isec -O u $refs/$chr_id.sorted.vcf.gz  $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants.sorted.vcf.gz -p $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants_ovl &> /dev/null

        pred_snp="$(grep -E "\t[ACGT]\t[ACGT]\t" $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants.vcf | wc -l)"
        pred_indel="$(grep "INDEL" $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants.vcf  | wc -l)"
        ovl_total="$(cat $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants_ovl/sites.txt | wc -l)"
        ovl_snp="$(grep -E "\t[ACGT]\t[ACGT]\t" $outroot/$chr_id/$read_length.strobealign_PE_0.5_variants_ovl/sites.txt| wc -l)"
        ovl_indel=$(( ovl_total - ovl_snp))
        echo $chr_id,$read_length,strobealign-0.5_PE,align,$pred_snp,$pred_indel,$ovl_snp,$ovl_indel

        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -A 2 -B 8 -O 12 -N 5 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_PE_MM.sam 2>  $outroot/$chr_id/$read_length.strobealign_PE_MM.stderr
        # echo -n $chr_id,$read_length,strobealign-0.4.1_MM,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_PE_MM.sam --time_mem $outroot/$chr_id/$read_length.strobealign_PE_MM.stderr

        ################################
        ################################

        # echo 

        #####################################
        # ########### SINGLE END ###########

        # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.minimap2_SE.sam 2>  $outroot/$chr_id/$read_length.minimap2_SE.stderr
        # echo -n $chr_id,$read_length,minimap2_SE,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.minimap2_SE.sam --time_mem $outroot/$chr_id/$read_length.minimap2_SE.stderr

        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign-0.4 -t 1 -r $read_length -o $outroot/$chr_id/$read_length.strobealign_SE.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq &>  $outroot/$chr_id/$read_length.strobealign_SE.stderr
        # echo -n $chr_id,$read_length,strobealign-0.4_SE,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_SE.sam --time_mem $outroot/$chr_id/$read_length.strobealign_SE.stderr

        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -A 2 -B 8 -O 12 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.strobealign_SE_PM.sam 2>  $outroot/$chr_id/$read_length.strobealign_SE_PM.stderr
        # echo -n $chr_id,$read_length,strobealign-0.5-SE_PM,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_SE_PM.sam --time_mem $outroot/$chr_id/$read_length.strobealign_SE_PM.stderr

        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -N 5 -A 2 -B 8 -O 12 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.strobealign_SE_MM.sam 2>  $outroot/$chr_id/$read_length.strobealign_SE_MM.stderr
        # echo -n $chr_id,$read_length,strobealign-0.5-SE_MM,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_SE_MM.sam --time_mem $outroot/$chr_id/$read_length.strobealign_SE_MM.stderr

        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -N 5 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq 1> $outroot/$chr_id/$read_length.strobealign_SE_MM.sam 2>  $outroot/$chr_id/$read_length.strobealign_SE_MM.stderr
        # echo -n $chr_id,$read_length,strobealign-0.4.1-SE_MM,align,
        # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.strobealign_SE_MM.sam --time_mem $outroot/$chr_id/$read_length.strobealign_SE_MM.stderr

        #################################
        #################################
        
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


