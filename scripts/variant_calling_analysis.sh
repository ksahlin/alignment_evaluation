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

# Split TRUTH into SNP and INDELS
touch $refs/$chr_id.variants.SNV.vcf
grep "#"  $refs/$chr_id.vcf  > $refs/$chr_id.variants.SNV.vcf
grep -E "sim_snp" $refs/$chr_id.vcf  >> $refs/$chr_id.variants.SNV.vcf

touch $refs/$chr_id.variants.INDEL.vcf
grep "#"  $refs/$chr_id.vcf  > $refs/$chr_id.variants.INDEL.vcf
grep -E -e "sim_small_indel" $refs/$chr_id.vcf  >> $refs/$chr_id.variants.INDEL.vcf

true_indel=$(grep "sim_small_indel" $refs/$chr_id.vcf | wc -l)
true_snv=$(grep "sim_snp" $refs/$chr_id.vcf | wc -l)
echo "TRUE INDELS" "${true_indel}"
echo "TRUE SNPS" "${true_snv}"

bcftools sort -Oz $refs/$chr_id.variants.SNV.vcf -o $refs/$chr_id.sorted.SNV.vcf.gz
bcftools index $refs/$chr_id.sorted.SNV.vcf.gz

bcftools sort -Oz $refs/$chr_id.variants.INDEL.vcf -o $refs/$chr_id.sorted.INDEL.vcf.gz
bcftools index $refs/$chr_id.sorted.INDEL.vcf.gz


if [[ ! -f $refs/$chr_id.fa.bwt ]]; then
    bwa index $refs/$chr_id.fa
fi



for read_length in 100 150 200 250 300
do 
    for chr_id in sim_repeat_genome500 # sim_repeat_genome5 #hg38_chr21 hg38_chr18 hg38_chrX hg38_chr1 hg38_chr15 # sim_repeat_genome500 # sim_50contigs # hg38_chr1_2 hg38_chr6_15_18_X_Y  #  # hg38_chr1_2 hg38_chr6_15_18_X_Y  #
    do
        # Simulate reads if not already there
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


        ###### run the tools ######

        # /usr/bin/time  -l $strobealign_dev_dir/./strobealign -t 1 -r $read_length $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.strobealign_PE.sam 2>  $outroot/$chr_id/$read_length.strobealign_PE.stderr
        # /usr/bin/time -l minimap2 -t 1 --eqx -ax sr $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 1> $outroot/$chr_id/$read_length.minimap2_PE.sam 2>  $outroot/$chr_id/$read_length.minimap2_PE.stderr
        # /usr/bin/time -l  bwa mem -t 1 -o $outroot/$chr_id/$read_length.bwa_mem_PE.sam $refs/$chr_id.fa $outroot/$chr_id/$read_length.L.fq $outroot/$chr_id/$read_length.R.fq 2>  $outroot/$chr_id/$read_length.bwa_mem_PE.stderr

        for tool in strobealign_PE minimap2_PE bwa_mem_PE 
        do

            ################################
            ######### PAIRED END ###########

            # echo -n $chr_id,$read_length,$tool,align,
            # python $eval_script_dir/get_accuracy.py --truth $outroot/$chr_id/$read_length.sam --predicted_sam $outroot/$chr_id/$read_length.$tool.sam --time_mem $outroot/$chr_id/$read_length.$tool.stderr

            # # sam to sorted indexed bam
            # samtools view -u $outroot/$chr_id/$read_length.$tool.sam | samtools sort -o $outroot/$chr_id/$read_length.$tool.bam &> /dev/null
            # rm $outroot/$chr_id/$read_length.$tool.sam
            # samtools index $outroot/$chr_id/$read_length.$tool.bam &> /dev/null

            # # call variants
            # # echo bcftools mpileup -O v --fasta-ref $refs/$chr_id.fa $outroot/$chr_id/$read_length.$tool.bam > $outroot/$chr_id/$read_length.$tool.vcf
            # bcftools mpileup -O z --fasta-ref $refs/$chr_id.fa $outroot/$chr_id/$read_length.$tool.bam > $outroot/$chr_id/$read_length.$tool.vcf.gz 2> /dev/null
            # bcftools call -v -c -O v $outroot/$chr_id/$read_length.$tool.vcf.gz > $outroot/$chr_id/$read_length.$tool.variants.vcf 2> /dev/null

            # # Split into SNP and INDELS
            # grep -v -E -e "INDEL;"  $outroot/$chr_id/$read_length.$tool.variants.vcf > $outroot/$chr_id/$read_length.$tool.variants.SNV.vcf
            # touch $outroot/$chr_id/$read_length.$tool.variants.INDEL.vcf
            # grep "#"  $outroot/$chr_id/$read_length.$tool.variants.vcf > $outroot/$chr_id/$read_length.$tool.variants.INDEL.vcf
            # grep -E -e "INDEL;"  $outroot/$chr_id/$read_length.$tool.variants.vcf >> $outroot/$chr_id/$read_length.$tool.variants.INDEL.vcf

            # # Get true SNV calls
            # bcftools sort -Oz $outroot/$chr_id/$read_length.$tool.variants.SNV.vcf -o $outroot/$chr_id/$read_length.$tool.variants.sorted.SNV.vcf.gz #&> /dev/null
            # bcftools index $outroot/$chr_id/$read_length.$tool.variants.sorted.SNV.vcf.gz #&> /dev/null
            # bcftools isec --nfiles 2  -O u $refs/$chr_id.sorted.SNV.vcf.gz  $outroot/$chr_id/$read_length.$tool.variants.sorted.SNV.vcf.gz -p $outroot/$chr_id/$read_length.$tool.variants.SNV.ovl #&> /dev/null

            # # Get true INDEL calls
            # bcftools sort -Oz $outroot/$chr_id/$read_length.$tool.variants.INDEL.vcf -o $outroot/$chr_id/$read_length.$tool.variants.sorted.INDEL.vcf.gz &> /dev/null
            # bcftools index $outroot/$chr_id/$read_length.$tool.variants.sorted.INDEL.vcf.gz &> /dev/null
            # bcftools isec --nfiles 2 -O u $refs/$chr_id.sorted.INDEL.vcf.gz  $outroot/$chr_id/$read_length.$tool.variants.sorted.INDEL.vcf.gz -p $outroot/$chr_id/$read_length.$tool.variants.INDEL.ovl &> /dev/null

            echo -n \| $read_length  \| $tool \|  
            python $eval_script_dir/get_sv_prediction_stats.py -T $true_snv -A $tool -R $read_length -V SNV  $outroot/$chr_id/$read_length.$tool.variants.SNV.vcf $outroot/$chr_id/$read_length.$tool.variants.SNV.ovl/sites.txt 
            python $eval_script_dir/get_sv_prediction_stats.py -T $true_indel -A $tool -R $read_length -V INDEL $outroot/$chr_id/$read_length.$tool.variants.INDEL.vcf $outroot/$chr_id/$read_length.$tool.variants.INDEL.ovl/sites.txt 

            # pred_snp="$(grep -E "\t[ACGT]\t[ACGT]\t" $outroot/$chr_id/$read_length.$tool.variants.SNV.vcf | wc -l)"
            # pred_indel="$(grep "INDEL" $outroot/$chr_id/$read_length.$tool.variants.INDEL.vcf  | wc -l)"

            # ovl_snv="$(cat $outroot/$chr_id/$read_length.$tool.variants.SNV.ovl/sites.txt | wc -l)"
            # ovl_indel="$(cat $outroot/$chr_id/$read_length.$tool.variants.INDEL.ovl/sites.txt | wc -l)"
            # echo $chr_id,$read_length,$tool,align,$pred_snp,$pred_indel,$ovl_snv,$ovl_indel


            
            echo
        done
 
    done
    echo
    echo
done


