#! /bin/zsh
# zsh option to split on whitespace (otherwise for-loops don't work)
setopt shwordsplit

# cat /restricted/projectnb/crem-bioinfo/reference_data/aux/genome_without_reporters.fa \
#   /restricted/projectnb/crem-bioinfo/reference_data/tdtomato.fa > \
#   /restricted/projectnb/crem-bioinfo/reference_data/hg19_tdtomato_genome/genome.fa

# source environmental variables and load tools
source /restricted/projectnb/crem-bioinfo/project_code/18_03_21_kostas/.source_me

# Export variables used in the script
export CELLRANGER_PATH="/restricted/projectnb/crem-bioinfo/project_code/00_pan_project/executables/cellranger-2.0.1/"
export script_name="cellranger_mkref"
export num_cores=16

export REF=/restricted/projectnb/crem-bioinfo/reference_data

# Set up housekeeping variables
group="crem-bioinfo"
scripts_dir=$PW/qsub_scripts
script=$scripts_dir/${script_name}
log_dir=$PW/logs/${script_name}
log=$log_dir/${script_name}.log
err=$log_dir/${script_name}.err
mkdir -p $scripts_dir $log_dir $genome_dir
echo > $log > $err # reset, because otherwise it appends

# cat $REF/refdata-GRCh38-2.1.0/fasta/genome.fa \
#   $REF/reporters.fasta > $REF/10x_GRCh38_reporters.fasta

# Create the script
echo '#!/bin/zsh
    printf "Start: %s\n" "$(date)"

      cd $REF
      $CELLRANGER_PATH/cellranger mkref \
        --genome=hg19_tdtomato_genome_cellranger \
        --fasta=$REF/hg19_tdtomato_genome/genome.fa \
        --genes=$REF/hg19_tdtomato_genome/genes.gtf
        
    printf "End: %s\n" "$(date)\n"
' > $script
chmod +x $script

# Submit
qsub -P $group -N $script_name -o $log -e $err -V \
  -pe omp $num_cores \
  -l mem_total=32G \
  $script