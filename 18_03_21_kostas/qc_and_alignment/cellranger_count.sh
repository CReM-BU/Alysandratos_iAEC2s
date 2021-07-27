#! /bin/zsh
# zsh option to split on whitespace (otherwise for-loops don't work)
setopt shwordsplit

# source environmental variables and load tools
source /restricted/projectnb/crem-bioinfo/project_code/18_03_21_kostas/.source_me

export CELLRANGER_PATH=/restricted/projectnb/crem-bioinfo/project_code/00_pan_project/executables/cellranger-2.0.1
export script_name="cellranger_count_with_tdTomato"
# export num_cores=32
export num_cores=16
export localmem=250 # 32*6 = 192G is the 10x recommendation (6G per core)
export GENOME_DIR=/restricted/projectnb/crem-bioinfo/reference_data/hg19_tdtomato_genome_cellranger
export output_dir=$PW/calculations/$script_name
export input_dir=/restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations/HM75JBGX5/outs/fastq_path/HM75JBGX5
## Run cellranger

# TODO: create hg19 with NKX2-1

# This takes 4 hours (with 24 cores, I couldn't get 32)

for output_prefix in B2 C11 ; do
  # Export variables used in the script
  export output_name=$output_prefix # It doesn't work without this

  # Set up housekeeping variables
  export name=${output_prefix}_${script_name}
  group="crem-seq"
  scripts_dir=$PW/qsub_scripts/$script_name
  script=$scripts_dir/${name}.qsub
  log_dir=$PW/logs/$script_name
  log=$log_dir/${name}.log
  err=$log_dir/${name}.err
  mkdir -p $scripts_dir $log_dir $output_dir
  echo > $log > $err # reset, because otherwise it appends
  
  # Create the script
  echo '#!/usr3/bustaff/nacho/bin/zsh
      source /restricted/projectnb/crem-bioinfo/project_code/18_03_21_kostas/.source_me
      printf "Start: %s\n" "$(date)"
      
      echo "job id - job name:" $JOB_ID " - " $JOB_NAME
      mkdir -p $output_dir
      cd $output_dir
      $CELLRANGER_PATH/cellranger count \
        --id=$output_name \
        --transcriptome=$GENOME_DIR \
        --fastqs=$input_dir/DK-${output_name}-Live-2000 \
        --localcores=$num_cores \
        --localmem=$localmem

      printf "End: %s\n\n" "$(date)"
  ' > $script
  chmod +x $script

  # Submit
  qsub -P $group -N $name -o $log -e $err -V \
    -pe omp $num_cores \
    -l mem_total="${localmem}G" \
    $script
done

# After this is done, we have to run cellranger_aggr