#!/bin/zsh

project_name="18_09_25_kostas"
BASE="/restricted/projectnb/crem-bioinfo"
PW=$BASE/project_workspace/${project_name}
CELLRANGER_PATH=$BASE/project_code/00_pan_project/executables/cellranger-2.0.1
script_name="cellranger_count_with_tdTomato_index_corrected" # previous run had wrong index attribution file

num_cores=16
localmem=504 # 16 cores * 6G = 96G is the minimum RAM recommendation (6G per core). Min cores: 8
GENOME_DIR=$BASE/reference_data/hg19_tdtomato_genome_cellranger
output_dir=$PW/calculations/$script_name
input_dir=$PW/calculations/HMNMWBGX7/outs/fastq_path/HMNMWBGX7
## Run cellranger

for sample_dir in ${input_dir}/* ; do
  # Set up housekeeping variables
  output_prefix=`basename ${sample_dir}`  # get B2 from /restricted/projectnb/crem-bioinfo/project_workspace/18_07_11_kostas/calculations/HKJGJBGX7/outs/fastq_path/HKJGJBGX7/B2 
  name=${output_prefix}_${script_name} 
  group="crem-seq"
  scripts_dir=$PW/qsub_scripts/$script_name
  script=$scripts_dir/${name}.qsub
  log_dir=$PW/logs/$script_name
  log=$log_dir/${name}.log
  err=$log_dir/${name}.err
  mkdir -p $scripts_dir $log_dir $output_dir
  echo > $log > $err # reset, because otherwise it appends
  
  # Create the script
  echo '#!/bin/zsh
      date
      echo "job id - job name:" '$JOB_ID' " - " '$JOB_NAME'
      cd '$output_dir'
      time '$CELLRANGER_PATH'/cellranger count \
        --id='$output_prefix' \
        --transcriptome='$GENOME_DIR' \
        --fastqs='${sample_dir}' \
        --localcores='$num_cores' \
        --localmem='$localmem'
      date
  ' > $script
  chmod +x $script

  # Submit
  qsub -P $group -N $name -o $log -e $err -V \
    -pe omp $num_cores \
    -l mem_total="${localmem}G" \
    $script
done

# After this is done, we have to run cellranger_aggr