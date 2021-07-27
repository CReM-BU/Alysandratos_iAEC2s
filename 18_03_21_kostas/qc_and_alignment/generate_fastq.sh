#! /bin/zsh
# zsh option to split on whitespace (otherwise for-loops don't work)
setopt shwordsplit

# source environmental variables and load tools
source /restricted/projectnb/crem-bioinfo/project_code/18_03_21_kostas/.source_me
module load gcc/4.9.2
module load bcl2fastq/2.18

export input_folder="/restricted/projectnb/crem-bioinfo/project_data/${project_name}/180314_NB500996_0140_AHM75JBGX5_Kotton"
export reference_folder="/restricted/projectnb/crem-bioinfo/reference_data"
export task_name="mkfastq"
export executables_path="/restricted/projectnb/crem-bioinfo/project_code/00_pan_project/executables"

# loop variables
script_name="${task_name}"

# housekeeping variables
group="crem-seq"
scripts_dir=$PW/qsub_scripts
script=$scripts_dir/${script_name}
log_dir=$PW/logs/$task_name/
log=$log_dir/${script_name}.log
err=$log_dir/${script_name}.err
mkdir -p $scripts_dir $log_dir $output_dir
echo > $log > $err # reset, because otherwise it appends

export num_cores="16"
export localmem="500"

# Create the script
echo '#!/usr3/bustaff/nacho/bin/zsh
  module load gcc/4.9.2
  module load bcl2fastq/2.18
  # initialize variables
  source /restricted/projectnb/crem-bioinfo/project_code/18_03_21_kostas/.source_me
  printf "Start: %s\n" "$(date)"
    echo $project_name
    cd /restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations
    rm -rf HM75JBGX5
    $executables_path/cellranger-2.0.1/cellranger mkfastq \
      --run="/restricted/projectnb/crem-bioinfo/project_data/18_03_21_kostas/180314_NB500996_0140_AHM75JBGX5_Kotton" \
      --csv="/restricted/projectnb/crem-bioinfo/project_data/18_03_21_kostas/180314_NB500996_0140_AHM75JBGX5_Kotton/2018-03-14_Kotton.csv" \
      --localcores=1
      --localmem=300
  printf "End: %s\n" "$(date)\n"
' > $script
chmod +x $script

# Submit
qsub -P $group -N $script_name -o $log -e $err -V \
  -pe omp $num_cores \
  -l mem_total="${localmem}G" \
  $script
