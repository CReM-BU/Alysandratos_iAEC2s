#!/bin/zsh
project_name="18_09_25_kostas"
BASE="/restricted/projectnb/crem-bioinfo"
PW=$BASE/project_workspace/${project_name}
PD=$BASE/project_data/${project_name}/2018-09-24_Kotton-10X
PC=$BASE/project_code/${project_name}
REF=$BASE/reference_data/REFERENCES
task_name="mkfastq_correctIndex"
executables_path="/restricted/projectnb/crem-bioinfo/project_code/00_pan_project/executables"

# loop variables
script_name="${task_name}"

# housekeeping variables
group="crem-seq"
scripts_dir=$PW/qsub_scripts
script=$scripts_dir/${script_name}.qsub
output_dir=$PW/calculations
log_dir=$PW/logs/$task_name/
log=$log_dir/${script_name}.log
err=$log_dir/${script_name}.err
mkdir -p $scripts_dir $log_dir $output_dir
echo > $log > $err # reset, because otherwise it appends

num_cores="16"
localmem="504" # current possible choices for SCC are 94, 125, 252 and 504

# ~55 minutes
# Create the script
echo '#!/bin/zsh
  module load gcc/4.9.2
  module load bcl2fastq/2.18
  date
  echo '$project_name'
  cd '$output_dir'
  time '$executables_path'/cellranger-2.0.1/cellranger mkfastq \
      --run='$PD' \
      --csv='$PD'/2018-09-24_Kotton.csv \
      --localcores='$num_cores' \
      --localmem='$localmem'
  date
' > $script
chmod +x $script

# Submit
qsub -P $group -N $script_name -o $log -e $err -V \
  -pe omp $num_cores \
  -l mem_total="${localmem}G" \
  $script
