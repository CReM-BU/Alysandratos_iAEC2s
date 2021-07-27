#!/bin/zsh
project_name="19_10_29_kostas"
BASE="/restricted/projectnb/crem-bioinfo"
PW=$BASE/project_workspace/${project_name}
PD=$BASE/project_data/${project_name}/2019-10-25_DKotton-10x
PC=$BASE/project_code/${project_name}
task_name="mkfastq"

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
localmem="256" # current possible choices for SCC are 94, 125, 252 and 504

# ~55 minutes
# Create the script
echo '#!/bin/zsh
  module load bcl2fastq/2.20
  module load cellranger/3.0.2
  date
  echo '$project_name'
  cd '$output_dir'
  time cellranger mkfastq \
      --run='$PD' \
      --csv='$PD'/../2019_10_28_DKotton_10x_sample_sheet.csv \
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
