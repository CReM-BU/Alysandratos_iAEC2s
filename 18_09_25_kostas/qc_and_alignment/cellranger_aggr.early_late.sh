#!/bin/zsh

project_name="18_09_25_kostas"
BASE="/restricted/projectnb/crem-bioinfo"
PW=$BASE/project_workspace/${project_name}
CELLRANGER_PATH=$BASE/project_code/00_pan_project/executables/cellranger-2.0.1
script_name="cellranger_aggr.early_late"

echo "library_id,molecule_h5
B2_early,/restricted/projectnb/crem-bioinfo/project_workspace/18_09_25_kostas/calculations/cellranger_count_with_tdTomato_index_corrected/B2/outs/molecule_info.h5
C11_early,/restricted/projectnb/crem-bioinfo/project_workspace/18_09_25_kostas/calculations/cellranger_count_with_tdTomato_index_corrected/C11/outs/molecule_info.h5
B2_late,/restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations/cellranger_count_with_tdTomato/B2/outs/molecule_info.h5
C11_late,/restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations/cellranger_count_with_tdTomato/C11/outs/molecule_info.h5
C11_late_b,/restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations/cellranger_count_with_tdTomato_2nd/C11/outs/molecule_info.h5" > \
  /restricted/projectnb/crem-bioinfo/project_data/18_09_25_kostas/aggr_early_late.csv


CELLRANGER_PATH=/restricted/projectnb/crem-bioinfo/project_code/00_pan_project/executables/cellranger-2.0.1
script_name="cellranger_aggr.early_late"
name="cellranger_aggr.early_late"
output_dir=$PW/calculations/$script_name/
# takes 40 minutes to run
num_cores=16
localmem=252
combined_metadata=/restricted/projectnb/crem-bioinfo/project_data/18_09_25_kostas/aggr_early_late.csv

# Set up housekeeping variables
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
    time '$CELLRANGER_PATH'/cellranger aggr \
      --id=combined_data \
      --csv='$combined_metadata' \
      --normalize=mapped \
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
