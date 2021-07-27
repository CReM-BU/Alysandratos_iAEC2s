#! /bin/zsh
# zsh option to split on whitespace (otherwise for-loops don't work)
setopt shwordsplit

# source environmental variables and load tools
source /restricted/projectnb/crem-bioinfo/project_code/18_03_21_kostas/.source_me

# echo "library_id,molecule_h5
# B2,/restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations/cellranger_count_with_tdTomato/B2/outs/molecule_info.h5
# C11,/restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations/cellranger_count_with_tdTomato/C11/outs/molecule_info.h5" > \
#   /restricted/projectnb/crem-bioinfo/project_data/18_03_21_kostas/180314_NB500996_0140_AHM75JBGX5_Kotton/combined_metadata.csv


export CELLRANGER_PATH=/restricted/projectnb/crem-bioinfo/project_code/00_pan_project/executables/cellranger-2.0.1
export script_name="cellranger_aggr"
export name="aggregate_cellranger"
export output_dir=$PW/calculations/$script_name/
# this takes 40 minutes to run
export num_cores=16
export localmem=96
export combined_metadata=$PD/180314_NB500996_0140_AHM75JBGX5_Kotton/combined_metadata.csv

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
echo '#!/usr3/bustaff/nacho/bin/zsh
    source /restricted/projectnb/crem-bioinfo/project_code/18_03_21_kostas/.source_me
    printf "Start: %s\n" "$(date)"
    
    echo "job id - job name:" $JOB_ID " - " $JOB_NAME
    cd $output_dir
    # rm -rf combined_data # otherwise it will fail
    $CELLRANGER_PATH/cellranger aggr \
      --id=combined_data \
      --csv=$combined_metadata \
      --normalize=mapped \
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
