#!/bin/zsh
# only last two samples
# aggregate samples of tracheal mesenchyme (previous run, 18_08_27), tracheal epithelium and lung epithelium (run 18_09_27) in order to compare the macrophages from each of them
project_name="19_10_29_kostas"
BASE="/restricted/projectnb/crem-bioinfo"
PW=$BASE/project_workspace/${project_name}

echo "library_id,molecule_h5,batch
B2_nov,/restricted/projectnb/crem-bioinfo/project_workspace/19_10_29_kostas/calculations/cellranger_count_old_genome/KA_B2/outs/molecule_info.h5,v3_lib
C11_nov,/restricted/projectnb/crem-bioinfo/project_workspace/19_10_29_kostas/calculations/cellranger_count_old_genome/KA_C11/outs/molecule_info.h5,v3_lib" > \
  /restricted/projectnb/crem-bioinfo/project_data/19_10_29_kostas/aggr_twosamples.csv


script_name="cellranger_aggr_twosamples"
name="cellranger_aggr_twosamples"
output_dir=$PW/calculations/$script_name/
# takes 40 minutes to run
num_cores=16
localmem=252
combined_metadata=/restricted/projectnb/crem-bioinfo/project_data/19_10_29_kostas/aggr_twosamples.csv

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
  module load bcl2fastq/2.20
  module load cellranger/3.0.2
    date    
    cd '$output_dir'
    time cellranger aggr \
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
