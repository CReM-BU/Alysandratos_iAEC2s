module load R/3.6.0
module load hdf5
module load pandoc/2.5 # knitr
R

mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/19_10_29_kostas/calculations/analysis
for pref in KA_B2  KA_C11; do
Rscript -e "require('rmarkdown'); render( \
input='/restricted/projectnb/crem-bioinfo/project_code/19_10_29_kostas/ka.10x.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_10_29_kostas/calculations/analysis/$pref.html',\
params= list( prefix= '$pref', resDEG= 'RNA_snn_res.0.25', percent.mito= 0.25 ))"  &
sleep 20
done 

R
params <- list()
params$prefix <- "KA_B2"
params$resDEG <- "RNA_snn_res.0.25"
params$percent.mito <- 0.25



# combine sctransform

Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_10_29_kostas/10x.merge.sctransform.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_10_29_kostas/calculations/analysis/merged.sctransform.html',\
params= list( prefix= 'merged.sctransform', resDEG= 'SCT_snn_res.0.25', percent.mito= 0.25, regress.batch='FALSE', sc.transform ='TRUE'))" &

# aggr:  19_10_29_kostas 19_07_09_kostas 18_09_25_kostas 18_03_21_kostas
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_10_29_kostas/ka.10x.aggr.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_10_29_kostas/calculations/analysis/aggr.C11_B2.all.html',\
params= list( prefix= 'combined_data', prefix.out = 'aggr.C11_B2.all', resDEG= 'RNA_snn_res.0.25', percent.mito= 0.25, regress.batch='FALSE' ))"

# aggr only two samples: 19_10_29_kostas B2 and C11 
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/19_10_29_kostas/ka.10x.aggr_twosamples.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/19_10_29_kostas/calculations/analysis/aggr.C11_B2.twosamples.html',\
params= list( prefix= 'combined_data', prefix.out = 'aggr.C11_B2.twosamples', resDEG= 'RNA_snn_res.0.25', percent.mito= 25, regress.batch='FALSE' ))"


params <- list()
params$prefix <- "merged.sctransform"
# params$prefix <- "combined_data"
# params$prefix.out  <- "aggr.C11_B2.all"
# params$prefix.out  <- "aggr.C11_B2.twosamples"
params$resDEG <- "SCT_snn_res.0.25"
params$percent.mito <- 25
params$regress.batch <- FALSE
params$sc.transform <- TRUE
