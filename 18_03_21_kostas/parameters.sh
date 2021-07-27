

module load R/3.6.0
module load pandoc/2.5 # knitr
module load hdf5
# individual analysis of each sample
R

Main analyses for paper: fgsea.Rmd and 10x.combined.Rmd (for report rendering)

pref="analysis.seurat.report"
mkdir -p /restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations/analysis/$pref
Rscript -e "require('rmarkdown'); render(input='/restricted/projectnb/crem-bioinfo/project_code/18_03_21_kostas/10x.combined.Rmd',\
output_file='/restricted/projectnb/crem-bioinfo/project_workspace/18_03_21_kostas/calculations/analysis/$pref/$pref.html',\
params= list( prefix= '$pref', resDEG= 'SCT_snn_res.0.2', percent.mito= 15, cell.source= 'human', harmonize = 'FALSE', sample.names = NULL ))" &

params <- list()
params$prefix <- "analysis.seurat.report"
# params$resDEG <- "Harmony_res.0.2"
params$resDEG <- "SCT_snn_res.0.2"
params$percent.mito <- 15 # range 0 to 100
params$cell.source <- 'human'
params$harmonize <- FALSE
params$sample.names <- NULL


