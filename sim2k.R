# run like Rscript sim2k.R /home/au/code/SCRABBLE/juliaEM/sim2k /home/au/code/SCRABBLE/analysis_library.R 2 && Rscript emjulia_cmd.R /home/au/code/SCRABBLE/juliaEM/sim2k /home/au/code/SCRABBLE/analysis_library.R 2

args <- commandArgs(trailingOnly=TRUE)

### arguments
datadir = args[1]
libpath = args[2]
tempcores = args[3]

# ensure that R.home() is the same as ENV["R_HOME"] when RCall.jl was built, otherwise segfaults will occur
Sys.setenv(JULIA_PROJECT = "/home/au/code/DTMwork")
Sys.setenv(LD_LIBRARY_PATH = "/usr/lib/R/lib")

# load the libraries

library(dplyr)
library(vsn)
library(Seurat)
library(scater)
library(edgeR)
library(gridExtra)
library(R.matlab)
library(cowplot)
library(biomaRt)
library(data.table)
library(lattice)
library(scImpute)
library(SCRABBLE)
library(VennDiagram)
library(Rtsne)
library(DT)
library(ggpubr)
library(ggsignif)
library(scatterplot3d)
library(ggplot2)
library(reshape2)
library(ggfortify)
library(refGenome)
library(pheatmap)
library(RColorBrewer)
library(dendsort)
library(entropy)
library(DrImpute)
library(splatter)
library(RColorBrewer)
library(mcriPalettes)
library(plotly)
library(factoextra)
library(cluster)
library(NbClust)
library(fpc)
library(class)
library(VIPER)
library(SC3)
library(doParallel)

source(libpath)
dir.create(datadir,recursive=TRUE)
setwd(datadir)

############################################
# begin cluster setup
############################################
if(exists("cl")){
  stopCluster(cl)
}

ncores = min(tempcores,detectCores(logical=FALSE)/2 - 1)
registerDoParallel(cores=ncores)
cl <- makeCluster(ncores, type="FORK")
############################################
# end cluster setup
############################################



dropout_mid = c(4, 6.5, 8)
ngenes = 20000 # default 800
nbulk = 8 # default 10
ldacores = min(15,nbulk)
ncellsperbulk = 2000 # default 100, if this is too low, svd has convergence issues
ndrop = length(dropout_mid)
nseed = 100
seed_drop = expand.grid(1:ndrop,1:nseed)
colnames(seed_drop) = c("dropout_index","seed_value")

# scrabble parameters
alpha_p <- 1 # default 1
beta_p <- 1e-06 # default 1e-06
gamma_p <- 1e-04 # default 1e-04

# the following script is to generate the simulation data. Here we use
# HPC to generate the simulation which could reduce the running time
foreach(i = 1:dim(seed_drop)[1], .combine = 'c',
        .export = c("generate_simulation_splatter","generate_save_data"),
        .packages = c("splatter", "doParallel","edgeR")) %dopar% {
          dropout_index = seed_drop$dropout_index[i]
          seed_value = seed_drop$seed_value[i]
          generate_save_data(dropout_index, seed_value,nGenes = ngenes, nbulk = nbulk, ncellsperbulk = ncellsperbulk, dropout_mid=dropout_mid)
        }

save(dropout_mid,ngenes,nbulk,ldacores,ncellsperbulk,ndrop,nseed,seed_drop, file = "dataparams.RData")
