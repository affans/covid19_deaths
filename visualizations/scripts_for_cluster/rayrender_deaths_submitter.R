#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# inputs <- matrix(nrow=length(args)/2,ncol=1,data=args[seq(2,length(args),2)])
# rownames(inputs) <- args[seq(1,length(args),2)]
# inputs[,1] <- trimws(inputs[,1])


# Create batch job and submit to the cluster

load("output_data/states_geom_and_data.RData")

for(index in 1:length(unique(states_geom_and_data$date))){
  
  file.create(paste("cluster_jobs/rayrender_deaths_",index,".sh",sep=""))
  sink(file = paste("cluster_jobs/rayrender_deaths_",index,".sh",sep=""))
  cat("#!/bin/bash\n")
  cat(paste("#SBATCH -n 1\n#SBATCH -c ","5","\n#SBATCH -J ",index,"_deaths\n#SBATCH --partition=","general","\n#SBATCH -t ","1:00:00","\n#SBATCH --mem-per-cpu=","5G","\n#SBATCH --mail-user=vincent.cannataro@yale.edu\n#SBATCH --mail-type=FAIL\n",sep=""),append=T)
  cat("module load GDAL\nmodule load GEOS\nmodule load PROJ\nmodule load UDUNITS/2.1.24-foss-2016b\n",append = T)
  cat(paste("~/../../pi/townsend/vlc24/R_source/R-4.0.2/bin/Rscript R/death_rayrender_render.R",index,sep=" "),append = T)
  sink()
  system(paste("sbatch cluster_jobs/rayrender_deaths_",index,".sh",sep=""))
  
}