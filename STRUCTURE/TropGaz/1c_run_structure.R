#install.packages("ParallelStructure", repos="http://R-Forge.R-project.org")
library(ParallelStructure)
library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run ParallelStructure                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Usage: nohup Rscript 2_run_structure.R &

options(scipen=999)

# Hier

gen <- fread("structure_trop_red.stru")

#~~ Specify in and out files for structure

infile <- "structure_trop_red.stru"
outpath <- "results/"

#~~ construct job matrix and write to job file

nrep <- 10
burnin <- 500000
niter <- 1000000
up_to_k <- 3


#~~ define variables for job matrix

k_var <- rep(1:up_to_k, each = nrep)
ID_var <- as.character(sapply(c(1:up_to_k), function(k) sapply(c(1:nrep), function(x) paste0("T",k, "_", x))))

#~~ make the job matrix
pop <- "1,2,3,4,5,6,7,8,9" #number of pops in the file

seal_jobs <- matrix(c(ID_var, rep(pop, nrep * up_to_k), k_var, rep(burnin, nrep * up_to_k),
                      rep(niter, nrep * up_to_k)), nrow = nrep * up_to_k)

write(t(seal_jobs), ncol = length(seal_jobs[1,]), file = "seal_jobs.txt")

#~~ file path to structure

STR_path='/usr/local/bin/'


#~~ Run Parallel Structure

system("mkdir results")

#~~ Run structure (from terminal, do not run this last part in Rstudio)
# change n_cpu depending on how busy server is (was 20)


ParallelStructure::parallel_structure(structure_path=STR_path, joblist='seal_jobs.txt', n_cpu=10, 
              infile='structure_trop_red.stru', 
              outpath='results/', 
              numinds = nrow(gen)/2, 
              numloci=ncol(gen)-2,
              noadmix = 0, # no admixture 1=Y (ie no admixture), 0=N (ie assumes admixture, default)
              inferalpha = 1, # 1 = Yes (default); 0 = No. (This option is ignored under the NOADMIX model)
              alpha = 1.0, #default = 1.0
              freqscorr=1, #freqscorr = 1 for correlated allele frequencies model, default = 1
              lambda = 1, #default = 1
              #popdata = 1, (default = 1)
              #popflag = 1, (default = 0)
              #usepopinfo = 1, (default = 0) 
              printqhat=1, 
              plot_output=0, 
              onerowperind=0)
