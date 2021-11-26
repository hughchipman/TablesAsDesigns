# This code carries out the simulation of the "statistical learning experiment"
# in the paper "Let's practice what we preach: Planning and interpreting simulation studies with design and analysis of experiments", 
# written by Hugh Chipman and Derek Bingham.

# First, load the function that will execute 1 replicate of the experiment.  
# This corresponds to 2^7 = 128 runs.
source("run_replicate.R")

# Use the parallel library, running 1 replicate on each core
library(parallel)
detectCores() #give the number of cores on your machine.
numrep <- 2 # The study actually uses 2 replicates

cl <- makeCluster(numrep, outfile="parallel_out.txt")
system.time(parLapply(cl, 1:numrep, run_replicate))

# For 5 replicates (not what is given above), we get on 2 occasions:
#    user  system elapsed <--- for 5 x 128 runs
#    0.68    0.40 5003.46 <--- about 83 minutes, or ~17 minutes per replicate

# the next day a second run of essentially the same code gave
#    user  system elapsed 
#    0.10    0.01 4369.47 


stopCluster(cl)


#system.time(run_replicate(1))
