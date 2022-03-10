# TablesAsDesigns
R Code associated with the paper "Let's practice what we preach: Planning and interpreting simulation studies with design and analysis of experiments" by Hugh Chipman and Derek Bingham.  

The code implements the 2 examples from the paper: 
  * A reanalysis of a simulation study from Krishnamoorthy,  Mallick and  Mathew (2011), "Inference for the Lognormal Mean and Quantiles Based on Samples With Left and Right Type I Censoring" *Technometrics*, 53, 72--83. (KMM.R, KMM-logistic.R)
  * An experiment developed for the paper, comparing the predictive performance of 2 statistical learning algorithms. (ff-analysis.R analyzes the runs contained in the 2 .RData files; run_cluster.R will re-run the experiment and create the .RData files; "replicate_12_combined_results.csv" contains the same results as the .RData files in CSV format.)