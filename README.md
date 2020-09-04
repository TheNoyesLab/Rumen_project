# NCBA Rumen project data



# Important components of this repository:
* rumen_staging_script.R
  * this is the main script to re-create the analysis in this repository.
* the "graphs/" directory
  * contains exploratory for the microbiome results and the resistome data (MEGARes, ACLAME, ICEBerg)
* the "network_counts/" directory
  * use the count table files in this directory to create the bayesian networks. 
  * If you can run the commands in the rumen_staging_script.R script up to line #206, you can open up the file, "scripts/1_network_counts_by_project.R" to see how these counts were made.
  * you should be able to use the file, "scripts/2_network_run_bnlearn.R", as an example for how to create the bayesian networks from the count tables available.
