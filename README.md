# Code for genomic simulations and Approximate Bayesian Computation for the paper: Substructured population growth in the Ashkenazi Jews inferred with Approximate Bayesian Computation
Created by Ariella Gladstein  
algladstein@gmail.com

## Documentation
Please go to the [wiki](https://bitbucket.org/agladstein/macsswig_simsaj/wiki/Home)!

## Usage

The main script run_sims_AJmodel1_chr1_all.py runs genome simulations and computes statistics on the simulation output, and is run with run_sims_AJmodel1_chr1_all.py  
Version 1


Run as  
`python run_sims_AJmodel1_chr1_all.py jobID inputfile simsize seed param_distribution germline output_dir`
  
jobID = can be any unique value to identify the output  
simsize = full or the length of the locus to be simulated in base pairs  
seed = a seed value, or 0 if no seed is to be used  
param_distribution = prior, min, or max  
* prior = simulations with parameter values given by prior distribution  
* min = simulations with predetermined parameter values that will produce   simulations with shorter run times and less memory - only for testing and profiling purposes  
* max = simulations with predetermined parameter values that will produce simulations with longer run times and more memory - only for testing and profiling purposes  
germline = 0 to run GERMLINE, 1 to not run GERMLINE (will try to read GERMLINE output from file)  

output_dir = path to the directory to output to. No argument will use the default of current dir "."

e.g.:  
``python run_sims_AJmodel1_chr1_all.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed full 0 prior 0 output_dir``

Test simulation:  
``python run_sims_AJmodel1_chr1_all.py 1 ill_650_test.bed 1000000 1278 prior 1 output_dir``

Uses c++ programs macs and GERMLINE. For more information on these programs, see:  
https://github.com/gchen98/macs  
https://github.com/sgusev/GERMLINE  

Intermediate files go to sim_data_AJ_M1 and germline_out_AJ_M1 and are NOT rm in python script.  
Output files go to sim_values_AJ_M1 and results_sims_AJ_M1.

## Citation
The manuscript has been submitted. If you make use of this work in your research, we would appreciate the following citation:

Gladstein, A.L. and Hammer, M.F. (2018). Substructured population growth in the Ashkenazi Jews inferred with Approximate Bayesian Computation. Manuscript submitted for publication.


This work is based on code from Consuelo Quinto-Cortes and Krishna Veeramah, used in the papers:  

Quinto-Cortés, C. D., Woerner, A. E., Watkins, J. C., & Hammer, M. F. (2018). Modeling SNP array ascertainment with Approximate Bayesian Computation for demographic inference. Scientific Reports, 8, 10209. http://doi.org/10.1038/s41598-018-28539-y  

Veeramah, K. R., Wegmann, D., Woerner, A., Mendez, F. L., Watkins, J. C., Destro-Bisol, G., … Hammer, M. F. (2012). An Early Divergence of KhoeSan Ancestors from Those of Other Modern Humans Is Supported by an ABC-Based Analysis of Autosomal Resequencing Data. Molecular Biology and Evolution, 29(2), 617–630. http://doi.org/10.1093/molbev/msr212
