# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* The main script run_sims_AJmodel1_chr1.py runs genome simulations and computes statistics on the simulation output.
* Version 1

### How do I get set up? ###

* Must create directories:
* sim_values_AJ_M1
* results_sims_AJ_M1
* sim_data_AJ_M1
* germline_out_AJ_M1

* Intermediate files go to sim_data_AJ_M1 and germline_out_AJ_M1 and are rm in python script.
* Output files go to sim_values_AJ_M1 and results_sims_AJ_M1.


### Contribution guidelines ###

* Run as
* python run_sims_AJmodel1_chr1.py jobID inputfile
* e.g.:
* python run_sims_AJmodel1_chr1.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed

* Uses c++ programs macs and germline. For more information on these programs, see:
* https://github.com/gchen98/macs
* https://github.com/sgusev/GERMLINE

### Who do I talk to? ###

* Created by Ariella Gladstein, based on code from Consuelo Quinto Cortes and Krishna Veeramah.
* agladstein@email.arizona.edu

* Other community or team contact