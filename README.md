Created by Ariella Gladstein, based on code from Consuelo Quinto Cortes and Krishna Veeramah.

agladstein@email.arizona.edu

The main script run_sims_AJmodel1_chr1_all.py runs genome simulations and computes statistics on the simulation output, and is run with run_sims_AJmodel1_chr1_all.py
Version 1


Run as
python run_sims_AJmodel1_chr1_all.py jobID inputfile simsize seed param_distribution germline
jobID = can be any unique value to identify the output
simsize = full or the length of the locus to be simulated in base pairs
seed = a seed value, or 0 if no seed is to be used
param_distribution = prior, min, or max
* prior = simulations with parameter values given by prior distribution
* min = simulations with predetermined parameter values that will produce simulations with shorter run times and less memory - only for testing and profiling purposes
* max = simulations with predetermined parameter values that will produce simulations with longer run times and more memory - only for testing and profiling purposes
germline = 0 to run GERMLINE, 1 to not run GERMLINE (will try to read GERMLINE output from file)
e.g.:
``python run_sims_AJmodel1_chr1_all.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed full 0 prior 0``
Test simulation:
``python run_sims_AJmodel1_chr1_all.py 1 ill_650_test.bed 1000000 1278 prior 1``

Uses c++ programs macs and GERMLINE. For more information on these programs, see:
https://github.com/gchen98/macs
https://github.com/sgusev/GERMLINE

Intermediate files go to sim_data_AJ_M1 and germline_out_AJ_M1 and are NOT rm in python script.
Output files go to sim_values_AJ_M1 and results_sims_AJ_M1.


If using Vagrant:

```bash
vagrant up
vagrant ssh
```

```bash
sudo apt-get update
sudo apt-get install python-virtualenv git python-dev
sudo easy_install -U distribute
virtualenv macss_env
source macss_env/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
macss_env/bin/python run_sims_AJmodel1_chr1_all.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed full 0 prior 0
```