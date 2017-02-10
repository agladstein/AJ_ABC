* Created by Ariella Gladstein, based on code from Consuelo Quinto Cortes and Krishna Veeramah.
* agladstein@email.arizona.edu

* The main script run_sims_AJmodel1_chr1.py runs genome simulations and computes statistics on the simulation output.
* Version 1


* Run as
* python run_sims_AJmodel1_chr1.py jobID inputfile simsize
* e.g.:
``python run_sims_AJmodel1_chr1.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed full``
* Test simulation:
``python run_sims_AJmodel1_chr1.py 1 ill_650_test.bed 1000000``

* Uses c++ programs macs and germline. For more information on these programs, see:
* https://github.com/gchen98/macs
* https://github.com/sgusev/GERMLINE

* Intermediate files go to sim_data_AJ_M1 and germline_out_AJ_M1 and are rm in python script.
* Output files go to sim_values_AJ_M1 and results_sims_AJ_M1.


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
python run_sims_AJmodel1_chr1.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed full
```