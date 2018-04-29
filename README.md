Created by Ariella Gladstein, based on code from Consuelo Quinto Cortes and Krishna Veeramah.
agladstein@email.arizona.edu

The main script run_sims_AJmodel2_chr1_all.py runs genome simulations and computes statistics on the simulation output, and is run with run_sims_AJmodel2_chr1_all.py  
Version 1


Run as  
`python run_sims_AJmodel2_chr1_all.py jobID arrayfile chr simsize seed param_distribution germline output_dir`
  
jobID = can be any unique value to identify the output   
arrayfile = array template  
chr = chromosome number  
simsize = full or the length of the locus to be simulated in base pairs    
seed = a seed value, or 0 if no seed is to be used  
param_distribution = ABCtoolbox output file \*MarginalPosteriorDensities\*   
germline = 0 to run GERMLINE, 1 to not run GERMLINE (will try to read GERMLINE output from file)  
output_dir = path to the directory to output to. No argument will use the default of current dir "."

e.g.:  
``python run_sims_AJmodel2_chr1_all.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed 2 full 0 ABC_M2_estimate_1446125_10pls_1000ret_model0_MarginalPosteriorDensities_Obs0.txt 0 output_dir``

Test simulation:  
``python run_sims_AJmodel2_chr1_all.py 1 ill_650_test.bed 2 1000000 0 ABC_M2_estimate_1446125_10pls_1000ret_model0_MarginalPosteriorDensities_Obs0.txt 1 output_dir``

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
macss_env/bin/python run_sims_AJmodel2_chr1_all.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed full 0 prior 0 output_dir
```

Alternatively,
```bash
sudo apt-get update
sudo apt-get install python-virtualenv git python-dev
sudo easy_install -U distribute
virtualenv macss_env
source macss_env/bin/activate
pip install pip-tools
pip-compile
pip-sync
macss_env/bin/python run_sims_AJmodel2_chr1_all.py 1 ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed full 0 prior 0 output_dir
```

Having trouble getting Vagrant started? Make sure you are in the correct directory.
* start up vagrant in working directory
* once in vagrant, go to `/vagrant` to find original working directory

-------------------------

## Getting the posteriors

We used the posterior densities from chr1 with 10 PLS components, and reduced parameters to estimate. 
See https://github.com/agladstein/ABC_results_AJ/blob/master/chr1_AJmodels/RunABC_pls.ipynb

These results are in the Data Store in `/iplant/home/agladstein/ABC_AJmodels_instant_newreal`

```bash
cd Posteriors/chr1
iget -K /iplant/home/agladstein/ABC_AJmodels_instant_newreal/ReduceParams/PLS/ABC_M2_genome_estimate_1446124_10pls_1000ret_model0_MarginalPosteriorDensities_Obs0.txt
```

-------------------------

## Running on University of Arizona HPC
There are four University of Arizona HPC systems - Ocelote, HTC, SMP, and Cluster. All four systems shared the same storage space. Log onto any any of the HPC systems with
``ssh name@hpc.arizona.edu``
Then enter ``ocelote`` for Ocelote or ``ice`` for HTC, SMP, or Cluster.

### Setting up virtualenv on HPC

#### Ocelote
```
module avail
module show python/2/2.7.14
module load python/2/2.7.14
virtualenv -p /cm/shared/uaapps/python/2.7.14/bin/python ~/env/macss_env_ocelote_2.7.14
source ~/env/macss_env_ocelote_2.7.14/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

If requirements.txt on Ocelote based on architecture and packages available, this may not work on other environments. You may have to delete requirements.txt and recompile it using `pip-compile`.

#### ICE
```
module avail
module show python/2.7.9
module load python/2.7.9
virtualenv -p /uaopt/python/2.7.9/bin/python macss_env_ICE_2.7.9
source macss_env_ICE_2.7.9/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

### Submit Update job to Ocelote
```bash
HPC_workflow/run_sims_update.sh HPC_workflow/PBS/run_sims_update_chr2.pbs
```

### Basic Commands on HPC
`qstat` shows all of the jobs currently in the que or running.  
`qstat -t` shows the status of all subjobs.  
`qstat -f` shows details of job.  
`qsub` submits a pbs script.  
`qdel` stops a job.  
`va` shows status of hours.  

-------------------------

## Running as a Workflow on Open Science Grid and Wisconsin CHTC

The workflow creates a Python virtual environment according to the requirements.in file, and then tars up this whole
directory. That tarball is then shipped with the jobs, with the result being that wherever the jobs start up, they 
still have access to any files currently in this directory. However, note that paths will be different, so only
reference files using relative paths.


### Submitting the Workflow

Figure out what the start job id and end job id you want to have in the run. These are provided by the user as a
convenience so that outputs from multiple runs uses unique names. To submit a new workflow, run:

    ./submit [modell] [startid] [endid]

For exaple, to run a small test workflow:

    ./submit run_sims_AJmodel2_chr1_all.py 1 100

The output of the command is something similar to:

    An outputs directory will be created within the base of the workflow directory.
    Directory: /local-scratch/rynge/workflows/macsswig_simsaj_1489465382/outputs

    2017.03.13 23:23:39.930 CDT:    
    2017.03.13 23:23:39.936 CDT:   ----------------------------------------------------------------------- 
    2017.03.13 23:23:39.941 CDT:   File for submitting this DAG to Condor           : macsswig_simsaj-0.dag.condor.sub 
    2017.03.13 23:23:39.947 CDT:   Log of DAGMan debugging messages                 : macsswig_simsaj-0.dag.dagman.out 
    2017.03.13 23:23:39.952 CDT:   Log of Condor library output                     : macsswig_simsaj-0.dag.lib.out 
    2017.03.13 23:23:39.958 CDT:   Log of Condor library error messages             : macsswig_simsaj-0.dag.lib.err 
    2017.03.13 23:23:39.963 CDT:   Log of the life of condor_dagman itself          : macsswig_simsaj-0.dag.dagman.log 
    2017.03.13 23:23:39.968 CDT:    
    2017.03.13 23:23:39.984 CDT:   ----------------------------------------------------------------------- 
    2017.03.13 23:23:41.049 CDT:   Your database is compatible with Pegasus version: 4.7.3 
    2017.03.13 23:23:41.178 CDT:   Submitting to condor macsswig_simsaj-0.dag.condor.sub 
    2017.03.13 23:23:41.360 CDT:   Submitting job(s). 
    2017.03.13 23:23:41.365 CDT:   1 job(s) submitted to cluster 3891204. 
    2017.03.13 23:23:41.371 CDT:    
    2017.03.13 23:23:41.376 CDT:   Your workflow has been started and is running in the base directory: 
    2017.03.13 23:23:41.381 CDT:    
    2017.03.13 23:23:41.387 CDT:     /local-scratch/rynge/workflows/macsswig_simsaj_1489465382/workflow/macsswig_simsaj_1489465382 
    2017.03.13 23:23:41.392 CDT:    
    2017.03.13 23:23:41.397 CDT:   *** To monitor the workflow you can run *** 
    2017.03.13 23:23:41.403 CDT:    
    2017.03.13 23:23:41.408 CDT:     pegasus-status -l /local-scratch/rynge/workflows/macsswig_simsaj_1489465382/workflow/macsswig_simsaj_1489465382 
    2017.03.13 23:23:41.413 CDT:    
    2017.03.13 23:23:41.419 CDT:   *** To remove your workflow run *** 
    2017.03.13 23:23:41.424 CDT:    
    2017.03.13 23:23:41.429 CDT:     pegasus-remove /local-scratch/rynge/workflows/macsswig_simsaj_1489465382/workflow/macsswig_simsaj_1489465382 
    2017.03.13 23:23:41.435 CDT:    
    2017.03.13 23:23:41.624 CDT:   Time taken to execute is 3.636 seconds 

Note how Pegasus uses a directory as "handle" to the workflow. This directory path can be used with various Pegasus commands.

### Starting up Jetstream Instances
When the fourth argument of `./submit` is `TRUE`, OSG will automatically fill active Jetstream instances.

Jetstream instances can be controlled via openstack.

1. find openstack api credentials
2. log into horizon
3. download openrc file
4. create virtualenv, install python openstack tools
5. source openrc file
6. run openstack command


user: aglad  
secret key: 8b2a9758bef242e3f8837725837eeac532ecb8e50463d2482f3ce584a86920f2  
horizon url: https://tacc.jetstream-cloud.org  

Set up openstack
```bash
pyenv install 2.7.14
pyenv shell 2.7.14 #need to do this every time log in
pip install python-openstackclient
source aglad-openrc3.sh #need to do this every time log in
```

show all instances on TACC
```
openstack server list
```

Start instances
```
openstack server create --image 6fc85a08-aad6-40f8-bd5c-ff953c731a09 --flavor 4 --min 2 --max 4 --wait --network aglad-net aj_update_worker
```

Delete instances
```
openstack server list --format value --column Name | grep aj_update_worker | xargs -L 1 openstack server delete
```

### Monitoring on OSG

The system will send email notifications when the workflow changes state, but if you want to see the current state, use the
`pegasus-status` command. For example:

    $ pegasus-status -l /local-scratch/rynge/workflows/macsswig_simsaj_1489465382/workflow/macsswig_simsaj_1489465382
    STAT  IN_STATE  JOB                                                                                                                
    Run      03:12  macsswig_simsaj-0 ( /local-scratch/rynge/workflows/macsswig_simsaj_1489465382/workflow/macsswig_simsaj_1489465382 )
    Idle     02:31   ┣━run-sim.sh_ID0000090                                                                                            
    ...
    Idle     00:40   ┗━run-sim.sh_ID0000077                                                                                            
    Summary: 101 Condor jobs total (I:74 R:27)
    
    UNRDY READY   PRE  IN_Q  POST  DONE  FAIL %DONE STATE   DAGNAME                                 
        6     0     0   100     0     3     0   2.8 Running *macsswig_simsaj-0.dag  

The last couple of lines will tell you the overall state. Note that jobs can be "READY" but not yet submitted to the queue. The
reason for this is that the workflow is configured to keep at most 1,000 idle jobs in the queue, in order to not overwhelm the
scheduler.

To see status of all jobs  
`~/bin/workflow-history`  

To see list of currently running jobs  
`cd ; pegasus-status | grep local-scratch`  

To get the run time of recently completed jobs  
`condor_history -w agladstein`

To get the run time of currently running jobs  
`condor_q agladstein`

To see if workflows are configured to run with Comet/Jetstream/Bridges resources  
`condor_q -const 'Owner == "agladstein" && JobUniverse == 5' -af Iwd -af WantsComet | sort | uniq`

To see how many jobs are currently running on Comet/Jetstream/Bridges  
`condor_q -const 'Owner == "agladstein" && JobUniverse == 5' -af Iwd -af RemoteHost | egrep -i 'jetstream|bridges|comet' | awk '{print $1;}' | sort | uniq -c`

If workflow fails, use the following commands to rescue results that haven't been merged
```bash
find . -type f -path \*.txt | head -1 | xargs head -1 >final.out
find . -type f -path \*.txt | xargs -L 1 tail -n +2 >>final.out
mv final.out final_results.txt
```

### Monitoring on CHTC

Get status of all jobs

```bash
cd /home/nu_agladstein/macsswig_simsaj/workflow/runs
ls | xargs -I % pegasus-status --noqueue /home/nu_agladstein/macsswig_simsaj/workflow/runs/%/workflow/%
```

Find dagid of jobs

```bash
pegasus-status -v | grep macsswig_simsaj-0
```

### Statistics / Debugging

For successful workflows, you can generate statistics such as cumulative runtimes using the `pegasus-statistics -s all [dir]`
command. For failed workflows, `pegasus-analyzer [dir]` can help pinpoint the failures.

look for jobs which have been retried multiple times (in the workflow directory):
```bash
find . -name \*.002
```

### Stopping / Restarting

If you want to stop a workflow, use the `pegasus-remove [dir]` command. Workflows which have stopped for some reason (removed
by the user or maybe from a longer system outage), can be restarted from where they left of with the `pegasus-run [dir]`
command.


### Links to documentation

 * [Pegasus Command Line Interface](https://pegasus.isi.edu/documentation/cli.php)
 * [Pegasus User Guide](https://pegasus.isi.edu/documentation/)
 * [OSG Connect Documentation](https://support.opensciencegrid.org/solution/categories)


-------------------------

## Post Processing

### Rsync OSG output to Atmosphere
After a workflow has completed, use the script `OSG_rsync_atmo.sh` to transfer OSG output to Atmosphere.
This should be run on OSG. It takes 1 argument - the workflow version (e.g. `instant`, `mfloat`, `rscale`).  
e.g.  
`~/macsswig_simsaj/OSG_rsync_atmo.sh instant`

### Combining OSG output
The Pegasus workflow outputs concatenated results_sims and sim_values for all the simulations in the workflow. 
The number of lines in the final output equals the number of simulations plus one for the header.  
To combine results_sims and sim_values across multiple workflows to create the input for ABCtoolbox use the shell script `combine_OSG_final.sh`  
e.g.  
`/vol_c/src/macsswig_simsaj/combine_OSG_final.sh /vol_c/results_macsSwig_AJmodels_instant/OSG`

### Fixing incorrect Headers
The following scripts use the Python package `multiprocessing` and should be run with all the cores of a node.  

First use `find_broken_headers.py` to find any results output files with incorrect headers. This will specifically look for files that have a duplicate of IBD_var_EE in the header. This will create a list of the files with incorrect headers.  
`find_broken_headers.py dir model >>files_to_fix.txt`
Then use `correct_header.py` to fix the results files with incorrect headers. This will create new files with the correct headers in `dir/results_AJ_M${model}_fixed`. Once these files are double checked, they should be moved to the original directory and overwrite incorrect files.     
`correct_header.py files_to_fix.txt`

### Tarring, backing up, transfering, and removing output files

1. tar output directories on HPC
2. upload tar files to google drive
3. transfer tar files to Atmosphere
4. remove files that have been tarred and transferred from HPC

Use crontab to automatically run `tar_rsync_rm.sh` script every hour.
```
MAILTO="agladstein@email.arizona.edu"
0 * * * * /rsgrps/mfh4/Ariella/macsSwig_AJmodels/tar_rsync_rm.sh 1
0 * * * * /rsgrps/mfh4/Ariella/macsSwig_AJmodels/tar_rsync_rm.sh 2
0 * * * * /rsgrps/mfh4/Ariella/macsSwig_AJmodels/tar_rsync_rm.sh 3
```
#### Google drive
We are using the program `drive` from  
https://github.com/odeke-em/drive

There is an executable version in my bin on HPC, Atmosphere, and OSG.
Initialize drive by mounting your Google Drive directory on your local file system.  Must be perfomed the directory you are backing up (does not apply to subdirectories).  
```
drive init
```

To back up all of OSG local-scratch     
```
tar cf - /local-scratch/agladstein | drive push -exclude-ops "delete,update" -no-prompt -piped backup_OSG_local-scratch/backup_OSG_local-scratch_$(date +"%m-%d-%Y-"%T"").tar.gz
```

Back up complete Atmosphere volume
```bash
tar cf - /vol_c | drive push -exclude-ops "delete,update" -no-prompt -piped backup_atmo_vol_1T/backup_atmo_vol_1T_$(date +%m%d%Y%T).tar
```

### Combining HPC output files
This should be run on Atmosphere, but can also be run with a pbs script on HPC, as the I/O operations are slow on HPC.

To combine the simulation output files, sim_values.txt and ms_output.summary, to make the input file for ABCtoolbox use the script
`post_process.py`  
This is will combine the files for one PBS_ID bucket (each bucket should contain about 2000 simulations).  
It uses multiprocessing, and will automatically detect the number of cores available to use.  
The arguments are:  
* path = the path to the sim_values_AJ_M and results_sims_AJ_M directories. On my atmosphere volume this is currently  
`/vol_c/results_macsSwig_AJmodels_instant` on HPC, this is `/rsgrps/mfh4/Ariella/results_macsSwig_AJmodels_instant`
* out_path = the path to write the ABCtoolbox input file to
* model = `1`, `2`, or `3`
* bucket_id = The PBS_ID that was used to make the bucket
* header_file_name = file with one line containing the desired header for the model. Should be tab delimited in the form    
Sim parameters  statistics
* combine_function = `original`, `same`, or `duplicate`  
For now, only use the `same` option. This will make a ABCtoolbox input file with exactly the same parameters and statsitics as in the header file.  

e.g.  
`python /vol_c/src/macsswig_simsaj/post_process.py /vol_c/results_macsSwig_AJmodels_instant /vol_c/results_macsSwig_AJmodels_instant/intermediate 1 691009 header_M1_222.txt same`

To run on all buckets use the script `run_post_process.sh`    
This will run `post_process.py` on all the buckets in the path, and combined files for each bucket in an intermediary directory.  
e.g.  
`/vol_c/src/macsswig_simsaj/run_post_process.sh /vol_c/results_macsSwig_AJmodels_instant /vol_c/results_macsSwig_AJmodels_instant/intermediate`

To combine all post processed bucket files for ABCtoolbox input use `combine_HPC_final.sh`,    
Takes one arguement, `output_path`  
e.g.  
`combine_HPC_final.sh /vol_c/results_macsSwig_AJmodels_instant/intermediate`
_____________________________


# Using ABC with ABCestimator

## Obtaining and compiling the code
ABCtoolbox is available from  
https://bitbucket.org/phaentu/abctoolbox-public/overview  

To download  
`git clone https://bitbucket.org/phaentu/abctoolbox-public.git`  

If openMP is installed, some functions inside ABCtoolbox, including findStatsModelchoice can be parallelized. Compile as follows  
`g++ -O3 -o ABCtoolbox *.cpp -DUSE_OMP -fopenmp`  
See [openMp forum](http://forum.openmp.org/forum/viewtopic.php?f=3&t=1993&p=7809#p7809) 

If openMP is not installed, compile as follows:  
`g++ -O3 -o ABCtoolbox *.cpp`


## Verifying simulation results
To plot the distribution of the parameters, summary statistics, and PCA of summary statsitics, use the R script dist_plot_stats.R  
There are three argument:  
* The simulation input file (created from post processing scripts)
* The summary statistics from the real data
* The header that contains the columns you want to use

e.g.
`Rscript /vol_c/src/macsswig_simsaj/dist_plot_stats.R /vol_c/results_macsSwig_AJmodels_instant/input_ABCtoolbox_M1_HPC.txt /vol_c/ABC_AJmodels/real_output_M23.summary /vol_c/results_macsSwig_AJmodels_instant/header_M1_132.txt`

## Removing and keeping summary statistics   
To remove summary statistics or keep summary statistics from ABCtoolbox input use the scripts  
`main_subset_sim.py` and
`main_subset_real.py`  

There are two options:    
1. Remove highly correlated statistics from ABCtoolbox input.   
The first input file should be the log file after running ABCestimator with the option pruneCorrelatedStats.  
2. Keep parameter values and statistics with high power from ABCtoolbox input.   
The first input file should be the output file *greedySearchForBestStatisticsForModelChoice.txt after running ABCestimator with the option findStatsModelChoice  

The arguments are    
- ABCtoolbox output to get summary staistics from  
- ABCtoolbox input (id parameters statistics)  
- keep or remove  
- if using keep option, number of sets of summary statistics  

Run as  
`subset_stats/main_subset_sim.py ABC_searchStatsForModelChoice_OSG_50000_100greedySearchForBestStatisticsForModelChoice.txt input_ABCtoolbox_M1_8.txt keep 4`  
or  
`subset_stats/main_subset_sim.py ABC_estimate_OSG_100_50000_100.log input_ABCtoolbox_M1_8.txt remove`

## Plotting the posterior distribution
To plot the prior, marginal, and posterior distribution use the script `plot_posterior_ABtoolbox_new.R`.

e.g.
```
Rscript /vol_c/src/macsswig_simsaj/plot_posterior_ABtoolbox_new.R keepPowerStats_input_ABCtoolbox_M2_HPC_OSG_2.txt ABC_correlatedstats6_1446125_pruneCorStats_90_model0_MarginalPosteriorDensities_Obs0.txt ABC_correlatedstats6_1446125_pruneCorStats_90_model0_BestSimsParamStats_Obs0.txt
```