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
* min = simulations with predetermined parameter values that will produce   simulations with shorter run times and less memory - only for testing and profiling purposes  
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
-------------------------

## Running on University of Arizona HPC
There are four University of Arizona HPC systems - Ocelote, HTC, SMP, and Cluster. All four systems shared the same storage space. Log onto any any of the HPC systems with
``ssh name@hpc.arizona.edu``
Then enter ``ocelote`` for Ocelote or ``ice`` for HTC, SMP, or Cluster.

### Submitting PBS from the command line
There are separate PBS files for each model and for each system:
main_function_AJmodel1_chr1.pbs - is run on Ocelote
main_function_AJmodel1_chr1_cluster.pbs, main_function_AJmodel1_chr1_htc.pbs, and main_function_AJmodel1_chr1_smp.pbs are run on ice.
You need to change the line ``#PBS -M agladstein@email.arizona.edu`` in all .pbs scripts to your email.
Submit a pbs script by:
`qsub main_function_AJmodel1_chr1.pbs`

### Automatically Submit PBS with crontab
`crontab -e` to edit the crontab file.
`crontab -l` to view the crontab file.
`crontab -r` to remove the crontab file.

You should use two seperate crontab files - one on Ocelote and one on ice.
```*/30 * * * * /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/checkque.sh 50000 /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/results_sims_AJ_M3 500 /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/main_function_AJmodel3_chr1.pbs >/rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/crontab_ocelote.log```

```*/30 * * * * /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/checkque_ice.sh 50000 /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/results_sims_AJ_M3 500 /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/main_function_AJmodel3_chr1_cluster.pbs clu; /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/checkque_ice.sh 50000 /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/results_sims_AJ_M3 500 /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/main_function_AJmodel3_chr1_smp.pbs smp; /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/checkque_ice.sh 50000 /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/results_sims_AJ_M3 500 /rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/main_function_AJmodel3_chr1_htc.pbs htc >/rsgrps/mfh/agladstein/Simulations/macsSwig_AJmodels/crontab_ice.log```
Use, your own absolute paths.

The crontab files run the shell scripts checkque.sh and checkque_ice.sh check the number of completed runs in the designated directory, the number of CPU hrs left to use, and the number of jobs currently in the que.
The shell scripts currently allow for a minimum of 350 hrs a day to be left for the group.
checkque.sh runs as
`checkque.sh sim_goal results_dir que_max pbs`
and checkque_ice.sh runs as
`checkque.sh sim_goal results_dir que_max pbs system`

### Basic Commands
`qstat` shows all of the jobs currently in the que or running.
`qstat -t` shows the status of all subjobs.
`qstat -f` shows details of job.
`qsub` submits a pbs script.
`qdel` stops a job.
`va` shows status of hours.



-------------------------

## Running as a Workflow on Open Science Grid

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


### Monitoring

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


### Statistics / Debugging

For successful workflows, you can generate statistics such as cumulative runtimes using the `pegasus-statistics -s all [dir]`
command. For failed workflows, `pegasus-analyzer [dir]` can help pinpoint the failures.


### Stopping / Restarting

If you want to stop a workflow, use the `pegasus-remove [dir]` command. Workflows which have stopped for some reason (removed
by the user or maybe from a longer system outage), can be restarted from where they left of with the `pegasus-run [dir]`
command.


### Links to documentation

 * [Pegasus Command Line Interface](https://pegasus.isi.edu/documentation/cli.php)
 * [Pegasus User Guide](https://pegasus.isi.edu/documentation/)
 * [OSG Connect Documentation](https://support.opensciencegrid.org/solution/categories)

