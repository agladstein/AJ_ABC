Created by Ariella Gladstein, based on code from Consuelo Quinto Cortes and Krishna Veeramah.
agladstein@email.arizona.edu

There are three scripts that 
1. generates parameters for the simulation for 22 human chromosomes (gen_macsargs_AJmodel1.py) 
2. runs one chromosome simulation and computes chromosome summary statistics (run_sims_AJmodel1_chr_all.py)
3. combines summary statistics across chromosomes (calc_genome_stats_AJmodel1.py)

This allows for all 22 chromosomes to be simulated in parallel with the same parameter values and combine their summary statistics to one results file.

Version 2


## Usage

### gen_macsargs_AJmodel1.py

Run as  
```bash
python gen_macsargs_AJmodel1.py jobID simsize param_distribution seed output_dir
```
  
`jobID` = can be any unique value to identify the output  
`simsize` = `full` or the length of the locus to be simulated in base pairs. If `full` is used, the full human chromosomes lengths will be used.     
`param_distribution` = `prior`, `min`, or `max`  
* `prior` = simulations with parameter values given by prior distribution  
* `min` = simulations with predetermined parameter values that will produce simulations with shorter run times and less memory - only for testing and profiling purposes  
* `max` = simulations with predetermined parameter values that will produce simulations with longer run times and more memory - only for testing and profiling purposes  

`seed` = a seed value, or 0 if no seed is to be used.  
`output_dir` = path to the directory to output to. No argument will use the default of current dir `.`

e.g.:
```bash
python gen_macsargs_AJmodel1.py 1 full prior 0 output_dir
```

Test simulation:
```bash
python gen_macsargs_AJmodel1.py 1 100000 prior 352121 output_dir
```

prints output file to `${output_dir}/macsargs_${jobID}.txt`

### run_sims_AJmodel1_chr_all.py

Run as  
```bash
python run_sims_AJmodel1_chr_all.py chr_number macsargs_file snp_file run_germline output_dir
```
  
`chr_number` = chromosome number (1-22)  
`macsargs_file` = output file from gen_macsargs_AJmodel1.py with macs_args.   
`snp_file` = SNP array template (bed format)    
`run_germline` = 0 to run GERMLINE, 1 to not run GERMLINE (will try to read GERMLINE output from file)  
`output_dir` = path to the directory to output to. No argument will use the default of current dir `.`

e.g.:  
```bash
python run_sims_AJmodel1_chr_all.py 1 output_dir/macsargs_1.txt ftDNA_hg18_auto_all_uniqSNPS_rmbadsites_pruned_chr1.bed 0 output_dir
```

Test simulation:  
``python run_sims_AJmodel1_chr_all.py 1 output_dir/macsargs_1.txt ill_650_test.bed 1000000 0 output_dir``

Uses c++ programs macs and GERMLINE. For more information on these programs, see:  
https://github.com/gchen98/macs  
https://github.com/sgusev/GERMLINE  

Intermediate files go to sim_data_AJ_M1 and germline_out_AJ_M1 and are NOT rm in python script.  
Output files go to `${output_dir}/results_AJ_M1/results_${jobID}_chr${chr_number}.txt`.

### calc_genome_stats_AJmodel1.py

Run as  
```bash
python calc_genome_stats_AJmodel1.py jobID output_dir
```
  
`jobID` = the same `jobID` used in the original `gen_macsargs_AJmodel1.py`.
`output_dir` = path to the directory to output to. No argument will use the default of current dir `.`

e.g.:  
```bash
python calc_genome_stats_AJmodel1.py 1 output_dir
```

prints output file to `${output_dir}/results_${jobID}.txt`.

## Environment Setup

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
```

Having trouble getting Vagrant started? Make sure you are in the correct directory.
* start up vagrant in working directory
* once in vagrant, go to `/vagrant` to find original working directory

## Calculating stats on real data

Use the script `main_function_realdata_M23.py` to calculate summary stats on one chromosome.
It takes 6 arguments:
1. chromosome number
2. directory with real data
3. Genome (CGI or 1000 Genomes) of HapMap pops file name (.tped format)
4. Pseudo array data of HapMap pops file name (.tped format)
5. Array data of AJ, Jews, Middle Eastern pops file name (.tped format)
6. Germline option (0 = run germline, 1 == don't run germline)

e.g.:
```bash
python main_function_real_data_M23.py 1 tests/test_data YRI9.CEU9.CHB4.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_snpsonly_maf0.005.tped YRI9.CEU9.CHB4.chr1.atDNA.biAllelicSNPnoDI.genotypes_hg18_Behar_HGDP_FtDNA.tped Behar_HGDP_FtDNA_Jews_MidEast_chr1_subset_21509.tped 1
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
module show python/2/2.7.11
module load python/2/2.7.11
virtualenv -p /cm/shared/uaapps/python/2.7.11/bin/python macss_env_ocelote_2.7.11
source macss_env_ocelote_2.7.11/bin/activate
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

### Automatically Submit PBS with crontab
Contab will run commands at timed intervals. See http://crontab-generator.org/

`crontab -e` to edit the crontab file.  
`crontab -l` to view the crontab file.  
`crontab -r` to remove the crontab file.  

You should use two seperate crontab files.    
Ocelote:
```
*/30 * * * * /home/u15/agladstein/ABC/macsSwig_AJmodels/checkque_ice.sh 500000 2000 /rsgrps/mfh4/Ariella/macsSwig_AJmodels_rscale4Trel100 2 ocelote >>/home/u15/agladstein/ABC/macsSwig_AJmodels/crontab_ocelote2.log 2>&1
*/30 * * * * /home/u15/agladstein/ABC/macsSwig_AJmodels/checkque_ice.sh 500000 2000 /rsgrps/mfh4/Ariella/macsSwig_AJmodels_rscale4Trel100 1 ocelote >>/home/u15/agladstein/ABC/macsSwig_AJmodels/crontab_ocelote1.log 2>&1
*/30 * * * * /home/u15/agladstein/ABC/macsSwig_AJmodels/checkque_ice.sh 500000 2000 /rsgrps/mfh4/Ariella/macsSwig_AJmodels_rscale4Trel100 3 ocelote >>/home/u15/agladstein/ABC/macsSwig_AJmodels/crontab_ocelote3.log 2>&1
```

ICE
```
*/30 * * * * /home/u15/agladstein/ABC/macsSwig_AJmodels/checkque_ice.sh 500000 2000 /rsgrps/mfh4/Ariella/macsSwig_AJmodels_rscale4Trel100 2 cluster >>/home/u15/agladstein/ABC/macsSwig_AJmodels/crontab_clu.log 2>&1
*/30 * * * * /home/u15/agladstein/ABC/macsSwig_AJmodels/checkque_ice.sh 500000 2000 /rsgrps/mfh4/Ariella/macsSwig_AJmodels_rscale4Trel100 1 smp >>/home/u15/agladstein/ABC/macsSwig_AJmodels/crontab_smp.log 2>&1
*/30 * * * * /home/u15/agladstein/ABC/macsSwig_AJmodels/checkque_ice.sh 500000 2000 /rsgrps/mfh4/Ariella/macsSwig_AJmodels_rscale4Trel100 3 htc >>/home/u15/agladstein/ABC/macsSwig_AJmodels/crontab_htc.log 2>&1
```

Use, your own absolute paths.  
If the file switch.txt exists in /home/u15/agladstein/ABC/macsSwig_AJmodels, checkque_ice.sh will submit PBS scripts. Once the goal is reached, switch.txt will be removed.

#### Checking the que and remaining hrs
The crontab files run the shell scripts checkque.sh and checkque_ice.sh check the number of completed runs in the designated directory, the number of CPU hrs left to use, and the number of jobs currently in the que.  
The shell scripts currently allow for a minimum of 350 hrs a day to be left for the group.  
If there are no more available standard hours, it will submit jobs to qualified (on smp and Ocelote) and windfall.  
checkque.sh runs as  
`checkque.sh sim_goal results_dir que_max pbs`
and checkque_ice.sh runs as  
`checkque.sh sim_goal que_max output_dir model system`
 

#### Generating PBS with jinja
To automatically create PBS scripts for all models and HPC systerms use the shell script main_function_AJmodel_j2.sh  
This should be run from the working directory.  

On ICE:  
`./main_function_AJmodel_j2.sh htc output_dir model`  
`./main_function_AJmodel_j2.sh smp output_dir model`  
`./main_function_AJmodel_j2.sh cluster output_dir model`

On Ocelote:  
`./main_function_AJmodel_j2.sh ocelote output_dir model`  

where model = 1, 2, or 3  
This will use the template template.pbs.j2 to create pbs files.  
*Note: the virtual env is specified in main_function_AJmodel_j2 - and must already be created (with requirements installed) to use jinja.*

##### jinja Documentation
http://jinja.pocoo.org/docs/2.9/  
https://github.com/kolypto/j2cli

### Submitting PBS from the command line
You will need to edit the jinja template or pbs scripts created from the jinja template.
You need to change the line ``#PBS -M agladstein@email.arizona.edu`` in all pbs scripts to your email and change the line ``#PBS -W group_list=mfh`` to your group.
You can find your group with ``groups``.

Submit a pbs script by:
`qsub main_function_AJmodel1_chr1.pbs`

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
```
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

Find bad macsargs.txt files (run in workflow directory) 
```bash
BASENAME=`basename $PWD`; for SH_FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err.*/.sh/'`; do cat $SH_FILE | grep stash | grep macsargs_ | perl -p -i -e 's/.*\/(\w+\/\w+\/macsargs_[0-9]*.txt).*/$1/'; done | sort | uniq | while read -r line; do echo "~/stash/public/$BASENAME/$line"; done
```

Find bad macsargs.txt files with timestamp (run in workflow directory)
```bash
BASENAME=`basename $PWD`; for SH_FILE in `pegasus-analyzer | grep "error file:" | sed 's/.*: //' | sed 's/\.err.*/.sh/'`; do cat $SH_FILE | grep stash | grep macsargs_ | perl -p -i -e 's/.*\/(\w+\/\w+\/macsargs_[0-9]*.txt).*/$1/'; done | sort | uniq | while read -r line; do ls -l $HOME/stash/public/$BASENAME/$line; done
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
Use a crontab on Atmosphere to rsync completed workflows from Ocelote, OSG, and CHTC.

```bash
0 */1 * * * /vol_c/src/macsswig_simsaj/backup_transfer/rsync_HPC.sh genome
0 */6 * * * /vol_c/src/macsswig_simsaj/backup_transfer/rsync_OSG.sh genome
0 */6 * * * /vol_c/src/macsswig_simsaj/backup_transfer/rsync_CHTC.sh genome
```
#### CHTC Tikal proxy
In order to rsync on Atmosphere from CHTC, we use Tikal as a proxy.

Set up `~/.ssh/config`:
```bash
host chtc
  HostName submit-4.chtc.wisc.edu
  User nu_agladstein
  ProxyCommand ssh -W %h:%p tikal
  ControlMaster auto
  ControlPath /tmp/ssh_mux_%h_%p_%r
  ControlPersist 4h

host tikal
  HostName tikal.arl.arizona.edu
  User agladstein
  ControlMaster auto
  ControlPath /tmp/ssh_mux_%h_%p_%r
  ControlPersist 4h
```


### Combining output

#### OSG and CHTC 
The Pegasus workflow outputs `final_results.txt` for all the simulations in the workflow. 
The number of lines in the final output equals the number of simulations plus one for the header.  
To combine `final_results.txt` across multiple workflows to create the input for ABCtoolbox use the shell script `combine_Pegasus_final.sh`  
e.g.  
`/vol_c/src/macsswig_simsaj/combine_Pegasus_final.sh genome 2`

*This script will overwrite the file `/vol_c/ABC_AJmodels_${VERSION}` if it already exists.*  
*DO NOT RUN THIS SCRIPT IF THE ORIGINAL RESULTS FILES ARE NOT IN THE DIRECTORY `/vol_c/results_macsSwig_AJmodels_${VERSION}`*

#### HPC
On Atmosphere use the script `combine_HPC_final.sh` to combine all the simulations from HPC into one file.

```bash
/vol_c/src/macsswig_simsaj/combine_HPC_final.sh genome 2
```

#### OSG, CHTC, and HPC 
```bash
cd /vol_c/ABC_AJmodels_genome
cp input_ABC_OSG_CHTC_2.txt input_ABC_OSG_CHTC_HPC_2.txt
tail -n +2 input_ABC_HPC_2.txt >>input_ABC_OSG_CHTC_HPC_2.txt
```

### Backing up Atmosphere files

- backup entire /vol_c to google drive
- backup results directory to google drive


Use crontab to automatically run backup scripts every day.
```
0 1 * * * /vol_c/src/macsswig_simsaj/backup_transfer/backup_drive_volc.sh
0 1 * * * /vol_c/src/macsswig_simsaj/backup_transfer/backup_drive_results.sh genome
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
cd /vol_c
tar cf - /vol_c | /vol_c/bin/drive push -exclude-ops delete,update -no-prompt -piped backup_atmo_vol_1T/backup_atmo_vol_1T_$(date +%m%d%Y%T).tar
```

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

# Analyzing ABCtoolbox results

Environment Setup

If using Vagrant:

```bash
vagrant up
vagrant ssh
```

```bash
virtualenv -p python3 pweave_env
source pweave_env/bin/activate
sudo apt-get install python3-dev
sudo apt-get install python3-tk
pip install --upgrade pip
pip install --upgrade Pweave
pip install --upgrade pandas
pip install --user ggplot
```

## Reformat ABCtoolbox output of posterior density characteristics
```bash
pweave_env/bin/python assess_ABC_results.py test_ABC_estimate_PLS.txt
```

## Plotting the posterior distribution
To plot the prior, marginal, and posterior distribution use the script `plot_posterior_ABtoolbox_new.R`.

e.g.
```
Rscript /vol_c/src/macsswig_simsaj/plot_posterior_ABtoolbox_new.R keepPowerStats_input_ABCtoolbox_M2_HPC_OSG_2.txt ABC_correlatedstats6_1446125_pruneCorStats_90_model0_MarginalPosteriorDensities_Obs0.txt ABC_correlatedstats6_1446125_pruneCorStats_90_model0_BestSimsParamStats_Obs0.txt
```

