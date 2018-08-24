# Removing and keeping summary statistics   
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
