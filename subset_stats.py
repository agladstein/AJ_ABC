from sys import argv
import pandas as pd
import itertools

def get_corr_stats(stats_file):
    '''This option should be used to remove highly correlated statistics from ABCtoolbox input.
    stats_file should be the log file after running ABCestimator with the option pruneCorrelatedStats'''

    print 'Reading summary statistics found in '+stats_file+'\n'
    stats = []
    with open(stats_file) as prune_file:
        for line in prune_file:
            if 'Removing' in line:
                columns = line.split('\'')
                stats.append(columns[1])
    return stats

def get_power_stats(stats_file, n_sets):
    ''' This option should be used to only keep parameter values and statistics with high power from ABCtoolbox input. 
    stats_file should be the output file *greedySearchForBestStatisticsForModelChoice.txt after running ABCestimator with the option findStatsModelChoice'''

    print 'Reading summary statistics found in ' + stats_file + '\n'
    keep_file_df = pd.read_csv(stats_file, sep = '\t')
    stats = list(itertools.chain(*keep_file_df.head(n_sets)['Statistics(Names)'].str.split(',').tolist()))
    return stats

def get_params(input_file_df):
    '''This assumes the first summary statistic is SegS_Af_CGI. 
    Make dataframe containing all columns up to the first summary statistic.'''

    print 'Extracting simulation id and parameters'
    n_params = input_file_df.columns.get_loc('SegS_Af_CGI')
    params_df = input_file_df.ix[:,:n_params]
    return params_df


def main():
    stats_file = argv[1]
    input_file = argv[2]
    option = argv[3] # remove or keep
    n_sets = int(argv[4]) # number of sets of statistics from model choice power analysis
    if len(argv) == 5:
        sim_id = 'sim'
    else:
        sim_id = argv[5] # the header for sim id

    input_file_df = pd.read_csv(input_file, sep="\t", index_col=sim_id)
    params_df = get_params(input_file_df)

    if option == "remove":
        stats = get_corr_stats(stats_file)
        print 'Creating new file with parameters and summary statistics with correlations passing filter'
        out_file_name = 'pruneCorStats_'+input_file
    elif option == "keep":
        if n_sets < 1:
            print 'Tell me how many sets of statistics to use!'
            return
        stats = get_power_stats(stats_file, n_sets)
        out_file_name = 'keepPowerStats_' + input_file
        print 'Creating new file with parameters and summary statistics with high power'
        pd.concat([params_df, input_file_df[stats]], axis=1, join='inner').to_csv(out_file_name, sep='\t')
    else:
        print "You must specify remove or keep!"
        return

if __name__ == '__main__':
    main()