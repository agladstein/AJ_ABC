from sys import argv
import pandas as pd

from subset_stats import get_params, get_corr_stats, get_power_stats


def main():
    stats_file = argv[1]
    input_file = argv[2]
    option = argv[3]  # remove or keep
    n_sets = int(argv[4])  # number of sets of statistics from model choice power analysis
    if len(argv) == 5:
        sim_id = 'sim'
    else:
        sim_id = argv[5]  # the header for sim id

    input_file_df = pd.read_csv(input_file, sep="\t", index_col=sim_id)
    params_df = get_params(input_file_df)

    if option == "remove":
        stats = get_corr_stats(stats_file)
        print 'Creating new file with parameters and summary statistics with correlations passing filter'
        out_file_name = 'pruneCorStats_' + input_file

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
