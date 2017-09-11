from sys import argv
import pandas as pd

from subset_stats import get_params, get_corr_stats, get_power_stats


def main():
    stats_file = argv[1]
    input_file = argv[2]
    option = argv[3]  # remove or keep

    chunksize = 10 ** 4

    print 'Reading input file'
    for input_file_df in pd.read_csv(input_file, sep="\t", chunksize=chunksize):

        if option == "remove":
            stats = get_corr_stats(stats_file)
            print 'Creating new file with parameters and summary statistics with correlations passing filter'
            out_file_name = 'pruneCorStats_' + input_file
            input_file_df.drop(input_file_df[stats], axis=1).to_csv(out_file_name, mode = 'a', sep='\t', index=False)

        elif option == "keep":
            n_sets = int(argv[4])  # number of sets of statistics from model choice power analysis

            params_df = get_params(input_file_df)
            if n_sets < 1:
                print 'Tell me how many sets of statistics to use!'
                return
            stats = get_power_stats(stats_file, n_sets)
            out_file_name = 'keepPowerStats_' + input_file
            print 'Creating new file with parameters and summary statistics with high power'
            pd.concat([params_df, input_file_df[stats]], axis=1, join='inner').to_csv(out_file_name, mode = 'a', sep='\t', index=False)
        else:
            print "You must specify remove or keep!"
            return


if __name__ == '__main__':
    main()
