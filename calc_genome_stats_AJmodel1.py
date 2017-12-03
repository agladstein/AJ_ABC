from sys import argv
import os
import pandas as pd

'''
This script combines the summary stats from all the chromosomes into genome summary stats
'''

def create_df(path, job):
    '''
    Create a pandas dataframe containing the stats from all chromosomes
    :param path: path to output directory
    :param job: job id (not including chr)
    :return: dataframe containing stats from all chromosomes
    '''
    print 'creating dataframe of stats from all chromosomes'

    chr_stats_df = pd.DataFrame()
    list_ = []
    for chr in range(1, 22+1):
        file = '{}/results_AJ_M1/results_{}_chr{}.txt'.format(path, job, chr)
        if os.path.isfile(file):
            print(file)
            df = pd.read_csv(file, sep = '\t')
            list_.append(df)
    chr_stats_df = pd.concat(list_)
    return chr_stats_df

def standardize_stats(chr_stats_df):
    '''
    standardize stats by chromosome length
    :param chr_stats_df: 
    :return: 
    '''

    # return chr_stats_stand_df

def combine_IBD_stats(path, job):

    print 'combining chromosome IBD'
    # return genome_IBD


def main():

    job = argv[1]  # must be a number
    path = argv[2]

    print 'JOB', job

    genome_results_file = '{}/results_{}.txt'.format(path, job)
    try:
        os.remove(genome_results_file)
    except OSError:
        pass
    out_file = open(genome_results_file, 'a')

    chr_stats_df = create_df(path, job)


if __name__ == '__main__':
    main()