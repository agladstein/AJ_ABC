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
    standardize stats by relative chromosome length
    :param chr_stats_df: 
    :return: 
    '''

    chr_stats_stand_df = pd.DataFrame()
    scalar = int(chr_stats_df.length[:1])
    chr_stats_stand_df['chr_ratio'] = chr_stats_df['length']/scalar
    for column in chr_stats_df:
        if 'Seg' in column or 'Sing' in column or 'Dupl' in column or 'Pi' in column :
            chr_stats_stand_df[str(column)+'_stand'] = chr_stats_df[column]/chr_stats_stand_df['chr_ratio']
        else:
            chr_stats_stand_df[column] = chr_stats_df[column]
    return chr_stats_stand_df

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
    chr_stats_stand_df = standardize_stats(chr_stats_df)

if __name__ == '__main__':
    main()