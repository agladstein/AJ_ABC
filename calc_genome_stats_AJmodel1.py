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
    :param chr_stats_df: dataframe with statistics for each chromosome. 
    :return: chr_stats_stand_df: dataframe with statistics for each chromosome, where number of segregating sites,
     singletons, doubletons, and pi are scaled by chromosome size relative to chromosome 1.
    '''

    chr_stats_stand_df = pd.DataFrame()
    scalar = int(chr_stats_df.length[:1])
    chr_stats_stand_df['chr_ratio'] = chr_stats_df['length']/scalar
    for column in chr_stats_df:
        if 'Seg' in column or 'Sing' in column or 'Dupl' in column or 'Pi' in column:
            chr_stats_stand_df[str(column)+'_stand'] = chr_stats_df[column] / chr_stats_stand_df['chr_ratio']
        else:
            chr_stats_stand_df[column] = chr_stats_df[column]
    return chr_stats_stand_df

def summarize_stats(chr_stats_stand_df):
    '''
    summarize standardize stats across chromosomes with mean and standard deviation
    :param chr_stats_stand_df: dataframe with statistics for each chromosome, where number of segregating sites,
     singletons, doubletons, and pi are scaled by chromosome size relative to chromosome 1.
    :return: names: list of headers
    :return: values: list of parameter values and mean and standard deviation of summary stats
    '''

    names = []
    values = []
    for column in chr_stats_stand_df.drop(['chr','chr_ratio','length'], axis=1):
        if 'Seg' in column or 'Sing' in column or 'Dupl' in column or 'Pi' in column or 'TajD' in column or 'Pi' in column:
            mean = chr_stats_stand_df[column].mean(axis=0)
            names.append('{}_mean'.format(column))
            values.append(mean)

            std = chr_stats_stand_df[column].std(axis=0)
            names.append('{}_std'.format(column))
            values.append(std)
        else:
            value = chr_stats_stand_df[column].iloc[0]
            names.append(column)
            values.append(value)

    return [names, values]

def combine_IBD_stats(path, job):

    print 'combining chromosome IBD'
    # return genome_IBD


def main():

    job = argv[1]  # must be a number
    path = argv[2]

    print 'JOB', job

    chr_stats_df = create_df(path, job)
    chr_stats_stand_df = standardize_stats(chr_stats_df)
    [names, values] = summarize_stats(chr_stats_stand_df)

    genome_results_file = '{}/results_{}.txt'.format(path, job)
    try:
        os.remove(genome_results_file)
    except OSError:
        pass
    out_file = open(genome_results_file, 'a')

if __name__ == '__main__':
    main()