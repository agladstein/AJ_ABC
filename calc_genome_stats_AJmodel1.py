from sys import argv
import os
import pandas as pd
import numpy as np
import shutil

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

def combine_IBD_lengths(path, job):
    '''
    gather IBD lengths from .match germline output for all chromosomes
    :param path: path to output directory
    :param job: job id (not including chr)
    :return: pairs: list of lists of lengths of IBD segments shared between populations
    '''

    IBDlengths_eAeA = []
    IBDlengths_wAwA = []
    IBDlengths_AA = []
    IBDlengths_JJ = []
    IBDlengths_MM = []
    IBDlengths_EE = []
    IBDlengths_eAwA = []
    IBDlengths_eAE = []
    IBDlengths_wAE = []
    IBDlengths_eAJ = []
    IBDlengths_wAJ = []
    IBDlengths_eAM = []
    IBDlengths_wAM = []
    IBDlengths_AE = []
    IBDlengths_AJ = []
    IBDlengths_AM = []
    IBDlengths_JM = []
    IBDlengths_JE = []
    IBDlengths_ME = []

    for chr in range(1, 22+1):
        germline_file_name = '{}/germline_out_AJ_M1/macs_asc_{}_chr{}.match'.format(path, job, chr)
        if os.path.isfile(germline_file_name):
            print 'reading Germline IBD output chr{}'.format(chr)
            filegermline = open(germline_file_name, 'r')
            for line in filegermline:
                pop1 = line.split()[0]
                pop2 = line.split()[2]
                segment = float(line.split()[10]) / 1000000
                pair = str(pop1) + '_' + str(pop2)
                if pair == 'EA_EA':
                    IBDlengths_eAeA.append(segment)
                if pair == 'WA_WA':
                    IBDlengths_wAwA.append(segment)
                if pair == 'A_A':
                    IBDlengths_AA.append(segment)
                if pair == 'J_J':
                    IBDlengths_JJ.append(segment)
                if pair == 'M_M':
                    IBDlengths_MM.append(segment)
                if pair == 'E_E':
                    IBDlengths_EE.append(segment)
                if pair == 'EA_WA' or pair == 'WA_EA':
                    IBDlengths_eAwA.append(segment)
                if pair == 'EA_E' or pair == 'E_EA':
                    IBDlengths_eAE.append(segment)
                if pair == 'WA_E' or pair == 'E_WA':
                    IBDlengths_wAE.append(segment)
                if pair == 'EA_J' or pair == 'J_EA':
                    IBDlengths_eAJ.append(segment)
                if pair == 'WA_J' or pair == 'J_WA':
                    IBDlengths_wAJ.append(segment)
                if pair == 'EA_M' or pair == 'M_EA':
                    IBDlengths_eAM.append(segment)
                if pair == 'WA_M' or pair == 'M_WA':
                    IBDlengths_wAM.append(segment)
                if pair == 'A_E' or pair == 'E_A':
                    IBDlengths_AE.append(segment)
                if pair == 'A_J' or pair == 'J_A':
                    IBDlengths_AJ.append(segment)
                if pair == 'A_M' or pair == 'M_A':
                    IBDlengths_AM.append(segment)
                if pair == 'J_M' or pair == 'M_J':
                    IBDlengths_JM.append(segment)
                if pair == 'J_E' or pair == 'E_J':
                    IBDlengths_JE.append(segment)
                if pair == 'M_E' or pair == 'E_M':
                    IBDlengths_ME.append(segment)
            filegermline.close()

    pairs = [IBDlengths_eAeA, IBDlengths_wAwA, IBDlengths_AA, IBDlengths_JJ, IBDlengths_MM, IBDlengths_EE,
             IBDlengths_eAwA,
             IBDlengths_eAE, IBDlengths_wAE, IBDlengths_eAJ, IBDlengths_wAJ, IBDlengths_eAM, IBDlengths_wAM,
             IBDlengths_AE, IBDlengths_AJ,
             IBDlengths_AM, IBDlengths_JM, IBDlengths_JE, IBDlengths_ME]
    return pairs

def calc_IBS_stats(pairs):
    '''
    Calculate summary statistics on lengths of IBD segments shared between populations
    :param pairs: list of lists of lengths of IBD segments shared between populations
    :return: IBD_stats: list of IBD statistics
    :return: IBD_head: list of names of IBD statistics, in the same order as the IBD_stats list.
    '''

    print 'calculating IBD summary stats'

    IBD_stats = []
    IBD_head = []
    IBDlengths_mean = []
    IBDlengths_median = []
    IBDlengths_num = []
    IBDlengths_var = []
    IBDlengths_mean30 = []
    IBDlengths_median30 = []
    IBDlengths_num30 = []
    IBDlengths_var30 = []

    for p in pairs:
        IBDlengths_num.append(len(p))
        if len(p) < 1:
            p.append(0)
        IBDlengths_mean.append(np.mean(p))
        IBDlengths_median.append(np.median(p))
        IBDlengths_var.append(np.var(p))
        #### Get IBD greater than 30 Mb
        IBDlengths30 = []
        for l in p:
            if l > 30:
                IBDlengths30.append(l)
        IBDlengths_num30.append(len(IBDlengths30))
        if len(IBDlengths30) == 0:
            IBDlengths30.append(0)
        IBDlengths_mean30.append(np.mean(IBDlengths30))
        IBDlengths_median30.append(np.median(IBDlengths30))
        IBDlengths_var30.append(np.var(IBDlengths30))

    IBD_stats.extend(IBDlengths_mean)
    IBD_head.extend(['IBD_mean_eAeA', 'IBD_mean_wAwA', 'IBD_mean_AA', 'IBD_mean_JJ', 'IBD_mean_MM', 'IBD_mean_EE', 'IBD_mean_eAwA', 'IBD_mean_eAE', 'IBD_mean_wAE', 'IBD_mean_eAJ', 'IBD_mean_wAJ', 'IBD_mean_eAM', 'IBD_mean_wAM', 'IBD_mean_AE', 'IBD_mean_AJ', 'IBD_mean_AM', 'IBD_mean_JM', 'IBD_mean_JE', 'IBD_mean_ME'])
    IBD_stats.extend(IBDlengths_median)
    IBD_head.extend(['IBD_median_eAeA', 'IBD_median_wAwA', 'IBD_median_AA', 'IBD_median_JJ', 'IBD_median_MM', 'IBD_median_EE', 'IBD_median_eAwA', 'IBD_median_eAE', 'IBD_median_wAE', 'IBD_median_eAJ', 'IBD_median_wAJ', 'IBD_median_eAM', 'IBD_median_wAM', 'IBD_median_AE', 'IBD_median_AJ', 'IBD_median_AM', 'IBD_median_JM', 'IBD_median_JE', 'IBD_median_ME'])
    IBD_stats.extend(IBDlengths_num)
    IBD_head.extend(['IBD_num_eAeA', 'IBD_num_wAwA', 'IBD_num_AA', 'IBD_num_JJ', 'IBD_num_MM', 'IBD_num_EE', 'IBD_num_eAwA', 'IBD_num_eAE', 'IBD_num_wAE', 'IBD_num_eAJ', 'IBD_num_wAJ', 'IBD_num_eAM', 'IBD_num_wAM', 'IBD_num_AE', 'IBD_num_AJ', 'IBD_num_AM', 'IBD_num_JM', 'IBD_num_JE', 'IBD_num_ME'])
    IBD_stats.extend(IBDlengths_var)
    IBD_head.extend(['IBD_var_eAeA', 'IBD_var_wAwA', 'IBD_var_AA', 'IBD_var_JJ', 'IBD_var_MM', 'IBD_var_EE', 'IBD_var_eAwA', 'IBD_var_eAE', 'IBD_var_wAE', 'IBD_var_eAJ', 'IBD_var_wAJ', 'IBD_var_eAM', 'IBD_var_wAM', 'IBD_var_AE', 'IBD_var_AJ', 'IBD_var_AM', 'IBD_var_JM', 'IBD_var_JE', 'IBD_var_ME'])

    IBD_stats.extend(IBDlengths_mean30)
    IBD_head.extend(['IBD30_mean_eAeA', 'IBD30_mean_wAwA', 'IBD30_mean_AA', 'IBD30_mean_JJ', 'IBD30_mean_MM', 'IBD30_mean_EE', 'IBD30_mean_eAwA', 'IBD30_mean_eAE', 'IBD30_mean_wAE', 'IBD30_mean_eAJ', 'IBD30_mean_wAJ', 'IBD30_mean_eAM', 'IBD30_mean_wAM', 'IBD30_mean_AE', 'IBD30_mean_AJ', 'IBD30_mean_AM', 'IBD30_mean_JM', 'IBD30_mean_JE', 'IBD30_mean_ME'])
    IBD_stats.extend(IBDlengths_median30)
    IBD_head.extend(['IBD30_median_eAeA', 'IBD30_median_wAwA', 'IBD30_median_AA', 'IBD30_median_JJ', 'IBD30_median_MM', 'IBD30_median_EE', 'IBD30_median_eAwA', 'IBD30_median_eAE', 'IBD30_median_wAE', 'IBD30_median_eAJ', 'IBD30_median_wAJ', 'IBD30_median_eAM', 'IBD30_median_wAM', 'IBD30_median_AE', 'IBD30_median_AJ', 'IBD30_median_AM', 'IBD30_median_JM', 'IBD30_median_JE', 'IBD30_median_ME'])
    IBD_stats.extend(IBDlengths_num30)
    IBD_head.extend(['IBD30_num_eAeA', 'IBD30_num_wAwA', 'IBD30_num_AA', 'IBD30_num_JJ', 'IBD30_num_MM', 'IBD30_num_EE', 'IBD30_num_eAwA', 'IBD30_num_eAE', 'IBD30_num_wAE', 'IBD30_num_eAJ', 'IBD30_num_wAJ', 'IBD30_num_eAM', 'IBD30_num_wAM', 'IBD30_num_AE', 'IBD30_num_AJ', 'IBD30_num_AM', 'IBD30_num_JM', 'IBD30_num_JE', 'IBD30_num_ME'])
    IBD_stats.extend(IBDlengths_var30)
    IBD_head.extend(['IBD30_var_eAeA', 'IBD30_var_wAwA', 'IBD30_var_AA', 'IBD30_var_JJ', 'IBD30_var_MM', 'IBD30_var_EE', 'IBD30_var_eAwA', 'IBD30_var_eAE', 'IBD30_var_wAE', 'IBD30_var_eAJ', 'IBD30_var_wAJ', 'IBD30_var_eAM', 'IBD30_var_wAM', 'IBD30_var_AE', 'IBD30_var_AJ', 'IBD30_var_AM', 'IBD30_var_JM', 'IBD30_var_JE', 'IBD30_var_ME'])

    return [IBD_stats, IBD_head]


def main():

    job = argv[1]  # must be a number
    path = argv[2]

    print 'JOB', job

    if len(os.listdir('{}/results_AJ_M1'.format(path))) == 0:
        print '{}/results_AJ_M1 is empty. You probably need to run run_sims_AJmodel1_chr_all.py'.format(path)
        exit()

    chr_stats_df = create_df(path, job)
    chr_stats_stand_df = standardize_stats(chr_stats_df)
    [names, values] = summarize_stats(chr_stats_stand_df)

    if len(os.listdir('{}/germline_out_AJ_M1'.format(path))) > 0:
        pairs = combine_IBD_lengths(path, job)
        [IBD_stats, IBD_head] = calc_IBS_stats(pairs)
    else:
        print '{}/results_AJ_M1 is empty. If you want IBD rerun run_sims_AJmodel1_chr_all.py with the germline option'.format(path)

    genome_results_file = '{}/results_{}.txt'.format(path, job)
    try:
        os.remove(genome_results_file)
    except OSError:
        pass
    out_file = open(genome_results_file, 'a')
    out_file.write('sim\t{}\t{}\n'.format('\t'.join(names), '\t'.join(IBD_head)))
    out_file.write('{}\t{}\t{}'.format(job, '\t'.join(map(str,values)), '\t'.join(map(str,IBD_stats))))
    out_file.close()

    if os.path.getsize(genome_results_file) > 0:
        shutil.rmtree('{}/germline_out_AJ_M1'.format(path))
        shutil.rmtree('{}/results_AJ_M1'.format(path))
        shutil.rmtree('{}/sim_data_AJ_M1'.format(path))
        os.remove('{}/macsargs_{}.txt'.format(path, job))

if __name__ == '__main__':
    main()