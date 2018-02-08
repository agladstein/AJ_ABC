from sys import argv
import os
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

# import math
# import random
# from random import randint
# from ggplot import *
# import functools

def read_abc_config(abc_config_name):
    """
    Get input files used for ABC and results files output by ABC.
    :param abc_config_name: The configuration file used to run ABCtoolbox 
    :return: simName:
    :return: obsName:
    :return: outputPrefix:
    """

    simName = ""
    obsName = ""
    outputPrefix = ""

    if os.path.isfile(abc_config_name):
        print(abc_config_name)
        abc_config = open(abc_config_name, 'r')
        for line in abc_config:
            line_lst = line.split()
            arg = line_lst[0]
            if arg == "simName":
                simName = line.split()[1]
            if arg == "obsName":
                obsName = line.split()[1]
            if arg == "outputPrefix":
                outputPrefix = line.split()[1]
        abc_config.close()

    if not simName:
        print("simName not included in ABCtoolbox config file")
        exit()
    if not obsName:
        print("obsName not included in ABCtoolbox config file")
        exit()
    if not outputPrefix:
        print("outputPrefix not included in ABCtoolbox config file")
        exit()

    return [simName, obsName, outputPrefix]

def get_results_files(outputPrefix):
    """
    Define names of ABCtoolbox results files
    :param outputPrefix: the output prefix provided in the ABCtoolbox config file
    :return: names of ABCtoolbox output files
    """

    BestSimsParamStats_name = '{}model0_BestSimsParamStats_Obs0.txt'.format(outputPrefix)
    MarginalPosteriorDensities_name = '{}model0_MarginalPosteriorDensities_Obs0.txt'.format(outputPrefix)
    MarginalPosteriorCharacteristics_name = '{}model0_MarginalPosteriorCharacteristics.txt'.format(outputPrefix)
    jointPosterior_name = '{}model0_jointPosterior_8_9_Obs0.txt'.format(outputPrefix)
    MarginalPosteriorCharacteristics_reformat_name = '{}model0_MarginalPosteriorCharacteristicsReformat.txt'.format(outputPrefix)

    return [BestSimsParamStats_name,
            MarginalPosteriorDensities_name,
            MarginalPosteriorCharacteristics_name,
            jointPosterior_name,
            MarginalPosteriorCharacteristics_reformat_name]

def reformat_Characteristics(MarginalPosteriorCharacteristics_name):
    """
    reformat the ABCtoolbox output file MarginalPosteriorCharacteristics to a table with parameter as the rows and
     posterior density characteristics as columns.
    :param MarginalPosteriorCharacteristics_name: file name of ABCtoolbox output file with characteristics of posterior
    density.
    :return: df_table: pandas dataframe with parameters as rows and posterior density characteristics as columns
    """

    characteristics = ['mode', 'mean', 'median', 'q50_lower', 'q50_upper', 'q90_lower', 'q90_upper', 'q95_lower', 'q95_upper',
                  'q99_lower', 'q99_upper', 'HDI50_lower', 'HDI50_upper', 'HDI90_lower', 'HDI90_upper', 'HDI95_lower',
                  'HDI95_upper', 'HDI99_lower', 'HDI99_upper']
    n_chars = len(characteristics)

    if os.path.isfile(MarginalPosteriorCharacteristics_name):
        print('parsing {}'.format(MarginalPosteriorCharacteristics_name))

        df = pd.read_csv(MarginalPosteriorCharacteristics_name, sep = '\t').drop('dataSet', 1)
        header = list(df)

        df_list = []
        start = 0
        for i in range(1, int(len(df.columns)/n_chars)):
            param = header[start].split(characteristics[0])[0].strip('_')
            df_param = df.loc[:, header[start]:header[start + n_chars - 1]]
            df_param.columns = characteristics
            df_param['param'] = param
            df_param.set_index('param')
            df_list.append(df_param)
            start = n_chars * i
        df_table = pd.concat(df_list)

    else:
        print('{} does not exist'.format(MarginalPosteriorCharacteristics_name))
        print('Did you run ABCtoolbox in this directory?')
        exit()

    return df_table

def create_joint_df(jointPosterior_name):

    if os.path.isfile(jointPosterior_name):
        joint_NEA_NWA_df = pd.read_csv(jointPosterior_name, sep = '\t')
    else:
        print('{} does not exist'.format(jointPosterior_name))
        print('Did you run ABCtoolbox in this directory?')
        exit()

    return joint_NEA_NWA_df

def plot_joint_mtpltlb(jointPosterior_name, df_chrs_reformat):

    tbl = np.genfromtxt(jointPosterior_name, names=True)

    # density map
    NWA, NEA, z = tbl['Log10_NWA'], tbl['Log10_NEA'], tbl['density']
    NWA = np.unique(NWA)
    NEA = np.unique(NEA)
    X, Y = np.meshgrid(NWA, NEA)
    Z = z.reshape(len(NEA), len(NWA))
    plt.pcolormesh(X, Y, Z, cmap='viridis')
    colorbar = plt.colorbar()
    colorbar.set_label('Density')

    # y = x line
    plt.plot(NWA, NWA, color='black')

    # Scatterplot point
    NEA_mode = df_chrs_reformat.loc[df_chrs_reformat['param'] == 'Log10_NEA']['mode']
    NWA_mode = df_chrs_reformat.loc[df_chrs_reformat['param'] == 'Log10_NWA']['mode']
    plt.scatter(NWA_mode, NEA_mode, marker='*', facecolor='black', edgecolor='none')

    # Axes limits and labels
    plt.xlim(np.min(NWA), np.max(NWA))
    plt.xlabel('$\log_{10}$ NWA')

    plt.ylabel('$\log_{10}$ NEA')
    plt.ylim(np.min(NEA), np.max(NEA))

    plot_name = '{}.png'.format(jointPosterior_name.strip(''))
    plt.savefig(plot_name)

    return

# def plot_joint_gg(joint_NEA_NWA_df):
#
#     plot = ggplot(aes(x=Log10_NWA, y=Log10_NEA), data = joint_NEA_NWA_df) + \
#            geom_raster(aes(fill=density)) + \
#            geom_abline(intercept=0) + \
#            geom_point(x=4.12121, y=6.4, shape=8) + \
#            scale_fill_viridis() + \
#            labs(y="Log10_NEA", x="Log10_NWA") + \
#            theme(legend.title = element_blank()) + \
#         theme_bw()
#
#     return


def main():
    abc_config_name = argv[1]

    [simName, obsName, outputPrefix] = read_abc_config(abc_config_name)

    [BestSimsParamStats_name,
     MarginalPosteriorDensities_name,
     MarginalPosteriorCharacteristics_name,
     jointPosterior_name,
     MarginalPosteriorCharacteristics_reformat_name] = get_results_files(outputPrefix)

    df_chrs_reformat = reformat_Characteristics(MarginalPosteriorCharacteristics_name)
    df_chrs_reformat.to_csv(MarginalPosteriorCharacteristics_reformat_name, sep='\t')

    joint_NEA_NWA_df = create_joint_df(jointPosterior_name)

    plot_joint_mtpltlb(jointPosterior_name, df_chrs_reformat)

if __name__ == '__main__':
    main()