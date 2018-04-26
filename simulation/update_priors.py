from __future__ import division
import pandas as pd
import numpy as np
from def_params_M2_inst import param_sim_asc_rand, choose_case


def create_hist_df(histogram_name):
    histogram_df = pd.read_csv(histogram_name, sep='\t')
    return histogram_df


def sample_from_hist(histogram_df, param):
    bins = histogram_df[param]
    density_index = histogram_df.columns.get_loc(param) + 1
    density = histogram_df.iloc[:,density_index]
    probability = density / sum(density)
    value = np.random.choice(bins, p=probability)
    return value


def posterior_to_param_value(histogram_name, param):
    histogram_df = create_hist_df(histogram_name)
    return_value = str(sample_from_hist(histogram_df, param))
    return return_value


def assign_param_value(histogram_df):

    parameters_update = {}

    [parameters, para_out, daf] = param_sim_asc_rand()

    case = None
    while case == None:
        for column in range(1, len(histogram_df.columns) - 1, 2):
            param = histogram_df.columns[column]
            value = sample_from_hist(histogram_df, param)
            if 'Log10' in param:
                parameters_update[param.split('Log10_')[1]] = float(round(10 ** value))
            elif 'Asc' in param:
                parameters_update['ASC_{}'.format(param.split('_')[1])] = value
            elif 'TAF' in param:
                parameters_update['Taf'] = value
            elif param == 'daf':
                daf = value
            else:
                parameters_update[param] = value

        case, modified_Tgrowth_Af = choose_case(parameters)
        parameters['Tgrowth_Af'] = modified_Tgrowth_Af

    for param in parameters:
        if param not in parameters_update:
            parameters_update[param] = parameters[param]

    return [parameters_update, para_out, case, daf]