from __future__ import division
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde


def kde(x, x_grid, bandwidth=0.2, **kwargs):
    """Kernel Density Estimation with Scipy"""
    kde = gaussian_kde(x, bw_method=bandwidth / x.std(ddof=1), **kwargs)
    return kde.evaluate(x_grid)

def generate_rand_from_pdf(pdf, x_grid):
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    values = np.random.rand(1000)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = x_grid[value_bins]
    return random_from_cdf

def sample_from_kernal1(histogram_name):

    histogram = open(histogram_name, 'r')
    header = histogram.readline()
    values = []
    densities = []
    for bin_line in histogram:
        print(bin_line)
        columns = bin_line.strip('\n').split('\t')
        values.append(float(columns[1]))
        densities.append(float(columns[2]))
    histogram.close()

    values_np = np.array(values)
    densities_np = np.array(densities)

    return [values_np, densities_np]

def sample_from_kernal2(histogram_df):

    data = np.random.normal(size=1000)
    hist, bins = np.histogram(data, bins=50)

    x_grid = np.linspace(min(data), max(data), 1000)
    kdepdf = kde(data, x_grid, bandwidth=0.1)
    random_from_kde = generate_rand_from_pdf(kdepdf, x_grid)

    bin_midpoints = bins[:-1] + np.diff(bins) / 2
    random_from_cdf = generate_rand_from_pdf(hist, bin_midpoints)

    return


def create_hist_df(histogram_name):
    histogram_df = pd.read_csv(histogram_name, sep='\t')
    return histogram_df

def create_param_hist(histogram_df, param):
    bins = histogram_df[param]
    density_head = str(param) + '.density'
    density = histogram_df[density_head]
    probability = density / sum(density)
    return [bins, probability]

def sample_from_hist(histogram_df, param):
    bins = histogram_df[param]
    density_head = str(param) + '.density'
    density = histogram_df[density_head]
    probability = density / sum(density)
    value = np.random.choice(bins, p=probability)
    return value

def assign_param_value(histogram_df):

    para_out = []
    parameters = {}

    for column in range(1, len(histogram_df.columns) - 1, 2):
        param = histogram_df.columns[column]
        value = sample_from_hist(histogram_df, param)

    return
