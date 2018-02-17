import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

from simulation import update_priors

    
def plot_hist(bins, probability):
    plt.plot(bins, probability, color='black')
    return 

def plot_samples_from_hist(bins, probability):
    
    data = np.random.choice(bins, 10000, p=probability)
    # plt.hist(data, 50, normed=True, alpha=0.5, label='hist')
    
    #normalize so sums to 1.
    results, edges = np.histogram(data, 50, normed=True)
    binWidth = edges[1] - edges[0]
    plt.bar(edges[:-1], results*binWidth, binWidth, alpha=0.5)
    return

def main():

    histogram_name = "ABC_M2_estimate_1446125_10pls_1000ret_model0_MarginalPosteriorDensities_Obs0.txt"
    histogram_df = update_priors.create_hist_df(histogram_name)
    
    for column in range(1,len(histogram_df.columns) - 1, 2):
        param = histogram_df.columns[column]
        [bins, probability] = update_priors.create_param_hist(histogram_df, param)
        plot_samples_from_hist(bins, probability)
        plot_hist(bins, probability)
        plt.xlabel(param)
        plt.show()

if __name__ == '__main__':
    main()