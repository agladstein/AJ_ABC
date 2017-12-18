import random
from sys import argv
from simulation import def_params_M1_inst, run_sim_M1_inst
import os
import macsSwig

'''
This script picks parameter values from priors for 22 chromosomes
and prints a file with macs_args for each chromosome.
'''

def main():

    #####set up simulations#####################################
    ############################################################

    job = argv[1]  # must be a number
    size = argv[2] # full or an integer (if it is an integer, the same integer will be used 22 times)
    param_type = argv[3] # prior, min, or max
    seed_option = int(argv[4]) # 0 (random, no seed) or integer > 0
    path = argv[5]

    print 'JOB', job

    ## Make output directory if it doesn't exist
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise

    para_macsargs_file = '{}/macsargs_{}.txt'.format(path, job)
    try:
        os.remove(para_macsargs_file)
    except OSError:
        pass
    out_file = open(para_macsargs_file, 'a')
    out_file.write('job:{}\n'.format(job))

    if seed_option > int(0):
        print 'seed: ', seed_option
        random.seed(seed_option)
    else:
        print 'seed: None'

    ####Get parameter values from priors
    if param_type == 'prior':
        param_model = def_params_M1_inst.param_sim_asc_rand()
    if param_type == 'min':
        param_model = def_params_M1_inst.param_sim_asc_min()
    if param_type == 'max':
        param_model = def_params_M1_inst.param_sim_asc_max()
    ###parameters is a dictionary with the parameter values
    parameters = param_model[0]
    para_out = param_model[1]
    daf = param_model[3]

    out_file.write('para_out:{}\n'.format(para_out))
    out_file.write('daf:{}\n'.format(daf))

    ####Samples to be simulated

    naf_CGI = 18
    neu_CGI = 18
    nas_CGI = 8
    nA = 76  # 528
    nJ = 28
    nM = 28  # 114

    print 'naf_CGI ' + str(naf_CGI)
    print 'neu_CGI ' + str(neu_CGI)
    print 'nas_CGI ' + str(nas_CGI)
    print 'nA ' + str(nA)
    print 'nJ ' + str(nJ)
    print 'nM ' + str(nM)

    total_CGI = naf_CGI + neu_CGI + nas_CGI + nA + nJ + nM
    print 'total samples ' + str(total_CGI)

    ###Discovery panel
    asc_nb_af = para_out[0]
    asc_nb_eu = para_out[1]
    asc_nb_as = para_out[2]

    total_naf = naf_CGI + asc_nb_af
    total_neu = neu_CGI + asc_nb_eu
    total_nas = nas_CGI + asc_nb_as

    ###Total number of chromosomes
    total_asc = asc_nb_af + asc_nb_eu + asc_nb_as
    total = total_CGI + total_asc


    ###Test parameters for simulation before printing file
    macs_args_test = run_sim_M1_inst.run_sim(parameters, 100, 1, total, total_naf, total_nas, total_neu, nJ, nM, nA,seed_option)
    print 'running test simulation'
    print macs_args_test
    sim = macsSwig.swigMain(len(macs_args_test), macs_args_test)
    print 'finished test simulation'


    ##Length of chromosomes
    lengths = [249163442, 243078003, 197813415, 191015739, 180695227, 170959304, 159091448, 146137372, 141069069,
               135430928, 134747460, 133630123, 96085774, 87668527, 82491127, 90079543, 81032226, 78003657, 58843222,
               62887650, 37234222, 35178458]
    for chr, length in enumerate(lengths, start=1):

        ###define simulation size
        if size != "full":
            length = int(size)

        macs_args = run_sim_M1_inst.run_sim(parameters, length, chr, total, total_naf, total_nas, total_neu, nJ, nM, nA, seed_option)
        out_file.write('macs_args_{}:{}\n'.format(chr,macs_args))

    print 'gen_macsargs_AJmodel1.py Completed'

if __name__ == '__main__':
    main()
