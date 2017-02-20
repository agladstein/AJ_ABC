from sys import argv
import macsSwig
from alleles_generator.macs_swig_alleles import AllelesMacsSwig


# python test_macs_swig_ag_ROH_asc.py job num_chr na nb T

def main():
    #####set up simulations#####################################
    ############################################################

    ###parameters of the simulations
    para_out = []
    mu = 2.5e-8
    rho = 1e-8
    na = int(argv[3])  # number of chromosomes
    nb = int(argv[4])
    nsamples = na + nb
    # effective population size
    NA = float(10000)
    # Time of split between A and B ##Split_time (assume 25 years per generation)
    T = int(argv[5])
    para_out.extend([T])  # you are never printing this
    # scaled parameters
    macs_theta = float(mu * 4 * NA)
    macs_rho = float(rho * 4 * NA)
    scaled_T = float(T / (4 * NA))
    lengths = [249163442, 243078003, 197813415, 191015739, 180695227, 170959304, 159091448, 146137372, 141069069,
               135430928, 134747460, 133630123, 96085774, 87668527, 82491127, 90079543, 81032226, 78003657, 58843222,
               62887650, 37234222, 35178458]

    lengths_sum = 0
    for i in lengths:
        lengths_sum = lengths_sum + i

    ##output of the simulations
    job = str(argv[1])  # id of the simulations

    chr = int(argv[2])
    for chr_number in xrange(1, chr + 1):
        print 'running chr' + str(chr_number)

        ###macs_commands
        ##macs N chr_length -t 0.001 -r 0.0004 -h 1e5 -R geneticmap.txt -I 2 na nab -ej t 1 2
        # macs_args = ['macs',str(na+nb),str(lengths[chr_number-1]),'-s','3239187','-t',str(macs_theta),'-r',str(macs_rho),'-h','1e5','-R','genetic_map_b37/genetic_map_GRCh37_chr'+str(chr_number)+'.txt.macshs','-I','2',str(na),str(nb),'-ej',str(scaled_T),'1','2']
        macs_args = ['macs', str(na + nb), '1000', '-s', '3239188', '-t', str(macs_theta), '-r', str(macs_rho), '-h',
                     '1e5', '-R', 'genetic_map_b37/genetic_map_GRCh37_chr' + str(chr_number) + '.txt.macshs', '-I', '2',
                     str(na), str(nb), '-ej', str(scaled_T), '1', '2']

        print macs_args

        #######run simulations
        print 'running simulation'
        sim = macsSwig.swigMain(len(macs_args), macs_args)

        # number of segregating sites
        nbss = sim.getNumSites()

        pos = []
        for i in xrange(nbss):
            position = round(sim.getPosition(i) * (lengths[chr_number - 1]))  # you have to give the number of the snp
            pos.append(position)

        alleles_macsswig = AllelesMacsSwig(nbss,sim,nsamples)
        alleles = alleles_macsswig.make_lists()
        print 'number of sites:', len(alleles)  # number of elements in alleles
        print 'number of chromosomes:', len(alleles[0])
        print(alleles)

        del sim

        Talleles = zip(*alleles)  # transposes alleles

        # identify chromsomes from each population by rows
        seqA = Talleles[0:na]
        # print len(seqA) #number of chromosomes from pop A
        seqB = Talleles[na:na + nb]
        # print len(seqB) #number of chromosomes from pop B

if __name__ == '__main__':
    main()
