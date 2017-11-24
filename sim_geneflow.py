import os
from sys import argv
from bitarray import bitarray
import random
from alleles_generator.real_file import AllelesReal

def sample_size(file):
    """
    get the haploid number of indviduals from a tped file. 
    :param file: tped file
    :return size: haploid sample size 
    """

    with open(file, 'r') as f:
        first_line = f.readline()
    size = len(first_line.split(' ')) - 4
    return size

def geneflow(receiving_bits, source_bits, n_receiving, n_source, p):
    """
    for each site in receiving individuals switch allele with probability p to random allele from source population.
    :param receiving_bits: bitarray of receiving samples
    :param source_bits: bitarray of source samples
    :param p: probability of switching allele
    :return admixed_bits: bitarray of generated admixed samples 
    """

    admixed_bits = bitarray(receiving_bits)

    site = 0
    for index in range(0, receiving_bits.length(), n_receiving):
        freq_source = float(source_bits[site*n_source:site*n_source+n_source].count(1))/float(source_bits[site*n_source:site*n_source+n_source].length())
        for indiv in range(n_receiving):
            if random.random() < p:
                if random.random() < freq_source:
                    admixed_bits[index + indiv] = 1
                else:
                    admixed_bits[index + indiv] = 0
        site+=1

    return admixed_bits


def print_admixed(admixed_bits, receiving_file, admixed_file, n_receiving):
    """
    Print tped file with gene flow from source population to receiving population
    :param admixed_bits: bitarray with admixed sequence
    :param receiving_file: original tped file of receiving population
    :param admixed_file: new admixed tped file
    :param n_receiving: number of haploid individuals in receiving population
    :return: 
    """

    file = open(receiving_file, 'r')
    site = 0
    for line in file:
        info = ' '.join(line.split(' ')[0:4])
        alleles = ' '.join(admixed_bits[site*n_receiving : site*n_receiving+n_receiving].to01())
        admixed_line = info + ' ' + alleles + '\n'
        admixed_file.write(admixed_line)
        site+=1

    file.close()
    admixed_file.close()


def main():
    """
    This script simulates gene flow from a source population to a receiving population.
    The two populations are made from real data (PLINK tped).
    For each site with probability p switch allele in receiving population to random allele from source population.
    
    Recombination and LD are not taken into account. 
    """

    receiving_file = argv[1] # Western AJ
    source_file = argv[2] # Khazar
    p = float(argv[3]) # probability of gene flow at site

    print 'file with receiving population data: ' + str(receiving_file)
    print 'file with source population data: ' + str(source_file)
    print 'probability of gene flow at each site: ' + str(p)

    receiving_file_type = receiving_file.split('.')[-1]
    source_file_type = source_file.split('.')[-1]

    if receiving_file_type != 'tped' or source_file_type != 'tped':
        print '***** Error: Must be a tped file ****'
        quit()
    if p > 1 or p < 0:
        print '***** Error: p must be between 0 and 1 *****'
        quit()

    # name output file
    n = (len(receiving_file.split('.')) - 1)
    base_name = '.'.join(receiving_file.split('.')[:-n])
    admixed_file_name = base_name + '_admixed' + str(p) + '.tped'

    if os.path.isfile(admixed_file_name):
        print admixed_file_name + ' already exists'
        quit()

    n_receiving = sample_size(receiving_file)
    receiving_bits = AllelesReal(receiving_file).make_bitarray_seq(0, n_receiving)

    n_source = sample_size(source_file)
    source_bits = AllelesReal(source_file).make_bitarray_seq(0, n_source)

    admixed_bits = geneflow(receiving_bits, source_bits, n_receiving, n_source, p)


    admixed_file = open(admixed_file_name, 'a')
    print_admixed(admixed_bits, receiving_file, admixed_file, n_receiving)

if __name__ == '__main__':
    main()