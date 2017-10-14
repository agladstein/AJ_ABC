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

    for site in range(0, receiving_bits.length(), n_receiving):
        freq_source = source_bits[site:site+n_source].count(1)/source_bits[site:site+n_source].length()
        for indiv in range(n_receiving):
            if random.random() < p:
                if random.random() < freq_source:
                    admixed_bits[site + indiv] = 1
                else:
                    admixed_bits[site + indiv] = 0

    return admixed_bits


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

    n_receiving = sample_size(receiving_file)
    receiving_bits = AllelesReal(receiving_file).make_bitarray_seq(0, n_receiving)

    n_source = sample_size(source_file)
    source_bits = AllelesReal(source_file).make_bitarray_seq(0, n_source)

    geneflow(receiving_bits, source_bits, n_receiving, n_source, p)


if __name__ == '__main__':
    main()