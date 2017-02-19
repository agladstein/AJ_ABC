import re


class AllelesMacsFile(object):
    """Make list of lists containing alleles from macs sites output file.
    macs sites output file is made by
    macs 20 1000 -t 0.001 -r 0.0004 -s 100000 -h 1e5 -I 2 10 10 -ej 0.0025 1 2 >trees.txt 1> sites.txt"""

    def __init__(self, macs_file_name):
        self.macs_file_name = macs_file_name

    def make_lists(self):
        alleles = []
        macs_file = open(self.macs_file_name, 'r')
        for line in macs_file:
            if re.match('SITE', line):
                columns = line.split('\t')
                site_alleles = list(columns[4].strip())
                alleles.append(site_alleles)
        macs_file.close()
        return alleles