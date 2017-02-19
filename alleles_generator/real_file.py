class AllelesReal(object):
    """Make list of lists containing alleles from real data modified .tped file (from Consuelo's bash code)
    haploid individuals in rows, sites in columns
    no whitespace"""

    def __init__(self, real_file_name):
        self.real_file_name = real_file_name

    def make_lists(self):
        talleles = []
        real_file = open(self.real_file_name, 'r')
        for line in real_file:
            talleles.append(line)
        real_file.close()
        alleles = zip(*talleles)
        return alleles