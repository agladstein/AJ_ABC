class AllelesMacsSwig(object):
    """Make list of lists containing alleles from macsSwig."""

    def __init__(self, nbss, sim, total):
        """

        :param nbss: Number of simulated sites
        :param sim: Result of macsSwig simulation
        :param total: Total number of haploid simulated individuals
        """
        self.nbss = nbss
        self.sim = sim
        self.total = total

    def make_lists(self):
        alleles = []
        for x in xrange(0, self.nbss):
            loc = []
            for m in xrange(0, self.total):
                loc.append(self.sim.getSite(x, m))
            alleles.append(loc)
        return alleles