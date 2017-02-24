from math import sqrt

def base_S_ss(seq_bits,indivs):
    """Finds the number of segregating sites, number of singletons, number of doubletons, and allele frequency spectrum from bitarray"""
    spec_zero = []
    for g in range(indivs - 1):
        spec_zero.append(0)
    var_ss = 0  # Segregating sites
    for site in xrange(0, seq_bits.length(), indivs):
        print seq_bits[site:site+indivs]
        if seq_bits[site:site+indivs].any() and not seq_bits[site:site+indivs].all(): ##this ignores sites that have all zeros, or all ones
            var_ss = var_ss + 1
            spec_zero[seq_bits[site:site+indivs].count(1) - 1] = spec_zero[seq_bits[site:site+indivs].count(1) - 1] + 1
    if var_ss > 0:
        Ns = spec_zero[0] + spec_zero[-1]  ##number of singletons
        Nd = spec_zero[1] + spec_zero[-2]  ##number of dupletons
    else:
        Ns = 0
        Nd = 0
    return [var_ss, Ns, Nd, spec_zero]