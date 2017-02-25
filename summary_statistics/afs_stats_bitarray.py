from math import sqrt

def base_S_ss(seq_bits,n):
    """Finds the number of segregating sites, number of singletons, number of doubletons, and allele frequency spectrum from bitarray
    n = sample size (in chromosomes)"""
    spec_zero = []
    for g in range(n - 1):
        spec_zero.append(0)
    var_ss = 0  # Segregating sites
    for site in xrange(0, seq_bits.length(), n):
        if seq_bits[site:site+n].any() and not seq_bits[site:site+n].all(): ##this ignores sites that have all zeros, or all ones
            var_ss = var_ss + 1
            spec_zero[seq_bits[site:site+n].count(1) - 1] = spec_zero[seq_bits[site:site+n].count(1) - 1] + 1
    if var_ss > 0:
        Ns = spec_zero[0] + spec_zero[-1]  ##number of singletons
        Nd = spec_zero[1] + spec_zero[-2]  ##number of dupletons
    else:
        Ns = 0
        Nd = 0
    return [var_ss, Ns, Nd, spec_zero]

def Pi2(spec,n):
    """Standard pi, calculated from allele frequency spectrum from bitarray
    n= sample size (in chromosomes)"""
    theta_pi=0.0
    for g in range(len(spec)):
        theta_pi=theta_pi+(2.0*float(spec[g])*(g+1.0)*(n-(g+1.0)))/(n*(n-1.0))
    return theta_pi

def Tajimas(p,S,n):
    """Calculates Tajima's D from bitarray.
    Takes in pi, number of segregating sites, and number of chromosomes"""
    if (S==0):
        return 0
    else:
        a1=0.0
        for g in range(n-1):
            a1=a1+(1.0/(g+1.0))
        a2=0.0
        for g in range(n-1):
            a2=a2+(1.0/((g+1.0)**2))
        b1=(n+1.0)/(3.0*(n-1.0))
        b2=(2*((n**2.0)+n+3))/((9*n)*(n-1))
        c1=b1-(1.0/a1)
        c2=b2-((n+2.0)/(a1*n))+(a2/(a1**2.0))
        e1=c1/a1
        e2=c2/((a1**2.0)+a2)
        TajD=(p-(S/a1))/(sqrt((e1*S)+((e2*S)*(S-1.0))))
        return TajD




def FST2(seq1_bits,pi1,nseq1,seq2_bits,pi2,nseq2):
    """FST based on pi within populations and between populations"""
    k3=0
    #Pi within populations
    pw=(pi1+pi2)/2
    for i in xrange(0, seq1_bits.length(), seq1_bits.length()/nseq1):
        for j in xrange(0, seq2_bits.length(), seq2_bits.length()/nseq2):
            #k3 = k3 + (~(~seq1_bits[i:i + seq1_bits.length()/nseq1] ^ seq2_bits[j:j + seq2_bits.length()/nseq2])).count(True)
            k3 = k3 + hamming_distance(seq1_bits[i:i + seq1_bits.length()/nseq1], seq2_bits[j:j + seq2_bits.length()/nseq2])
    pb=k3/(float(nseq1)*float(nseq2))
    if (pb==0):
        return '0'
    else:
        fst=float(1-(pw/pb))
        return fst

def hamming_distance(s1, s2):
    """Hamming distance between two strings of equal length is the number of positions at which the corresponding symbols are different"""
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))
