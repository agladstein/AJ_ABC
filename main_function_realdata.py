from alleles_generator.real_file import AllelesReal
from bitarray import bitarray
from summary_statistics import afs_stats_bitarray
import os
import itertools


# from memory_profiler import profile

#@profile()
def main():
    """This takes in 3 real PLINK .tped files.
    PLINK .tped with CGI or 1000 genomes HapMap populations.
    PLINK .tped snp array data made from same HapMap individuals with PLINK exclude.
    PLINK .tped snp array data with populations of interest."""

    naf_CGI = 18
    neu_CGI = 18
    nas_CGI = 8
    nJ = 28
    nM = 28
    nA = 76


    print 'naf_CGI ' + str(naf_CGI)
    print 'neu_CGI ' + str(neu_CGI)
    print 'nas_CGI ' + str(nas_CGI)
    print 'nA ' + str(nA)
    print 'nJ ' + str(nJ)
    print 'nM ' + str(nM)


    seq_real_CGI_file = AllelesReal('real_data/YRI9.CEU9.CHB4.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_snpsonly_maf0.005.tped')
    seqAF_CGI_bits = seq_real_CGI_file.make_bitarray_seq(0, naf_CGI)
    seqEu_CGI_bits = seq_real_CGI_file.make_bitarray_seq(naf_CGI, naf_CGI + neu_CGI)
    seqAs_CGI_bits = seq_real_CGI_file.make_bitarray_seq(naf_CGI + neu_CGI, naf_CGI + neu_CGI + nas_CGI)

    seq_real_CGIarray_file = AllelesReal('real_data/YRI9.CEU9.CHB4.chr1.atDNA.biAllelicSNPnoDI.genotypes_hg18_Behar_HGDP_FtDNA.tped')
    seqAf_asc_bits = seq_real_CGIarray_file.make_bitarray_seq(0, naf_CGI)
    seqEu_asc_bits = seq_real_CGIarray_file.make_bitarray_seq(naf_CGI, naf_CGI + neu_CGI)
    seqAs_asc_bits = seq_real_CGIarray_file.make_bitarray_seq(naf_CGI + neu_CGI, naf_CGI + neu_CGI + nas_CGI)

    seq_real_array_file = AllelesReal('real_data/Behar_HGDP_FtDNA_Jews_MidEast_chr1_subset.tped')
    seqJ_asc_bits = seq_real_array_file.make_bitarray_seq(naf_CGI + neu_CGI + nas_CGI, naf_CGI + neu_CGI + nas_CGI + nJ)
    seqM_asc_bits = seq_real_array_file.make_bitarray_seq(naf_CGI + neu_CGI + nas_CGI + nJ, naf_CGI + neu_CGI + nas_CGI + nJ + nM)
    seqA_asc_bits = seq_real_array_file.make_bitarray_seq(naf_CGI + neu_CGI + nas_CGI + nJ + nM, naf_CGI + neu_CGI + nas_CGI + nJ + nM + nA)

    res = []

    res = []
    Af_res = []
    Af_res.extend(afs_stats_bitarray.base_S_ss(seqAF_CGI_bits, naf_CGI))
    pi_AfCGI = afs_stats_bitarray.Pi2(Af_res[3], naf_CGI)
    Af_res.append(afs_stats_bitarray.Tajimas(pi_AfCGI, Af_res[0], naf_CGI))
    del (Af_res[3])
    res.extend(Af_res)
    head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'

    Eu_res = []
    Eu_res.extend(afs_stats_bitarray.base_S_ss(seqEu_CGI_bits, neu_CGI))
    pi_EuCGI = afs_stats_bitarray.Pi2(Eu_res[3], neu_CGI)
    Eu_res.append(afs_stats_bitarray.Tajimas(pi_EuCGI, Eu_res[0], neu_CGI))
    del (Eu_res[3])
    res.extend(Eu_res)
    head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'

    As_res = []
    As_res.extend(afs_stats_bitarray.base_S_ss(seqAs_CGI_bits, nas_CGI))
    pi_AsCGI = afs_stats_bitarray.Pi2(As_res[3], nas_CGI)
    As_res.append(afs_stats_bitarray.Tajimas(pi_AsCGI, As_res[0], nas_CGI))
    del (As_res[3])
    res.extend(As_res)
    head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'

    ##fst between populations
    res.append(afs_stats_bitarray.FST2(seqAF_CGI_bits, pi_AfCGI, naf_CGI, seqEu_CGI_bits, pi_EuCGI, neu_CGI))
    res.append(afs_stats_bitarray.FST2(seqAF_CGI_bits, pi_AfCGI, naf_CGI, seqAs_CGI_bits, pi_AsCGI, nas_CGI))
    res.append(afs_stats_bitarray.FST2(seqEu_CGI_bits, pi_EuCGI, neu_CGI, seqAs_CGI_bits, pi_AsCGI, nas_CGI))
    head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'

    Af_asc = []
    ss_Af_asc = afs_stats_bitarray.base_S_ss(seqAf_asc_bits, naf_CGI)
    if (ss_Af_asc[0] == 0):
        for i in xrange(5):
            Af_asc.append(0)
        pi_Af_asc = 0
    else:
        Af_asc.extend(afs_stats_bitarray.base_S_ss(seqAf_asc_bits, naf_CGI))
        pi_Af_asc = afs_stats_bitarray.Pi2(Af_asc[3], naf_CGI)
        Af_asc.append(pi_Af_asc)
        Af_asc.append(afs_stats_bitarray.Tajimas(pi_Af_asc, Af_asc[0], naf_CGI))
        del (Af_asc[3])
    res.extend(Af_asc)
    head = 'SegS_Af_ASC\tSing_Af_ASC\tDupl_Af_ASC\tPi_Af_ASC\tTajD_Af_ASC\t'

    Eu_asc = []
    ss_Eu_asc = afs_stats_bitarray.base_S_ss(seqEu_asc_bits, neu_CGI)
    if (ss_Eu_asc[0] == 0):
        for i in xrange(5):
            Eu_asc.append(0)
        pi_Eu_asc = 0
    else:
        Eu_asc.extend(afs_stats_bitarray.base_S_ss(seqEu_asc_bits, neu_CGI))
        pi_Eu_asc = afs_stats_bitarray.Pi2(Eu_asc[3], neu_CGI)
        Eu_asc.append(pi_Eu_asc)
        Eu_asc.append(afs_stats_bitarray.Tajimas(pi_Eu_asc, Eu_asc[0], neu_CGI))
        del (Eu_asc[3])
    res.extend(Eu_asc)
    head = head + 'SegS_Eu_ASC\tSing_Eu_ASC\tDupl_Eu_ASC\tPi_Eu_ASC\tTajD_Eu_ASC\t'

    As_asc = []
    ss_As_asc = afs_stats_bitarray.base_S_ss(seqAs_asc_bits, nas_CGI)
    if (ss_As_asc[0] == 0):
        for i in xrange(5):
            As_asc.append(0)
        pi_As_asc = 0
    else:
        As_asc.extend(afs_stats_bitarray.base_S_ss(seqAs_asc_bits, nas_CGI))
        pi_As_asc = afs_stats_bitarray.Pi2(As_asc[3], nas_CGI)
        As_asc.append(pi_As_asc)
        As_asc.append(afs_stats_bitarray.Tajimas(pi_As_asc, As_asc[0], nas_CGI))
        del (As_asc[3])
    res.extend(As_asc)
    head = head + 'SegS_As_ASC\tSing_As_ASC\tDupl_As_ASC\tPi_As_ASC\tTajD_As_ASC\t'

    J_asc = []
    ss_J_asc = afs_stats_bitarray.base_S_ss(seqJ_asc_bits, nJ)
    if (ss_J_asc[0] == 0):
        for i in xrange(5):
            J_asc.append(0)
        pi_J_asc = 0
    else:
        J_asc.extend(afs_stats_bitarray.base_S_ss(seqJ_asc_bits, nJ))
        pi_J_asc = afs_stats_bitarray.Pi2(J_asc[3], nJ)
        J_asc.append(pi_J_asc)
        J_asc.append(afs_stats_bitarray.Tajimas(pi_J_asc, J_asc[0], nJ))
        del (J_asc[3])
    res.extend(J_asc)
    head = head + 'SegS_J_ASC\tSing_J_ASC\tDupl_J_ASC\tPi_J_ASC\tTajD_J_ASC\t'

    M_asc = []
    ss_M_asc = afs_stats_bitarray.base_S_ss(seqM_asc_bits, nM)
    if (ss_M_asc[0] == 0):
        for i in xrange(5):
            M_asc.append(0)
        pi_M_asc = 0
    else:
        M_asc.extend(afs_stats_bitarray.base_S_ss(seqM_asc_bits, nM))
        pi_M_asc = afs_stats_bitarray.Pi2(M_asc[3], nM)
        M_asc.append(pi_M_asc)
        M_asc.append(afs_stats_bitarray.Tajimas(pi_M_asc, M_asc[0], nM))
        del (M_asc[3])
    res.extend(M_asc)
    head = head + 'SegS_M_ASC\tSing_M_ASC\tDupl_M_ASC\tPi_M_ASC\tTajD_M_ASC\t'

    A_asc = []
    ss_A_asc = afs_stats_bitarray.base_S_ss(seqA_asc_bits, nA)
    if (ss_A_asc[0] == 0):
        for i in xrange(5):
            A_asc.append(0)
        pi_A_asc = 0
    else:
        A_asc.extend(afs_stats_bitarray.base_S_ss(seqA_asc_bits, nA))
        pi_A_asc = afs_stats_bitarray.Pi2(A_asc[3], nA)
        A_asc.append(pi_A_asc)
        A_asc.append(afs_stats_bitarray.Tajimas(pi_A_asc, A_asc[0], nA))
        del (A_asc[3])
    res.extend(A_asc)
    head = head + 'SegS_A_ASC\tSing_A_ASC\tDupl_A_ASC\tPi_A_ASC\tTajD_A_ASC\t'

    res.append(afs_stats_bitarray.FST2(seqAf_asc_bits, pi_Af_asc, naf_CGI, seqEu_asc_bits, pi_Eu_asc, neu_CGI))
    res.append(afs_stats_bitarray.FST2(seqAf_asc_bits, pi_Af_asc, naf_CGI, seqAs_asc_bits, pi_As_asc, nas_CGI))
    res.append(afs_stats_bitarray.FST2(seqEu_asc_bits, pi_Eu_asc, neu_CGI, seqAs_asc_bits, pi_As_asc, nas_CGI))
    head = head + 'FST_AfEu_ASC\tFST_AfAs_ASC_m\tFST_EuAs_ASC\t'

    res.append(afs_stats_bitarray.FST2(seqA_asc_bits, pi_A_asc, nA, seqEu_asc_bits, pi_Eu_asc, neu_CGI))
    res.append(afs_stats_bitarray.FST2(seqA_asc_bits, pi_A_asc, nA, seqJ_asc_bits, pi_J_asc, nJ))
    res.append(afs_stats_bitarray.FST2(seqA_asc_bits, pi_A_asc, nA, seqM_asc_bits, pi_M_asc, nM))
    res.append(afs_stats_bitarray.FST2(seqM_asc_bits, pi_M_asc, nM, seqJ_asc_bits, pi_J_asc, nJ))
    head = head + 'FST_AEu_ASC\tFST_AJ_ASC\tFST_AM_ASC\tFST_MJ_ASC\n'

    filesummary='real_output.summary'
    filesumm=open(filesummary,'w')
    filesumm.write(head)

    out=''
    for g in range(len(res)):
        out=out+str(res[g])+'\t'
    out=out[:-1]+'\n'

    filesumm.write(out)
    filesumm.close()

    return res


if __name__ == '__main__':
    main()
