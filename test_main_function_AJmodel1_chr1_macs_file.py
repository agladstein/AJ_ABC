from sys import argv
from alleles_generator.macs_file import AllelesMacsFile
from bitarray import bitarray
from summary_statistics import afs_stats
import itertools


# from memory_profiler import profile

#@profile()
def main():
    """This is for memory_profiling purposes of the differences between lists and bitarray.
    Must give list or bitarray as argument"""

    option = argv[1]

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
    asc_nb_af = 10
    asc_nb_eu = 16
    asc_nb_as = 8

    total_naf = naf_CGI + asc_nb_af
    total_neu = neu_CGI + asc_nb_eu
    total_nas = nas_CGI + asc_nb_as

    ###Total number of chromosomes
    total_asc = asc_nb_af + asc_nb_eu + asc_nb_as
    total = total_CGI + total_asc

    if option == 'list':
        alleles_macs_file = AllelesMacsFile('tests/test_data/sites1000000.txt')
        alleles = alleles_macs_file.make_lists()
        nbss = len(alleles)

        #########get data from the simulations
        Talleles = zip(*alleles)

        ###Get data from the simulations
        seqAf = Talleles[0:total_naf]
        seqEu = Talleles[total_naf:total_naf + total_neu]
        seqAs = Talleles[total_naf + total_neu:total_naf + total_neu + total_nas]

        ####CGI data
        seqAfCGI = seqAf[:naf_CGI]
        seqEuCGI = seqEu[:neu_CGI]
        seqAsCGI = seqAs[:nas_CGI]

        ####Discovery subset
        seqAf_ds = seqAf[naf_CGI:total_naf]
        seqEu_ds = seqEu[neu_CGI:total_neu]
        seqAs_ds = seqAs[nas_CGI:total_nas]

        #####put all the samples together to calculate the daf and select SNPs (matching distance as the array)
        asc_panel = []
        asc_panel.extend(seqAf_ds)
        asc_panel.extend(seqEu_ds)
        asc_panel.extend(seqAs_ds)

        res = []
        Af_res = []
        Af_res.extend(afs_stats.base_S_ss(seqAfCGI, nbss))
        pi_AfCGI = afs_stats.Pi2(Af_res[3], len(seqAfCGI))
        Af_res.append(afs_stats.Tajimas(pi_AfCGI, Af_res[0], len(seqAfCGI)))
        del (Af_res[3])
        res.extend(Af_res)
        head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'

        Eu_res = []
        Eu_res.extend(afs_stats.base_S_ss(seqEuCGI, nbss))
        pi_EuCGI = afs_stats.Pi2(Eu_res[3], len(seqEuCGI))
        Eu_res.append(afs_stats.Tajimas(pi_EuCGI, Eu_res[0], len(seqEuCGI)))
        del (Eu_res[3])
        res.extend(Eu_res)
        head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'

        As_res = []
        As_res.extend(afs_stats.base_S_ss(seqAsCGI, nbss))
        pi_AsCGI = afs_stats.Pi2(As_res[3], len(seqAsCGI))
        As_res.append(afs_stats.Tajimas(pi_AsCGI, As_res[0], len(seqAsCGI)))
        del (As_res[3])
        res.extend(As_res)
        head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'

        ##fst between populations
        res.append(afs_stats.FST2(seqAfCGI, pi_AfCGI, naf_CGI, seqEuCGI, pi_EuCGI, neu_CGI))
        res.append(afs_stats.FST2(seqAfCGI, pi_AfCGI, naf_CGI, seqAsCGI, pi_AsCGI, nas_CGI))
        res.append(afs_stats.FST2(seqEuCGI, pi_EuCGI, neu_CGI, seqAsCGI, pi_AsCGI, nas_CGI))
        head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'

        return res

    elif option == 'bitarray':
        seq_macs_file = AllelesMacsFile('tests/test_data/sites1000000.txt')
        print 'seqAF_bits'
        seqAF_bits = seq_macs_file.make_bitarray_seq(0, total_naf)
        print 'seqEu_bits'
        seqEu_bits = seq_macs_file.make_bitarray_seq(total_naf, total_naf + total_neu)
        print 'seqAs_bits'
        seqAs_bits = seq_macs_file.make_bitarray_seq(total_naf + total_neu, total_naf + total_neu + total_nas)


        nbss = len(seqAF_bits) / total_naf

        ####CGI data
        seqAfCGI_bits = bitarray()
        for first_index in xrange(0, len(seqAF_bits), total_naf):
            seqAfCGI_bits.extend(seqAF_bits[first_index:first_index + naf_CGI])

        seqEuCGI_bits = bitarray()
        for first_index in xrange(0, len(seqEu_bits), total_neu):
            seqEuCGI_bits.extend(seqEu_bits[first_index:first_index + neu_CGI])

        seqAsCGI_bits = bitarray()
        for first_index in xrange(0, len(seqAs_bits), total_nas):
            seqAsCGI_bits.extend(seqAs_bits[first_index:first_index + nas_CGI])



        ####Discovery subset
        seqAf_ds_bits = bitarray()
        seqAf_ds_length = total_naf - naf_CGI
        for first_index in xrange(naf_CGI, len(seqAF_bits), total_naf):
            seqAf_ds_bits.extend(seqAF_bits[first_index:first_index + seqAf_ds_length])

        seqEu_ds_bits = bitarray()
        seqEu_ds_length = total_neu - neu_CGI
        for first_index in xrange(neu_CGI, len(seqEu_bits), total_neu):
            seqEu_ds_bits.extend(seqEu_bits[first_index:first_index + seqEu_ds_length])


        seqAs_ds_bits = bitarray()
        seqAs_ds_length = total_nas - nas_CGI
        for first_index in xrange(nas_CGI, len(seqAs_bits), total_nas):
            seqAs_ds_bits.extend(seqAs_bits[first_index:first_index + seqAs_ds_length])


        #####put all the samples together to calculate the daf and select SNPs (matching distance as the array)
        asc_panel_bits = bitarray()
        print 'asc_panel_bits'
        asc_panel_bits.extend(seqAf_ds_bits)
        print 'asc_panel_bits'
        asc_panel_bits.extend(seqEu_ds_bits)
        print 'asc_panel_bits'
        asc_panel_bits.extend(seqAs_ds_bits)


        # res = []
        # Af_res = []
        # Af_res.extend(afs_stats.base_S_ss(seqAfCGI_bits, nbss))
        # pi_AfCGI = afs_stats.Pi2(Af_res[3], len(seqAfCGI_bits))
        # Af_res.append(afs_stats.Tajimas(pi_AfCGI, Af_res[0], len(seqAfCGI_bits)))
        # del (Af_res[3])
        # res.extend(Af_res)
        # head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'
        #
        # return res


if __name__ == '__main__':
    main()
