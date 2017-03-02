from sys import argv
from alleles_generator.macs_file import AllelesMacsFile
from bitarray import bitarray
from summary_statistics import afs_stats, afs_stats_bitarray
from ascertainment.pseudo_array import pseudo_array, pseudo_array_bits
import itertools


# from memory_profiler import profile

#@profile()
def main():
    """This is for memory_profiling purposes of the differences between lists and bitarray.
    Must give list or bitarray as argument"""

    option = argv[1]

    length = 249163442

    naf_CGI = 18
    neu_CGI = 18
    nas_CGI = 8
    nA = 66 #76  # 528
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
    daf = 0.05

    total_naf = naf_CGI + asc_nb_af
    total_neu = neu_CGI + asc_nb_eu
    total_nas = nas_CGI + asc_nb_as

    ###Total number of chromosomes
    total_asc = asc_nb_af + asc_nb_eu + asc_nb_as
    total = total_CGI + total_asc

    ####Get the positions of the SNPs that are on the chip
    snp_file = 'ill_650_test.bed'  # SNP file
    fileSNP = open(snp_file, 'r')
    SNP = []
    for line in fileSNP:
        SNP.append(line)
    fileSNP.close()

    ###get sites from snp array
    snps = []
    for line_snp in SNP:
        columns = line_snp.split('\t')
        snps.append(int(columns[2]))

    if option == 'list':
        alleles_macs_file = AllelesMacsFile('tests/test_data/sites1000000.txt')
        alleles = alleles_macs_file.make_lists()[0]
        nbss = len(alleles)

        pos = []
        sim_position = alleles_macs_file.make_lists()[1]
        for i in xrange(nbss):
            position = round(float(sim_position[i]) * (float(length)))
            pos.append(position)

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

        pos_asc, nbss_asc, index_avail_sites = pseudo_array(asc_panel, daf, pos, snps)

        ############
        ##Transpose the data
        TseqAf = zip(*seqAfCGI)
        TseqEu = zip(*seqEuCGI)
        TseqAs = zip(*seqAsCGI)
        seqJ = Talleles[total_naf + total_neu + total_nas:total_naf + total_neu + total_nas + nJ]
        TseqJ = zip(*seqJ)
        seqM = Talleles[total_naf + total_neu + total_nas + nJ:total_naf + total_neu + total_nas + nJ + nM]
        TseqM = zip(*seqM)
        seqA = Talleles[total_naf + total_neu + total_nas + nJ + nM:total_naf + total_neu + total_nas + nJ + nM + nA]
        TseqA = zip(*seqA)

        ###get genotypes for pseudo array
        allelesAf_asc = []
        allelesEu_asc = []
        allelesAs_asc = []
        allelesJ_asc = []
        allelesM_asc = []
        allelesA_asc = []

        ##get the ascertained SNPs in the populations of interest
        if (nbss_asc == len(index_avail_sites)):
            for x in xrange(nbss_asc):
                allelesAf_asc.append(TseqAf[pos_asc[x]])
                allelesEu_asc.append(TseqEu[pos_asc[x]])
                allelesAs_asc.append(TseqAs[pos_asc[x]])
                allelesJ_asc.append(TseqJ[pos_asc[x]])
                allelesM_asc.append(TseqM[pos_asc[x]])
                allelesA_asc.append(TseqA[pos_asc[x]])

        elif (len(index_avail_sites) > nbss_asc):
            for x in xrange(len(pos_asc)):
                allelesAf_asc.append(TseqAf[index_avail_sites[pos_asc[x]]])
                allelesEu_asc.append(TseqEu[index_avail_sites[pos_asc[x]]])
                allelesAs_asc.append(TseqAs[index_avail_sites[pos_asc[x]]])
                allelesJ_asc.append(TseqJ[index_avail_sites[pos_asc[x]]])
                allelesM_asc.append(TseqM[index_avail_sites[pos_asc[x]]])
                allelesA_asc.append(TseqA[index_avail_sites[pos_asc[x]]])

        ###Genotypes for the ascertained SNPs
        seqAf_asc = zip(*allelesAf_asc)
        seqEu_asc = zip(*allelesEu_asc)
        seqAs_asc = zip(*allelesAs_asc)
        seqJ_asc = zip(*allelesJ_asc)
        seqM_asc = zip(*allelesM_asc)
        seqA_asc = zip(*allelesA_asc)

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

        seqJ = Talleles[total_naf + total_neu + total_nas:total_naf + total_neu + total_nas + nJ]
        TseqJ = zip(*seqJ)
        seqM = Talleles[total_naf + total_neu + total_nas + nJ:total_naf + total_neu + total_nas + nJ + nM]
        TseqM = zip(*seqM)
        seqA = Talleles[total_naf + total_neu + total_nas + nJ + nM:total_naf + total_neu + total_nas + nJ + nM + nA]
        TseqA = zip(*seqA)


    elif option == 'bitarray':
        seq_macs_file = AllelesMacsFile('tests/test_data/sites1000000.txt')
        seqAF_bits = seq_macs_file.make_bitarray_seq(0, total_naf)
        seqEu_bits = seq_macs_file.make_bitarray_seq(total_naf, total_naf + total_neu)
        seqAs_bits = seq_macs_file.make_bitarray_seq(total_naf + total_neu, total_naf + total_neu + total_nas)
        seqJ_bits = seq_macs_file.make_bitarray_seq(total_naf + total_neu + total_nas, total_naf + total_neu + total_nas + nJ)
        seqM_bits = seq_macs_file.make_bitarray_seq(total_naf + total_neu + total_nas + nJ, total_naf + total_neu + total_nas + nJ + nM)
        seqA_bits = seq_macs_file.make_bitarray_seq(total_naf + total_neu + total_nas + nJ + nM, total_naf + total_neu + total_nas + nJ + nM + nA)

        nbss = len(seqAF_bits) / total_naf

        pos = []
        sim_position = seq_macs_file.make_bitarray()[1]
        for i in xrange(nbss):
            position = round(float(sim_position[i]) * (float(length)))
            pos.append(position)

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

        asc_panel_bits = bitarray()
        for site in xrange(0, nbss):
            asc_panel_bits.extend(seqAF_bits[site*total_naf+naf_CGI:site*total_naf+total_naf])
            asc_panel_bits.extend(seqEu_bits[site * total_neu + neu_CGI:site * total_neu + total_neu])
            asc_panel_bits.extend(seqAs_bits[site * total_nas + nas_CGI:site * total_nas + total_nas])

        pos_asc, nbss_asc, index_avail_sites = pseudo_array_bits(asc_panel_bits, daf, pos, snps)

        seqAf_asc_bits = bitarray()
        seqEu_asc_bits = bitarray()
        seqAs_asc_bits = bitarray()
        seqJ_asc_bits = bitarray()
        seqM_asc_bits = bitarray()
        seqA_asc_bits = bitarray()
        if (nbss_asc == len(index_avail_sites)):
            for x in xrange(0,nbss_asc):
                seqAf_asc_bits.extend(seqAfCGI_bits[pos_asc[x]*naf_CGI:pos_asc[x]*naf_CGI+naf_CGI])
                seqEu_asc_bits.extend(seqEuCGI_bits[pos_asc[x] * neu_CGI:pos_asc[x] * neu_CGI + neu_CGI])
                seqAs_asc_bits.extend(seqAsCGI_bits[pos_asc[x] * nas_CGI:pos_asc[x] * nas_CGI + nas_CGI])
                seqJ_asc_bits.extend(seqJ_bits[pos_asc[x] * nJ:pos_asc[x] * nJ + nJ])
                seqM_asc_bits.extend(seqM_bits[pos_asc[x] * nM:pos_asc[x] * nM + nM])
                seqA_asc_bits.extend(seqA_bits[pos_asc[x] * nA:pos_asc[x] * nA + nA])
        elif (len(index_avail_sites) > nbss_asc):
            for x in xrange(0, len(pos_asc)):
                seqAf_asc_bits.extend(seqAfCGI_bits[index_avail_sites[pos_asc[x]] * naf_CGI:index_avail_sites[pos_asc[x]] * naf_CGI + naf_CGI])
                seqEu_asc_bits.extend(seqEuCGI_bits[index_avail_sites[pos_asc[x]] * neu_CGI:index_avail_sites[pos_asc[x]] * neu_CGI + neu_CGI])
                seqAs_asc_bits.extend(seqAsCGI_bits[index_avail_sites[pos_asc[x]] * nas_CGI:index_avail_sites[pos_asc[x]] * nas_CGI + nas_CGI])
                seqJ_asc_bits.extend(seqJ_bits[index_avail_sites[pos_asc[x]] * nJ:index_avail_sites[pos_asc[x]] * nJ + nJ])
                seqM_asc_bits.extend(seqM_bits[index_avail_sites[pos_asc[x]] * nM:index_avail_sites[pos_asc[x]] * nM + nM])
                seqA_asc_bits.extend(seqA_bits[index_avail_sites[pos_asc[x]] * nA:index_avail_sites[pos_asc[x]] * nA + nA])

        res = []
        Af_res = []
        Af_res.extend(afs_stats_bitarray.base_S_ss(seqAfCGI_bits, naf_CGI))
        pi_AfCGI = afs_stats_bitarray.Pi2(Af_res[3], naf_CGI)
        Af_res.append(afs_stats_bitarray.Tajimas(pi_AfCGI, Af_res[0], naf_CGI))
        del (Af_res[3])
        res.extend(Af_res)
        head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'

        Eu_res = []
        Eu_res.extend(afs_stats_bitarray.base_S_ss(seqEuCGI_bits, neu_CGI))
        pi_EuCGI = afs_stats_bitarray.Pi2(Eu_res[3], neu_CGI)
        Eu_res.append(afs_stats_bitarray.Tajimas(pi_EuCGI, Eu_res[0], neu_CGI))
        del (Eu_res[3])
        res.extend(Eu_res)
        head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'

        As_res = []
        As_res.extend(afs_stats_bitarray.base_S_ss(seqAsCGI_bits, nas_CGI))
        pi_AsCGI = afs_stats_bitarray.Pi2(As_res[3], nas_CGI)
        As_res.append(afs_stats_bitarray.Tajimas(pi_AsCGI, As_res[0], nas_CGI))
        del (As_res[3])
        res.extend(As_res)
        head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'

        ##fst between populations
        res.append(afs_stats_bitarray.FST2(seqAfCGI_bits, pi_AfCGI, naf_CGI, seqEuCGI_bits, pi_EuCGI, neu_CGI))
        res.append(afs_stats_bitarray.FST2(seqAfCGI_bits, pi_AfCGI, naf_CGI, seqAsCGI_bits, pi_AsCGI, nas_CGI))
        res.append(afs_stats_bitarray.FST2(seqEuCGI_bits, pi_EuCGI, neu_CGI, seqAsCGI_bits, pi_AsCGI, nas_CGI))
        head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'

        return res


if __name__ == '__main__':
    main()
