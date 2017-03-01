from sys import argv
from alleles_generator.macs_file import AllelesMacsFile
from bitarray import bitarray
from summary_statistics import afs_stats
from summary_statistics import afs_stats_bitarray
from ascertainment.pseudo_array import find2, add_snps
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
    daf = 0.05

    total_naf = naf_CGI + asc_nb_af
    total_neu = neu_CGI + asc_nb_eu
    total_nas = nas_CGI + asc_nb_as

    ###Total number of chromosomes
    total_asc = asc_nb_af + asc_nb_eu + asc_nb_as
    total = total_CGI + total_asc

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

        Tasc_panel = zip(*asc_panel)
        print 'number of sites in Tasc_panel:', len(Tasc_panel)
        print 'number of chromosomes in Tasc_panel:', len(Tasc_panel[0])

        #######Array with the available sites given the frequency cut off
        ##array with the frequency of all the simulated snps
        sites_freq = []
        ##array with the available sites, that pass the frequency cut-off
        avail_sites = []  ##this one has the positions of the snps
        index_avail_sites = []  ##this one has the indexes of the snps

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

        for n in xrange(len(Tasc_panel)):
            freq_site = float(Tasc_panel[n][0:len(asc_panel)].count('1')) / float(len(asc_panel))
            if freq_site >= daf and freq_site <= 1 - daf:
                sites_freq.append(freq_site)
                avail_sites.append(pos[n])
                index_avail_sites.append(n)

        nb_avail_sites = len(avail_sites)

        if (len(avail_sites) == len(snps)):
            print "number of avail_sites is equal to the number of Array snps"
            pos_asc = []
            pos_asc = index_avail_sites
            nbss_asc = len(pos_asc)
            flag_nb_asc_snps = 1


        elif (len(avail_sites) > len(snps)):
            print "number of avail_sites greater than number of Array snps"
            pos_asc = [None] * int(len(snps))  ##indexes of the SNPs that pass the frequency cut-off and position
            for i in xrange(len(snps)):  # each snp on the snp array on a chromosome
                ## y is the position of the SNPs in the array
                y = snps[i]
                ##find the closest SNP in the array
                closestleft = find2(avail_sites, y)

                if (i > 0 and pos_asc[i - 1] == closestleft and closestleft + 1 < len(avail_sites)):  ##avoid duplicates
                    closestleft = closestleft + 1  ##move one position to the right
                    pos_asc[i] = closestleft

                elif (i > 0 and pos_asc[i - 1] > closestleft and pos_asc[i - 1] + 1 < len(avail_sites)):
                    closestleft = pos_asc[i - 1] + 1
                    pos_asc[i] = closestleft

                else:
                    pos_asc[i] = closestleft
                    ###if I have duplicates at this point, it means that there were not anyt more snps to choose from
                    ###closestleft+1 or pos_asc[i-1]+1 == len(avail_sites)

            #####smoothing
            ##last index of the pos_asc
            i = len(pos_asc) - 1

            ##check if there is another position that might work better
            for j in xrange(0, i):
                if (j == i - 1 and pos_asc[j] + 1 < pos_asc[j + 1] and pos_asc[j] < (len(avail_sites) - 1) and (
                            j + 1) < len(avail_sites)):
                    d1 = abs(snps[j] - avail_sites[pos_asc[j]])
                    d2 = abs(snps[j] - avail_sites[pos_asc[j] + 1])
                    if (d2 < d1):
                        pos_asc[j] = pos_asc[j] + 1

            ##removes duplicates
            pos_asc = (list(set(pos_asc)))
            pos_asc.sort()

            nbss_asc = len(pos_asc)

            if (len(snps) == nbss_asc):
                flag_nb_asc_snps = 1
                print 'nb of asc snps equal to nb array snps'

            if (len(snps) != len(pos_asc)):
                flag_nb_asc_snps = 0
                print 'nb of asc snps not equal to nb array snps'

                diff = int(len(snps) - len(pos_asc))

                for m in xrange(1, diff + 1):
                    pos_asc2 = []
                    pos_asc2 = add_snps(avail_sites, nb_avail_sites, pos_asc, nbss_asc, nb_array_snps)
                    pos_asc = pos_asc2
                    nbss_asc = len(pos_asc)

                    if nbss_asc == len(snps):
                        flag_nb_asc_snps = 1
                        break

                    else:
                        flag_nb_asc_snps = 0

            if (
                flag_nb_asc_snps == 0):  ##it means that the 1st index in pos_asc is 0; and the last is len(avail_sites)-1
                diff = int(len(snps) - len(pos_asc))
                while (len(pos_asc) != len(snps)):
                    rand_numb = randint(0, len(avail_sites) - 1)
                    # print 'random',rand_numb
                    if rand_numb not in pos_asc:
                        pos_asc.append(rand_numb)
                pos_asc.sort()
                nbss_asc = len(pos_asc)

        print 'finished making pseudo array'

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
        asc_panel_bits.extend(seqAf_ds_bits)
        asc_panel_bits.extend(seqEu_ds_bits)
        asc_panel_bits.extend(seqAs_ds_bits)


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
