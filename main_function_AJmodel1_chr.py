from subprocess import Popen
import os
from bitarray import bitarray
import itertools
from alleles_generator.macs_swig_alleles import AllelesMacsSwig
from ascertainment.pseudo_array import pseudo_array_bits
import macsSwig
from summary_statistics import afs_stats_bitarray

# from memory_profiler import profile

#@profile()
def main(arguments):

    chr_number = int(arguments[1]) # chromosome number
    macsargs_name = arguments[2] #Path to file created by gen_macsargs_AJmodel1.py
    snp_file = arguments[3] # SNP array template file
    run_germline = int(arguments[4]) # 0 to run GERMLINE, 1 to not run GERMLINE
    path = arguments[5] # Path for output files

    if chr_number < 1 or chr_number > 22:
        print 'Error: chromosome number must be an integer 1-22'
        exit()

    macs_args = []
    if os.path.isfile(macsargs_name):
        macsargs_file = open(macsargs_name, 'r')
        for line in macsargs_file:
            strip_line = line.strip('\n')
            if 'job' in strip_line:
                ID = strip_line.split(':')[1]
            if 'para_out' in strip_line:
                para_out = eval(strip_line.split(':')[1])
            if 'daf' in strip_line:
                daf = float(strip_line.split(':')[1])
            if 'macs_args_'+str(chr_number) == strip_line.split(':')[0]:
                macs_args = eval(strip_line.split(':')[1])
        macsargs_file.close()

        job = str(ID)+'_chr'+str(chr_number)
        length = macs_args[2]

        naf_CGI = 18
        neu_CGI = 18
        nas_CGI = 8
        nA = 76  # 528
        nWA = 38
        nEA = 38
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
        asc_nb_af = para_out[0]
        asc_nb_eu = para_out[1]
        asc_nb_as = para_out[2]

        total_naf = naf_CGI + asc_nb_af
        total_neu = neu_CGI + asc_nb_eu
        total_nas = nas_CGI + asc_nb_as

        ###Total number of chromosomes
        total_asc = asc_nb_af + asc_nb_eu + asc_nb_as
        total = total_CGI + total_asc

        # flag to check if the nb of asc SNPs is the same as the nb of Array SNPs
        flag_nb_asc_snps = 0

        ####Get the positions of the SNPs that are on the chip
        fileSNP = open(snp_file, 'r')
        # print "read SNP file"
        SNP = []
        for line in fileSNP:
            if line.split()[0] == 'chr{}'.format(chr_number):
                SNP.append(line)
        fileSNP.close()

        ###get sites from snp array
        # print "get sites from snp array"
        snps = []
        for line_snp in SNP:
            columns = line_snp.split('\t')
            snps.append(int(columns[2]))

        ##############START SIMULATIONS##############

        ##flag to check if the simulation work (generate the number of file
        flag_sim = False
        rep = 1
        while flag_sim == False:

            #####Run simulations
            print macs_args
            print 'running simulation'
            sim = macsSwig.swigMain(len(macs_args), macs_args)
            print 'finished simulation'

            ##number of segregating sites
            nbss = sim.getNumSites()
            print 'number sites in simulation', nbss
            print 'number of chromosomes', total

            ##get position of the simulated sites and scale it to the "real" position in the SNP chip
            pos = []
            for i in xrange(nbss):
                position = round(sim.getPosition(i) * (float(length)))  # you have to give the number of the snp
                pos.append(position)

            ###Get data from the simulations
            seq_macsswig = AllelesMacsSwig(nbss, sim, total)
            seqAF_bits = seq_macsswig.make_bitarray_seq(0, total_naf)
            seqEu_bits = seq_macsswig.make_bitarray_seq(total_naf, total_naf + total_neu)
            seqAs_bits = seq_macsswig.make_bitarray_seq(total_naf + total_neu, total_naf + total_neu + total_nas)
            seqJ_bits = seq_macsswig.make_bitarray_seq(total_naf + total_neu + total_nas, total_naf + total_neu + total_nas + nJ)
            seqM_bits = seq_macsswig.make_bitarray_seq(total_naf + total_neu + total_nas + nJ, total_naf + total_neu + total_nas + nJ + nM)
            seqA_bits = seq_macsswig.make_bitarray_seq(total_naf + total_neu + total_nas + nJ + nM, total_naf + total_neu + total_nas + nJ + nM + nA)
            seqEA_bits = seq_macsswig.make_bitarray_seq(total_naf + total_neu + total_nas + nJ + nM, total_naf + total_neu + total_nas + nJ + nM + nEA)
            seqWA_bits = seq_macsswig.make_bitarray_seq(total_naf + total_neu + total_nas + nJ + nM + nEA, total_naf + total_neu + total_nas + nJ + nM + nEA + nWA)

            del sim

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

            ####Discovery subset. Put all the samples together to calculate the daf and select SNPs (matching distance as the array)
            asc_panel_bits = bitarray()
            for site in xrange(0, nbss):
                asc_panel_bits.extend(seqAF_bits[site*total_naf+naf_CGI:site*total_naf+total_naf])
                asc_panel_bits.extend(seqEu_bits[site * total_neu + neu_CGI:site * total_neu + total_neu])
                asc_panel_bits.extend(seqAs_bits[site * total_nas + nas_CGI:site * total_nas + total_nas])

            print 'number of chromosomes in asc_panel:', asc_panel_bits.length()/nbss

            ####Get pseudo array sites
            print 'Make pseudo array'
            print 'nb Array snps', len(snps)
            pos_asc, nbss_asc, index_avail_sites, avail_sites = pseudo_array_bits(asc_panel_bits, daf, pos, snps)

            nb_avail_sites = len(avail_sites)
            if (nb_avail_sites >= len(snps)):
                flag_sim = True
            else:
                flag_sim = False
                rep = rep + 1

        #######################
        #####Calculate summary statistics from the regions for the CGI data
        res = []

        if nbss > 0:  ###no segregating sites in the simulations which is not possible

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

        (seqAf_asc_bits, seqEu_asc_bits) = (bitarray(), bitarray())
        seqAs_asc_bits = bitarray()
        seqJ_asc_bits = bitarray()
        seqM_asc_bits = bitarray()
        seqA_asc_bits = bitarray()
        seqEA_asc_bits = bitarray()
        seqWA_asc_bits = bitarray()
        if (nbss_asc == len(index_avail_sites)):
            for x in xrange(0, nbss_asc):
                seqAf_asc_bits.extend(seqAfCGI_bits[pos_asc[x] * naf_CGI:pos_asc[x] * naf_CGI + naf_CGI])
                seqEu_asc_bits.extend(seqEuCGI_bits[pos_asc[x] * neu_CGI:pos_asc[x] * neu_CGI + neu_CGI])
                seqAs_asc_bits.extend(seqAsCGI_bits[pos_asc[x] * nas_CGI:pos_asc[x] * nas_CGI + nas_CGI])
                seqJ_asc_bits.extend(seqJ_bits[pos_asc[x] * nJ:pos_asc[x] * nJ + nJ])
                seqM_asc_bits.extend(seqM_bits[pos_asc[x] * nM:pos_asc[x] * nM + nM])
                seqA_asc_bits.extend(seqA_bits[pos_asc[x] * nA:pos_asc[x] * nA + nA])
                seqEA_asc_bits.extend(seqEA_bits[pos_asc[x] * nEA:pos_asc[x] * nEA + nEA])
                seqWA_asc_bits.extend(seqWA_bits[pos_asc[x] * nWA:pos_asc[x] * nWA + nWA])
        elif (len(index_avail_sites) > nbss_asc):
            for x in xrange(0, len(pos_asc)):
                seqAf_asc_bits.extend(seqAfCGI_bits[index_avail_sites[pos_asc[x]] * naf_CGI:index_avail_sites[pos_asc[
                    x]] * naf_CGI + naf_CGI])
                seqEu_asc_bits.extend(seqEuCGI_bits[index_avail_sites[pos_asc[x]] * neu_CGI:index_avail_sites[pos_asc[
                    x]] * neu_CGI + neu_CGI])
                seqAs_asc_bits.extend(seqAsCGI_bits[index_avail_sites[pos_asc[x]] * nas_CGI:index_avail_sites[pos_asc[
                    x]] * nas_CGI + nas_CGI])
                seqJ_asc_bits.extend(seqJ_bits[index_avail_sites[pos_asc[x]] * nJ:index_avail_sites[pos_asc[x]] * nJ + nJ])
                seqM_asc_bits.extend(seqM_bits[index_avail_sites[pos_asc[x]] * nM:index_avail_sites[pos_asc[x]] * nM + nM])
                seqA_asc_bits.extend(seqA_bits[index_avail_sites[pos_asc[x]] * nA:index_avail_sites[pos_asc[x]] * nA + nA])
                seqEA_asc_bits.extend(seqEA_bits[index_avail_sites[pos_asc[x]] * nEA:index_avail_sites[pos_asc[x]] * nEA + nEA])
                seqWA_asc_bits.extend(seqWA_bits[index_avail_sites[pos_asc[x]] * nWA:index_avail_sites[pos_asc[x]] * nWA + nWA])

        ##Make ped file
        print 'Make ped and map files'
        filenameped = str(path) + '/macs_asc_' + str(job) + '.ped'
        try:
            os.remove(filenameped)
        except OSError:
            pass
        fileped = open(filenameped, 'w')
        for indiv in xrange(0, neu_CGI, 2):
            fileped.write('E ' + str(indiv / 2 + 1) + '_E 0 0 1 -9 ')
            for bit in itertools.chain.from_iterable(
                    [seqEu_asc_bits[i:i + 2] for i in xrange(indiv, seqEu_asc_bits.length(), neu_CGI)]):
                if bit:
                    fileped.write('2 ')
                else:
                    fileped.write('1 ')
            fileped.write('\n')
        for indiv in xrange(0, nJ, 2):
            fileped.write('J ' + str(indiv / 2 + 1) + '_J 0 0 1 -9 ')
            for bit in itertools.chain.from_iterable(
                    [seqJ_asc_bits[i:i + 2] for i in xrange(indiv, seqJ_asc_bits.length(), nJ)]):
                if bit:
                    fileped.write('2 ')
                else:
                    fileped.write('1 ')
            fileped.write('\n')
        for indiv in xrange(0, nM, 2):
            fileped.write('M ' + str(indiv / 2 + 1) + '_M 0 0 1 -9 ')
            for bit in itertools.chain.from_iterable(
                    [seqM_asc_bits[i:i + 2] for i in xrange(indiv, seqM_asc_bits.length(), nM)]):
                if bit:
                    fileped.write('2 ')
                else:
                    fileped.write('1 ')
            fileped.write('\n')
        for indiv in xrange(0, nEA, 2):
            fileped.write('EA ' + str(indiv / 2 + 1) + '_EA 0 0 1 -9 ')
            for bit in itertools.chain.from_iterable(
                    [seqEA_asc_bits[i:i + 2] for i in xrange(indiv, seqEA_asc_bits.length(), nEA)]):
                if bit:
                    fileped.write('2 ')
                else:
                    fileped.write('1 ')
            fileped.write('\n')
        for indiv in xrange(0, nWA, 2):
            fileped.write('WA ' + str(indiv / 2 + 1) + '_WA 0 0 1 -9 ')
            for bit in itertools.chain.from_iterable(
                    [seqWA_asc_bits[i:i + 2] for i in xrange(indiv, seqWA_asc_bits.length(), nWA)]):
                if bit:
                    fileped.write('2 ')
                else:
                    fileped.write('1 ')
        for indiv in xrange(0, nA, 2):
            fileped.write('A ' + str(indiv / 2 + 1) + '_A 0 0 1 -9 ')
            for bit in itertools.chain.from_iterable(
                    [seqA_asc_bits[i:i + 2] for i in xrange(indiv, seqA_asc_bits.length(), nA)]):
                if bit:
                    fileped.write('2 ')
                else:
                    fileped.write('1 ')
            fileped.write('\n')
        fileped.close()

        ##Make map file
        filenamemap = str(path) + '/macs_asc_' + str(job) + '.map'
        try:
            os.remove(filenamemap)
        except OSError:
            pass
        filemap = open(filenamemap, 'a')
        map = ''
        for g in range(0, len(pos_asc)):
            map = str(chr_number) + ' ' + 'chr' + str(chr_number) + '_' + str(pos_asc[g]) + ' ' + str(
                int(avail_sites[pos_asc[g]] - 1)) + ' ' + str(int(avail_sites[pos_asc[g]])) + '\n'
            filemap.write(map)
        filemap.close()

        ########Use Germline to find IBD on pseduo array ped and map files
        filenameout = str(path) + '/macs_asc_' + str(job)

        print 'run germline? '+str(run_germline)
        if (run_germline == 0):
            print 'Running Germline on ' + str(filenameped) + ' ' + str(filenamemap)

            ###Germline seems to be outputting in the wrong unit - so I am putting the min at 3000000 so that it is 3Mb, but should be the default.
            if int(length) < 35178458:
                print 'bash ./bin/phasing_pipeline/gline.sh ./bin/germline-1-5-1/germline  ' + str(filenameped) + ' ' + str(filenamemap) + ' ' + str(filenameout) + ' "-bits 10 -min_m 10000'
                germline = Popen.wait(Popen('bash ./bin/phasing_pipeline/gline.sh ./bin/germline-1-5-1/germline  ' + str(filenameped) + ' ' + str(filenamemap) + ' ' + str(filenameout) + ' "-bits 10 -min_m 10000"', shell=True))
            else:
                print 'bash ./bin/phasing_pipeline/gline.sh ./bin/germline-1-5-1/germline  ' + str(filenameped) + ' ' + str(filenamemap) + ' ' + str(filenameout) + ' "-bits 10 -min_m 3000000"'
                germline = Popen.wait(Popen('bash ./bin/phasing_pipeline/gline.sh ./bin/germline-1-5-1/germline  ' + str(filenameped) + ' ' + str(filenamemap) + ' ' + str(filenameout) + ' "-bits 10 -min_m 3000000"',shell=True))

            print 'finished running germline'

        #########calculate summary stats from the ascertained SNPs

        if nbss_asc > 0:
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
            head = head + 'SegS_Af_ASC\tSing_Af_ASC\tDupl_Af_ASC\tPi_Af_ASC\tTajD_Af_ASC\t'

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

            EA_asc = []
            ss_EA_asc = afs_stats_bitarray.base_S_ss(seqEA_asc_bits, nEA)
            if (ss_EA_asc[0] == 0):
                for i in xrange(5):
                    EA_asc.append(0)
                pi_EA_asc = 0
            else:
                EA_asc.extend(afs_stats_bitarray.base_S_ss(seqEA_asc_bits, nEA))
                pi_EA_asc = afs_stats_bitarray.Pi2(EA_asc[3], nEA)
                EA_asc.append(pi_EA_asc)
                EA_asc.append(afs_stats_bitarray.Tajimas(pi_EA_asc, EA_asc[0], nEA))
                del (EA_asc[3])
            res.extend(EA_asc)
            head = head + 'SegS_EA_ASC\tSing_EA_ASC\tDupl_EA_ASC\tPi_EA_ASC\tTajD_EA_ASC\t'

            WA_asc = []
            ss_WA_asc = afs_stats_bitarray.base_S_ss(seqWA_asc_bits, nWA)
            if (ss_WA_asc[0] == 0):
                for i in xrange(5):
                    WA_asc.append(0)
                pi_WA_asc = 0
            else:
                WA_asc.extend(afs_stats_bitarray.base_S_ss(seqWA_asc_bits, nWA))
                pi_WA_asc = afs_stats_bitarray.Pi2(WA_asc[3], nWA)
                WA_asc.append(pi_WA_asc)
                WA_asc.append(afs_stats_bitarray.Tajimas(pi_WA_asc, WA_asc[0], nWA))
                del (WA_asc[3])
            res.extend(WA_asc)
            head = head + 'SegS_WA_ASC\tSing_WA_ASC\tDupl_WA_ASC\tPi_WA_ASC\tTajD_WA_ASC\t'


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
            res.append(afs_stats_bitarray.FST2(seqEA_asc_bits, pi_EA_asc, nEA, seqWA_asc_bits, pi_WA_asc, nWA))
            res.append(afs_stats_bitarray.FST2(seqEA_asc_bits, pi_EA_asc, nEA, seqEu_asc_bits, pi_Eu_asc, neu_CGI))
            res.append(afs_stats_bitarray.FST2(seqEA_asc_bits, pi_EA_asc, nEA, seqJ_asc_bits, pi_J_asc, nJ))
            res.append(afs_stats_bitarray.FST2(seqEA_asc_bits, pi_EA_asc, nEA, seqM_asc_bits, pi_M_asc, nM))
            res.append(afs_stats_bitarray.FST2(seqM_asc_bits, pi_M_asc, nM, seqJ_asc_bits, pi_J_asc, nJ))
            res.append(afs_stats_bitarray.FST2(seqWA_asc_bits, pi_WA_asc, nWA, seqEu_asc_bits, pi_Eu_asc, neu_CGI))
            res.append(afs_stats_bitarray.FST2(seqWA_asc_bits, pi_WA_asc, nWA, seqJ_asc_bits, pi_J_asc, nJ))
            res.append(afs_stats_bitarray.FST2(seqWA_asc_bits, pi_WA_asc, nWA, seqM_asc_bits, pi_M_asc, nM))
            head = head + 'FST_AEu_ASC\tFST_AJ_ASC\tFST_AM_ASC\tFST_eAwA_ASC\tFST_eAEu_ASC\tFST_eAJ_ASC\tFST_eAM_ASC\tFST_MJ_ASC\tFST_wAEu_ASC\tFST_wAJ_ASC\tFST_wAM_ASC\n'

        print 'finished calculating ss'

        ################
        #####write parameter values to file
        head_param='Asc_NAF\tAsc_NEU\tAsc_NCHB\tdaf\tLog10_NAF\tLog10_NANC\tLog10_NCEU\tLog10_NCHB\tLog10_NA\tLog10_NAg\tLog10_NJ\tLog10_NM\tm\tTgrowth_Af\tTAF\tTEM\tTeu_as\tTA\tTMJ\tTm\tTAg'

        result = '{}/results_{}.txt'.format(path, job)
        out_file = open(result, 'w')

        header = 'chr\tlength\t'+str(head_param)+'\t'+str(head)
        params = '\t'.join([str(r) for r in para_out])
        stats = '\t'.join([str(r) for r in res])

        out_file.write(str(header)+str(chr_number)+'\t'+str(length)+'\t'+str(params)+'\t'+str(stats))
        out_file.close()

        return [res,para_out]

    else:
        print '{} does not exist'.format(macsargs_name)
        print 'run the script gen_macsargs_AJmodel1.py first'
        exit()