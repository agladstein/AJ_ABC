from subprocess import Popen
import numpy as np
import os
from sys import argv
from alleles_generator.real_file import AllelesReal
from summary_statistics import afs_stats_bitarray

def main():
    """This takes in 3 real PLINK .tped files.
    PLINK .tped with CGI or 1000 genomes HapMap populations.
    PLINK .tped snp array data made from same HapMap individuals with PLINK exclude.
    PLINK .tped snp array data with populations of interest."""

    test = argv[1]
    if test == 'test':
        dir_data = 'tests/test_data/'
    elif test == 'real':
        dir_data = 'real_data/'

    naf_CGI = 18
    neu_CGI = 18
    nas_CGI = 8
    nJ = 28
    nM = 28
    nWA = 38
    nEA = 38

    print 'naf_CGI ' + str(naf_CGI)
    print 'neu_CGI ' + str(neu_CGI)
    print 'nas_CGI ' + str(nas_CGI)
    print 'nWA ' + str(nWA)
    print 'nEA ' + str(nEA)
    print 'nJ ' + str(nJ)
    print 'nM ' + str(nM)

    CGI_file = str(dir_data)+'YRI9.CEU9.CHB4.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_snpsonly_maf0.005'
    CGIarray_file = str(dir_data)+'YRI9.CEU9.CHB4.chr1.atDNA.biAllelicSNPnoDI.genotypes_hg18_Behar_HGDP_FtDNA'
    array_file = str(dir_data)+'Behar_HGDP_FtDNA_Jews_MidEast_chr1_subset_21509'
    print CGI_file
    print CGIarray_file
    print array_file

    seq_real_CGI_file = AllelesReal(str(CGI_file)+'.tped')
    seqAF_CGI_bits = seq_real_CGI_file.make_bitarray_seq(0, naf_CGI)
    seqEu_CGI_bits = seq_real_CGI_file.make_bitarray_seq(naf_CGI, naf_CGI + neu_CGI)
    seqAs_CGI_bits = seq_real_CGI_file.make_bitarray_seq(naf_CGI + neu_CGI, naf_CGI + neu_CGI + nas_CGI)

    seq_real_CGIarray_file = AllelesReal(str(CGIarray_file)+'.tped')
    seqAf_asc_bits = seq_real_CGIarray_file.make_bitarray_seq(0, naf_CGI)
    seqEu_asc_bits = seq_real_CGIarray_file.make_bitarray_seq(naf_CGI, naf_CGI + neu_CGI)
    seqAs_asc_bits = seq_real_CGIarray_file.make_bitarray_seq(naf_CGI + neu_CGI, naf_CGI + neu_CGI + nas_CGI)

    seq_real_array_file = AllelesReal(str(array_file)+'.tped')
    seqJ_asc_bits = seq_real_array_file.make_bitarray_seq(0, nJ)
    seqM_asc_bits = seq_real_array_file.make_bitarray_seq(nJ, nJ + nM)
    seqEA_asc_bits = seq_real_array_file.make_bitarray_seq(nJ + nM, nJ + nM + nEA)
    seqWA_asc_bits = seq_real_array_file.make_bitarray_seq(nJ + nM + nEA, nJ + nM + nEA + nWA)


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

    ########Use Germline to find IBD on pseduo array ped and map files
    run_germline = int(argv[2])
    filenameped = str(dir_data)+'Behar_HGDP_FtDNA_Jews_MidEast_YRI9.CEU9.CHB4.chr1.ped'
    filenamemap = str(dir_data)+'Behar_HGDP_FtDNA_Jews_MidEast_YRI9.CEU9.CHB4.chr1.map'
    filenameout = str(dir_data)+'Behar_HGDP_FtDNA_Jews_MidEast_YRI9.CEU9.CHB4.chr1'

    print 'run germline? '+str(run_germline)
    if (run_germline == 0):
        print 'Running Germline on ' + str(filenameped) + ' ' + str(filenamemap)
        print 'p  ' + str(filenameped) + ' ' + str(filenamemap) + ' ' + str(filenameout) + ' "-bits 10"'
        germline = Popen.wait(Popen('bash ./bin/phasing_pipeline/gline.sh ./bin/germline-1-5-1/germline  ' + str(filenameped) + ' ' + str(filenamemap) + ' ' + str(filenameout) + ' "-bits 10"', shell=True))

        print 'finished running germline'

    ########Get IBD stats from Germline output
    if os.path.isfile(str(filenameout) + '.match'):
        print 'reading Germline IBD output'
        filegermline = open(str(filenameout) + '.match', 'r')
        IBDlengths_eAeA = []
        IBDlengths_wAwA = []
        IBDlengths_JJ = []
        IBDlengths_MM = []
        IBDlengths_EE = []
        IBDlengths_eAwA = []
        IBDlengths_eAE = []
        IBDlengths_wAE = []
        IBDlengths_eAJ = []
        IBDlengths_wAJ = []
        IBDlengths_eAM = []
        IBDlengths_wAM = []
        IBDlengths_JM = []
        IBDlengths_JE = []
        IBDlengths_ME = []
        for line in filegermline:
            pop1 = line.split()[0]
            pop2 = line.split()[2]
            segment = float(line.split()[10])
            pair = str(pop1) + '_' + str(pop2)
            if pair == 'EA_EA':
                IBDlengths_eAeA.append(segment)
            if pair == 'WA_WA':
                IBDlengths_wAwA.append(segment)
            if pair == 'J_J':
                IBDlengths_JJ.append(segment)
            if pair == 'M_M':
                IBDlengths_MM.append(segment)
            if pair == 'E_E':
                IBDlengths_EE.append(segment)
            if pair == 'EA_WA' or pair == 'WA_EA':
                IBDlengths_eAwA.append(segment)
            if pair == 'EA_E' or pair == 'E_EA':
                IBDlengths_eAE.append(segment)
            if pair == 'WA_E' or pair == 'E_WA':
                IBDlengths_wAE.append(segment)
            if pair == 'EA_J' or pair == 'J_EA':
                IBDlengths_eAJ.append(segment)
            if pair == 'WA_J' or pair == 'J_WA':
                IBDlengths_wAJ.append(segment)
            if pair == 'EA_M' or pair == 'M_EA':
                IBDlengths_eAM.append(segment)
            if pair == 'WA_M' or pair == 'M_WA':
                IBDlengths_eAM.append(segment)
            if pair == 'J_M' or pair == 'M_J':
                IBDlengths_JM.append(segment)
            if pair == 'J_E' or pair == 'E_J':
                IBDlengths_JE.append(segment)
            if pair == 'M_E' or pair == 'E_M':
                IBDlengths_ME.append(segment)
        filegermline.close()

        print 'calculating summary stats'

        IBDlengths_mean = []
        IBDlengths_median = []
        IBDlengths_num = []
        IBDlengths_var = []
        IBDlengths_mean30 = []
        IBDlengths_median30 = []
        IBDlengths_num30 = []
        IBDlengths_var30 = []

        pairs = [IBDlengths_eAeA, IBDlengths_wAwA, IBDlengths_JJ, IBDlengths_MM, IBDlengths_EE, IBDlengths_eAwA,
                 IBDlengths_eAE, IBDlengths_wAE, IBDlengths_eAJ, IBDlengths_wAJ, IBDlengths_eAM, IBDlengths_wAM,
                 IBDlengths_JM, IBDlengths_JE, IBDlengths_ME]
        for p in pairs:
            IBDlengths_num.append(len(p))
            if len(p) < 1:
                p.append(0)
            IBDlengths_mean.append(np.mean(p))
            IBDlengths_median.append(np.median(p))
            IBDlengths_var.append(np.var(p))
            #### Get IBD greater than 30 Mb
            IBDlengths30 = []
            for l in p:
                if l > 30:
                    IBDlengths30.append(l)
            IBDlengths_num30.append(len(IBDlengths30))
            if len(IBDlengths30) == 0:
                IBDlengths30.append(0)
            IBDlengths_mean30.append(np.mean(IBDlengths30))
            IBDlengths_median30.append(np.median(IBDlengths30))
            IBDlengths_var30.append(np.var(IBDlengths30))


        res.extend(IBDlengths_mean)
        head = head + 'IBD_mean_eAeA\tIBD_mean_wAwA\tIBD_mean_JJ\tIBD_mean_MM\tIBD_mean_EE\tIBD_mean_eAwA\tIBD_mean_eAE\tIBD_mean_wAE\tIBD_mean_eAJ\tIBD_mean_wAJ\tIBD_mean_eAM\tIBD_mean_wAM\tIBD_mean_JM\tIBD_mean_JE\tIBD_mean_ME\t'
        res.extend(IBDlengths_median)
        head = head + 'IBD_median_eAeA\tIBD_median_wAwA\tIBD_median_JJ\tIBD_median_MM\tIBD_median_EE\tIBD_median_eAwA\tIBD_median_eAE\tIBD_median_wAE\tIBD_median_eAJ\tIBD_median_wAJ\tIBD_median_eAM\tIBD_median_wAM\tIBD_median_JM\tIBD_median_JE\tIBD_median_ME\t'
        res.extend(IBDlengths_num)
        head = head + 'IBD_num_eAeA\tIBD_num_wAwA\tIBD_num_JJ\tIBD_num_MM\tIBD_num_EE\tIBD_num_eAwA\tIBD_num_eAE\tIBD_num_wAE\tIBD_num_eAJ\tIBD_num_wAJ\tIBD_num_eAM\tIBD_num_wAM\tIBD_num_JM\tIBD_num_JE\tIBD_num_ME\t'
        res.extend(IBDlengths_var)
        head = head + 'IBD_var_eAeA\tIBD_var_wAwA\tIBD_var_JJ\tIBD_var_MM\tIBD_var_EE\tIBD_var_eAwA\tIBD_var_eAE\tIBD_var_wAE\tIBD_var_eAJ\tIBD_var_wAJ\tIBD_var_eAM\tIBD_var_wAM\tIBD_var_JM\tIBD_var_JE\tIBD_var_ME\t'

        res.extend(IBDlengths_mean30)
        head = head + 'IBD30_mean_AA\tIBD30_mean_JJ\tIBD30_mean_MM\tIBD30_mean_EE\tIBD30_mean_AE\tIBD30_mean_AJ\tIBD30_mean_AM\tIBD30_mean_JM\tIBD30_mean_JE\tIBD30_mean_ME\t'
        res.extend(IBDlengths_median30)
        head = head + 'IBD30_median_AA\tIBD30_median_JJ\tIBD30_median_MM\tIBD30_median_EE\tIBD30_median_AE\tIBD30_median_AJ\tIBD30_median_AM\tIBD30_median_JM\tIBD30_median_JE\tIBD30_median_ME\t'
        res.extend(IBDlengths_num30)
        head = head + 'IBD30_num_AA\tIBD30_num_JJ\tIBD30_num_MM\tIBD_num_EE\tIBD30_num_AE\tIBD30_num_AJ\tIBD30_num_AM\tIBD30_num_JM\tIBD30_num_JE\tIBD30_num_ME\t'
        res.extend(IBDlengths_var30)
        head = head + 'IBD30_var_AA\tIBD30_var_JJ\tIBD30_var_MM\tIBD_var_EE\tIBD30_var_AE\tIBD30_var_AJ\tIBD30_var_AM\tIBD30_var_JM\tIBD30_var_JE\tIBD30_var_ME\t'


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

    res.append(afs_stats_bitarray.FST2(seqAf_asc_bits, pi_Af_asc, naf_CGI, seqEu_asc_bits, pi_Eu_asc, neu_CGI))
    res.append(afs_stats_bitarray.FST2(seqAf_asc_bits, pi_Af_asc, naf_CGI, seqAs_asc_bits, pi_As_asc, nas_CGI))
    res.append(afs_stats_bitarray.FST2(seqEu_asc_bits, pi_Eu_asc, neu_CGI, seqAs_asc_bits, pi_As_asc, nas_CGI))
    head = head + 'FST_AfEu_ASC\tFST_AfAs_ASC_m\tFST_EuAs_ASC\t'

    res.append(afs_stats_bitarray.FST2(seqEA_asc_bits, pi_EA_asc, nEA, seqWA_asc_bits, pi_WA_asc, nWA))
    res.append(afs_stats_bitarray.FST2(seqEA_asc_bits, pi_EA_asc, nEA, seqEu_asc_bits, pi_Eu_asc, neu_CGI))
    res.append(afs_stats_bitarray.FST2(seqEA_asc_bits, pi_EA_asc, nEA, seqJ_asc_bits, pi_J_asc, nJ))
    res.append(afs_stats_bitarray.FST2(seqEA_asc_bits, pi_EA_asc, nEA, seqM_asc_bits, pi_M_asc, nM))
    res.append(afs_stats_bitarray.FST2(seqM_asc_bits, pi_M_asc, nM, seqJ_asc_bits, pi_J_asc, nJ))
    res.append(afs_stats_bitarray.FST2(seqWA_asc_bits, pi_WA_asc, nWA, seqEu_asc_bits, pi_Eu_asc, neu_CGI))
    res.append(afs_stats_bitarray.FST2(seqWA_asc_bits, pi_WA_asc, nWA, seqJ_asc_bits, pi_J_asc, nJ))
    res.append(afs_stats_bitarray.FST2(seqWA_asc_bits, pi_WA_asc, nWA, seqM_asc_bits, pi_M_asc, nM))
    head = head + 'FST_eAwA_ASC\tFST_eAEu_ASC\tFST_eAJ_ASC\tFST_eAM_ASC\tFST_MJ_ASC\tFST_wAEu_ASC\tFST_wAJ_ASC\tFST_wAM_ASC\n'

    filesummary='real_output_M23.summary'
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
