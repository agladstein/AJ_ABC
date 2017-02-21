from sys import argv
#from main_function import main

import random
import time
from random import randint
from subprocess import Popen

import numpy as np
import os

from alleles_generator.macs_swig_alleles import AllelesMacsSwig
from ascertainment.pseudo_array import find2, add_snps
from simulation import def_params, run_sim
from summary_statistics import afs_stats

#def main():
def main(arguments):

    ###### START

    #####set up simulations#####################################
    ############################################################

    start_time = time.time()
    elapsed_time = time.time() - start_time

    ##Length of chromosomes
    lengths = [249163442, 243078003, 197813415, 191015739, 180695227, 170959304, 159091448, 146137372, 141069069,
               135430928, 134747460, 133630123, 96085774, 87668527, 82491127, 90079543, 81032226, 78003657, 58843222,
               62887650, 37234222, 35178458]

    job = int(arguments[1])  # must be a number
    print 'JOB', job

    ####Get parameter values from priors
    if arguments[4] > int(0):
        random.seed(arguments[4])
    if arguments[5] == 'rand':
        param_model = def_params.param_sim_asc_rand()
    if arguments[5] == 'min':
        param_model = def_params.param_sim_asc_min()
    if arguments[5] == 'max':
        param_model = def_params.param_sim_asc_max()
    ###parameters is a dictionary with the parameter values
    parameters = param_model[0]
    ###case is an integer that indicates which topology/model is going to be simulated
    para_out = param_model[1]
    case = param_model[2]
    daf = param_model[3]

    ####Samples to be simulated

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
    asc_nb_af = para_out[0]
    asc_nb_eu = para_out[1]
    asc_nb_as = para_out[2]

    total_naf = naf_CGI + asc_nb_af
    total_neu = neu_CGI + asc_nb_eu
    total_nas = nas_CGI + asc_nb_as

    ###Total number of chromosomes
    total_asc = asc_nb_af + asc_nb_eu + asc_nb_as
    total = total_CGI + total_asc

    ###### END

    chr=1

    #### Check if necessary directories exist.
    sim_data_dir = './sim_data_AJ_M1'
    germline_out_dir='./germline_out_AJ_M1'
    sim_values_dir='./sim_values_AJ_M1'
    results_sims_dir='./results_sims_AJ_M1'

    if not os.path.exists(sim_data_dir):
        os.makedirs(sim_data_dir)
    if not os.path.exists(germline_out_dir):
        os.makedirs(germline_out_dir)
    if not os.path.exists(sim_values_dir):
        os.makedirs(sim_values_dir)
    if not os.path.exists(results_sims_dir):
        os.makedirs(results_sims_dir)

    ##############START SIMULATIONS
    ##############

    cont=0
    reg_use=0

    ####

    for chr_number in range(1,chr+1):
        elapsed_time=time.time()-start_time
        print '***********'+str(elapsed_time)+'***********'
        print 'running chr'+str(chr_number)

        res=[]

        #flag to check if the nb of asc SNPs is the same as the nb of Array SNPs
        flag_nb_asc_snps=0

        ####Get the positions of the SNPs that are on the chip
        snp_file=arguments[2] #SNP file
        fileSNP=open(snp_file,'r')
        #print "read SNP file"
        SNP=[]
        for line in fileSNP:
            SNP.append(line)
        fileSNP.close()

        ###get sites from snp array
        #print "get sites from snp array"
        snps=[]
        for line_snp in SNP:
            columns=line_snp.split('\t')
            snps.append(int(columns[2]))

        print 'nb Array snps', len(snps)
        nb_array_snps=len(snps)
        #print snps

        ###define simulation size
        size=arguments[3]
        if size == "full":
            length=lengths[chr_number-1]
        else:
            length=size


        ##flag to check if the simulation work (generate the number of file
        flag_sim=False
        rep=1

        while flag_sim==False:


            #####Run simulations
            print 'running simulation'
            #sim=run_sim(parameters,case,length,chr_number)
            sim = run_sim.run_sim(parameters,case,length,chr_number,total,total_naf,total_nas,total_neu,nJ,nM,nA)
            print 'finished simulation'

            elapsed_time=time.time()-start_time
            print '***********'+str(elapsed_time)+'***********'

            ##number of segregating sites
            nbss=sim.getNumSites()
            #print 'number sites in simulation', nbss


            ##get position of the simulated sites and scale it to the "real" position in the SNP chip
            pos=[]
            for i in xrange(nbss):
                position=round(sim.getPosition(i)*(float(length)))#you have to give the number of the snp
                pos.append(position)

            ##alleles will have all the information of all the simulated sites for all the pops
            #sites in each row
            alleles_macsswig = AllelesMacsSwig(nbss, sim, total)
            alleles = alleles_macsswig.make_lists()
            print 'total number of sites:',len(alleles) #number of elements in alleles
            print 'total number of chromosomes:',len(alleles[0])


            del sim


            #########get data from the simulations
            Talleles=zip(*alleles)

            ###Get data from the simulations
            seqAf=Talleles[0:total_naf]
            seqEu=Talleles[total_naf:total_naf+total_neu]
            seqAs=Talleles[total_naf+total_neu:total_naf+total_neu+total_nas]


            ####CGI data
            seqAfCGI=seqAf[:naf_CGI]
            seqEuCGI=seqEu[:neu_CGI]
            seqAsCGI=seqAs[:nas_CGI]

            ####Discovery subset
            seqAf_ds=seqAf[naf_CGI:total_naf]
            seqEu_ds=seqEu[neu_CGI:total_neu]
            seqAs_ds=seqAs[nas_CGI:total_nas]

            #####put all the samples together to calculate the daf and select SNPs (matching distance as the array)
            asc_panel=[]
            asc_panel.extend(seqAf_ds)
            asc_panel.extend(seqEu_ds)
            asc_panel.extend(seqAs_ds)


            #print asc_panel
            Tasc_panel=zip(*asc_panel)
            print 'number of sites in Tasc_panel:',len(Tasc_panel)
            print 'number of chromosomes in Tasc_panel:',len(Tasc_panel[0])
            #print Tasc_panel


            #######Array with the available sites given the frequency cut off
            ##array with the frequency of all the simulated snps
            sites_freq=[]
            ##array with the available sites, that pass the frequency cut-off
            avail_sites=[] ##this one has the positions of the snps
            index_avail_sites=[] ##this one has the indexes of the snps

            for n in xrange(len(Tasc_panel)):
                freq_site=float(Tasc_panel[n][0:len(asc_panel)].count('1'))/float(len(asc_panel))
                if freq_site>=daf and freq_site<=1-daf:
                    sites_freq.append(freq_site)
                    avail_sites.append(pos[n])
                    index_avail_sites.append(n)

            #print sites_freq
            #print 'nb avail sites', len(avail_sites)
            #print index_avail_sites

            nb_avail_sites=len(avail_sites)
            #print 'nb avail and seg sites', nb_avail_sites
            #print 'index_avail_sites', index_avail_sites


            if (nb_avail_sites>=len(snps)):
                flag_sim=True

            else:
                flag_sim=False
                rep=rep+1


        #######################
        #####Calculate summary statistics from the regions for the CGI data
        print 'calculating summary stats for CGI data'
        elapsed_time=time.time()-start_time
        print '***********'+str(elapsed_time)+'***********'

        if nbss>0: ###no segregating sites in the simulations which is not possible

            Af_res=[]
            Af_res.extend(afs_stats.base_S_ss(seqAfCGI,nbss))
            pi_AfCGI=afs_stats.Pi2(Af_res[3],len(seqAfCGI))
            Af_res.append(afs_stats.Tajimas(pi_AfCGI,Af_res[0],len(seqAfCGI)))
            del(Af_res[3])
            res.extend(Af_res)
            head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'

            Eu_res=[]
            Eu_res.extend(afs_stats.base_S_ss(seqEuCGI,nbss))
            pi_EuCGI=afs_stats.Pi2(Eu_res[3],len(seqEuCGI))
            Eu_res.append(afs_stats.Tajimas(pi_EuCGI,Eu_res[0],len(seqEuCGI)))
            del(Eu_res[3])
            res.extend(Eu_res)
            head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'

            As_res=[]
            As_res.extend(afs_stats.base_S_ss(seqAsCGI,nbss))
            pi_AsCGI=afs_stats.Pi2(As_res[3],len(seqAsCGI))
            As_res.append(afs_stats.Tajimas(pi_AsCGI,As_res[0],len(seqAsCGI)))
            del(As_res[3])
            res.extend(As_res)
            head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'


            ##fst between populations
            res.append(afs_stats.FST2(seqAfCGI,pi_AfCGI,naf_CGI,seqEuCGI,pi_EuCGI,neu_CGI))
            res.append(afs_stats.FST2(seqAfCGI,pi_AfCGI,naf_CGI,seqAsCGI,pi_AsCGI,nas_CGI))
            res.append(afs_stats.FST2(seqEuCGI,pi_EuCGI,neu_CGI,seqAsCGI,pi_AsCGI,nas_CGI))
            head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'
            #print 'len(res)', len(res)


        print 'Done calculating ss from chomosomes'
        elapsed_time=time.time()-start_time
        print '***********'+str(elapsed_time)+'***********'
        ############


        print 'Make pseudo array'

        if flag_sim==False:
            print 'the sims did not work'
            pos_asc=[]
            pos_asc=index_avail_sites
            nbss_asc=len(pos_asc)
            #print 'nbss_asc', nbss_asc

        if (len(avail_sites)==len(snps)):
            print "number of avail_sites is equal to the number of Array snps"
            pos_asc=[]
            pos_asc=index_avail_sites
            nbss_asc=len(pos_asc)
            #print pos_asc
            flag_nb_asc_snps=1


        elif (len(avail_sites)>len(snps)):
            print "number of avail_sites greater than number of Array snps"
            pos_asc=[None]*int(len(snps)) ##indexes of the SNPs that pass the frequency cut-off and position
            for i in xrange(len(snps)): #each snp on the snp array on a chromosome
                ## y is the position of the SNPs in the array
                y=snps[i]
                ##find the closest SNP in the array
                closestleft=find2(avail_sites,y)
                #print 'closestleft', closestleft, avail_sites[closestleft]

                if(i>0 and pos_asc[i-1]==closestleft and closestleft+1<len(avail_sites)): ##avoid duplicates
                    #print 'duplicate 1'
                    closestleft=closestleft+1 ##move one position to the right
                    #print 'new closestleft', closestleft, avail_sites[closestleft]
                    pos_asc[i]=closestleft

                elif(i>0 and pos_asc[i-1]>closestleft and pos_asc[i-1]+1<len(avail_sites)):
                    #print 'duplicate 2'
                    closestleft=pos_asc[i-1]+1
                    #print 'new closestleft', closestleft, avail_sites[closestleft]
                    pos_asc[i]=closestleft

                else:
                    pos_asc[i]=closestleft

                #print 'after checking for dupl', pos_asc[i]
                ###if I have duplicates at this point, it means that there were not anyt more snps to choose from
                ###closestleft+1 or pos_asc[i-1]+1 == len(avail_sites)

            #print 'before smoothing'
            #print pos_asc

            #####smoothing
            ##last index of the pos_asc
            i=len(pos_asc)-1
            #print 'last i', i

            ##check if there is another position that might work better
            for j in xrange(0,i):
                if(j==i-1 and pos_asc[j]+1<pos_asc[j+1] and pos_asc[j]<(len(avail_sites)-1) and (j+1)<len(avail_sites)):
                    d1=abs(snps[j]-avail_sites[pos_asc[j]])
                    d2=abs(snps[j]-avail_sites[pos_asc[j]+1])
                    if(d2<d1):
                        pos_asc[j]=pos_asc[j]+1

            #print 'after smoothing'
            #print pos_asc
            ##removes duplicates
            pos_asc=(list(set(pos_asc)))
            pos_asc.sort()
            #print 'check for duplicates ', pos_asc


            ##might need to put this at the end of this part again
            nbss_asc=len(pos_asc)
            #print 'number of ascertained SNPs:',nbss_asc
            #print 'position of ascertained SNPs:', pos_asc


            if (len(snps)==nbss_asc):
                flag_nb_asc_snps=1
                print 'nb of asc snps equal to nb array snps'

            if (len(snps)!=len(pos_asc)):
                flag_nb_asc_snps=0
                print 'nb of asc snps not equal to nb array snps'

                diff=int(len(snps)-len(pos_asc))
                #print 'diff', diff

                for m in xrange(1,diff+1):
                    #print 'm', m

                    pos_asc2=[]

                    #print 'avail_sites', avail_sites
                    #print 'nb_avail_sites', nb_avail_sites

                    pos_asc2=add_snps(avail_sites, nb_avail_sites, pos_asc, nbss_asc, nb_array_snps)
                    pos_asc=pos_asc2

                    #print 'new len of pos_asc', len(pos_asc)
                    #print pos_asc
                    nbss_asc=len(pos_asc)

                    if nbss_asc==len(snps):
                        flag_nb_asc_snps=1
                        break

                    else:
                        flag_nb_asc_snps=0

            if (flag_nb_asc_snps==0): ##it means that the 1st index in pos_asc is 0; and the last is len(avail_sites)-1

                #print 'asc snps still missing'

                diff=int(len(snps)-len(pos_asc))
                #print 'diff', diff

                while (len(pos_asc)!=len(snps)):
                    rand_numb=randint(0,len(avail_sites)-1)
                    #print 'random',rand_numb
                    if rand_numb not in pos_asc:
                        pos_asc.append(rand_numb)

                #print 'pos_asc', pos_asc
                pos_asc.sort()
                #print 'pos_asc', pos_asc
                nbss_asc=len(pos_asc)

        print 'finished making pseudo array'
        elapsed_time=time.time()-start_time
        print '***********'+str(elapsed_time)+'***********'

        ############
        ##Transpose the data
        TseqAf=zip(*seqAfCGI)
        #print 'len(TseqAf)', len(TseqAf)
        TseqEu=zip(*seqEuCGI)
        #print 'len(TseqEu)', len(TseqEu)
        TseqAs=zip(*seqAsCGI)
        #print 'len(TseqAs)', len(TseqAs)

        seqJ=Talleles[total_naf+total_neu+total_nas:total_naf+total_neu+total_nas+nJ]
        #print 'len(seqJ)', len(seqJ)
        TseqJ=zip(*seqJ)
        #print 'len(TseqJ)', len(TseqJ)

        seqM=Talleles[total_naf+total_neu+total_nas+nJ:total_naf+total_neu+total_nas+nJ+nM]
        #print 'len(seqM)', len(seqM)
        TseqM=zip(*seqM)
        #print 'len(TseqM)', len(TseqM)

        seqA=Talleles[total_naf+total_neu+total_nas+nJ+nM:total_naf+total_neu+total_nas+nJ+nM+nA]
        #print 'len(seqA)', len(seqA)
        TseqA=zip(*seqA)
        #print 'len(TseqA)', len(TseqA)



        ###get genotypes for pseudo array
        allelesAf_asc=[]
        allelesEu_asc=[]
        allelesAs_asc=[]

        allelesJ_asc=[]
        allelesM_asc=[]
        allelesA_asc=[]


        ##get the ascertained SNPs in the populations of interest
        if (nbss_asc==len(index_avail_sites)):
            for x in xrange(nbss_asc):
                allelesAf_asc.append(TseqAf[pos_asc[x]])
                allelesEu_asc.append(TseqEu[pos_asc[x]])
                allelesAs_asc.append(TseqAs[pos_asc[x]])
                allelesJ_asc.append(TseqJ[pos_asc[x]])
                allelesM_asc.append(TseqM[pos_asc[x]])
                allelesA_asc.append(TseqA[pos_asc[x]])

        elif (len(index_avail_sites)>nbss_asc):
            for x in xrange(len(pos_asc)):
                allelesAf_asc.append(TseqAf[index_avail_sites[pos_asc[x]]])
                allelesEu_asc.append(TseqEu[index_avail_sites[pos_asc[x]]])
                allelesAs_asc.append(TseqAs[index_avail_sites[pos_asc[x]]])
                allelesJ_asc.append(TseqJ[index_avail_sites[pos_asc[x]]])
                allelesM_asc.append(TseqM[index_avail_sites[pos_asc[x]]])
                allelesA_asc.append(TseqA[index_avail_sites[pos_asc[x]]])

        print 'Make ped and map files'
        elapsed_time=time.time()-start_time
        print '***********'+str(elapsed_time)+'***********'

        ##Make ped file
        filenameped=str(sim_data_dir)+'/macs_asc_'+str(job)+'_chr'+str(chr_number)+'.ped'
        fileped=open(filenameped,'w')
        ped=''
        for e in range(2,int(neu_CGI)+2,2):
            ped=ped+'E '+str(e/2)+'_E 0 0 1 -9 '
            for g in range(0,len(pos_asc)):
                ped=ped+str(int(allelesEu_asc[g-1][e-2])+1)+' '+str(int(allelesEu_asc[g-1][e-1])+1)
                if g<len(pos_asc)-1:
                    ped=ped+' '
            ped=ped+'\n'
        for j in range(2,int(nJ)+2,2):
            ped=ped+'J '+str(j/2)+'_J 0 0 1 -9 '
            for g in range(0,len(pos_asc)):
                ped=ped+str(int(allelesJ_asc[g-1][j-2])+1)+' '+str(int(allelesJ_asc[g-1][j-1])+1)
                if g<len(pos_asc)-1:
                    ped=ped+' '
            ped=ped+'\n'
        for m in range(2,int(nM)+2,2):
            ped=ped+'M '+str(m/2)+'_M 0 0 1 -9 '
            for g in range(0,len(pos_asc)):
                ped=ped+str(int(allelesM_asc[g-1][m-2])+1)+' '+str(int(allelesM_asc[g-1][m-1])+1)
                if g<len(pos_asc)-1:
                    ped=ped+' '
            ped=ped+'\n'
        for a in range(2,int(nA)+2,2):
            ped=ped+'A '+str(a/2)+'_A 0 0 1 -9 '
            for g in range(0,len(pos_asc)):
                ped=ped+str(int(allelesA_asc[g-1][a-2])+1)+' '+str(int(allelesA_asc[g-1][a-1])+1)
                if g<len(pos_asc)-1:
                    ped=ped+' '
            ped=ped+'\n'
        fileped.write(ped)
        fileped.close()


        ##Make map file
        filenamemap=str(sim_data_dir)+'/macs_asc_'+str(job)+'_chr'+str(chr_number)+'.map'
        filemap=open(filenamemap,'a')
        map=''
        for g in range(0,len(pos_asc)):
            map=str(chr_number)+' '+'chr'+str(chr_number)+'_'+str(pos_asc[g])+' '+str(int(avail_sites[pos_asc[g]]-1))+' '+str(int(avail_sites[pos_asc[g]]))+'\n'
            filemap.write(map)
        filemap.close()


        ###Genotypes for the ascertained SNPs
        seqAf_asc=zip(*allelesAf_asc)
        seqEu_asc=zip(*allelesEu_asc)
        seqAs_asc=zip(*allelesAs_asc)

        seqJ_asc=zip(*allelesJ_asc)
        seqM_asc=zip(*allelesM_asc)
        seqA_asc=zip(*allelesA_asc)


        ########Use Germline to find IBD on pseduo array ped and map files
        run_germline = int(arguments[6])
        filenameout = str(germline_out_dir) + '/macs_asc_' + str(job) + '_chr' + str(chr_number)

        print run_germline

        if (run_germline == 0):
            print 'Running Germline on '+str(filenameped)+' '+str(filenamemap)
            elapsed_time=time.time()-start_time
            print '***********'+str(elapsed_time)+'***********'

            ###Germline seems to be outputting in the wrong unit - so I am putting the min at 3000000 so that it is 3Mb, but should be the default.
            print 'bash ./bin/phasing_pipeline/gline.sh ./bin/germline-1-5-1/germline  '+str(filenameped)+' '+str(filenamemap)+' '+str(filenameout)+' "-bits 10 -min_m 3000000"'
            germline=Popen.wait(Popen('bash ./bin/phasing_pipeline/gline.sh ./bin/germline-1-5-1/germline  '+str(filenameped)+' '+str(filenamemap)+' '+str(filenameout)+' "-bits 10 -min_m 3000000"',shell=True))

            print 'finished running germline'
            elapsed_time=time.time()-start_time
            print '***********'+str(elapsed_time)+'***********'


        ########Get IBD stats from Germline output
        if os.path.isfile(str(filenameout) + '.match'):
#			rmped=Popen.wait(Popen('rm '+str(filenameped),shell=True))
#			rmmap=Popen.wait(Popen('rm '+str(filenamemap),shell=True))
#			rmlog = Popen.wait(Popen('rm ' + str(filenameout) + '.log', shell=True))

            print 'reading Germline IBD output'
            filegermline = open(str(filenameout) + '.match', 'r')
            IBDlengths_AA=[]
            IBDlengths_JJ=[]
            IBDlengths_MM=[]
            IBDlengths_EE=[]
            IBDlengths_AE=[]
            IBDlengths_AJ=[]
            IBDlengths_AM=[]
            IBDlengths_JM=[]
            IBDlengths_JE=[]
            IBDlengths_ME=[]
            for line in filegermline:
                pop1 = line.split()[0]
                pop2 = line.split()[2]
                segment = float(line.split()[10])/1000000
                pair = str(pop1)+'_'+str(pop2)
                if pair=='A_A':
                    IBDlengths_AA.append(segment)
                if pair=='J_J':
                    IBDlengths_JJ.append(segment)
                if pair=='M_M':
                    IBDlengths_MM.append(segment)
                if pair=='E_E':
                    IBDlengths_EE.append(segment)
                if pair=='A_E' or pair=='E_A':
                    IBDlengths_AE.append(segment)
                if pair=='A_J' or pair=='J_A':
                    IBDlengths_AJ.append(segment)
                if pair=='A_M' or pair=='M_A':
                    IBDlengths_AM.append(segment)
                if pair=='J_M' or pair=='M_J':
                    IBDlengths_JM.append(segment)
                if pair=='J_E' or pair=='E_J':
                    IBDlengths_JE.append(segment)
                if pair=='M_E' or pair=='E_M':
                    IBDlengths_ME.append(segment)
            filegermline.close()
#			rmmatch=Popen.wait(Popen('rm '+str(filenameout)+'.match',shell=True))

            elapsed_time=time.time()-start_time
            print '***********'+str(elapsed_time)+'***********'
            print 'calculating summary stats'

            IBDlengths_mean=[]
            IBDlengths_median=[]
            IBDlengths_num=[]
            IBDlengths_var=[]
            IBDlengths_mean30 = []
            IBDlengths_median30 = []
            IBDlengths_num30 = []
            IBDlengths_var30 = []

            pairs=[IBDlengths_AA,IBDlengths_JJ,IBDlengths_MM,IBDlengths_EE,IBDlengths_AE,IBDlengths_AJ,IBDlengths_AM,IBDlengths_JM,IBDlengths_JE,IBDlengths_ME]
            for p in pairs:
                IBDlengths_num.append(len(p))
                if len(p)<1:
                    p.append(0)
                IBDlengths_mean.append(np.mean(p))
                IBDlengths_median.append(np.median(p))
                IBDlengths_var.append(np.var(p))
                #### Get IBD greater than 30 Mb
                for l in p:
                    n=0
                    IBDlengths30=[]
                    if l>30:
                        n=n+1
                        IBDlengths30.append(l)
                    else:
                        IBDlengths30.append(0)
                IBDlengths_mean30.append(np.mean(IBDlengths30))
                IBDlengths_median30.append(np.median(IBDlengths30))
                IBDlengths_var30.append(np.var(IBDlengths30))
                IBDlengths_num30.append(n)


            res.extend(IBDlengths_mean)
            head = head + 'IBD_mean_AA\tIBD_mean_JJ\tIBD_mean_MM\tIBD_mean_EE\tIBD_mean_AE\tIBD_mean_AJ\tIBD_mean_AM\tIBD_mean_JM\tIBD_mean_JE\tIBD_mean_ME\t'
            res.extend(IBDlengths_median)
            head = head + 'IBD_median_AA\tIBD_median_JJ\tIBD_median_MM\tIBD_median_EE\tIBD_median_AE\tIBD_median_AJ\tIBD_median_AM\tIBD_median_JM\tIBD_median_JE\tIBD_median_ME\t'
            res.extend(IBDlengths_num)
            head = head + 'IBD_num_AA\tIBD_num_JJ\tIBD_num_MM\tIBD_num_EE\tIBD_num_AE\tIBD_num_AJ\tIBD_num_AM\tIBD_num_JM\tIBD_num_JE\tIBD_num_ME\t'
            res.extend(IBDlengths_var)
            head = head + 'IBD_var_AA\tIBD_var_JJ\tIBD_var_MM\tIBD_var_EE\tIBD_var_AE\tIBD_var_AJ\tIBD_var_AM\tIBD_var_JM\tIBD_var_JE\tIBD_var_ME\t'

            res.extend(IBDlengths_mean30)
            head = head + 'IBD30_mean_AA\tIBD30_mean_JJ\tIBD30_mean_MM\tIBD30_mean_EE\tIBD30_mean_AE\tIBD30_mean_AJ\tIBD30_mean_AM\tIBD30_mean_JM\tIBD30_mean_JE\tIBD30_mean_ME\t'
            res.extend(IBDlengths_median30)
            head = head + 'IBD30_median_AA\tIBD30_median_JJ\tIBD30_median_MM\tIBD30_median_EE\tIBD30_median_AE\tIBD30_median_AJ\tIBD30_median_AM\tIBD30_median_JM\tIBD30_median_JE\tIBD30_median_ME\t'
            res.extend(IBDlengths_num30)
            head = head + 'IBD30_num_AA\tIBD30_num_JJ\tIBD30_num_MM\tIBD_num_EE\tIBD30_num_AE\tIBD30_num_AJ\tIBD30_num_AM\tIBD30_num_JM\tIBD30_num_JE\tIBD30_num_ME\t'
            res.extend(IBDlengths_var30)
            head = head + 'IBD30_var_AA\tIBD30_var_JJ\tIBD30_var_MM\tIBD_var_EE\tIBD30_var_AE\tIBD30_var_AJ\tIBD30_var_AM\tIBD30_var_JM\tIBD30_var_JE\tIBD30_var_ME\t'

        #######
        #########calculate summary stats from the ascertained SNPs
        if nbss_asc>0:

            Af_asc=[]
            ss_Af_asc=afs_stats.base_S_ss(seqAf_asc,nbss_asc)
            if (ss_Af_asc[0]==0):
                for i in xrange(5):
                    Af_asc.append(0)
                pi_Af_asc=0
            else:
                Af_asc.extend(afs_stats.base_S_ss(seqAf_asc,nbss_asc))
                pi_Af_asc=afs_stats.Pi2(Af_asc[3],len(seqAf_asc))
                Af_asc.append(pi_Af_asc)
                Af_asc.append(afs_stats.Tajimas(pi_Af_asc,Af_asc[0],len(seqAf_asc)))
                del(Af_asc[3])

            res.extend(Af_asc)
            head = head + 'SegS_Af_ASC\tSing_Af_ASC\tDupl_Af_ASC\tPi_Af_ASC\tTajD_Af_ASC\t'

            ############

            Eu_asc=[]
            ss_Eu_asc=afs_stats.base_S_ss(seqEu_asc,nbss_asc)
            if (ss_Eu_asc[0]==0):
                for i in xrange(5):
                    Eu_asc.append(0)
                pi_Eu_asc=0
            else:
                Eu_asc.extend(afs_stats.base_S_ss(seqEu_asc,nbss_asc))
                pi_Eu_asc=afs_stats.Pi2(Eu_asc[3],len(seqEu_asc))
                Eu_asc.append(pi_Eu_asc)
                Eu_asc.append(afs_stats.Tajimas(pi_Eu_asc,Eu_asc[0],len(seqEu_asc)))
                del(Eu_asc[3])

            res.extend(Eu_asc)
            head = head + 'SegS_Eu_ASC\tSing_Eu_ASC\tDupl_Eu_ASC\tPi_Eu_ASC\tTajD_Eu_ASC\t'
            ###########

            As_asc=[]
            ss_As_asc=afs_stats.base_S_ss(seqAs_asc,nbss_asc)
            if (ss_As_asc[0]==0):
                for i in xrange(5):
                    As_asc.append(0)
                pi_As_asc=0
            else:
                As_asc.extend(afs_stats.base_S_ss(seqAs_asc,nbss_asc))
                pi_As_asc=afs_stats.Pi2(As_asc[3],len(seqAs_asc))
                As_asc.append(pi_As_asc)
                As_asc.append(afs_stats.Tajimas(pi_As_asc,As_asc[0],len(seqAs_asc)))
                del(As_asc[3])

            res.extend(As_asc)
            head = head + 'SegS_As_ASC\tSing_As_ASC\tDupl_As_ASC\tPi_As_ASC\tTajD_As_ASC\t'
            ############

            J_asc=[]
            ss_J_asc=afs_stats.base_S_ss(seqJ_asc,nbss_asc)
            if (ss_J_asc[0]==0):
                for i in xrange(5):
                    J_asc.append(0)
                pi_J_asc=0
            else:
                J_asc.extend(afs_stats.base_S_ss(seqJ_asc,nbss_asc))
                pi_J_asc=afs_stats.Pi2(J_asc[3],len(seqJ_asc))
                J_asc.append(pi_J_asc)
                J_asc.append(afs_stats.Tajimas(pi_J_asc,J_asc[0],len(seqJ_asc)))
                del(J_asc[3])

            res.extend(J_asc)
            head = head + 'SegS_J_ASC\tSing_J_ASC\tDupl_J_ASC\tPi_J_ASC\tTajD_J_ASC\t'
            #############

            M_asc=[]
            ss_M_asc=afs_stats.base_S_ss(seqM_asc,nbss_asc)
            if (ss_M_asc[0]==0):
                for i in xrange(5):
                    M_asc.append(0)
                pi_M_asc=0
            else:
                M_asc.extend(afs_stats.base_S_ss(seqM_asc,nbss_asc))
                pi_M_asc=afs_stats.Pi2(M_asc[3],len(seqM_asc))
                M_asc.append(pi_M_asc)
                M_asc.append(afs_stats.Tajimas(pi_M_asc,M_asc[0],len(seqM_asc)))
                del(M_asc[3])

            res.extend(M_asc)
            head = head + 'SegS_M_ASC\tSing_M_ASC\tDupl_M_ASC\tPi_M_ASC\tTajD_M_ASC\t'
            #############

            A_asc=[]
            ss_A_asc=afs_stats.base_S_ss(seqA_asc,nbss_asc)
            if (ss_A_asc[0]==0):
                for i in xrange(5):
                    A_asc.append(0)
                pi_A_asc=0
            else:
                A_asc.extend(afs_stats.base_S_ss(seqA_asc,nbss_asc))
                pi_A_asc=afs_stats.Pi2(A_asc[3],len(seqA_asc))
                A_asc.append(pi_A_asc)
                A_asc.append(afs_stats.Tajimas(pi_A_asc,A_asc[0],len(seqA_asc)))
                del(A_asc[3])

            res.extend(A_asc)
            head = head + 'SegS_A_ASC\tSing_A_ASC\tDupl_A_ASC\tPi_A_ASC\tTajD_A_ASC\t'
            #############


            ##fst between populations
            res.append(afs_stats.FST2(seqAf_asc,pi_Af_asc,naf_CGI,seqEu_asc,pi_Eu_asc,neu_CGI))
            res.append(afs_stats.FST2(seqAf_asc,pi_Af_asc,naf_CGI,seqAs_asc,pi_As_asc,nas_CGI))
            res.append(afs_stats.FST2(seqEu_asc,pi_Eu_asc,neu_CGI,seqAs_asc,pi_As_asc,nas_CGI))
            head = head + 'FST_AfEu_ASC\tFST_AfAs_ASC_m\tFST_EuAs_ASC\t'


            res.append(afs_stats.FST2(seqA_asc,pi_A_asc,nA,seqEu_asc,pi_Eu_asc,neu_CGI))
            res.append(afs_stats.FST2(seqA_asc,pi_A_asc,nA,seqJ_asc,pi_J_asc,nJ))
            res.append(afs_stats.FST2(seqA_asc,pi_A_asc,nA,seqM_asc,pi_M_asc,nM))
            res.append(afs_stats.FST2(seqM_asc,pi_M_asc,nM,seqJ_asc,pi_J_asc,nJ))
            head = head + 'FST_AEu_ASC\tFST_AJ_ASC\tFST_AM_ASC\tFST_MJ_ASC\n'

            cont=cont+1


    print 'finished calculating ss'
    elapsed_time=time.time()-start_time
    print '***********'+str(elapsed_time)+'***********'

    ################
    #####write parameter values to file

    param_file=str(sim_values_dir)+'/sim_'+str(job)+'_values.txt'
    fileoutparam=open(param_file,'w')

    ##write parameter values
    head_param='Asc_NAF\tAsc_NEU\tAsc_NCHB\tdaf\tLog10_NAF\tLog10_NANC\tLog10_NCEU\tLog10_NCHB\tLog10_NA\tLog10_NJ\tLog10_NM\trA\trMJ\tm\tTgrowth_Af\tTAF\tTEM\tTeu_as\tTA\tTMJ\tTm\n'
    fileoutparam.write(head_param)



    for z in range(len(para_out)):
        if z==(len(para_out)-1):
            fileoutparam.write("%s\n" % para_out[z])
        else:
            fileoutparam.write("%s\t" % para_out[z])

    fileoutparam.close()

    #####
    #####

    filesummary=str(results_sims_dir)+'/ms_output_'+str(job)+'.summary'
    filesumm=open(filesummary,'w')
    filesumm.write(head)

    out=''

    for g in range(len(res)):
        out=out+str(res[g])+'\t'
    out=out[:-1]+'\n'

    filesumm.write(out)
    filesumm.close()


if __name__ == '__main__':
    main(argv)
    #main()
