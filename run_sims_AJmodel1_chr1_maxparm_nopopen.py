import time
import subprocess 
from subprocess import Popen,PIPE,call
import os
import string
from string import join
import math
import random
from random import randint
from random import uniform
from random import gauss
from random import gammavariate
from random import betavariate
from math import sqrt
import sys
from sys import argv
import datetime
import macsSwig
from bisect import bisect_left
from bisect import bisect_right
import numpy as np
import re
#from memory_profiler import LogFile,profile
#import cProfile
#import re
import psutil


###summary statistics####################################

def hamming_distance(s1, s2): 
    #Hamming distance between two strings of equal length is the number of positions at which the corresponding symbols are different
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def base_S_ss(seq,nbsites):

    #print 'base_S_ss'
    #print 'nbsites', nbsites

    spec_zero=[]
    for g in range(len(seq)-1):
        spec_zero.append(0)

    var_ss=0 #Segregating sites
    #het_sum=0.0

    alleles=zip(*seq)
    for g in range(nbsites):
        #print 'g', g
        if 0<list(alleles[g]).count('1')<(list(alleles[g]).count('1')+list(alleles[g]).count('0')): ##this ignores sites that have all zeros, or all ones
            var_ss=var_ss+1
            #het_sum=het_sum+((float(list(alleles[g]).count('1'))/float((list(alleles[g]).count('1')+list(alleles[g]).count('0'))))**2)
            spec_zero[list(alleles[g]).count('1')-1]=spec_zero[list(alleles[g]).count('1')-1]+1

    if var_ss>0:
        #het=1.0-((1/float(var_ss))*het_sum)
        Ns=spec_zero[0]+spec_zero[-1] ##number of singletons
        Nd=spec_zero[1]+spec_zero[-2] ##number of dupletons
    else:
        #het='NA'
        Ns=0
        Nd=0


    ##Nds=spec_zero[0] ##number of derived singletons

    return [var_ss,Ns,Nd,spec_zero]

def base_h_ss(seq):
    #print seq
    #A set is an unordered collection with no duplicate elements
    haps=list(set(seq)) #makes a list of the different haplotypes found in seq
    #print haps
    cnt=[]
    hi=0 #hi is the mode of the haplotypes
    for g in range(len(haps)):
        cnt.append(seq.count(haps[g])) #count() returns the number of occurrences of substring sub (in this case, haplotypes) in seq
        #print cnt[g]
        if cnt[g]>hi: #you are keeping the highest occurrence of all the haplotypes, which the mode
            hi=cnt[g]


    #The average pairwise difference within a population can be calculated as the sum of the pairwise differences divided by the number of pairs.
    pwd_dist=0 #parwise differences between haplotypes??
    for g in range(len(cnt)-1): #count will have the same size as haps
        m=g+1
        while m<len(cnt):
            pwd_dist=pwd_dist+(hamming_distance(haps[g], haps[m])*(cnt[g]*cnt[m]))
            m=m+1

    #p is the total sum of the pairwise differences
    #p=pwd_dist*(2.0/(float(len(seq))*(float(len(seq))-1)))
    #print p

    #this function returns the number of different haplotypes and the mode
    return [len(haps),hi]


def base_h_ss_nosing(seq): #calculates the number of different haplotypes and the mode, but ignoring singleton haplotypes

    d={}
    d=dict((i,seq.count(i)) for i in seq)

    d2=[] #d2 has the list of the haplotypes that occur more than once in the population

    for k in d.keys():
        #print k  #this gives you the values/haplotypes of the dictionary
        #print d[k]
        if d[k]!=1:
            #print k
            d2.append(k)

    cnt=[]
    hi=0 #hi is the mode of the haplotypes
    for g in range(len(d2)):
        cnt.append(seq.count(d2[g])) #count() returns the number of occurrences of substring sub (in this case, haplotypes) in seq
        #print cnt[g]
        if cnt[g]>hi: #you are keeping the highest occurrence of all the haplotypes, which the mode
            hi=cnt[g]

    #this function returns the number of different haplotypes and the mode
    return [len(d2),hi]


def pri_sha_h_nosing(seqs1,seqs2):

    d1={}
    d1=dict((i,seqs1.count(i)) for i in seqs1)

    d1_2=[] #d2 has the list of the haplotypes that occur more than once in the population

    for k in d1.keys():
        if d1[k]!=1:
            d1_2.append(k)

    d2={}
    d2=dict((i,seqs2.count(i)) for i in seqs2)

    d2_2=[] #d2 has the list of the haplotypes that occur more than once in the population

    for l in d2.keys():
        if d2[l]!=1:
            d2_2.append(l)

    priA=0
    priB=0
    sha=0

    #this puts the two sequences together in the same array
    seqs=d1_2[:]
    seqs.extend(d2_2)
    haps=list(set(seqs)) #makes a list of the different haplotypes found in all populations

    for g in range(len(haps)):
        pop1_cnt=d1_2.count(haps[g]) #count() returns the number of occurrences of substring sub (in this case, haplotypes) in seqs1, which population 1
        pop2_cnt=d2_2.count(haps[g])

        if pop1_cnt>0 and pop2_cnt>0:
            sha=sha+1
        elif pop1_cnt>0 and pop2_cnt==0:
            priA=priA+1
        elif pop1_cnt==0 and pop2_cnt>0:
            priB=priB+1

    return [sha,priA,priB]


def foldedAFS(array):

    delta=0

    n=len(array)
    #print n

    unfolded=[]
    for g in range((n/2)+1):
        unfolded.append(0)

    #print len(unfolded)

    for g in xrange((n/2)+1):

        #print g

        if g != (n-1-g):
            delta=0
        else:
            delta=1

        #print 'array[g]',array[g]
        #print 'array[n]',array[(n-1)]
        #print 'array[n-g]', array[(n-1)-g]

        #print 'delta', delta
        #print g
        #print (n-1)-g

        unfolded[g]=(array[g]+array[(n-1)-g])/(1+delta)
        #print 'unfolded[g]', unfolded[g]

    return unfolded

def Pi2(spec,n): #standard pi, n = sample size (in chromosomes)
    theta_pi=0.0

    for g in range(len(spec)):
        theta_pi=theta_pi+(2.0*float(spec[g])*(g+1.0)*(n-(g+1.0)))/(n*(n-1.0))

    return theta_pi


def Tajimas(p,S,n):
    ### pi, number of segregating sites, and number of chromosomes

    if (S==0):
        #return 'NA'
        return 0

    else:

        a1=0.0
        for g in range(n-1):
            a1=a1+(1.0/(g+1.0))

        #print a1
        a2=0.0
        for g in range(n-1):
            a2=a2+(1.0/((g+1.0)**2))

        b1=(n+1.0)/(3.0*(n-1.0))

        b2=(2*((n**2.0)+n+3))/((9*n)*(n-1))
        c1=b1-(1.0/a1)
        #print 'c1', c1
        c2=b2-((n+2.0)/(a1*n))+(a2/(a1**2.0))
        e1=c1/a1
        #print 'e1', e1
        e2=c2/((a1**2.0)+a2)
        #print 'e2', e2
        TajD=(p-(S/a1))/(sqrt((e1*S)+((e2*S)*(S-1.0))))

        return TajD


def FST(geno1,geno2,nbsites):
    r=2.0
    n1=float(len(geno1))#how many individuals
    n2=float(len(geno2))
    n_bar=(n1/r)+(n2/r)#average sample size
    nc=((r*n_bar)-(((n1**2)/(r*n_bar))+((n2**2)/(r*n_bar))))/(r-1.0)

    a_sum=0.0
    abc_sum=0.0

    nsit_us=0
    for g in range(nbsites):
        #print g
        aa1=0 #how many homozygotes 0 in 1
        ab1=0 #how many heterozygotes in 1
        bb1=0 #how many homozygotes 1 in 1
        aa2=0
        ab2=0
        bb2=0
        for n in range(len(geno1)):
            ball=geno1[n][0][g].count('1')+geno1[n][1][g].count('1')#to get the information from geno, count how many ones in one pair
            if ball==0:
                aa1=aa1+1
            elif ball==1:
                ab1=ab1+1
            elif ball==2:
                bb1=bb1+1
        for n in range(len(geno2)):
            ball=geno2[n][0][g].count('1')+geno2[n][1][g].count('1')
            if ball==0:
                aa2=aa2+1
            elif ball==1:
                ab2=ab2+1
            elif ball==2:
                bb2=bb2+1

        if 0<aa1+aa2<(len(geno1)+len(geno2)):
            #print 'if'

            p1=float((bb1*2.0)+ab1)/(n1*2.0)#get frequency of the derived allele in the two populations
            #print 'p1'
            #print p1
            p2=float((bb2*2.0)+ab2)/(n2*2.0)
            #print 'p2'
            #print p2
            p_bar=((n1*p1)/(r*n_bar))+((n2*p2)/(r*n_bar)) #average allele frequency for that site
            #print 'p_bar'
            #print p_bar
            s_sq=(n1*((p1-p_bar)**2.0))/((r-1.0)*n_bar)+(n2*((p2-p_bar)**2.0))/((r-1.0)*n_bar)
            #print 's_sq'
            #print s_sq
            h1=float(ab1)/n1 #frequency of the heterozygotes in population 1
            #print 'h1'
            #print h1
            h2=float(ab2)/n2
            ##print 'h2'
            #print h2
            h_bar=((n1*h1)/(r*n_bar))+((n2*h2)/(r*n_bar))
            #print 'h_bar'
            #print h_bar

            a=(n_bar/nc)*((s_sq)-(1.0/(n_bar-1.0))*((p_bar*(1.0-p_bar))-(((r-1.0)/r)*s_sq)-(h_bar/4.0)))
            #print 'a'
            #print a
            b=(n_bar/(n_bar-1))*(((p_bar)*(1.0-p_bar))-(((r-1.0)/r)*s_sq)-((((2.0*n_bar)-1)/(4*n_bar))*h_bar))
            #print 'b'
            #print b
            c=(1.0/2.0)*h_bar
            #print 'c'
            #print c

            a_sum=a_sum+a
            abc_sum=abc_sum+(a+b+c)
            #print 'abc_sum'
            #print abc_sum

            nsit_us=nsit_us+1
        else:
            nsit_us=nsit_us

    if abc_sum==0.0:
        theta='NA'
    else:
        theta=a_sum/abc_sum

    return theta

def Pi(seq1,nseq1):
    k1=0
    for i in xrange(0,nseq1):
        for j in xrange(i+1,nseq1):

            k1=k1+hamming_distance(seq1[i],seq1[j])

    p1=(2/(float(nseq1)*(float(nseq1)-1)))*k1

    return p1


def FST2(seq1,pi1,nseq1,seq2,pi2,nseq2): ###FST based on pi within populations and between populations
###number of chromosomes

    k3=0

    ##Pi within populations
    pw=(pi1+pi2)/2
    #print 'pw', pw

    for i in xrange(len(seq1)):
        for j in xrange(len(seq2)):

            k3=k3+hamming_distance(seq1[i],seq2[j])

    pb=k3/(float(nseq1)*float(nseq2))
    #print 'pb', pb

    if (pb==0):
        #return 'NA'
        return '0'

    else:
        fst=float(1-(pw/pb))
        return fst


def pri_sha_h(seqs1,seqs2):
#how many haplotypes are shared and how many are private to the populations	

    priA=0
    priB=0
    sha=0

    #this puts the two sequences together in the same array
    seqs=seqs1[:]
    seqs.extend(seqs2)
    haps=list(set(seqs)) #makes a list of the different haplotypes found in all populations

    for g in range(len(haps)):
        pop1_cnt=seqs1.count(haps[g]) #count() returns the number of occurrences of substring sub (in this case, haplotypes) in seqs1, which population 1
        pop2_cnt=seqs2.count(haps[g])

        #print pop1_cnt
        #print pop2_cnt

        if pop1_cnt>0 and pop2_cnt>0:
            sha=sha+1
        elif pop1_cnt>0 and pop2_cnt==0:
            priA=priA+1
        elif pop1_cnt==0 and pop2_cnt>0:
            priB=priB+1

    return [sha,priA,priB]


#find simulated snps closest to snps from chip in pop A with maf>0.05 in n indivs
def find(a, x, alleles, na_samp, daf, sim):
    if x <= a[0]: #x is smaller than first value in a
        k=0
        while float(alleles[k][0:na_samp].count('1')/float(na_samp))<daf or float(alleles[k][0:na_samp].count('1')/float(na_samp))>1-daf:
            k=k+1
        return k
    if x >= a[len(a)-1]: #x is larger than last value in a
        j=len(a)-1
        while float(alleles[j][0:na_samp].count('1')/float(na_samp))<daf or float(alleles[j][0:na_samp].count('1')/float(na_samp))>1-daf:
            j=j-1
        return j
    #Find leftmost item greater than or equal to x
    else:
        i = bisect_left(a, x)
        k=i
        j=i
        #walk up simulated sites
        while k<int(len(a)-1) and float(alleles[k][0:na_samp].count('1')/float(na_samp))<daf or float(alleles[k][0:na_samp].count('1')/float(na_samp))>1-daf:
            k=k+1
        #walk down simulated sites
        while int(j)>0 and float(alleles[j][0:na_samp].count('1')/float(na_samp))<daf or float(alleles[j][0:na_samp].count('1')/float(na_samp))>1-daf:
            j=j-1
        if float(alleles[j][0:na_samp].count('1')/float(na_samp))<daf or float(alleles[j][0:na_samp].count('1')/float(na_samp))>1-daf:
            if float(alleles[k][0:na_samp].count('1')/float(na_samp))<daf or float(alleles[k][0:na_samp].count('1')/float(na_samp))>1-daf: #no simulated sites meet daf cutoff
                print 'poopjk!'
            else: #all simulated sites to the left of i do not meet daf cutoff
                return k
        elif float(alleles[k][0:na_samp].count('1')/float(na_samp))<daf or float(alleles[k][0:na_samp].count('1')/float(na_samp))>1-daf: #all simulated sites to the right of i do not meet daf cutoff
            return j
        elif int(a[k]-x)<int(x-a[j]): #the position greater than x is closer
            return k
        else: #the position less than x is closer
            return j


##This function will receive the array with available sites (sites that passed the frequency cut-off)
def find2(a, x):
    #print a

    if x <= a[0]: #x is smaller than first value in a
        #print a[0]
        return 0
    elif x >= a[len(a)-1]: #x is larger than last value in a
        #print a[len(a)-1]
        return len(a)-1
    #Find leftmost item greater than or equal to x
    else:
        #print 'entro al ese, bisect'
        i = bisect_right(a, x)
        #print 'i', i
        #print a[i]
        j = i-1
        #print 'j', j
        #print a[j]
        d1=abs(a[i]-x)
        d2=abs(a[j]-x)
        if d2<d1:
            #print 'return j'
            return j
        elif d1<d2:
            #print 'return i'
            return i
        elif d1==d2:
            return i


        #print a[i-1]
        #return i-1

def add_snps(avail_sites, nb_avail_sites, pos_asc, nbss_asc, nb_array_snps):

    first_index=pos_asc[0]
    last_index=pos_asc[-1]

    if(nb_avail_sites>nbss_asc): ###this should happen all the time

        if (last_index<nb_avail_sites-1):

            try:
                avail_sites[last_index+1]
            except:
                print "well, it WASN'T defined after all!"

                try:
                    avail_sites[first_index-1]
                except:
                    print "well, it WASN'T defined after all!"
                else:
                    for i in xrange(len(pos_asc)):
                        pos_asc[i]=pos_asc[i]-1

            else:
                print "sure, it was defined"
                pos_asc.append(last_index+1)

        elif (last_index==nb_avail_sites-1):

            if(first_index-1)>=0:
                pos_asc.insert(0,first_index-1)


    return pos_asc


def param_sim_asc(): ##get parameter values from the priors

    #print "choosing param values"

    para_out=[]
    parameters={}

    ###Discovery panel
    low=2
    high=20

    asc_nb_af=high
    asc_nb_eu=high
    asc_nb_as=high

    para_out.extend([asc_nb_af])
    para_out.extend([asc_nb_eu])
    para_out.extend([asc_nb_as])

    daf=random.uniform(0.05,0.10)
    para_out.extend([daf])


    ####Demographic model
    #population size in Africa
    NAF=float(round(10**5))
    para_out.extend([math.log10(NAF)])
    parameters['NAF']=NAF
    #print NAF

    #Ancestral population size, before population growth in AF
    #choose ancestral Ne based on being some value smaller than NAF between 1 or 0.1 x.
    #population growth
    Nrat_High=0.0   #Allow only growth for now
    Nrat_Low=-1.0
    NANC=float(round(10**Nrat_High*NAF))
    para_out.extend([math.log10(NANC)])
    parameters['NANC']=NANC

    NCEU=float(round(10**5.0))
    para_out.extend([math.log10(NCEU)])
    parameters['NCEU']=NCEU

    NCHB=float(round(10**5.0))
    para_out.extend([math.log10(NCHB)])
    parameters['NCHB']=NCHB

    #Ancestral population size of Eurasians
    #NEU_AS=float(randint(1500,5000))
    #para_out.extend([NEU_AS])
    #parameters['NEU_AS']=NEU_AS

    #Population size of AJ
    NA=float(round(10**6.7))
    para_out.extend([math.log10(NA)])
    parameters['NA']=NA

    #Population size of Jews
    NJ=float(round(10**6.0))
    para_out.extend([math.log10(NJ)])
    parameters['NJ']=NJ

    #Population size of Middle Easterns
    NM=float(round(10**6.0))
    para_out.extend([math.log10(NM)])
    parameters['NM']=NM

    #Growth rate in AJ
    rA=float(round(10**random.uniform(0,1))) #should it be -1?
    para_out.extend([math.log10(rA)])
    parameters['rA']=rA

    #Growth rate in Jews and Middle East
    rMJ=float(round(10**random.uniform(0,1))) #should it be -1?
    para_out.extend([math.log10(rMJ)])
    parameters['rMJ']=rMJ

    #migration rate from Europe to AJ
    m_High=1
    m_Low=0
    m=float(randint(m_Low,m_High))
    para_out.extend([m])
    parameters['m']=m

    #Time of the instantaneous growth in Africa, before the split between Africans and non Africans
    Tgrowth_Low=1
    Tgrowth_High=4100
    Tgrowth_Af=float(randint(Tgrowth_Low,Tgrowth_High))
    parameters['Tgrowth_Af']=Tgrowth_Af
    para_out.extend([Tgrowth_Af])

    #Time of split between YRI and CEU/CHB
    Taf_High=4100				  #102,500 years using 25 years per generation
    Taf_Low=1600				  #40,000 years using 25 years per generation.
    Taf=float(Taf_High)
    para_out.extend([Taf])
    parameters['Taf']=Taf

    #Time of split between Europe and Middle East
    TEM_High=1200
    TEM_Low=400
    TEM=float(TEM_High)
    para_out.extend([TEM])
    parameters['TEM']=TEM

    #Time of split between CEU and CHB
    Teu_as_High=int(Taf)-1
    Teu_as_Low=int(TEM)+1
    Teu_as=float(Teu_as_High)
    para_out.extend([Teu_as])
    parameters['Teu_as']=Teu_as

    #Time of split between Jews and AJ
    TA_High=36
    TA_Low=20
    TA=float(TA_High)
    para_out.extend([TA])
    parameters['TA']=TA

    #Time of split between Jews and Middle East
    TMJ_High=int(TEM)-1
    TMJ_Low=int(TA)+1
    TMJ=float(TMJ_High)
    para_out.extend([TMJ])
    parameters['TMJ']=TMJ

    #Time of migration
    Tm_High=int(TA)-1
    Tm_Low=16
    Tm=float(randint(Tm_Low,Tm_High))
    para_out.extend([Tm])
    parameters['Tm']=Tm



    ##choose model/topology
    print "choosing case"


    #################
    if(parameters['Tgrowth_Af'] > parameters['Taf']):
        case=1

    if(parameters['Tgrowth_Af'] == parameters['Taf']):
        Tgrowth_Af+=0.00001
        parameters['Tgrowth_Af']=Tgrowth_Af
        case=1

    ################

    if(parameters['Taf'] > parameters['Tgrowth_Af'] > parameters['Teu_as']):
        case=2

    if(parameters['Taf'] > parameters['Tgrowth_Af'] == parameters['Teu_as']):
        Tgrowth_Af+=0.00001
        parameters['Tgrowth_Af']=Tgrowth_Af
        case=2

    ################

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['Tgrowth_Af'] > parameters['TEM']):
        case=3

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['Tgrowth_Af'] == parameters['TEM']):
        Tgrowth_Af+=0.00001
        parameters['Tgrowth_Af']=Tgrowth_Af
        case=3

    ################

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['Tgrowth_Af'] > parameters['TMJ']):
        case=4

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['Tgrowth_Af'] == parameters['TMJ']):
        Tgrowth_Af+=0.00001
        parameters['Tgrowth_Af']=Tgrowth_Af
        case=4


    ################

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['Tgrowth_Af'] > parameters['TA']):
        case=5

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['Tgrowth_Af'] == parameters['TA']):
        Tgrowth_Af+=0.00001
        parameters['Tgrowth_Af']=Tgrowth_Af
        case=5


    ################

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] > parameters['Tgrowth_Af'] > parameters['Tm']):
        case=6

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] > parameters['Tgrowth_Af'] == parameters['Tm']):
        Tgrowth_Af+=0.00001
        parameters['Tgrowth_Af']=Tgrowth_Af
        case=6


    ################

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] > parameters['Tm'] > parameters['Tgrowth_Af']):
        case=7

    if(parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] > parameters['Tm'] == parameters['Tgrowth_Af']):
        Tgrowth_Af+=-0.00001
        parameters['Tgrowth_Af']=Tgrowth_Af
        case=7



    ################
    #################
    #print 'params:', parameters
    #print 'case:',case

    return [parameters, para_out, case, daf]

def run_sim(parameters,case,length,chr_number):

    mu=2.5e-8
    rho=1e-8

    NAF=float(parameters['NAF'])
    NANC=float(parameters['NANC'])
    NCEU=float(parameters['NCEU'])
    NCHB=float(parameters['NCHB'])
    #NEU_AS=float(parameters['NEU_AS'])
    NA=float(parameters['NA'])
    NJ=float(parameters['NJ'])
    NM=float(parameters['NM'])

    rA=float(parameters['rA'])
    rMJ=float(parameters['rMJ'])

    m=float(parameters['m'])

    Tgrowth_Af=float(parameters['Tgrowth_Af'])
    Taf=float(parameters['Taf'])
    TEM=float(parameters['TEM'])
    Teu_as=float(parameters['Teu_as'])
    TA=float(parameters['TA'])
    TMJ=float(parameters['TMJ'])
    Tm=float(parameters['Tm'])


    macs_theta=float(mu*4*NANC)
    macs_rho=float(rho*4*NANC)
    scaled_NAF=float(NAF/NANC)
    scaled_NANC=float(NANC/NANC)
    scaled_NCEU=float(NCEU/NANC)
    scaled_NCHB=float(NCHB/NANC)
    #scaled_NEU_AS=float(NEU_AS/NANC)
    scaled_NA=float(NA/NANC)
    scaled_NJ=float(NJ/NANC)
    scaled_NM=float(NM/NANC)

    scaled_m=float(4*m*NANC)

    scaled_Tgrowth_Af=float(Tgrowth_Af/(4*NANC))
    scaled_Taf=float(Taf/(4*NANC))
    scaled_TEM=float(TEM/(4*NANC))
    scaled_Teu_as=float(Teu_as/(4*NANC))
    scaled_TA=float(TA/(4*NANC))
    scaled_TMJ=float(TMJ/(4*NANC))
    scaled_Tm=float(Tm/(4*NANC))


    #############
    if case==1:

        macs_args = ['./bin/macs',str(total),str(length),'-t',str(macs_theta),'-r',str(macs_rho),'-h','1e5','-R','genetic_map_b37/genetic_map_GRCh37_chr'+str(chr_number)+'.txt.macshs','-I','6',str(total_naf),str(total_nas),str(total_neu),str(nJ),str(nM),str(nA),'-n','1',str(scaled_NAF),'-n','2',str(scaled_NCHB),'-n','3',str(scaled_NCEU),'-n','4',str(scaled_NJ),'-n','5',str(scaled_NM),'-n','6',str(scaled_NA),'-eg','0','6',str(rA),'-eg','0.000001','4',str(rMJ),'-eg','0.000002','5',str(rMJ),'-em',str(scaled_Tm),'6','3',str(scaled_m),'-em',str(scaled_Tm+0.000001),'6','3','0','-ej',str(scaled_TA),'6','4','-ej',str(scaled_TMJ),'5','4','-ej',str(scaled_TEM),'4','3','-ej',str(scaled_Teu_as),'3','2','-ej',str(scaled_Taf),'2','1','-en',str(scaled_Tgrowth_Af),'1',str(scaled_NANC)]


    if case==2:

        macs_args = ['./bin/macs',str(total),str(length),'-t',str(macs_theta),'-r',str(macs_rho),'-h','1e5','-R','genetic_map_b37/genetic_map_GRCh37_chr'+str(chr_number)+'.txt.macshs','-I','6',str(total_naf),str(total_nas),str(total_neu),str(nJ),str(nM),str(nA),'-n','1',str(scaled_NAF),'-n','2',str(scaled_NCHB),'-n','3',str(scaled_NCEU),'-n','4',str(scaled_NJ),'-n','5',str(scaled_NM),'-n','6',str(scaled_NA),'-eg','0','6',str(rA),'-eg','0.000001','4',str(rMJ),'-eg','0.000002','5',str(rMJ),'-em',str(scaled_Tm),'6','3',str(scaled_m),'-em',str(scaled_Tm+0.000001),'6','3','0','-ej',str(scaled_TA),'6','4','-ej',str(scaled_TMJ),'5','4','-ej',str(scaled_TEM),'4','3','-ej',str(scaled_Teu_as),'3','2','-en',str(scaled_Tgrowth_Af),'1',str(scaled_NANC),'-ej',str(scaled_Taf),'2','1']


    if case==3:

        macs_args = ['./bin/macs',str(total),str(length),'-t',str(macs_theta),'-r',str(macs_rho),'-h','1e5','-R','genetic_map_b37/genetic_map_GRCh37_chr'+str(chr_number)+'.txt.macshs','-I','6',str(total_naf),str(total_nas),str(total_neu),str(nJ),str(nM),str(nA),'-n','1',str(scaled_NAF),'-n','2',str(scaled_NCHB),'-n','3',str(scaled_NCEU),'-n','4',str(scaled_NJ),'-n','5',str(scaled_NM),'-n','6',str(scaled_NA),'-eg','0','6',str(rA),'-eg','0.000001','4',str(rMJ),'-eg','0.000002','5',str(rMJ),'-em',str(scaled_Tm),'6','3',str(scaled_m),'-em',str(scaled_Tm+0.000001),'6','3','0','-ej',str(scaled_TA),'6','4','-ej',str(scaled_TMJ),'5','4','-ej',str(scaled_TEM),'4','3','-en',str(scaled_Tgrowth_Af),'1',str(scaled_NANC),'-ej',str(scaled_Teu_as),'3','2','-ej',str(scaled_Taf),'2','1']

    if case==4:

        macs_args = ['./bin/macs',str(total),str(length),'-t',str(macs_theta),'-r',str(macs_rho),'-h','1e5','-R','genetic_map_b37/genetic_map_GRCh37_chr'+str(chr_number)+'.txt.macshs','-I','6',str(total_naf),str(total_nas),str(total_neu),str(nJ),str(nM),str(nA),'-n','1',str(scaled_NAF),'-n','2',str(scaled_NCHB),'-n','3',str(scaled_NCEU),'-n','4',str(scaled_NJ),'-n','5',str(scaled_NM),'-n','6',str(scaled_NA),'-eg','0','6',str(rA),'-eg','0.000001','4',str(rMJ),'-eg','0.000002','5',str(rMJ),'-em',str(scaled_Tm),'6','3',str(scaled_m),'-em',str(scaled_Tm+0.000001),'6','3','0','-ej',str(scaled_TA),'6','4','-ej',str(scaled_TMJ),'5','4','-en',str(scaled_Tgrowth_Af),'1',str(scaled_NANC),'-ej',str(scaled_TEM),'4','3','-ej',str(scaled_Teu_as),'3','2','-ej',str(scaled_Taf),'2','1']

    if case==5:
        macs_args = ['./bin/macs',str(total),str(length),'-t',str(macs_theta),'-r',str(macs_rho),'-h','1e5','-R','genetic_map_b37/genetic_map_GRCh37_chr'+str(chr_number)+'.txt.macshs','-I','6',str(total_naf),str(total_nas),str(total_neu),str(nJ),str(nM),str(nA),'-n','1',str(scaled_NAF),'-n','2',str(scaled_NCHB),'-n','3',str(scaled_NCEU),'-n','4',str(scaled_NJ),'-n','5',str(scaled_NM),'-n','6',str(scaled_NA),'-eg','0','6',str(rA),'-eg','0.000001','4',str(rMJ),'-eg','0.000002','5',str(rMJ),'-em',str(scaled_Tm),'6','3',str(scaled_m),'-em',str(scaled_Tm+0.000001),'6','3','0','-ej',str(scaled_TA),'6','4','-en',str(scaled_Tgrowth_Af),'1',str(scaled_NANC),'-ej',str(scaled_TMJ),'5','4','-ej',str(scaled_TEM),'4','3','-ej',str(scaled_Teu_as),'3','2','-ej',str(scaled_Taf),'2','1']

    if case==6:

        macs_args = ['./bin/macs',str(total),str(length),'-t',str(macs_theta),'-r',str(macs_rho),'-h','1e5','-R','genetic_map_b37/genetic_map_GRCh37_chr'+str(chr_number)+'.txt.macshs','-I','6',str(total_naf),str(total_nas),str(total_neu),str(nJ),str(nM),str(nA),'-n','1',str(scaled_NAF),'-n','2',str(scaled_NCHB),'-n','3',str(scaled_NCEU),'-n','4',str(scaled_NJ),'-n','5',str(scaled_NM),'-n','6',str(scaled_NA),'-eg','0','6',str(rA),'-eg','0.000001','4',str(rMJ),'-eg','0.000002','5',str(rMJ),'-em',str(scaled_Tm),'6','3',str(scaled_m),'-em',str(scaled_Tm+0.000001),'6','3','0','-en',str(scaled_Tgrowth_Af),'1',str(scaled_NANC),'-ej',str(scaled_TA),'6','4','-ej',str(scaled_TMJ),'5','4','-ej',str(scaled_TEM),'4','3','-ej',str(scaled_Teu_as),'3','2','-ej',str(scaled_Taf),'2','1']

    if case==7:

        macs_args = ['./bin/macs',str(total),str(length),'-t',str(macs_theta),'-r',str(macs_rho),'-h','1e5','-R','genetic_map_b37/genetic_map_GRCh37_chr'+str(chr_number)+'.txt.macshs','-I','6',str(total_naf),str(total_nas),str(total_neu),str(nJ),str(nM),str(nA),'-n','1',str(scaled_NAF),'-n','2',str(scaled_NCHB),'-n','3',str(scaled_NCEU),'-n','4',str(scaled_NJ),'-n','5',str(scaled_NM),'-n','6',str(scaled_NA),'-eg','0','6',str(rA),'-eg','0.000001','4',str(rMJ),'-eg','0.000002','5',str(rMJ),'-en',str(scaled_Tgrowth_Af),'1',str(scaled_NANC),'-em',str(scaled_Tm),'6','3',str(scaled_m),'-em',str(scaled_Tm+0.000001),'6','3','0','-ej',str(scaled_TA),'6','4','-ej',str(scaled_TMJ),'5','4','-ej',str(scaled_TEM),'4','3','-ej',str(scaled_Teu_as),'3','2','-ej',str(scaled_Taf),'2','1']


    ##########

    #print macs_args
    sim=macsSwig.swigMain(len(macs_args),macs_args)

    return sim



#####set up simulations#####################################
############################################################ 

start_time=time.time()
elapsed_time=time.time()-start_time


##Length of chromosomes
lengths = [249163442, 243078003, 197813415, 191015739, 180695227, 170959304, 159091448, 146137372, 141069069, 135430928, 134747460, 133630123, 96085774, 87668527, 82491127, 90079543, 81032226, 78003657, 58843222, 62887650, 37234222, 35178458]

job=int(argv[1]) #must be a number
print 'JOB', job 


####Get parameter values from priors

param_model=param_sim_asc()
###parameters is a dictionary with the parameter values  
parameters=param_model[0]
###case is an integer that indicates which topology/model is going to be simulated  
para_out=param_model[1]
case=param_model[2]
daf=param_model[3]

####Samples to be simulated

naf_CGI=18
neu_CGI=18
nas_CGI=8

nA=76#528
nJ=28
nM=28#114

print 'naf_CGI '+str(naf_CGI)
print 'neu_CGI '+str(neu_CGI)
print 'nas_CGI '+str(nas_CGI)
print 'nA '+str(nA)
print 'nJ '+str(nJ)
print 'nM '+str(nM)

total_CGI=naf_CGI+neu_CGI+nas_CGI+nA+nJ+nM
print 'total samples '+str(total_CGI)

###Discovery panel
asc_nb_af=para_out[0]	
asc_nb_eu=para_out[1]
asc_nb_as=para_out[2]

total_naf=naf_CGI+asc_nb_af
total_neu=neu_CGI+asc_nb_eu
total_nas=nas_CGI+asc_nb_as


###Total number of chromosomes
total_asc=asc_nb_af+asc_nb_eu+asc_nb_as
total=total_CGI+total_asc


#fp=open('memory_profiler'+str(job)+'.log','w+')
#@profile(stream=fp)
def main():

    fprof=open('profile'+str(job)+'.log','w')
    p = psutil.Process()
    with p.oneshot():
        p.name()  # execute internal routine once collecting multiple info
        p.cpu_times()  # return cached value
        p.cpu_percent()  # return cached value
        p.create_time()  # return cached value
        p.ppid()  # return cached value
        p.status()  # return cached value

    profile = p.cpu_times(),p.memory_full_info()
    fprof.write('main_start\t'+str(profile)+'\n')

    chr=1

    #### Check if necessary directories exist.
    sim_data_dir = './sim_data_AJ_M1'
    sim_values_dir='./sim_values_AJ_M1'
    results_sims_dir='./results_sims_AJ_M1'

    if not os.path.exists(sim_data_dir):
        os.makedirs(sim_data_dir)
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
        snp_file=argv[2] #SNP file
        fileSNP=open(snp_file,'r')
        #print "read SNP file"
        SNP=[]
        for line in fileSNP:
            SNP.append(line)
        fileSNP.close()

        profile = p.cpu_times(),p.memory_full_info()
        fprof.write('read_snpfile_list\t'+str(profile)+'\n')


        ###get sites from snp array
        #print "get sites from snp array"
        snps=[]
        for line_snp in SNP:
            columns=line_snp.split('\t')
            snps.append(int(columns[2]))

        profile = p.cpu_times(),p.memory_full_info()
        fprof.write('append_snplist\t'+str(profile)+'\n')

        print 'nb Array snps', len(snps)
        nb_array_snps=len(snps)
        #print snps

        ###define simulation size
        size=argv[3]
        if size == "full":
            length=lengths[chr_number-1]
        else:
            length=size


        ##flag to check if the simulation work (generate the number of file
        flag_sim=False
        rep=1

        while flag_sim==False:

            #print 'rep', rep

            #####Run simulations
            print 'running simulation'
            sim=run_sim(parameters,case,length,chr_number)
            print 'finished simulation'

            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('run_sim\t'+str(profile)+'\n')

            elapsed_time=time.time()-start_time
            print '***********'+str(elapsed_time)+'***********'


            ##number of segregating sites
            nbss=sim.getNumSites()
            #print 'number sites in simulation', nbss

            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('genNumSites\t'+str(profile)+'\n')


            ##get position of the simulated sites and scale it to the "real" position in the SNP chip
            pos=[]
            for i in xrange(nbss):
                position=round(sim.getPosition(i)*(float(length)))#you have to give the number of the snp
                pos.append(position)
            
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('getPosition\t'+str(profile)+'\n')

                
            ##alleles will have all the information of all the simulated sites for all the pops
            #sites in each row
            alleles=[]
            for x in xrange(0,nbss):
                loc=[]
                for m in xrange(0,total):
                    loc.append(sim.getSite(x,m))
                alleles.append(loc)
            print 'total number of sites:',len(alleles) #number of elements in alleles
            print 'total number of chromosomes:',len(alleles[0])

            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('append_sim.getSite\t'+str(profile)+'\n')

            del sim


            #########get data from the simulations
            Talleles=zip(*alleles)
            
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('zip_alleles\t'+str(profile)+'\n')

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

            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('seq_lists\t'+str(profile)+'\n')

            #print asc_panel
            Tasc_panel=zip(*asc_panel)
            print 'number of sites in Tasc_panel:',len(Tasc_panel)
            print 'number of chromosomes in Tasc_panel:',len(Tasc_panel[0])
            #print Tasc_panel

            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('zip_asc_panel\t'+str(profile)+'\n')

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

            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('append_lists\t'+str(profile)+'\n')

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
            
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('len_lists\t'+str(profile)+'\n')


        #######################
        #####Calculate summary statistics from the regions for the CGI data
        print 'calculating summary stats for CGI data'
        elapsed_time=time.time()-start_time
        print '***********'+str(elapsed_time)+'***********'

        if nbss>0: ###no segregating sites in the simulations which is not possible

            Af_res=[]
            Af_res.extend(base_S_ss(seqAfCGI,nbss))
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('bases_S_ss\t'+str(profile)+'\n')
            pi_AfCGI=Pi2(Af_res[3],len(seqAfCGI))
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('Pi2\t'+str(profile)+'\n')
            Af_res.append(Tajimas(pi_AfCGI,Af_res[0],len(seqAfCGI)))
            del(Af_res[3])
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('Tajimas\t'+str(profile)+'\n')
            res.extend(Af_res)
            head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'

            Eu_res=[]
            Eu_res.extend(base_S_ss(seqEuCGI,nbss))
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('base_S_ss\t'+str(profile)+'\n')
            pi_EuCGI=Pi2(Eu_res[3],len(seqEuCGI))
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('Pi2\t'+str(profile)+'\n')
            Eu_res.append(Tajimas(pi_EuCGI,Eu_res[0],len(seqEuCGI)))
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('Tajimas\t'+str(profile)+'\n')
            del(Eu_res[3])
            res.extend(Eu_res)
            head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'

            As_res=[]
            As_res.extend(base_S_ss(seqAsCGI,nbss))
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('base_S_ss\t'+str(profile)+'\n')
            pi_AsCGI=Pi2(As_res[3],len(seqAsCGI))
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('Pi2\t'+str(profile)+'\n')
            As_res.append(Tajimas(pi_AsCGI,As_res[0],len(seqAsCGI)))
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('Tajimas\t'+str(profile)+'\n')
            del(As_res[3])
            res.extend(As_res)
            head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'


            ##fst between populations
            res.append(FST2(seqAfCGI,pi_AfCGI,naf_CGI,seqEuCGI,pi_EuCGI,neu_CGI))
            res.append(FST2(seqAfCGI,pi_AfCGI,naf_CGI,seqAsCGI,pi_AsCGI,nas_CGI))
            res.append(FST2(seqEuCGI,pi_EuCGI,neu_CGI,seqAsCGI,pi_AsCGI,nas_CGI))
            head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'
            #print 'len(res)', len(res)
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('FST2\t'+str(profile)+'\n')

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

        profile = p.cpu_times(),p.memory_full_info()
        fprof.write('pseudo_array\t'+str(profile)+'\n')
        
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

        profile = p.cpu_times(),p.memory_full_info()
        fprof.write('zip_lists\t'+str(profile)+'\n')


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

        profile = p.cpu_times(),p.memory_full_info()
        fprof.write('alleles_asc\t'+str(profile)+'\n')
                
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

        profile = p.cpu_times(),p.memory_full_info()
        fprof.write('pedfile\t'+str(profile)+'\n')

        ##Make map file
        filenamemap=str(sim_data_dir)+'/macs_asc_'+str(job)+'_chr'+str(chr_number)+'.map'
        filemap=open(filenamemap,'a')
        map=''
        for g in range(0,len(pos_asc)):
            map=str(chr_number)+' '+'chr'+str(chr_number)+'_'+str(pos_asc[g])+' '+str(int(avail_sites[pos_asc[g]]-1))+' '+str(int(avail_sites[pos_asc[g]]))+'\n'
            filemap.write(map)
        filemap.close()

        profile = p.cpu_times(),p.memory_full_info()
        fprof.write('map_file\t'+str(profile)+'\n')


        ###Genotypes for the ascertained SNPs
        seqAf_asc=zip(*allelesAf_asc)
        seqEu_asc=zip(*allelesEu_asc)
        seqAs_asc=zip(*allelesAs_asc)

        seqJ_asc=zip(*allelesJ_asc)
        seqM_asc=zip(*allelesM_asc)
        seqA_asc=zip(*allelesA_asc)

        profile = p.cpu_times(),p.memory_full_info()
        fprof.write('zip_alleles_asc\t'+str(profile)+'\n')

        #######
        #########calculate summary stats from the ascertained SNPs
        if nbss_asc>0:

            Af_asc=[]
            ss_Af_asc=base_S_ss(seqAf_asc,nbss_asc)
            if (ss_Af_asc[0]==0):
                #print "zeros"
                for i in xrange(5):
                    Af_asc.append(0)
                pi_Af_asc=0
            else:
                Af_asc.extend(base_S_ss(seqAf_asc,nbss_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('base_S_ss\t'+str(profile)+'\n')
                pi_Af_asc=Pi2(Af_asc[3],len(seqAf_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Pi2\t'+str(profile)+'\n')
                Af_asc.append(pi_Af_asc)
                Af_asc.append(Tajimas(pi_Af_asc,Af_asc[0],len(seqAf_asc)))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Tajimas\t'+str(profile)+'\n')
                del(Af_asc[3])

            res.extend(Af_asc)
            head = head + 'SegS_Af_ASC\tSing_Af_ASC\tDupl_Af_ASC\tPi_Af_ASC\tTajD_Af_ASC\t'
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('extend_res\t'+str(profile)+'\n')
            #print 'len(res)', len(res)

            ############

            Eu_asc=[]
            ss_Eu_asc=base_S_ss(seqEu_asc,nbss_asc)
            if (ss_Eu_asc[0]==0):
                #print "zeros"
                for i in xrange(5):
                    Eu_asc.append(0)
                pi_Eu_asc=0
            else:
                Eu_asc.extend(base_S_ss(seqEu_asc,nbss_asc))
                #print Eu_asc
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('base_S_ss\t'+str(profile)+'\n')
                pi_Eu_asc=Pi2(Eu_asc[3],len(seqEu_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Pi2\t'+str(profile)+'\n')
                Eu_asc.append(pi_Eu_asc)
                Eu_asc.append(Tajimas(pi_Eu_asc,Eu_asc[0],len(seqEu_asc)))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Tajimas\t'+str(profile)+'\n')
                del(Eu_asc[3])

            res.extend(Eu_asc)
            head = head + 'SegS_Eu_ASC\tSing_Eu_ASC\tDupl_Eu_ASC\tPi_Eu_ASC\tTajD_Eu_ASC\t'

            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('extend_res\t'+str(profile)+'\n')
            #print 'len(res)', len(res)
            ###########

            As_asc=[]
            ss_As_asc=base_S_ss(seqAs_asc,nbss_asc)
            if (ss_As_asc[0]==0):
                #print "zeros"
                for i in xrange(5):
                    As_asc.append(0)
                pi_As_asc=0
            else:
                As_asc.extend(base_S_ss(seqAs_asc,nbss_asc))
                #print As_asc
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('base_S_ss\t'+str(profile)+'\n')
                pi_As_asc=Pi2(As_asc[3],len(seqAs_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Pi2\t'+str(profile)+'\n')
                As_asc.append(pi_As_asc)
                As_asc.append(Tajimas(pi_As_asc,As_asc[0],len(seqAs_asc)))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Tajimas\t'+str(profile)+'\n')
                del(As_asc[3])

            res.extend(As_asc)
            head = head + 'SegS_As_ASC\tSing_As_ASC\tDupl_As_ASC\tPi_As_ASC\tTajD_As_ASC\t'
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('extend_res\t'+str(profile)+'\n')
            #print 'len(res)', len(res)
            ############

            J_asc=[]
            ss_J_asc=base_S_ss(seqJ_asc,nbss_asc)
            if (ss_J_asc[0]==0):
                #print "zeros"
                for i in xrange(5):
                    J_asc.append(0)
                pi_J_asc=0
            else:
                J_asc.extend(base_S_ss(seqJ_asc,nbss_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('base_S_ss\t'+str(profile)+'\n')
                pi_J_asc=Pi2(J_asc[3],len(seqJ_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Pi2\t'+str(profile)+'\n')
                J_asc.append(pi_J_asc)
                J_asc.append(Tajimas(pi_J_asc,J_asc[0],len(seqJ_asc)))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Tajimas\t'+str(profile)+'\n')
                del(J_asc[3])

            res.extend(J_asc)
            head = head + 'SegS_J_ASC\tSing_J_ASC\tDupl_J_ASC\tPi_J_ASC\tTajD_J_ASC\t'
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('extend_res\t'+str(profile)+'\n')
            #print 'len(res)', len(res)
            #############

            M_asc=[]
            ss_M_asc=base_S_ss(seqM_asc,nbss_asc)
            if (ss_M_asc[0]==0):
                #print "zeros"
                for i in xrange(5):
                    M_asc.append(0)
                pi_M_asc=0
            else:
                M_asc.extend(base_S_ss(seqM_asc,nbss_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('base_S_ss\t'+str(profile)+'\n')
                pi_M_asc=Pi2(M_asc[3],len(seqM_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Pi2\t'+str(profile)+'\n')
                M_asc.append(pi_M_asc)
                M_asc.append(Tajimas(pi_M_asc,M_asc[0],len(seqM_asc)))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Tajimas\t'+str(profile)+'\n')
                del(M_asc[3])

            res.extend(M_asc)
            head = head + 'SegS_M_ASC\tSing_M_ASC\tDupl_M_ASC\tPi_M_ASC\tTajD_M_ASC\t'
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('extend_res\t'+str(profile)+'\n')
            #print 'len(res)', len(res)
            #############

            A_asc=[]
            ss_A_asc=base_S_ss(seqA_asc,nbss_asc)
            if (ss_A_asc[0]==0):
                #print "zeros"
                for i in xrange(5):
                    A_asc.append(0)
                pi_A_asc=0
            else:
                A_asc.extend(base_S_ss(seqA_asc,nbss_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('base_S_ss\t'+str(profile)+'\n')
                pi_A_asc=Pi2(A_asc[3],len(seqA_asc))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Pi2\t'+str(profile)+'\n')
                A_asc.append(pi_A_asc)
                A_asc.append(Tajimas(pi_A_asc,A_asc[0],len(seqA_asc)))
                profile = p.cpu_times(),p.memory_full_info()
                fprof.write('Tajimas\t'+str(profile)+'\n')
                del(A_asc[3])

            res.extend(A_asc)
            head = head + 'SegS_A_ASC\tSing_A_ASC\tDupl_A_ASC\tPi_A_ASC\tTajD_A_ASC\t'
            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('extend_res\t'+str(profile)+'\n')
            #print 'len(res)', len(res)
            #############


            ##fst between populations
            res.append(FST2(seqAf_asc,pi_Af_asc,naf_CGI,seqEu_asc,pi_Eu_asc,neu_CGI))
            res.append(FST2(seqAf_asc,pi_Af_asc,naf_CGI,seqAs_asc,pi_As_asc,nas_CGI))
            res.append(FST2(seqEu_asc,pi_Eu_asc,neu_CGI,seqAs_asc,pi_As_asc,nas_CGI))
            head = head + 'FST_AfEu_ASC\tFST_AfAs_ASC_m\tFST_EuAs_ASC\t'


            res.append(FST2(seqA_asc,pi_A_asc,nA,seqEu_asc,pi_Eu_asc,neu_CGI))
            res.append(FST2(seqA_asc,pi_A_asc,nA,seqJ_asc,pi_J_asc,nJ))
            res.append(FST2(seqA_asc,pi_A_asc,nA,seqM_asc,pi_M_asc,nM))
            res.append(FST2(seqM_asc,pi_M_asc,nM,seqJ_asc,pi_J_asc,nJ))
            head = head + 'FST_AEu_ASC\tFST_AJ_ASC\tFST_AM_ASC\tFST_MJ_ASC\n'
            #print 'len(res) FST', len(res)

            profile = p.cpu_times(),p.memory_full_info()
            fprof.write('FST2\t'+str(profile)+'\n')

            cont=cont+1
            #print 'AFS!!'

            #print res
            #print 'len(res) final '+str(len(res))



    #print 'cont'
    #print cont
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

    #sys.stdout = LogFile('memory_profiler'+str(job)+'.log')
    fprof.close()

if __name__ == '__main__':
    main()
