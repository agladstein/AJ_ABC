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
import numpy as np
import re

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

###end summary statistics###################################
############################################################


def main():

	chr=1

	naf_CGI = 18
	neu_CGI = 18
	nas_CGI = 8
	nA = 76
	nJ = 28
	nM = 28
	total = naf_CGI+neu_CGI+nas_CGI+nA+nJ+nM

	res=[]

	infile=argv[1] #real data
	print 'reading file '+str(infile)

	file=open(infile, 'r')
	# Talleles = []
	# for line in file:
	# 	columns = line.split(' ')
	# 	Talleles.append(columns[6:len(columns)-6])
	# file.close()


	Talleles=file.read()
	Talleles=string.split(Talleles,'\n')
	file.close()


	del(Talleles[-1]) ### why are you doing this?
	# del(Talleles[0:5])

	nbss=len(Talleles[1])
	#print nbss

	if len(Talleles)!=266: #int((naf_CGI+neu_CGI+nas_CGI+nA+nJ+nM)/2):
		print "something is wrong"
		return
	else:
		alleles=zip(*Talleles)
		print 'total number of sites:',len(alleles)
		print 'total number of individuals:',len(alleles[0])
		
        ###get genotypes for array
		seqAf=Talleles[0:naf]
		seqEu=Talleles[naf:naf+neu]
		seqAs=Talleles[naf+neu:naf+neu+nas]
		seqJ=Talleles[naf+neu+nas:naf+neu+nas+nJ]
		seqM=Talleles[naf+neu+nas+nJ:naf+neu+nas+nJ+nM]
        seqA=Talleles[naf+neu+nas+nJ+nM:naf+neu+nas+nJ+nM+nA]


		#####Calculate summary statistics from the CGI data
		Af_res=[]
		Af_res.extend(base_S_ss(seqAfCGI,nbss))
		pi_AfCGI=Pi2(Af_res[3],len(seqAfCGI))
		Af_res.append(Tajimas(pi_AfCGI,Af_res[0],len(seqAfCGI)))
		del(Af_res[3])
		res.extend(Af_res)

		Eu_res=[]
		Eu_res.extend(base_S_ss(seqEuCGI,nbss))
		pi_EuCGI=Pi2(Eu_res[3],len(seqEuCGI))
		Eu_res.append(Tajimas(pi_EuCGI,Eu_res[0],len(seqEuCGI)))
		del(Eu_res[3])
		res.extend(Eu_res)
		
		As_res=[]
		As_res.extend(base_S_ss(seqAsCGI,nbss))
		pi_AsCGI=Pi2(As_res[3],len(seqAsCGI))
		As_res.append(Tajimas(pi_AsCGI,As_res[0],len(seqAsCGI)))
		del(As_res[3])
		res.extend(As_res)

		###fst between populations    
		res.append(FST2(seqAfCGI,pi_AfCGI,naf_CGI,seqEuCGI,pi_EuCGI,neu_CGI))
		res.append(FST2(seqAfCGI,pi_AfCGI,naf_CGI,seqAsCGI,pi_AsCGI,nas_CGI))
		res.append(FST2(seqEuCGI,pi_EuCGI,neu_CGI,seqAsCGI,pi_AsCGI,nas_CGI))

	outfile = open('ss_allpops.txt', 'w')

	####head of the file with the summary statistics of the regions
	head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'
	head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'
	head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'

	head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'

	#################

	head = head + 'IBD_mean_AA\tIBD_mean_JJ\tIBD_mean_MM\tIBD_mean_EE\tIBD_mean_AE\tIBD_mean_AJ\tIBD_mean_AM\tIBD_mean_JM\tIBD_mean_JE\tIBD_mean_ME\t'
	head = head + 'IBD_median_AA\tIBD_median_JJ\tIBD_median_MM\tIBD_median_EE\tIBD_median_AE\tIBD_median_AJ\tIBD_median_AM\tIBD_median_JM\tIBD_median_JE\tIBD_median_ME\t'
	head = head + 'IBD_num_AA\tIBD_num_JJ\tIBD_num_MM\tIBD_num_EE\tIBD_num_AE\tIBD_num_AJ\tIBD_num_AM\tIBD_num_JM\tIBD_num_JE\tIBD_num_ME\t'
	head = head + 'IBD_var_AA\tIBD_var_JJ\tIBD_var_MM\tIBD_var_EE\tIBD_var_AE\tIBD_var_AJ\tIBD_var_AM\tIBD_var_JM\tIBD_var_JE\tIBD_var_ME\t'

	#################

	head = head + 'SegS_Af_ASC\tSing_Af_ASC\tDupl_Af_ASC\tPi_Af_ASC\tTajD_Af_ASC\t'
	head = head + 'SegS_Eu_ASC\tSing_Eu_ASC\tDupl_Eu_ASC\tPi_Eu_ASC\tTajD_Eu_ASC\t'
	head = head + 'SegS_As_ASC\tSing_As_ASC\tDupl_As_ASC\tPi_As_ASC\tTajD_As_ASC\t'

	head = head + 'SegS_AJ_ASC\tSing_AJ_ASC\tDupl_AJ_ASC\tPi_AJ_ASC\tTajD_AJ_ASC\t'

	head = head + 'FST_AfEu_ASC\tFST_AfAs_ASC_m\tFST_EuAs_ASC\t'
	head = head + 'FST_AJEu_ASC\n'

	#################
	#################

	outfile.write(head)
	outfile.close()


if __name__ == '__main__':
	main()
