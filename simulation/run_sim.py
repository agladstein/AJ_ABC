import macsSwig

def run_sim(parameters,case,length,chr_number,total,total_naf,total_nas,total_neu,nJ,nM,nA):

    mu=2.5e-8
    rho=1e-8

    NAF=float(parameters['NAF'])
    NANC=float(parameters['NANC'])
    NCEU=float(parameters['NCEU'])
    NCHB=float(parameters['NCHB'])
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

    sim=macsSwig.swigMain(len(macs_args),macs_args)

    return sim
