import math
import random
from random import randint


def param_sim_asc_rand():
    """Define the parameters for macs simulation Model 1.
    If running ABC, parameter values come from the priors."""

    para_out = []
    parameters = {}

    ###Discovery panel
    low = 2
    high = 20

    asc_nb_af = randint(low, high)
    asc_nb_eu = randint(low, high)
    asc_nb_as = randint(low, high)

    para_out.extend([asc_nb_af])
    para_out.extend([asc_nb_eu])
    para_out.extend([asc_nb_as])

    daf = random.uniform(0.05, 0.10)
    para_out.extend([daf])

    ####Demographic model
    # population size in Africa
    NAF = float(round(10 ** random.uniform(3.7, 5.0)))
    para_out.extend([math.log10(NAF)])
    parameters['NAF'] = NAF
    # print NAF

    # Ancestral population size, before population growth in AF
    # choose ancestral Ne based on being some value smaller than NAF between 1 or 0.1 x.
    # population growth
    Nrat_High = 0.0  # Allow only growth for now
    Nrat_Low = -1.0
    NANC = float(round((10 ** random.uniform(Nrat_Low, Nrat_High)) * NAF))
    para_out.extend([math.log10(NANC)])
    parameters['NANC'] = NANC

    NCEU = float(round(10 ** random.uniform(3.0, 5.0)))
    para_out.extend([math.log10(NCEU)])
    parameters['NCEU'] = NCEU

    NCHB = float(round(10 ** random.uniform(3.0, 5.0)))
    para_out.extend([math.log10(NCHB)])
    parameters['NCHB'] = NCHB

    # Ancestral population size of Eurasians
    # NEU_AS=float(randint(1500,5000))
    # para_out.extend([NEU_AS])
    # parameters['NEU_AS']=NEU_AS

    # Population size of WAJ
    NWA = float(round(10 ** random.uniform(3.0, 6.7)))
    para_out.extend([math.log10(NWA)])
    parameters['NWA'] = NWA

    # Population size of EAJ
    NEA = float(round(10 ** random.uniform(4.0, 6.7)))
    para_out.extend([math.log10(NEA)])
    parameters['NEA'] = NEA

    # Population size of Jews
    NJ = float(round(10 ** random.uniform(3.0, 6.0)))
    para_out.extend([math.log10(NJ)])
    parameters['NJ'] = NJ

    # Population size of Middle Easterns
    NM = float(round(10 ** random.uniform(3.0, 6.0)))
    para_out.extend([math.log10(NM)])
    parameters['NM'] = NM

    # Growth rate in WAJ
    rWA = float(round(10 ** random.uniform(0, 1)))
    para_out.extend([math.log10(rWA)])
    parameters['rWA'] = rWA

    # Growth rate in EAJ
    rEA = float(round(10 ** random.uniform(0, 1)))
    para_out.extend([math.log10(rEA)])
    parameters['rEA'] = rEA

    # Growth rate in Jews and Middle East
    rMJ = float(round(10 ** random.uniform(0, 1)))
    para_out.extend([math.log10(rMJ)])
    parameters['rMJ'] = rMJ

    # migration rate from Europe to AJ
    m_High = 1.0
    m_Low = 0.0
    m = random.uniform(m_Low,m_High)
    para_out.extend([m])
    parameters['m'] = m

    # Time of the instantaneous growth in Africa, before the split between Africans and non Africans
    Tgrowth_Low = 1
    Tgrowth_High = 4100
    Tgrowth_Af = float(randint(Tgrowth_Low, Tgrowth_High))
    parameters['Tgrowth_Af'] = Tgrowth_Af
    para_out.extend([Tgrowth_Af])

    # Time of split between YRI and CEU/CHB
    Taf_High = 4100  # 102,500 years using 25 years per generation
    Taf_Low = 1600  # 40,000 years using 25 years per generation.
    Taf = float(randint(Taf_Low, Taf_High))
    para_out.extend([Taf])
    parameters['Taf'] = Taf

    # Time of split between Europe and Middle East
    TEM_High = 1200
    TEM_Low = 400
    TEM = float(randint(TEM_Low, TEM_High))
    para_out.extend([TEM])
    parameters['TEM'] = TEM

    # Time of split between CEU and CHB
    Teu_as_High = int(Taf) - 1
    Teu_as_Low = int(TEM) + 1
    Teu_as = float(randint(Teu_as_Low, Teu_as_High))
    para_out.extend([Teu_as])
    parameters['Teu_as'] = Teu_as

    # Time of split between Jews and AJ
    TA_High = 36
    TA_Low = 20
    TA = float(randint(TA_Low, TA_High))
    para_out.extend([TA])
    parameters['TA'] = TA

    # Time of split between Jews and Middle East
    TMJ_High = int(TEM) - 1
    TMJ_Low = int(TA) + 1
    TMJ = float(randint(TMJ_Low, TMJ_High))
    para_out.extend([TMJ])
    parameters['TMJ'] = TMJ

    # Time of split between Eastern and Western AJ
    TAEW_High = (TA) - 2
    TAEW_Low = 1
    TAEW = float(randint(TAEW_Low, TAEW_High))
    para_out.extend([TAEW])
    parameters['TAEW'] = TAEW

    # Time of migration
    Tm_High = int(TA) - 1
    Tm_Low = int(TAEW) + 1
    Tm = float(randint(Tm_Low, Tm_High))
    para_out.extend([Tm])
    parameters['Tm'] = Tm

    case, modified_Tgrowth_Af = choose_case(parameters)

    parameters['Tgrowth_Af'] = modified_Tgrowth_Af

    return [parameters, para_out, case, daf]

def param_sim_asc_min():
    """Define the parameters for macs simulation Model 2.
    If running ABC, parameter values come from the priors.
    Some of the parameters are the minimum values from the prior distribution.
    The minimum values are chosen to minimize coalescent times.
    This function should only be used for testing purposes, not for running ABC.
    This option will cause the pipeline to run faster and use less memory."""

    para_out=[]
    parameters={}

    ###Discovery panel
    low=2

    asc_nb_af=low
    asc_nb_eu=low
    asc_nb_as=low

    para_out.extend([asc_nb_af])
    para_out.extend([asc_nb_eu])
    para_out.extend([asc_nb_as])

    daf=random.uniform(0.05,0.10)
    para_out.extend([daf])


    ####Demographic model
    #population size in Africa
    NAF=float(round(10**3.7))
    para_out.extend([math.log10(NAF)])
    parameters['NAF']=NAF
    #print NAF

    #Ancestral population size, before population growth in AF
    #choose ancestral Ne based on being some value smaller than NAF between 1 or 0.1 x.
    #population growth
    Nrat_Low=-1.0
    NANC=float(round(10**Nrat_Low*NAF))
    para_out.extend([math.log10(NANC)])
    parameters['NANC']=NANC

    NCEU=float(round(10**3.0))
    para_out.extend([math.log10(NCEU)])
    parameters['NCEU']=NCEU

    NCHB=float(round(10**3.0))
    para_out.extend([math.log10(NCHB)])
    parameters['NCHB']=NCHB

    # Population size of WAJ
    NWA = float(round(10 **3.0))
    para_out.extend([math.log10(NWA)])
    parameters['NWA'] = NWA

    # Population size of EAJ
    NEA = float(round(10 **4.0))
    para_out.extend([math.log10(NEA)])
    parameters['NEA'] = NEA

    #Population size of Jews
    NJ=float(round(10**3.0))
    para_out.extend([math.log10(NJ)])
    parameters['NJ']=NJ

    #Population size of Middle Easterns
    NM=float(round(10**3.0))
    para_out.extend([math.log10(NM)])
    parameters['NM']=NM

    # Growth rate in WAJ
    rWA = float(round(10 ** random.uniform(0, 1)))
    para_out.extend([math.log10(rWA)])
    parameters['rWA'] = rWA

    # Growth rate in EAJ
    rEA = float(round(10 ** random.uniform(0, 1)))
    para_out.extend([math.log10(rEA)])
    parameters['rEA'] = rEA

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
    Taf_Low=1600				  #40,000 years using 25 years per generation.
    Taf=float(Taf_Low)
    para_out.extend([Taf])
    parameters['Taf']=Taf

    #Time of split between Europe and Middle East
    TEM_Low=400
    TEM=float(TEM_Low)
    para_out.extend([TEM])
    parameters['TEM']=TEM

    #Time of split between CEU and CHB
    Teu_as_Low=int(TEM)+1
    Teu_as=float(Teu_as_Low)
    para_out.extend([Teu_as])
    parameters['Teu_as']=Teu_as

    #Time of split between Jews and AJ
    TA_Low=20
    TA=float(TA_Low)
    para_out.extend([TA])
    parameters['TA']=TA

    #Time of split between Jews and Middle East
    TMJ_Low=int(TA)+1
    TMJ=float(TMJ_Low)
    para_out.extend([TMJ])
    parameters['TMJ']=TMJ

    # Time of split between Eastern and Western AJ
    TAEW_Low = 1
    TAEW = float(TAEW_Low)
    para_out.extend([TAEW])
    parameters['TAEW'] = TAEW

    #Time of migration
    Tm_High=int(TA)-1
    Tm_Low=16
    Tm=float(randint(Tm_Low,Tm_High))
    para_out.extend([Tm])
    parameters['Tm']=Tm

    case, modified_Tgrowth_Af = choose_case(parameters)

    parameters['Tgrowth_Af'] = modified_Tgrowth_Af

    return [parameters, para_out, case, daf]


def param_sim_asc_max():
    """Define the parameters for macs simulation Model 2.
    If running ABC, parameter values come from the priors.
    Some of the parameters are the maximum values from the prior distribution.
    The maximum values are chosen to maximize coalescent times.
    This function should only be used for testing purposes, not for running ABC.
    This option will cause the pipeline to run slower and use more memory."""

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

    #Ancestral population size, before population growth in AF
    #choose ancestral Ne based on being some value smaller than NAF between 1 or 0.1 x.
    #population growth
    Nrat_High=0.0   #Allow only growth for now
    NANC=float(round(10**Nrat_High*NAF))
    para_out.extend([math.log10(NANC)])
    parameters['NANC']=NANC

    NCEU=float(round(10**5.0))
    para_out.extend([math.log10(NCEU)])
    parameters['NCEU']=NCEU

    NCHB=float(round(10**5.0))
    para_out.extend([math.log10(NCHB)])
    parameters['NCHB']=NCHB

    # Population size of WAJ
    NWA = float(round(10 **6.7))
    para_out.extend([math.log10(NWA)])
    parameters['NWA'] = NWA

    # Population size of EAJ
    NEA = float(round(10 **6.7))
    para_out.extend([math.log10(NEA)])
    parameters['NEA'] = NEA

    #Population size of Jews
    NJ=float(round(10**6.0))
    para_out.extend([math.log10(NJ)])
    parameters['NJ']=NJ

    #Population size of Middle Easterns
    NM=float(round(10**6.0))
    para_out.extend([math.log10(NM)])
    parameters['NM']=NM

    # Growth rate in WAJ
    rWA = float(round(10 ** random.uniform(0, 1)))
    para_out.extend([math.log10(rWA)])
    parameters['rWA'] = rWA

    # Growth rate in EAJ
    rEA = float(round(10 ** random.uniform(0, 1)))
    para_out.extend([math.log10(rEA)])
    parameters['rEA'] = rEA

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
    Taf=float(Taf_High)
    para_out.extend([Taf])
    parameters['Taf']=Taf

    #Time of split between Europe and Middle East
    TEM_High=1200
    TEM=float(TEM_High)
    para_out.extend([TEM])
    parameters['TEM']=TEM

    #Time of split between CEU and CHB
    Teu_as_High=int(Taf)-1
    Teu_as=float(Teu_as_High)
    para_out.extend([Teu_as])
    parameters['Teu_as']=Teu_as

    #Time of split between Jews and AJ
    TA_High=36
    TA=float(TA_High)
    para_out.extend([TA])
    parameters['TA']=TA

    #Time of split between Jews and Middle East
    TMJ_High=int(TEM)-1
    TMJ=float(TMJ_High)
    para_out.extend([TMJ])
    parameters['TMJ']=TMJ

    # Time of split between Eastern and Western AJ
    TAEW_High = (TA) - 2
    TAEW = float(TAEW_High)
    para_out.extend([TAEW])
    parameters['TAEW'] = TAEW

    #Time of migration
    Tm_High=int(TA)-1
    Tm_Low=16
    Tm=float(randint(Tm_Low,Tm_High))
    para_out.extend([Tm])
    parameters['Tm']=Tm


    case, modified_Tgrowth_Af = choose_case(parameters)

    parameters['Tgrowth_Af'] = modified_Tgrowth_Af

    return [parameters, para_out, case, daf]


def choose_case(parameters):
    ##choose model/topology
    print "choosing case"
    Tgrowth_Af = parameters['Tgrowth_Af']

    #################
    if (parameters['Tgrowth_Af'] > parameters['Taf']):
        case = 1

    if (parameters['Tgrowth_Af'] == parameters['Taf']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 1

    ################

    if (parameters['Taf'] > parameters['Tgrowth_Af'] > parameters['Teu_as']):
        case = 2

    if (parameters['Taf'] > parameters['Tgrowth_Af'] == parameters['Teu_as']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 2

    ################

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['Tgrowth_Af'] > parameters['TEM']):
        case = 3

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['Tgrowth_Af'] == parameters['TEM']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 3

    ################

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['Tgrowth_Af'] > parameters['TMJ']):
        case = 4

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['Tgrowth_Af'] == parameters['TMJ']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 4

    ################

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['Tgrowth_Af'] >
            parameters['TA']):
        case = 5

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['Tgrowth_Af'] ==
            parameters['TA']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 5

    ################

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] >
            parameters['Tgrowth_Af'] > parameters['Tm']):
        case = 6

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] >
            parameters['Tgrowth_Af'] == parameters['Tm']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 6

    ################

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] >
            parameters['Tm'] > parameters['Tgrowth_Af']):
        case = 7

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] >
            parameters['Tm'] == parameters['Tgrowth_Af']):
        Tgrowth_Af += -0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 7

    ################

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] >
            parameters['Tm'] > parameters['TAEW'] > parameters['Tgrowth_Af']):
        case = 8

    if (parameters['Taf'] > parameters['Teu_as'] > parameters['TEM'] > parameters['TMJ'] > parameters['TA'] >
            parameters['Tm'] > parameters['TAEW'] == parameters['Tgrowth_Af']):
        Tgrowth_Af += -0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 8

    return case, Tgrowth_Af