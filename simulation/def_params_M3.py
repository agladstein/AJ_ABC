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

    # migration rate from Europe to EAJ
    mE_High = 1.0
    mE_Low = 0.0
    mE = random.uniform(mE_Low,mE_High)
    para_out.extend([mE])
    parameters['mE'] = mE

    # migration rate from Europe to WAJ
    mW_High = 1.0
    mW_Low = 0.0
    mW = random.uniform(mW_Low, mW_High)
    para_out.extend([mW])
    parameters['mW'] = mW

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
    TAEW_Low = 3
    TAEW = float(randint(TAEW_Low, TAEW_High))
    para_out.extend([TAEW])
    parameters['TAEW'] = TAEW

    # Time of migration to EAJ
    TmE_High = int(TAEW) - 1
    TmE_Low = 1
    TmE = float(randint(TmE_Low, TmE_High))
    para_out.extend([TmE])
    parameters['TmE'] = TmE

    # Time of migration to WAJ
    TmW_High = int(TAEW) - 1
    TmW_Low = 1
    TmW = float(randint(TmW_Low, TmW_High))
    para_out.extend([TmW])
    parameters['TmW'] = TmW

    # Growth rate in WAJ
    rWA_High = -(1/TAEW) * math.log(10/NWA) # set max growth rate so the minimum number of individuals at East West split is 10
    rWA_Low = 0.0
    rWA = random.uniform(rWA_Low,rWA_High)
    para_out.extend([rWA])
    parameters['rWA'] = rWA

    # Growth rate in EAJ
    rEA_High = -(1/TAEW) * math.log(10/NEA) # set max growth rate so the minimum number of individuals at East West split is 10
    rEA_Low = 0.0
    rEA = random.uniform(rEA_Low,rEA_High)
    para_out.extend([rEA])
    parameters['rEA'] = rEA

    # Growth rate in Jews and Middle East
    if NM < NJ:
        rMJ_High = -(1/TMJ) * math.log(10/NM)
    else:
        rMJ_High = -(1 / TMJ) * math.log(10 / NJ)
    rMJ_Low = 0.0
    rMJ = random.uniform(rMJ_Low, rMJ_High)
    para_out.extend([rMJ])
    parameters['rMJ'] = rMJ

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

    # migration rate from Europe to EAJ
    mE_High = 1
    mE_Low = 0
    mE = float(randint(mE_Low, mE_High))
    para_out.extend([mE])
    parameters['mE'] = mE

    # migration rate from Europe to WAJ
    mW_High = 1
    mW_Low = 0
    mW = float(randint(mW_Low, mW_High))
    para_out.extend([mW])
    parameters['mW'] = mW

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
    TAEW_Low = 3
    TAEW = float(TAEW_Low)
    para_out.extend([TAEW])
    parameters['TAEW'] = TAEW

    # Time of migration to EAJ
    TmE_High = int(TAEW) - 1
    TmE_Low = 1
    TmE = float(randint(TmE_Low, TmE_High))
    para_out.extend([TmE])
    parameters['TmE'] = TmE

    # Time of migration to WAJ
    TmW_High = int(TAEW) - 1
    TmW_Low = 1
    TmW = float(randint(TmW_Low, TmW_High))
    para_out.extend([TmW])
    parameters['TmW'] = TmW

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

    # migration rate from Europe to EAJ
    mE_High = 1
    mE_Low = 0
    mE = float(randint(mE_Low, mE_High))
    para_out.extend([mE])
    parameters['mE'] = mE

    # migration rate from Europe to WAJ
    mW_High = 1
    mW_Low = 0
    mW = float(randint(mW_Low, mW_High))
    para_out.extend([mW])
    parameters['mW'] = mW

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

    # Time of migration to EAJ
    TmE_High = int(TAEW) - 1
    TmE_Low = 1
    TmE = float(randint(TmE_Low, TmE_High))
    para_out.extend([TmE])
    parameters['TmE'] = TmE

    # Time of migration to WAJ
    TmW_High = int(TAEW) - 1
    TmW_Low = 1
    TmW = float(randint(TmW_Low, TmW_High))
    para_out.extend([TmW])
    parameters['TmW'] = TmW


    case, modified_Tgrowth_Af = choose_case(parameters)

    parameters['Tgrowth_Af'] = modified_Tgrowth_Af

    return [parameters, para_out, case, daf]


def choose_case(parameters):
    ##choose model/topology
    print "choosing case"
    Tgrowth_Af = parameters['Tgrowth_Af']
    TmE = parameters['TmE']

    #################
    if (parameters['Tgrowth_Af'] > parameters['Taf'] > parameters['TmE'] > parameters['TmW']):
        case = 1

    if (parameters['Tgrowth_Af'] == parameters['Taf'] > parameters['TmE'] > parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 1

    if (parameters['Tgrowth_Af'] > parameters['Taf'] > parameters['TmE'] == parameters['TmW']):
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 1

    if (parameters['Tgrowth_Af'] == parameters['Taf'] > parameters['TmE'] == parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 1

    #################
    if (parameters['Tgrowth_Af'] > parameters['Taf'] > parameters['TmW'] > parameters['TmE']):
        case = 2

    if (parameters['Tgrowth_Af'] == parameters['Taf'] > parameters['TmW'] > parameters['TmE']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 2

    ################

    if (parameters['Taf'] > parameters['Tgrowth_Af'] > parameters['Teu_as'] > parameters['TmE'] > parameters['TmW']):
        case = 3

    if (parameters['Taf'] > parameters['Tgrowth_Af'] == parameters['Teu_as'] > parameters['TmE'] > parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 3

    if (parameters['Taf'] > parameters['Tgrowth_Af'] > parameters['Teu_as'] > parameters['TmE'] == parameters['TmW']):
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 3

    if (parameters['Taf'] > parameters['Tgrowth_Af'] == parameters['Teu_as'] > parameters['TmE'] == parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 3

    ################


    if (parameters['Taf'] > parameters['Tgrowth_Af'] > parameters['Teu_as'] > parameters['TmW'] > parameters['TmE']):
        case = 4

    if (parameters['Taf'] > parameters['Tgrowth_Af'] == parameters['Teu_as'] > parameters['TmW'] > parameters['TmE']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 4

    ################

    if (parameters['Teu_as'] > parameters['Tgrowth_Af'] > parameters['TEM'] > parameters['TmE'] > parameters['TmW']):
        case = 5

    if (parameters['Teu_as'] > parameters['Tgrowth_Af'] == parameters['TEM'] > parameters['TmE'] > parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 5

    if (parameters['Teu_as'] > parameters['Tgrowth_Af'] > parameters['TEM'] > parameters['TmE'] == parameters['TmW']):
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 5

    if (parameters['Teu_as'] > parameters['Tgrowth_Af'] == parameters['TEM'] > parameters['TmE'] == parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 5

    ################

    if (parameters['Teu_as'] > parameters['Tgrowth_Af'] > parameters['TEM'] > parameters['TmW'] > parameters['TmE']):
        case = 6

    if (parameters['Teu_as'] > parameters['Tgrowth_Af'] == parameters['TEM'] > parameters['TmW'] > parameters['TmE']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 6

    ################

    if (parameters['TEM'] > parameters['Tgrowth_Af'] > parameters['TMJ'] > parameters['TmE'] > parameters['TmW']):
        case = 7

    if (parameters['TEM'] > parameters['Tgrowth_Af'] == parameters['TMJ'] > parameters['TmE'] > parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 7


    if (parameters['TEM'] > parameters['Tgrowth_Af'] > parameters['TMJ'] > parameters['TmE'] == parameters['TmW']):
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 7

    if (parameters['TEM'] > parameters['Tgrowth_Af'] == parameters['TMJ'] > parameters['TmE'] == parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 7

    ################

    if (parameters['TEM'] > parameters['Tgrowth_Af'] > parameters['TMJ'] > parameters['TmW'] > parameters['TmE']):
        case = 8

    if (parameters['TEM'] > parameters['Tgrowth_Af'] == parameters['TMJ'] > parameters['TmW'] > parameters['TmE']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 8

    ################

    if (parameters['TMJ'] > parameters['Tgrowth_Af'] > parameters['TA'] > parameters['TmE'] > parameters['TmW']):
        case = 9

    if (parameters['TMJ'] > parameters['Tgrowth_Af'] == parameters['TA'] > parameters['TmE'] > parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 9

    if (parameters['TMJ'] > parameters['Tgrowth_Af'] > parameters['TA'] > parameters['TmE'] == parameters['TmW']):
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 9

    if (parameters['TMJ'] > parameters['Tgrowth_Af'] == parameters['TA'] > parameters['TmE'] == parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 9

    ################

    if (parameters['TMJ'] > parameters['Tgrowth_Af'] > parameters['TA'] > parameters['TmW'] > parameters['TmE']):
        case = 10

    if (parameters['TMJ'] > parameters['Tgrowth_Af'] == parameters['TA'] > parameters['TmW'] > parameters['TmE']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 10

    ################

    if (parameters['TA'] > parameters['Tgrowth_Af'] > parameters['TAEW'] > parameters['TmE'] > parameters['TmW']):
        case = 11

    if (parameters['TA'] > parameters['Tgrowth_Af'] == parameters['TAEW'] > parameters['TmE'] > parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 11

    if (parameters['TA'] > parameters['Tgrowth_Af'] > parameters['TAEW'] > parameters['TmE'] == parameters['TmW']):
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 11

    if (parameters['TA'] > parameters['Tgrowth_Af'] == parameters['TAEW'] > parameters['TmE'] == parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 11

    ################

    if (parameters['TA'] > parameters['Tgrowth_Af'] > parameters['TAEW'] > parameters['TmW'] > parameters['TmE']):
        case = 12

    if (parameters['TA'] > parameters['Tgrowth_Af'] == parameters['TAEW'] > parameters['TmW'] > parameters['TmE']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 12

    ################

    if (parameters['TAEW'] > parameters['Tgrowth_Af'] > parameters['TmE'] > parameters['TmW']):
        case = 13

    if (parameters['TAEW'] > parameters['Tgrowth_Af'] == parameters['TmE'] > parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 13

    if (parameters['TAEW'] > parameters['Tgrowth_Af'] > parameters['TmE'] == parameters['TmW']):
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 13

    if (parameters['TAEW'] > parameters['Tgrowth_Af'] == parameters['TmE'] == parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 13

    ################

    if (parameters['TAEW'] > parameters['Tgrowth_Af'] > parameters['TmW'] > parameters['TmE']):
        case = 14

    if (parameters['TAEW'] > parameters['Tgrowth_Af'] == parameters['TmW'] > parameters['TmE']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 14

    ################

    if (parameters['TmE'] > parameters['Tgrowth_Af'] > parameters['TmW']):
        case = 15

    if (parameters['TmE'] > parameters['Tgrowth_Af'] == parameters['TmW']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 15


    ################

    if (parameters['TmW'] > parameters['Tgrowth_Af'] > parameters['TmE']):
        case = 16

    if (parameters['TmW'] > parameters['Tgrowth_Af'] == parameters['TmE']):
        Tgrowth_Af += 0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 16

    ################

    if (parameters['TmE'] > parameters['TmW'] > parameters['Tgrowth_Af']):
        case = 17

    if (parameters['TmE'] > parameters['TmW'] == parameters['Tgrowth_Af']):
        Tgrowth_Af += -0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 17

    if (parameters['TmE'] == parameters['TmW'] > parameters['Tgrowth_Af']):
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 17

    if (parameters['TmE'] == parameters['TmW'] == parameters['Tgrowth_Af']):
        Tgrowth_Af += -0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        TmE += 0.00001
        parameters['TmE'] = TmE
        case = 17

    ################

    if (parameters['TmW'] > parameters['TmE'] > parameters['Tgrowth_Af']):
        case = 18

    if (parameters['TmW'] > parameters['TmE'] == parameters['Tgrowth_Af']):
        Tgrowth_Af += -0.00001
        parameters['Tgrowth_Af'] = Tgrowth_Af
        case = 18

    print case
    return case, Tgrowth_Af