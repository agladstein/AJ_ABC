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


    daf = random.uniform(0.05, 0.10)

    ####Demographic model
    # population size in Africa
    NAF = float(round(10 ** random.uniform(3.7, 5.0)))
    parameters['NAF'] = NAF
    # print NAF

    # Ancestral population size, before population growth in AF
    # choose ancestral Ne based on being some value smaller than NAF between 1 or 0.1 x.
    # population growth
    Nrat_High = 0.0  # Allow only growth for now
    Nrat_Low = -1.0
    NANC = float(round((10 ** random.uniform(Nrat_Low, Nrat_High)) * NAF))
    parameters['NANC'] = NANC

    NCEU = float(round(10 ** random.uniform(3.0, 5.0)))
    parameters['NCEU'] = NCEU

    NCHB = float(round(10 ** random.uniform(3.0, 5.0)))
    parameters['NCHB'] = NCHB

    # Population size of WAJ
    NWA = float(round(10 ** random.uniform(3.0, 6.7)))
    parameters['NWA'] = NWA

    # Population size of EAJ
    NEA = float(round(10 ** random.uniform(4.0, 6.7)))
    parameters['NEA'] = NEA

    # Population size of AJ before growth
    if (NWA < NEA):
        NAg_High = math.log10(NWA)
    else:
        NAg_High = math.log10(NEA)
    NAg_Low = 2.0
    NAg = float(round(10 ** random.uniform(NAg_Low, NAg_High)))
    parameters['NAg'] = NAg

    # Population size of Jews
    NJ = float(round(10 ** random.uniform(3.0, 6.0)))
    parameters['NJ'] = NJ

    # Population size of Middle Easterns
    NM = float(round(10 ** random.uniform(3.0, 6.0)))
    parameters['NM'] = NM

    # migration rate from Europe to EAJ
    mE_High = 1.0
    mE_Low = 0.0
    mE = random.uniform(mE_Low,mE_High)
    parameters['mE'] = mE

    # migration rate from Europe to WAJ
    mW_High = 1.0
    mW_Low = 0.0
    mW = random.uniform(mW_Low, mW_High)
    parameters['mW'] = mW

    # Time of the instantaneous growth in Africa, before the split between Africans and non Africans
    Tgrowth_Low = 1
    Tgrowth_High = 4100
    Tgrowth_Af = float(randint(Tgrowth_Low, Tgrowth_High))
    parameters['Tgrowth_Af'] = Tgrowth_Af

    # Time of split between YRI and CEU/CHB
    Taf_High = 4100  # 102,500 years using 25 years per generation
    Taf_Low = 1600  # 40,000 years using 25 years per generation.
    Taf = float(randint(Taf_Low, Taf_High))
    parameters['Taf'] = Taf

    # Time of split between Europe and Middle East
    TEM_High = 1200
    TEM_Low = 400
    TEM = float(randint(TEM_Low, TEM_High))
    parameters['TEM'] = TEM

    # Time of split between CEU and CHB
    Teu_as_High = int(Taf) - 1
    Teu_as_Low = int(TEM) + 1
    Teu_as = float(randint(Teu_as_Low, Teu_as_High))
    parameters['Teu_as'] = Teu_as

    # Time of split between Jews and AJ
    TA_High = 36
    TA_Low = 20
    TA = float(randint(TA_Low, TA_High))
    parameters['TA'] = TA

    # Time of split between Jews and Middle East
    TMJ_High = int(TEM) - 1
    TMJ_Low = int(TA) + 1
    TMJ = float(randint(TMJ_Low, TMJ_High))
    parameters['TMJ'] = TMJ

    # Time of split between Eastern and Western AJ
    TAEW_High = (TA) - 2
    TAEW_Low = 2
    TAEW = float(randint(TAEW_Low, TAEW_High))
    parameters['TAEW'] = TAEW

    # Time of migration to EAJ
    TmE_High = int(TAEW) - 1
    TmE_Low = 1
    TmE = float(randint(TmE_Low, TmE_High))
    parameters['TmE'] = TmE

    # Time of migration to WAJ
    TmW_High = int(TAEW) - 1
    TmW_Low = 1
    TmW = float(randint(TmW_Low, TmW_High))
    parameters['TmW'] = TmW

    # Time of growth in AJ
    TAg_High = int(TAEW) - 1
    TAg_Low = 1
    TAg = float(randint(TAg_Low, TAg_High))
    parameters['TAg'] = TAg


    para_out.extend([asc_nb_af])
    para_out.extend([asc_nb_eu])
    para_out.extend([asc_nb_as])
    para_out.extend([daf])
    para_out.extend([math.log10(NAF)])
    para_out.extend([math.log10(NANC)])
    para_out.extend([math.log10(NCEU)])
    para_out.extend([math.log10(NCHB)])
    para_out.extend([math.log10(NWA)])
    para_out.extend([math.log10(NEA)])
    para_out.extend([math.log10(NAg)])
    para_out.extend([math.log10(NJ)])
    para_out.extend([math.log10(NM)])
    para_out.extend([mE])
    para_out.extend([mW])
    para_out.extend([Tgrowth_Af])
    para_out.extend([Taf])
    para_out.extend([TEM])
    para_out.extend([Teu_as])
    para_out.extend([TA])
    para_out.extend([TMJ])
    para_out.extend([TAEW])
    para_out.extend([TmE])
    para_out.extend([TmW])
    para_out.extend([TAg])


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

    daf=random.uniform(0.05,0.10)


    ####Demographic model
    #population size in Africa
    NAF=float(round(10**3.7))
    parameters['NAF']=NAF
    #print NAF

    #Ancestral population size, before population growth in AF
    #choose ancestral Ne based on being some value smaller than NAF between 1 or 0.1 x.
    #population growth
    Nrat_Low=-1.0
    NANC=float(round(10**Nrat_Low*NAF))
    parameters['NANC']=NANC

    NCEU=float(round(10**3.0))
    parameters['NCEU']=NCEU

    NCHB=float(round(10**3.0))
    parameters['NCHB']=NCHB

    # Population size of WAJ
    NWA = float(round(10 **3.0))
    parameters['NWA'] = NWA

    # Population size of EAJ
    NEA = float(round(10 **4.0))
    parameters['NEA'] = NEA

    # Population size of AJ before growth
    if (NWA < NEA):
        NAg_High = math.log10(NWA)
    else:
        NAg_High = NEA
    NAg_Low = 2.0
    NAg = float(round(10 ** random.uniform(NAg_Low, NAg_High)))
    parameters['NAg'] = NAg

    #Population size of Jews
    NJ=float(round(10**3.0))
    parameters['NJ']=NJ

    #Population size of Middle Easterns
    NM=float(round(10**3.0))
    parameters['NM']=NM

    # migration rate from Europe to EAJ
    mE_High = 1
    mE_Low = 0
    mE = random.uniform(mE_Low, mE_High)
    parameters['mE'] = mE

    # migration rate from Europe to WAJ
    mW_High = 1
    mW_Low = 0
    mW = random.uniform(mW_Low, mW_High)
    parameters['mW'] = mW

    #Time of the instantaneous growth in Africa, before the split between Africans and non Africans
    Tgrowth_Low=1
    Tgrowth_High=4100
    Tgrowth_Af=float(randint(Tgrowth_Low,Tgrowth_High))
    parameters['Tgrowth_Af']=Tgrowth_Af

    #Time of split between YRI and CEU/CHB
    Taf_Low=1600				  #40,000 years using 25 years per generation.
    Taf=float(Taf_Low)
    parameters['Taf']=Taf

    #Time of split between Europe and Middle East
    TEM_Low=400
    TEM=float(TEM_Low)
    parameters['TEM']=TEM

    #Time of split between CEU and CHB
    Teu_as_Low=int(TEM)+1
    Teu_as=float(Teu_as_Low)
    parameters['Teu_as']=Teu_as

    #Time of split between Jews and AJ
    TA_Low=20
    TA=float(TA_Low)
    parameters['TA']=TA

    #Time of split between Jews and Middle East
    TMJ_Low=int(TA)+1
    TMJ=float(TMJ_Low)
    parameters['TMJ']=TMJ

    # Time of split between Eastern and Western AJ
    TAEW_Low = 2
    TAEW = float(TAEW_Low)
    parameters['TAEW'] = TAEW

    # Time of migration to EAJ
    TmE_High = int(TAEW) - 1
    TmE_Low = 1
    TmE = float(randint(TmE_Low, TmE_High))
    parameters['TmE'] = TmE

    # Time of migration to WAJ
    TmW_High = int(TAEW) - 1
    TmW_Low = 1
    TmW = float(randint(TmW_Low, TmW_High))
    parameters['TmW'] = TmW

    # Time of growth in AJ
    TAg_High = int(TAEW) - 1
    TAg_Low = 1
    TAg = float(randint(TAg_Low, TAg_High))
    parameters['TAg'] = TAg


    para_out.extend([asc_nb_af])
    para_out.extend([asc_nb_eu])
    para_out.extend([asc_nb_as])
    para_out.extend([daf])
    para_out.extend([math.log10(NAF)])
    para_out.extend([math.log10(NANC)])
    para_out.extend([math.log10(NCEU)])
    para_out.extend([math.log10(NCHB)])
    para_out.extend([math.log10(NWA)])
    para_out.extend([math.log10(NEA)])
    para_out.extend([math.log10(NAg)])
    para_out.extend([math.log10(NJ)])
    para_out.extend([math.log10(NM)])
    para_out.extend([mE])
    para_out.extend([mW])
    para_out.extend([Tgrowth_Af])
    para_out.extend([Taf])
    para_out.extend([TEM])
    para_out.extend([Teu_as])
    para_out.extend([TA])
    para_out.extend([TMJ])
    para_out.extend([TAEW])
    para_out.extend([TmE])
    para_out.extend([TmW])
    para_out.extend([TAg])


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


    daf=random.uniform(0.05,0.10)


    ####Demographic model
    #population size in Africa
    NAF=float(round(10**5))
    parameters['NAF']=NAF

    #Ancestral population size, before population growth in AF
    #choose ancestral Ne based on being some value smaller than NAF between 1 or 0.1 x.
    #population growth
    Nrat_High=0.0   #Allow only growth for now
    NANC=float(round(10**Nrat_High*NAF))
    parameters['NANC']=NANC

    NCEU=float(round(10**5.0))
    parameters['NCEU']=NCEU

    NCHB=float(round(10**5.0))
    parameters['NCHB']=NCHB

    # Population size of WAJ
    NWA = float(round(10 **6.7))
    parameters['NWA'] = NWA

    # Population size of EAJ
    NEA = float(round(10 **6.7))
    parameters['NEA'] = NEA

    # Population size of AJ before growth
    if (NWA < NEA):
        NAg_High = math.log10(NWA)
    else:
        NAg_High = NEA
    NAg_Low = 2.0
    NAg = float(round(10 ** random.uniform(NAg_Low, NAg_High)))
    parameters['NAg'] = NAg

    #Population size of Jews
    NJ=float(round(10**6.0))
    parameters['NJ']=NJ

    #Population size of Middle Easterns
    NM=float(round(10**6.0))
    parameters['NM']=NM

    # migration rate from Europe to EAJ
    mE_High = 1.0
    mE_Low = 0.0
    mE = random.uniform(mE_Low,mE_High)
    parameters['mE'] = mE

    # migration rate from Europe to WAJ
    mW_High = 1.0
    mW_Low = 0.0
    mW = random.uniform(mW_Low, mW_High)
    parameters['mW'] = mW

    #Time of the instantaneous growth in Africa, before the split between Africans and non Africans
    Tgrowth_Low=1
    Tgrowth_High=4100
    Tgrowth_Af=float(randint(Tgrowth_Low,Tgrowth_High))
    parameters['Tgrowth_Af']=Tgrowth_Af

    #Time of split between YRI and CEU/CHB
    Taf_High=4100				  #102,500 years using 25 years per generation
    Taf=float(Taf_High)
    parameters['Taf']=Taf

    #Time of split between Europe and Middle East
    TEM_High=1200
    TEM=float(TEM_High)
    parameters['TEM']=TEM

    #Time of split between CEU and CHB
    Teu_as_High=int(Taf)-1
    Teu_as=float(Teu_as_High)
    parameters['Teu_as']=Teu_as

    #Time of split between Jews and AJ
    TA_High=36
    TA=float(TA_High)
    parameters['TA']=TA

    #Time of split between Jews and Middle East
    TMJ_High=int(TEM)-1
    TMJ=float(TMJ_High)
    parameters['TMJ']=TMJ

    # Time of split between Eastern and Western AJ
    TAEW_High = (TA) - 2
    TAEW = float(TAEW_High)
    parameters['TAEW'] = TAEW

    # Time of migration to EAJ
    TmE_High = int(TAEW) - 1
    TmE_Low = 1
    TmE = float(randint(TmE_Low, TmE_High))
    parameters['TmE'] = TmE

    # Time of migration to WAJ
    TmW_High = int(TAEW) - 1
    TmW_Low = 1
    TmW = float(randint(TmW_Low, TmW_High))
    parameters['TmW'] = TmW

    # Time of growth in AJ
    TAg_High = int(TAEW) - 1
    TAg_Low = 1
    TAg = float(randint(TAg_Low, TAg_High))
    parameters['TAg'] = TAg


    para_out.extend([asc_nb_af])
    para_out.extend([asc_nb_eu])
    para_out.extend([asc_nb_as])
    para_out.extend([daf])
    para_out.extend([math.log10(NAF)])
    para_out.extend([math.log10(NANC)])
    para_out.extend([math.log10(NCEU)])
    para_out.extend([math.log10(NCHB)])
    para_out.extend([math.log10(NWA)])
    para_out.extend([math.log10(NEA)])
    para_out.extend([math.log10(NAg)])
    para_out.extend([math.log10(NJ)])
    para_out.extend([math.log10(NM)])
    para_out.extend([mE])
    para_out.extend([mW])
    para_out.extend([Tgrowth_Af])
    para_out.extend([Taf])
    para_out.extend([TEM])
    para_out.extend([Teu_as])
    para_out.extend([TA])
    para_out.extend([TMJ])
    para_out.extend([TAEW])
    para_out.extend([TmE])
    para_out.extend([TmW])
    para_out.extend([TAg])


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