import random

def run_sim(parameters,length,chr_number,total,total_naf,total_nas,total_neu,nJ,nM,nA,seed_option):

    mu=2.5e-8
    rho=1e-8

    NAF=float(parameters['NAF'])
    NANC=float(parameters['NANC'])
    NCEU=float(parameters['NCEU'])
    NCHB=float(parameters['NCHB'])
    NA=float(parameters['NA'])
    NAg = float(parameters['NAg'])
    NJ=float(parameters['NJ'])
    NM=float(parameters['NM'])

    m=float(parameters['m'])

    Tgrowth_Af=float(parameters['Tgrowth_Af'])
    Taf=float(parameters['Taf'])
    TEM=float(parameters['TEM'])
    Teu_as=float(parameters['Teu_as'])
    TA=float(parameters['TA'])
    TMJ=float(parameters['TMJ'])
    Tm=float(parameters['Tm'])
    TAg = float(parameters['TAg'])

    macs_theta=float(mu*4*NANC)
    macs_rho=float(rho*4*NANC)
    scaled_NAF=float(NAF/NANC)
    scaled_NANC=float(NANC/NANC)
    scaled_NCEU=float(NCEU/NANC)
    scaled_NCHB=float(NCHB/NANC)
    scaled_NA=float(NA/NANC)
    scaled_NAg = float(NAg / NANC)
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
    scaled_TAg = float(TAg / (4 * NANC))

    adjust = random.uniform(0.0000010, 0.0000099)
    if scaled_Tm == scaled_TAg:
        scaled_Tm = float(scaled_Tm + adjust)

    adjust = random.uniform(0.0000000010, 0.0000000099)
    scaled_Tm2 = float(scaled_Tm + adjust)


    if seed_option > int(0):
        macs_args = ['./bin/macs', str(total), str(length), '-t', str(macs_theta), '-s', str(seed_option), '-r',
                     str(macs_rho), '-h', '1e5',
                     '-R', 'genetic_map_b37/genetic_map_GRCh37_chr' + str(chr_number) + '.txt.macshs', '-I', '6',
                     str(total_naf), str(total_nas), str(total_neu), str(nJ), str(nM), str(nA), '-n', '1',
                     str(scaled_NAF), '-n', '2', str(scaled_NCHB), '-n', '3', str(scaled_NCEU), '-n', '4',
                     str(scaled_NJ), '-n', '5', str(scaled_NM), '-n', '6', str(scaled_NA)]
    else:
        macs_args = ['./bin/macs', str(total), str(length), '-t', str(macs_theta), '-r', str(macs_rho), '-h', '1e5',
                     '-R', 'genetic_map_b37/genetic_map_GRCh37_chr' + str(chr_number) + '.txt.macshs', '-I', '6',
                     str(total_naf), str(total_nas), str(total_neu), str(nJ), str(nM), str(nA), '-n', '1',
                     str(scaled_NAF), '-n', '2', str(scaled_NCHB), '-n', '3', str(scaled_NCEU), '-n', '4',
                     str(scaled_NJ), '-n', '5', str(scaled_NM), '-n', '6', str(scaled_NA)]

    ej = [[str(scaled_TA),'6', '4'], [str(scaled_TMJ), '5', '4'], [str(scaled_TEM), '4', '3'], [str(scaled_Teu_as), '3', '2'], [str(scaled_Taf), '2', '1']]
    em = [[str(scaled_Tm), '6', '3', str(scaled_m)], [str(scaled_Tm2), '6', '3', '0']]
    en = [[str(scaled_Tgrowth_Af), '1', str(scaled_NANC)], [str(scaled_TAg), '6', str(scaled_NAg)]]
    vars = {"-ej": ej, "-em": em, "-en": en}

    # pull out all the time sensitive vars '-e*' and order them
    seasons = []
    for thing in vars.keys():
        '''
        make a list of all the vars with a -e 
        sort them by -e such that the smallest -e value is first
        check out lamda sorts(?)
        '''

        if thing.startswith('-e'):  # makes the list of -t vars.
            if thing not in seasons:
                for i in range(len(vars[thing])):
                    into = [thing]  # thing is the key, I just need key to not be used later so we're going with thing.
                    for j in range(len(vars[thing][i])):
                        into.append(vars[thing][i][j])
                    if into[1] != '':
                        seasons.append(into)

    # sorts by time
    # extend the data to mac_args in the correct order.
    seasons = sorted(seasons, key=lambda x: float(x[1]))  # sorts the vars in order of time
    for i in range(len(seasons)):
        seasons[i][1] = str((seasons[i][1]))
        macs_args.extend(seasons[i])

    return macs_args

