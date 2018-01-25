import macsSwig
import random


def run_sim(parameters, length, chr_number, total, total_naf, total_nas, total_neu, nJ, nM, nEA, nWA, seed_option):

    mu = 2.5e-8
    rho = 1e-8

    NAF = float(parameters['NAF'])
    NANC = float(parameters['NANC'])
    NCEU = float(parameters['NCEU'])
    NCHB = float(parameters['NCHB'])
    NWA = float(parameters['NWA'])
    NEA = float(parameters['NEA'])
    NAg = float(parameters['NAg'])

    NJ = float(parameters['NJ'])
    NM = float(parameters['NM'])

    mE = float(parameters['mE'])
    mW = float(parameters['mW'])

    Tgrowth_Af = float(parameters['Tgrowth_Af'])
    Taf = float(parameters['Taf'])
    TEM = float(parameters['TEM'])
    Teu_as = float(parameters['Teu_as'])
    TAEW = float(parameters['TAEW'])
    TA = float(parameters['TA'])
    TMJ = float(parameters['TMJ'])
    TmE = float(parameters['TmE'])
    TmW = float(parameters['TmW'])
    TAg = float(parameters['TAg'])

    macs_theta = float(mu * 4 * NANC)
    macs_rho = float(rho * 4 * NANC)
    scaled_NAF = float(NAF / NANC)
    scaled_NANC = float(NANC / NANC)
    scaled_NCEU = float(NCEU / NANC)
    scaled_NCHB = float(NCHB / NANC)
    scaled_NWA = float(NWA / NANC)
    scaled_NEA = float(NEA / NANC)
    scaled_NAg = float(NAg / NANC)
    scaled_NJ = float(NJ / NANC)
    scaled_NM = float(NM / NANC)

    scaled_mE = float(4 * mE * NANC)
    scaled_mW = float(4 * mW * NANC)

    scaled_Tgrowth_Af = float(Tgrowth_Af / (4 * NANC))
    scaled_Taf = float(Taf / (4 * NANC))
    scaled_TEM = float(TEM / (4 * NANC))
    scaled_Teu_as = float(Teu_as / (4 * NANC))
    scaled_TAEW = float(TAEW / (4 * NANC))
    scaled_TA = float(TA / (4 * NANC))
    scaled_TMJ = float(TMJ / (4 * NANC))
    scaled_TmE = float(TmE / (4 * NANC))
    scaled_TmW = float(TmW / (4 * NANC))
    scaled_TAg = float(TAg / (4 * NANC))

    print "scaled_TmE:", scaled_TmE
    print "scaled_TmW:", scaled_TmW
    print "scaled_TAg:", scaled_TAg
    adjust = random.uniform(0.0000010, 0.0000099)
    if scaled_TmE == scaled_TmW:
        scaled_TmE = float(scaled_TmE + adjust)
        print "new scaled_TmE:", scaled_TmE
    adjust = random.uniform(0.0000010, 0.0000099)
    if scaled_TmE == scaled_TAg or scaled_TmW == scaled_TAg:
        scaled_TAg = float(scaled_TAg + adjust)
        print "new scaled_TAg:", scaled_TAg

    adjust = random.uniform(0.00000000010, 0.00000000099)
    scaled_TmE2 = float(scaled_TmE + adjust)
    adjust = random.uniform(0.00000000010, 0.00000000099)
    scaled_TmW2 = float(scaled_TmW + adjust)
    adjust = random.uniform(0.00000000010, 0.00000000099)
    scaled_TAg2 = float(scaled_TAg + adjust)


    if seed_option > int(0):
        macs_args = ['./bin/macs', str(total), str(length), '-t', str(macs_theta), '-s', str(seed_option), '-r',
                 str(macs_rho), '-h', '1e5', '-R', 'genetic_map_b37/genetic_map_GRCh37_chr' + str(chr_number) + '.txt.macshs', '-I', '7',
                 str(total_naf), str(total_nas), str(total_neu), str(nJ), str(nM), str(nEA), str(nWA), '-n',
                 '1', str(scaled_NAF), '-n', '2', str(scaled_NCHB), '-n', '3', str(scaled_NCEU), '-n', '4',
                 str(scaled_NJ), '-n', '5', str(scaled_NM), '-n', '6', str(scaled_NEA), '-n', '7',
                 str(scaled_NWA)]
    else:
        macs_args = ['./bin/macs', str(total), str(length), '-t', str(macs_theta), '-r',
                 str(macs_rho), '-h', '1e5', '-R', 'genetic_map_b37/genetic_map_GRCh37_chr' + str(chr_number) + '.txt.macshs', '-I', '7',
                 str(total_naf), str(total_nas), str(total_neu), str(nJ), str(nM), str(nEA), str(nWA), '-n',
                 '1', str(scaled_NAF), '-n', '2', str(scaled_NCHB), '-n', '3', str(scaled_NCEU), '-n', '4',
                 str(scaled_NJ), '-n', '5', str(scaled_NM), '-n', '6', str(scaled_NEA), '-n', '7',
                 str(scaled_NWA)]

    ej = [[str(scaled_TAEW), '6', '7'], [str(scaled_TA), '7', '4'], [str(scaled_TMJ), '5', '4'],
          [str(scaled_TEM), '4', '3'], [str(scaled_Teu_as), '3', '2'], [str(scaled_Taf), '2', '1']]
    em = [[str(scaled_TmW), '7', '3', str(scaled_mW)], [str(scaled_TmW2), '7', '3', '0'],
          [str(scaled_TmE), '6', '3', str(scaled_mE)], [str(scaled_TmE2), '6', '3', '0']]
    en = [[str(scaled_Tgrowth_Af), '1', str(scaled_NANC)], [str(scaled_TAg), '7', str(scaled_NAg)],
          [str(scaled_TAg2), '6', str(scaled_NAg)]]
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

