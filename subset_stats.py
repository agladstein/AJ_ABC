from sys import argv
import csv


def get_corr_stats(stats_file):
    '''This option should be used to remove highly correlated statistics from ABCtoolbox input.
    stats_file should be the log file after running ABCestimator with the option pruneCorrelatedStats'''

    print 'Reading summary statistics found in '+stats_file+'\n'
    stats = []
    with open(stats_file) as prune_file:
        for line in prune_file:
            if 'Removing' in line:
                columns = line.split('\'')
                stats.append(columns[1])
    return stats

def get_power_stats(stats_file):
    ''' This option should be used to only keep parameter values and statistics with high power from ABCtoolbox input. 
    stats_file should be the output file *greedySearchForBestStatisticsForModelChoice.txt after running ABCestimator with the option findStatsModelChoice'''

    print 'Reading summary statistics found in ' + stats_file + '\n'
    # with open(stats_file) as keep_file:
    #     for line in keep_file:
    #
    # return stats

def remove_stats(stats, input_file):
    print 'Pruning correlated summary statistics \n'
    with open(input_file) as prune_input_file:
        reader = csv.reader(prune_input_file, delimiter='\t')
        head = reader.next()
    rm_index = []
    for i, item in enumerate(head):
        if item in stats:
            continue
        else:
            rm_index.append(i)
    return rm_index


def keep_stats(stats, input_file):
    print 'Keeping higher power summary statistics \n'

    # return out_string


def get_string(input_file, rm_index):
    out_string = ''
    with open(input_file) as prune_input_file:
        for line in prune_input_file:
            line = line.strip('\n')
            columns = line.split('\t')
            for j in rm_index:
                out_string = out_string + columns[j] + '\t'
            out_string = out_string + '\n'
    return out_string


def main():
    stats_file = argv[1]
    input_file = argv[2]
    option = argv[3] # remove or keep

    if option == "remove":
        stats = get_corr_stats(stats_file)
        index = remove_stats(stats, input_file)
        out_file_name = 'pruneCorStats_'+input_file
    elif option == "keep":
        stats = get_power_stats(stats_file)
        index = keep_stats(stats, input_file)
        out_file_name = 'keepPowerStats_' + input_file
    else:
        print "You must specify remove or keep!"
        return

    out_string = get_string(input_file, index)

    fileout = open(out_file_name, 'w')
    fileout.write(out_string)
    fileout.close()

if __name__ == '__main__':
    main()