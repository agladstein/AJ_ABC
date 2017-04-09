import glob
import multiprocessing
import sys


def detect_broken_file(file_name):
    with open(file_name, 'r') as f:
        first_line = f.readline()
        count = first_line.count('IBD_var_EE')
        return file_name, count


def main():
    sim_path = sys.argv[1]
    model = sys.argv[2]
    pool = multiprocessing.Pool()

    listing = glob.glob('{}/results_sims_AJ_M{}/*.summary'.format(sim_path, model))
    results = pool.imap_unordered(detect_broken_file, listing, 1000)
    for file_name, count in results:
        if count == 1:
            # All good. Nothing to see here. Move along.
            sys.stderr.write('.')
        elif count == 0:
            # Not there. Weird.
            sys.stderr.write('\nGermline not run? {}\n'.format(file_name))
        elif count > 2:
            # Whoah Nelly. Throw a fit.
            sys.stderr.write('\nWhat the hell? {}\n'.format(file_name))
        elif count == 2:
            print(file_name)


if __name__ == '__main__':
    main()
