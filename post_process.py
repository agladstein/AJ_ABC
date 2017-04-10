from sys import stderr
from sys import argv
import os.path
import linecache
import multiprocessing
import functools
from os import listdir

def combine_files(n, sim_path_results, sim_path_values, f):
    if f.endswith(".summary"):
        if len(f.split("_")) == 4:
            job_id = str(f.split(".")[0].split("_")[2]) + "_" + str(f.split(".")[0].split("_")[3])
        elif len(f.split("_")) == 3:
            job_id = f.split("_")[2].split(".")[0]
        sim_values_file_name = str(sim_path_values) + "/sim_" + str(job_id) + "_values.txt"
        if os.path.exists(sim_values_file_name):
            line = str(job_id) + "\t" + linecache.getline(sim_values_file_name, 2).rstrip('\n') + "\t" + linecache.getline(str(sim_path_results) + "/" + str(f), 2)
            if len(line.split('\t')) == n:
                print str(sim_path_results) + "/" + str(f)
                print sim_values_file_name
                combined_string = str(job_id) + "\t" + linecache.getline(sim_values_file_name, 2).rstrip('\n') + "\t" + linecache.getline(str(sim_path_results) + "/" + str(f), 2)
                return combined_string
            else:
                stderr.write(str(sim_path_results) + " does not have the desired number of summary statistics\n")
        else:
            stderr.write(str(sim_values_file_name) + " does not exist\n")


def main():
    pool = multiprocessing.Pool()

    sim_path = argv[1]
    out_path = argv[2]
    model = argv[3]
    header_file_name = argv[4]

    header_file = open(header_file_name, 'r')
    for line in header_file:
        head = line
    header_file.close()
    n = len(head.split('\t'))

    sim_path_results = str(sim_path)+"/results_sims_AJ_M"+str(model)
    sim_path_values = str(sim_path)+"/sim_values_AJ_M"+str(model)

    files_results = listdir(sim_path_results)

    function_to_map = functools.partial(combine_files, n, sim_path_results, sim_path_values)
    results = pool.imap_unordered(function_to_map, files_results, 10)

    out_file_name = str(out_path)+"/input_ABCtoolbox_M"+str(model)+"_"+str(len(files_results))+".txt"
    fileout = open(out_file_name, 'w')
    fileout.write(head)
    for result in results:
        fileout.write(result)


if __name__ == '__main__':
    main()