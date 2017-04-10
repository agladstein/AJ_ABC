from sys import argv
import os.path
import linecache
import multiprocessing
import glob

def combine_files(f, fileout, n, sim_path_results, sim_path_values):
    print str(sim_path_results) + "/" + str(f)
    if len(f.split("_")) == 4:
        job_id = str(f.split(".")[0].split("_")[2]) + "_" + str(f.split(".")[0].split("_")[3])
    elif len(f.split("_")) == 3:
        job_id = f.split("_")[2].split(".")[0]
    sim_values_file_name = str(sim_path_values) + "/sim_" + str(job_id) + "_values.txt"
    if os.path.exists(sim_values_file_name):
        line = str(job_id) + "\t" + linecache.getline(sim_values_file_name, 2).rstrip('\n') + "\t" + linecache.getline(
            str(sim_path_results) + "/" + str(f), 2)
        if len(line.split('\t')) == n:
            print sim_values_file_name
            fileout.write(
                str(job_id) + "\t" + linecache.getline(sim_values_file_name, 2).rstrip('\n') + "\t" + linecache.getline(
                    str(sim_path_results) + "/" + str(f), 2))
    else:
        print str(sim_values_file_name) + " does not exist"


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

    # Find all .summary files in results directory
    files_results = glob.glob('{}/results_sims_AJ_M{}/*.summary'.format(sim_path, model))

    out_file_name = str(out_path)+"/input_ABCtoolbox_M"+str(model)+"_"+str(len(files_results))+".txt"
    fileout = open(out_file_name, 'w')
    fileout.write(head)


    results = pool.imap_unordered(combine_files(files_results, fileout, n, sim_path_results, sim_path_values), files_results, 10)



if __name__ == '__main__':
    main()