from sys import argv
from os import listdir
import os.path
import linecache

def main():

    sim_path = argv[1]
    out_path = argv[2]
    model = argv[3]
    header_file_name = argv[4]

    sim_path_results = str(sim_path)+"/results_sims_AJ_M"+str(model)
    sim_path_values = str(sim_path)+"/sim_values_AJ_M"+str(model)

    header_file = open(header_file_name, 'r')
    for line in header_file:
        head = line
    header_file.close()
    n = len(head.split('\t'))

    files_results = listdir(sim_path_results)

    out_file_name = str(out_path)+"/input_ABCtoolbox_M"+str(model)+"_"+str(len(files_results))+".txt"
    fileout = open(out_file_name, 'w')
    fileout.write(head)

    for f in files_results:
        if f.endswith(".summary"):
            print str(sim_path_results)+"/"+str(f)
            if len(f.split("_"))==4:
                job_id = str(f.split(".")[0].split("_")[2])+"_"+str(f.split(".")[0].split("_")[3])
            elif len(f.split("_"))==3:
                job_id = f.split("_")[2].split(".")[0]
            sim_values_file_name = str(sim_path_values)+"/sim_"+str(job_id)+"_values.txt"
            if os.path.exists(sim_values_file_name):
                line = str(job_id)+"\t"+linecache.getline(sim_values_file_name, 2).rstrip('\n')+"\t"+linecache.getline(str(sim_path_results)+"/"+str(f), 2)
                if len(line.split('\t')) == n:
                    print sim_values_file_name
                    fileout.write(str(job_id)+"\t"+linecache.getline(sim_values_file_name, 2).rstrip('\n')+"\t"+linecache.getline(str(sim_path_results)+"/"+str(f), 2))

                # string = job_id
                # sim_values_file = open(sim_values_file_name)
                # for i, line_values in enumerate(sim_values_file):
                #     if i == 1:
                #         string = str(string)+"\t"+line_values.rstrip('\n')+"\t"
                #         sim_results_file = open(str(sim_path_results)+"/"+str(f))
                #         for j, line_results in enumerate(sim_results_file):
                #             if j == 1:
                #                 string = str(string)+str(line_results)
                #         sim_results_file.close()
                # sim_values_file.close()
            else:
                print str(sim_values_file_name)+" does not exist"

if __name__ == '__main__':
    main()