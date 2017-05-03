from sys import stderr
from sys import argv
import os.path
import linecache
import multiprocessing
import functools
from os import listdir
import csv
import re


def get_file_name(f, sim_path_results, sim_path_values):
    # isolate the job id to get the file names in the sim_values directory
    if len(f.split("_")) == 4:
        job_id = str(f.split(".")[0].split("_")[2]) + "_" + str(f.split(".")[0].split("_")[3])
    elif len(f.split("_")) == 3:
        job_id = f.split("_")[2].split(".")[0]
    sim_values_file_name = str(sim_path_values) + "/sim_" + str(job_id) + "_values.txt"
    results_file_name = str(sim_path_results) + "/" + str(f)
    return job_id, results_file_name, sim_values_file_name


def combine_files(n, sim_path_results, sim_path_values, f):
    if f.endswith(".summary"):
        job_id, results_file_name, sim_values_file_name = get_file_name(f, sim_path_results, sim_path_values)
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


def combine_same_header(head, sim_path_results, sim_path_values, f):
    if f.endswith(".summary"):
        job_id, results_file_name, sim_values_file_name = get_file_name(f, sim_path_results, sim_path_values)

        # If the matching sim_values file exists get the headers and the lines
        if os.path.exists(sim_values_file_name):
            with open(sim_values_file_name) as sim_csvfile:
                reader = csv.reader(sim_csvfile, delimiter='\t')
                headers = reader.next()
                line = reader.next()
            with open(results_file_name) as results_csvfile:
                reader = csv.reader(results_csvfile, delimiter='\t')
                headers.extend(reader.next())
                line.extend(reader.next())

            # If the results file has IBD statistics (no files missing stats)
            if re.findall(r"(IBD)", str(headers)):
                line_same = str(job_id)
                for i, header in enumerate(headers):
                    # Only use stats that are included in the input header file
                    if header not in head:
                        continue
                    else:
                        line_same = line_same + '\t' + line[i]
                line_same = line_same + '\n'
                # check the resulting line has the same length as the input header
                if len(line_same.split('\t')) == (len(head)):
                    return line_same
                else:
                    stderr.write(str(sim_values_file_name) + " does not have the desired parameters or summary statistics\n")
                    return False
            else:
                stderr.write(str(sim_values_file_name) + " does not have IBD\n")
                return False
        else:
            stderr.write(str(sim_values_file_name) + " does not exist\n")
            return False


def combine_duplicate_header(head, sim_path_results, sim_path_values, to_duplicate, to_remove, f):
    # Only get results files that end in .summary
    if f.endswith(".summary"):
        job_id, results_file_name, sim_values_file_name = get_file_name(f, sim_path_results, sim_path_values)

        # If the matching sim_values file exists get the headers and the lines
        if os.path.exists(sim_values_file_name):
            with open(sim_values_file_name) as sim_csvfile:
                reader = csv.reader(sim_csvfile, delimiter='\t')
                headers = reader.next()
                line = reader.next()
            with open(results_file_name) as results_csvfile:
                reader = csv.reader(results_csvfile, delimiter='\t')
                headers.extend(reader.next())
                line.extend(reader.next())

            # If the results file has IBD statistics (no files missing stats)
            if re.findall(r"(IBD)", str(headers)):
                line_duplicate = str(job_id)
                for i, header in enumerate(headers):
                    # Only use stats that are included in the input header file
                    if header not in head:
                        continue
                    if header in to_remove:
                        continue
                    if header in to_duplicate:
                        assert isinstance(line[i], basestring)
                        if re.search(r"(AA)", header):
                            line_duplicate = line_duplicate + '\t' + line[i] + '\t' + line[i] + '\t' + line[i]
                        else:
                            line_duplicate = line_duplicate + '\t' + line[i] + '\t' + line[i]
                    else:
                        line_duplicate = line_duplicate + '\t' + line[i]
                line_duplicate = line_duplicate + '\n'
                return line_duplicate
            else:
                stderr.write(results_file_name + " does not have IBD\n")
                return False
        else:
            stderr.write(str(sim_values_file_name) + " does not exist\n")
            return False


def make_duplicate_header(head, to_duplicate, to_remove):
    duplicated_header = []
    for item in head:
        if item in to_remove:
            continue
        if item in to_duplicate:
            duplicated_header.extend(to_duplicate[item])
        else:
            duplicated_header.append(item)
    return duplicated_header


def main():
    pool = multiprocessing.Pool()

    sim_path = argv[1]
    results_path = argv[2]
    out_path = argv[3]
    model = argv[4]
    header_file_name = argv[5]
    combine_function = argv[6] # original, same, or duplicate

    with open(header_file_name) as f:
        reader = csv.reader(f, delimiter='\t')
        head = reader.next()
    n = len(head)

    sim_path_results = str(results_path) + "/results_sims_AJ_M" + str(model)
    sim_path_values = str(sim_path) + "/sim_values_AJ_M" + str(model)

    files_results = listdir(sim_path_results)
    out_file_name = str(out_path) + "/input_ABCtoolbox_M" + str(model) + "_" + str(len(files_results)) + ".txt"

    if combine_function == "original":
        function_to_map = functools.partial(combine_files, n, sim_path_results, sim_path_values)
        fileout = open(out_file_name, 'w')
        fileout.write("\t".join(head) + "\n")
    elif combine_function == "same":
        function_to_map = functools.partial(combine_same_header, head, sim_path_results, sim_path_values)
        fileout = open(out_file_name, 'w')
        fileout.write("\t".join(head) + "\n")
    elif combine_function == "duplicate":
        if model == str(2):
            to_duplicate = {
                'IBD_mean_AA': ['IBD_mean_eAeA', 'IBD_mean_wAwA', 'IBD_mean_eAwA'],
                'IBD_mean_AE': ['IBD_mean_eAE', 'IBD_mean_wAE'],
                'IBD_mean_AJ': ['IBD_mean_eAJ', 'IBD_mean_wAJ'],
                'IBD_mean_AM': ['IBD_mean_eAM', 'IBD_mean_wAM'],
                'IBD_median_AA': ['IBD_median_eAeA', 'IBD_median_wAwA', 'IBD_median_eAwA'],
                'IBD_median_AE': ['IBD_median_eAE', 'IBD_median_wAE'],
                'IBD_median_AJ': ['IBD_median_eAJ', 'IBD_median_wAJ'],
                'IBD_median_AM': ['IBD_median_eAM', 'IBD_median_wAM'],
                'IBD30_mean_AA': ['IBD30_mean_eAeA', 'IBD30_mean_wAwA', 'IBD30_mean_eAwA'],
                'IBD30_mean_AE': ['IBD30_mean_eAE', 'IBD30_mean_wAE'],
                'IBD30_mean_AJ': ['IBD30_mean_eAJ', 'IBD30_mean_wAJ'],
                'IBD30_mean_AM': ['IBD30_mean_eAM', 'IBD30_mean_wAM'],
                'IBD30_median_AA': ['IBD30_median_eAeA', 'IBD30_median_wAwA', 'IBD30_median_eAwA'],
                'IBD30_median_AE': ['IBD30_median_eAE', 'IBD30_median_wAE'],
                'IBD30_median_AJ': ['IBD30_median_eAJ', 'IBD30_median_wAJ'],
                'IBD30_median_AM': ['IBD30_median_eAM', 'IBD30_median_wAM'],
                'SegS_A_ASC': ['SegS_EA_ASC', 'SegS_WA_ASC'],
                'Sing_A_ASC': ['Sing_EA_ASC', 'Sing_WA_ASC'],
                'Dupl_A_ASC': ['Dupl_EA_ASC', 'Dupl_WA_ASC'],
                'Pi_A_ASC': ['Pi_EA_ASC', 'Pi_WA_ASC'],
                'TajD_A_ASC': ['TajD_EA_ASC', 'TajD_WA_ASC'],
                'FST_AEu_ASC': ['FST_EAEu_ASC', 'FST_WAEu_ASC'],
                'FST_AJ_ASC': ['FST_EAJ_ASC', 'FST_WAJ_ASC'],
                'FST_AM_ASC': ['FST_EAM_ASC', 'FST_WAM_ASC']
            }
            to_remove = ['IBD_num_AA', 'IBD_num_AE', 'IBD_num_AJ', 'IBD_num_AM', 'IBD_var_AA', 'IBD_var_AE',
                         'IBD_var_AJ', 'IBD_var_AM', 'IBD30_num_AA', 'IBD30_num_AE', 'IBD30_num_AJ', 'IBD30_num_AM',
                         'IBD30_var_AA', 'IBD30_var_AE', 'IBD30_var_AJ', 'IBD30_var_AM', 'FST_eAwA_ASC']
            duplicated_header = make_duplicate_header(head, to_duplicate, to_remove)
            fileout = open(out_file_name, 'w')
            fileout.write("\t".join(duplicated_header) + "\n")
            function_to_map = functools.partial(combine_duplicate_header, head, sim_path_results, sim_path_values,
                                                to_duplicate, to_remove)
        else:
            print "You should only use the duplicate option with Model 1"
            return False
    else:
        print "You didn't specificy how you want to combine the files."
        return False
    results = pool.imap_unordered(function_to_map, files_results, 10)

    for result in results:
        if result:
            fileout.write(result)


if __name__ == '__main__':
    main()
