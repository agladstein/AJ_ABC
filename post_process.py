from sys import argv
from os import listdir
import os.path
import linecache

def main():

    sim_path = argv[1]
    out_path = argv[2]
    model = argv[3]

    sim_path_results = str(sim_path)+"/results_sims_AJ_M"+str(model)
    sim_path_values = str(sim_path)+"/sim_values_AJ_M"+str(model)

    if model == "1":
        head_param = 'Sim\tAsc_NAF\tAsc_NEU\tAsc_NCHB\tdaf\tLog10_NAF\tLog10_NANC\tLog10_NCEU\tLog10_NCHB\tLog10_NA\tLog10_NJ\tLog10_NM\trA\trMJ\tm\tTgrowth_Af\tTAF\tTEM\tTeu_as\tTA\tTMJ\tTm\t'
        head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'
        head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'
        head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'
        head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'
        head = head + 'IBD_mean_AA\tIBD_mean_JJ\tIBD_mean_MM\tIBD_mean_EE\tIBD_mean_AE\tIBD_mean_AJ\tIBD_mean_AM\tIBD_mean_JM\tIBD_mean_JE\tIBD_mean_ME\t'
        head = head + 'IBD_median_AA\tIBD_median_JJ\tIBD_median_MM\tIBD_median_EE\tIBD_median_AE\tIBD_median_AJ\tIBD_median_AM\tIBD_median_JM\tIBD_median_JE\tIBD_median_ME\t'
        head = head + 'IBD_num_AA\tIBD_num_JJ\tIBD_num_MM\tIBD_num_EE\tIBD_num_AE\tIBD_num_AJ\tIBD_num_AM\tIBD_num_JM\tIBD_num_JE\tIBD_num_ME\t'
        head = head + 'IBD_var_AA\tIBD_var_JJ\tIBD_var_MM\tIBD_var_EE\tIBD_var_AE\tIBD_var_AJ\tIBD_var_AM\tIBD_var_JM\tIBD_var_JE\tIBD_var_ME\t'
        head = head + 'IBD30_mean_AA\tIBD30_mean_JJ\tIBD30_mean_MM\tIBD30_mean_EE\tIBD30_mean_AE\tIBD30_mean_AJ\tIBD30_mean_AM\tIBD30_mean_JM\tIBD30_mean_JE\tIBD30_mean_ME\t'
        head = head + 'IBD30_median_AA\tIBD30_median_JJ\tIBD30_median_MM\tIBD30_median_EE\tIBD30_median_AE\tIBD30_median_AJ\tIBD30_median_AM\tIBD30_median_JM\tIBD30_median_JE\tIBD30_median_ME\t'
        head = head + 'IBD30_num_AA\tIBD30_num_JJ\tIBD30_num_MM\tIBD_num_EE\tIBD30_num_AE\tIBD30_num_AJ\tIBD30_num_AM\tIBD30_num_JM\tIBD30_num_JE\tIBD30_num_ME\t'
        head = head + 'IBD30_var_AA\tIBD30_var_JJ\tIBD30_var_MM\tIBD_var_EE\tIBD30_var_AE\tIBD30_var_AJ\tIBD30_var_AM\tIBD30_var_JM\tIBD30_var_JE\tIBD30_var_ME\t'
        head = head + 'SegS_Af_ASC\tSing_Af_ASC\tDupl_Af_ASC\tPi_Af_ASC\tTajD_Af_ASC\t'
        head = head + 'SegS_Eu_ASC\tSing_Eu_ASC\tDupl_Eu_ASC\tPi_Eu_ASC\tTajD_Eu_ASC\t'
        head = head + 'SegS_As_ASC\tSing_As_ASC\tDupl_As_ASC\tPi_As_ASC\tTajD_As_ASC\t'
        head = head + 'SegS_J_ASC\tSing_J_ASC\tDupl_J_ASC\tPi_J_ASC\tTajD_J_ASC\t'
        head = head + 'SegS_M_ASC\tSing_M_ASC\tDupl_M_ASC\tPi_M_ASC\tTajD_M_ASC\t'
        head = head + 'SegS_A_ASC\tSing_A_ASC\tDupl_A_ASC\tPi_A_ASC\tTajD_A_ASC\t'
        head = head + 'FST_AfEu_ASC\tFST_AfAs_ASC_m\tFST_EuAs_ASC\t'
        head = head + 'FST_AEu_ASC\tFST_AJ_ASC\tFST_AM_ASC\tFST_MJ_ASC\n'

    elif model == "2":
        head_param = 'Sim\tAsc_NAF\tAsc_NEU\tAsc_NCHB\tdaf\tLog10_NAF\tLog10_NANC\tLog10_NCEU\tLog10_NCHB\tLog10_NWA\tLog10_NEA\tLog10_NJ\tLog10_NM\trWA\trEA\trMJ\tm\tTgrowth_Af\tTAF\tTEM\tTeu_as\tTA\tTMJ\tTAEW\tTm\t'
        head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'
        head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'
        head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'
        head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'
        head = head + 'IBD_mean_eAeA\tIBD_mean_wAwA\tIBD_mean_JJ\tIBD_mean_MM\tIBD_mean_EE\tIBD_mean_eAwA\tIBD_mean_eAE\tIBD_mean_wAE\tIBD_mean_eAJ\tIBD_mean_wAJ\tIBD_mean_eAM\tIBD_mean_wAM\tIBD_mean_JM\tIBD_mean_JE\tIBD_mean_ME\t'
        head = head + 'IBD_median_eAeA\tIBD_median_wAwA\tIBD_median_JJ\tIBD_median_MM\tIBD_median_EE\tIBD_median_eAwA\tIBD_median_eAE\tIBD_median_wAE\tIBD_median_eAJ\tIBD_median_wAJ\tIBD_median_eAM\tIBD_median_wAM\tIBD_median_JM\tIBD_median_JE\tIBD_median_ME\t'
        head = head + 'IBD_num_eAeA\tIBD_num_wAwA\tIBD_num_JJ\tIBD_num_MM\tIBD_num_EE\tIBD_num_eAwA\tIBD_num_eAE\tIBD_num_wAE\tIBD_num_eAJ\tIBD_num_wAJ\tIBD_num_eAM\tIBD_num_wAM\tIBD_num_JM\tIBD_num_JE\tIBD_num_ME\t'
        head = head + 'IBD_var_eAeA\tIBD_var_wAwA\tIBD_var_JJ\tIBD_var_MM\tIBD_var_EE\tIBD_var_eAwA\tIBD_var_eAE\tIBD_var_wAE\tIBD_var_eAJ\tIBD_var_wAJ\tIBD_var_eAM\tIBD_var_wAM\tIBD_var_JM\tIBD_var_JE\tIBD_var_ME\t'
        head = head + 'IBD30_mean_eAeA\tIBD30_mean_wAwA\tIBD30_mean_JJ\tIBD30_mean_MM\tIBD30_mean_EE\tIBD30_mean_eAwA\tIBD30_mean_eAE\tIBD30_mean_wAE\tIBD30_mean_eAJ\tIBD30_mean_wAJ\tIBD30_mean_eAM\tIBD30_mean_wAM\tIBD30_mean_JM\tIBD30_mean_JE\tIBD30_mean_ME\t'
        head = head + 'IBD30_median_eAeA\tIBD30_median_wAwA\tIBD30_median_JJ\tIBD30_median_MM\tIBD30_median_EE\tIBD30_median_eAwA\tIBD30_median_eAE\tIBD30_median_wAE\tIBD30_median_eAJ\tIBD30_median_wAJ\tIBD30_median_eAM\tIBD30_median_wAM\tIBD30_median_JM\tIBD30_median_JE\tIBD30_median_ME\t'
        head = head + 'IBD30_num_eAeA\tIBD30_num_wAwA\tIBD30_num_JJ\tIBD30_num_MM\tIBD30_num_EE\tIBD30_num_eAwA\tIBD30_num_eAE\tIBD30_num_wAE\tIBD30_num_eAJ\tIBD30_num_wAJ\tIBD30_num_eAM\tIBD30_num_wAM\tIBD30_num_JM\tIBD30_num_JE\tIBD30_num_ME\t'
        head = head + 'IBD30_var_eAeA\tIBD30_var_wAwA\tIBD30_var_JJ\tIBD30_var_MM\tIBD30_var_EE\tIBD30_var_eAwA\tIBD30_var_eAE\tIBD30_var_wAE\tIBD30_var_eAJ\tIBD30_var_wAJ\tIBD30_var_eAM\tIBD30_var_wAM\tIBD30_var_JM\tIBD30_var_JE\tIBD30_var_ME\t'
        head = head + 'SegS_Af_ASC\tSing_Af_ASC\tDupl_Af_ASC\tPi_Af_ASC\tTajD_Af_ASC\t'
        head = head + 'SegS_Eu_ASC\tSing_Eu_ASC\tDupl_Eu_ASC\tPi_Eu_ASC\tTajD_Eu_ASC\t'
        head = head + 'SegS_As_ASC\tSing_As_ASC\tDupl_As_ASC\tPi_As_ASC\tTajD_As_ASC\t'
        head = head + 'SegS_J_ASC\tSing_J_ASC\tDupl_J_ASC\tPi_J_ASC\tTajD_J_ASC\t'
        head = head + 'SegS_M_ASC\tSing_M_ASC\tDupl_M_ASC\tPi_M_ASC\tTajD_M_ASC\t'
        head = head + 'SegS_EA_ASC\tSing_EA_ASC\tDupl_EA_ASC\tPi_EA_ASC\tTajD_EA_ASC\t'
        head = head + 'SegS_WA_ASC\tSing_WA_ASC\tDupl_WA_ASC\tPi_WA_ASC\tTajD_WA_ASC\t'
        head = head + 'FST_AfEu_ASC\tFST_AfAs_ASC_m\tFST_EuAs_ASC\t'
        head = head + 'FST_eAwA_ASC\tFST_eAEu_ASC\tFST_eAJ_ASC\tFST_eAM_ASC\tFST_MJ_ASC\tFST_wAEu_ASC\tFST_wAJ_ASC\tFST_wAM_ASC\n'

    elif model == "3":
        head_param = 'Sim\tAsc_NAF\tAsc_NEU\tAsc_NCHB\tdaf\tLog10_NAF\tLog10_NANC\tLog10_NCEU\tLog10_NCHB\tLog10_NWA\tLog10_NEA\tLog10_NJ\tLog10_NM\trWA\trEA\trMJ\tmE\tmW\tTgrowth_Af\tTAF\tTEM\tTeu_as\tTA\tTMJ\tTAEW\tTmE\tTmW\t'
        head = 'SegS_Af_CGI\tSing_Af_CGI\tDupl_Af_CGI\tTajD_Af_CGI\t'
        head = head + 'SegS_Eu_CGI\tSing_Eu_CGI\tDupl_Eu_CGI\tTajD_Eu_CGI\t'
        head = head + 'SegS_As_CGI\tSing_As_CGI\tDupl_As_CGI\tTajD_As_CGI\t'
        head = head + 'FST_AfEu_CGI\tFST_AfAs_CGI\tFST_EuAs_CGI\t'
        head = head + 'IBD_mean_eAeA\tIBD_mean_wAwA\tIBD_mean_JJ\tIBD_mean_MM\tIBD_mean_EE\tIBD_mean_eAwA\tIBD_mean_eAE\tIBD_mean_wAE\tIBD_mean_eAJ\tIBD_mean_wAJ\tIBD_mean_eAM\tIBD_mean_wAM\tIBD_mean_JM\tIBD_mean_JE\tIBD_mean_ME\t'
        head = head + 'IBD_median_eAeA\tIBD_median_wAwA\tIBD_median_JJ\tIBD_median_MM\tIBD_median_EE\tIBD_median_eAwA\tIBD_median_eAE\tIBD_median_wAE\tIBD_median_eAJ\tIBD_median_wAJ\tIBD_median_eAM\tIBD_median_wAM\tIBD_median_JM\tIBD_median_JE\tIBD_median_ME\t'
        head = head + 'IBD_num_eAeA\tIBD_num_wAwA\tIBD_num_JJ\tIBD_num_MM\tIBD_num_EE\tIBD_num_eAwA\tIBD_num_eAE\tIBD_num_wAE\tIBD_num_eAJ\tIBD_num_wAJ\tIBD_num_eAM\tIBD_num_wAM\tIBD_num_JM\tIBD_num_JE\tIBD_num_ME\t'
        head = head + 'IBD_var_eAeA\tIBD_var_wAwA\tIBD_var_JJ\tIBD_var_MM\tIBD_var_EE\tIBD_var_eAwA\tIBD_var_eAE\tIBD_var_wAE\tIBD_var_eAJ\tIBD_var_wAJ\tIBD_var_eAM\tIBD_var_wAM\tIBD_var_JM\tIBD_var_JE\tIBD_var_ME\t'
        head = head + 'IBD30_mean_eAeA\tIBD30_mean_wAwA\tIBD30_mean_JJ\tIBD30_mean_MM\tIBD30_mean_EE\tIBD30_mean_eAwA\tIBD30_mean_eAE\tIBD30_mean_wAE\tIBD30_mean_eAJ\tIBD30_mean_wAJ\tIBD30_mean_eAM\tIBD30_mean_wAM\tIBD30_mean_JM\tIBD30_mean_JE\tIBD30_mean_ME\t'
        head = head + 'IBD30_median_eAeA\tIBD30_median_wAwA\tIBD30_median_JJ\tIBD30_median_MM\tIBD30_median_EE\tIBD30_median_eAwA\tIBD30_median_eAE\tIBD30_median_wAE\tIBD30_median_eAJ\tIBD30_median_wAJ\tIBD30_median_eAM\tIBD30_median_wAM\tIBD30_median_JM\tIBD30_median_JE\tIBD30_median_ME\t'
        head = head + 'IBD30_num_eAeA\tIBD30_num_wAwA\tIBD30_num_JJ\tIBD30_num_MM\tIBD30_num_EE\tIBD30_num_eAwA\tIBD30_num_eAE\tIBD30_num_wAE\tIBD30_num_eAJ\tIBD30_num_wAJ\tIBD30_num_eAM\tIBD30_num_wAM\tIBD30_num_JM\tIBD30_num_JE\tIBD30_num_ME\t'
        head = head + 'IBD30_var_eAeA\tIBD30_var_wAwA\tIBD30_var_JJ\tIBD30_var_MM\tIBD30_var_EE\tIBD30_var_eAwA\tIBD30_var_eAE\tIBD30_var_wAE\tIBD30_var_eAJ\tIBD30_var_wAJ\tIBD30_var_eAM\tIBD30_var_wAM\tIBD30_var_JM\tIBD30_var_JE\tIBD30_var_ME\t'
        head = head + 'SegS_Af_ASC\tSing_Af_ASC\tDupl_Af_ASC\tPi_Af_ASC\tTajD_Af_ASC\t'
        head = head + 'SegS_Eu_ASC\tSing_Eu_ASC\tDupl_Eu_ASC\tPi_Eu_ASC\tTajD_Eu_ASC\t'
        head = head + 'SegS_As_ASC\tSing_As_ASC\tDupl_As_ASC\tPi_As_ASC\tTajD_As_ASC\t'
        head = head + 'SegS_J_ASC\tSing_J_ASC\tDupl_J_ASC\tPi_J_ASC\tTajD_J_ASC\t'
        head = head + 'SegS_M_ASC\tSing_M_ASC\tDupl_M_ASC\tPi_M_ASC\tTajD_M_ASC\t'
        head = head + 'SegS_EA_ASC\tSing_EA_ASC\tDupl_EA_ASC\tPi_EA_ASC\tTajD_EA_ASC\t'
        head = head + 'SegS_WA_ASC\tSing_WA_ASC\tDupl_WA_ASC\tPi_WA_ASC\tTajD_WA_ASC\t'
        head = head + 'FST_AfEu_ASC\tFST_AfAs_ASC_m\tFST_EuAs_ASC\t'
        head = head + 'FST_eAwA_ASC\tFST_eAEu_ASC\tFST_eAJ_ASC\tFST_eAM_ASC\tFST_MJ_ASC\tFST_wAEu_ASC\tFST_wAJ_ASC\tFST_wAM_ASC\n'


    files_results = listdir(sim_path_results)

    out_file_name = str(out_path)+"/input_ABCtoolbox_M"+str(model)+"_"+str(len(files_results))+".txt"
    fileout = open(out_file_name, 'w')
    fileout.write(head_param+head)

    for f in files_results:
        if f.endswith(".summary"):
            print str(sim_path_results)+"/"+str(f)
            job_id = f.split("_")[2].split(".")[0]
            sim_values_file_name = str(sim_path_values)+"/sim_"+str(job_id)+"_values.txt"
            if os.path.exists(sim_values_file_name):
                print sim_values_file_name
                line = str(job_id)+"\t"+linecache.getline(sim_values_file_name, 2).rstrip('\n')+"\t"+linecache.getline(str(sim_path_results)+"/"+str(f), 2)
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

if __name__ == '__main__':
    main()