'''makes new file input_ABCtoolbox_M1_HPC_OSG_203.txt, with column headers listed header_M1_EW.txt from input_ABCtoolbox_M1_HPC_OSG.txt.
Generalize script later.
'''

import pandas as pd

input_file_df = pd.read_csv('input_ABCtoolbox_M1_HPC_OSG.txt', sep="\t")
h = open("header_M1_EW.txt")
line = h.readlines()
head=line[0].rstrip('\n').split('\t')
out_file=input_file_df[head]
out_file.to_csv('input_ABCtoolbox_M1_HPC_OSG_203.txt', sep='\t', index=False)