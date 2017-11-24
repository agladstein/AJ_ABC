from sys import argv
import pandas as pd

def main():
    input_file = argv[1]
    out_file_name = 'cleaned_'+str(input_file)

    chunksize = 10 ** 4

    n = 0
    for input_file_df in pd.read_csv(input_file, sep="\t", chunksize=chunksize, error_bad_lines=False):
        if n == 0:
            input_file_df.dropna().to_csv(out_file_name, sep='\t', index=False)
        else:
            input_file_df.dropna().to_csv(out_file_name, mode='a', sep='\t', index=False, header=False)
        n = n + 1

if __name__ == '__main__':
    main()