from sys import argv
from main_function_AJmodel1_chr1 import main


if __name__ == '__main__':
    if len(argv) == 7:
        argv.extend('.')
    main(argv)
