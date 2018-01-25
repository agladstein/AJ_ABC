from sys import argv
from main_function_AJmodel2_chr import main


if __name__ == '__main__':
    if len(argv) == 4:
        argv.extend('.')
    main(argv)
