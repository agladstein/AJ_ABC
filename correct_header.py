import fileinput
import multiprocessing
import sys
import os.path


def fix_broken_file(file_name):
    file_name = file_name.strip()
    headers_and_replacements = {'IBD_num_EE': 'IBD30_num_EE', 'IBD_var_EE': 'IBD30_var_EE'}
    replaced = False
    with open(file_name, 'r') as f:
        first_line = f.readline()
        second_line = f.readline()
    for header, replacement in headers_and_replacements.iteritems():
        if first_line.count(header) <= 1:
            continue
        parts = first_line.split(header)
        first_line = parts[0] + header + replacement.join(parts[1:])
        replaced = True

    path_parts = list(os.path.split(file_name))
    path_parts[0] = path_parts[0] + '_fixed'
    fixed_directory = '/'.join(path_parts[:-1])
    try:
        os.mkdir(fixed_directory)
    except:
        pass  # Ignore if it already exists
    new_file_name = '/'.join(path_parts)
    with open(new_file_name, 'w') as write_file:
        write_file.write(first_line)
        write_file.write(second_line)
    return file_name, replaced


def main():
    broken_files_file = sys.argv[1]
    pool = multiprocessing.Pool()

    broken_files_iterator = fileinput.input([broken_files_file])

    results = pool.imap_unordered(fix_broken_file, broken_files_iterator, 1000)
    for file_name, replaced in results:
        if replaced is True:
            # All good. Fixed file without complaining.
            sys.stderr.write('.')
        else:
            # Oops. Didn't replace things.
            sys.stderr.write('\n{}: {}\n'.format(file_name, replaced))

if __name__ == '__main__':
    main()
