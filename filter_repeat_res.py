#! /usr/bin/env python3

import argparse
import math
import statistics

def get_args(arg_list):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', required=True, help='file name processing.')
    parser.add_argument('-max_len', type=int, default=None, help='maximal length.')
    parser.add_argument('-min_len', type=int, default=None, help='minimal length.')
    if arg_list:
        args = parser.parse_args(arg_list)
    else:
        args = parser.parse_args()
    f_name, max_len, min_len = args.f, args.max_len, args.min_len
    return f_name, max_len, min_len


def structure_file(file_name):
    def structure_line(line):
        data_out = [ele.split(',') for ele in line.split(';')[:-1]]
        for ele in data_out:
            ele[0] = int(ele[0])
            ele[1] = int(ele[1])
        return data_out


    data_out = {}
    all_lines = [line.rstrip() for line in open(file_name)]
    for line in all_lines:
        if line[0] == '>':
            fasta_file = line[1:]
            data_out[fasta_file] = {}
        elif line[0] == '^':
            head_name = line[1:]
            data_out[fasta_file][head_name] = []
        else:
            line_struc = structure_line(line)
            data_out[fasta_file][head_name].append(line_struc)
    return data_out


def filter_data(data_in, max_len, min_len):
    data_out = {}
    if not max_len and not min_len:
        data_out = data_in
        return data_out
    elif not max_len and min_len:
        for fasta_file in data_in:
            data_out[fasta_file] = {}
            for head in data_in[fasta_file]:
                data_out[fasta_file][head] = []
                for record in data_in[fasta_file][head]:
                    if all([(ele[1] - ele[0] + 1) >= min_len for ele in record]):
                        data_out[fasta_file][head].append(record)
    elif max_len and not min_len:
         for fasta_file in data_in:
            data_out[fasta_file] = {}
            for head in data_in[fasta_file]:
                data_out[fasta_file][head] = []
                for record in data_in[fasta_file][head]:
                    if all([(ele[1] - ele[0] + 1) <= max_len for ele in record]):
                        data_out[fasta_file][head].append(record)
    else:
        for fasta_file in data_in:
            data_out[fasta_file] = {}
            for head in data_in[fasta_file]:
                data_out[fasta_file][head] = []
                for record in data_in[fasta_file][head]:
                    lengths = [ele[1] - ele[0] + 1 for ele in record]
                    if all([(ele >= min_len) and (ele <= max_len) for ele in lengths]):
                        data_out[fasta_file][head].append(record)
    return data_out


def output_data_to_file(file_data_filtered, f_out):
    for fasta_file in file_data_filtered:
        print('>' + fasta_file, file=f_out)
        for head in file_data_filtered[fasta_file]:
            print('^' + head, file=f_out)
            for record in file_data_filtered[fasta_file][head]:
                ava_len = math.floor(statistics.mean([ele[1] - ele[0] + 1 for ele in record]))
                string = ';'.join([','.join([str(ele[0]), str(ele[1]), ele[2]]) for ele in record]) + ';' + str(ava_len)
                print(string, file=f_out)
    return 0


def main(called_name='name', args=None):
    my_name = 'name'
    if my_name == called_name:
        f_name, max_len, min_len = get_args(args)
        f_out = open('repeat_filtered.fasta', 'w')
        file_data = structure_file(f_name)
        file_data_filtered = filter_data(file_data, max_len, min_len)
        fd = file_data_filtered
        output_data_to_file(file_data_filtered, f_out)
        f_out.close()
    return 0


if __name__ == '__main__':
    main()


