#! /usr/bin/env python3

import argparse


def get_args(arg_list):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', required=True, help='file name processing.')
    parser.add_argument('-max_len', type=int, default=None, help='maximal length.')
    parser.add_argument('-min_len', type=int, default=None, help='minimal length.')
    parser.add_argument('-dis_range_down', type=int, default=None, help='lower repeats distributing range length.')
    parser.add_argument('-dis_range_up', type=int, default=None, help='higher repeats distributing range length.')
    parser.add_argument('-gap_len', type=int, default=None, help='gap length.')
    if arg_list:
        args = parser.parse_args(arg_list)
    else:
        args = parser.parse_args()
    f_name, max_len, min_len, dis_range_down, dis_range_up, gap_len = args.f, args.max_len, args.min_len, args.dis_range_down, args.dis_range_up, args.gap_len
    return f_name, max_len, min_len, dis_range_down, dis_range_up, gap_len


def structure_file(file_name):
    def structure_line(line):
        data_out = []
        tmp = line.split('\t')
        for i in list(range(len(tmp) - 1))[::3]:
            data_out.append([int(tmp[i]), int(tmp[i + 1]), tmp[i + 2]])
        data_out.append(int(tmp[-1]))
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


def filter_data(data_in, max_len, min_len, dis_range_down, dis_range_up, gap_len):
    def juder(record, max_len, min_len, dis_range, gap_len):
        record = record[:-1]
        elements = [ele[1] - ele[0] + 1 for ele in record]
        elements.sort()
        max_piece = elements[-1]
        min_piece = elements[0]
        record.sort(key=lambda x: x[0])
        range_len = record[-1][1] - record[0][0]
        max_gaps = max([record[i + 1][0] - record[i][1] + 1 for i in range(len(record) - 1)])
        if max_len:
            if max_piece > max_len:
                return False
        if min_len:
            if min_piece < min_len:
                return False
        if dis_range_down:
            if range_len <= dis_range_down:
                return False
        if dis_range_up:
            if range_len > dis_range_up:
                return False
        if gap_len:
            if max_gaps > gap_len:
                return False
        return True


    data_out = {}
    for fasta_file in data_in:
        data_out[fasta_file] = {}
        for head in data_in[fasta_file]:
            data_out[fasta_file][head] = []
            for record in data_in[fasta_file][head]:
                if juder(record, max_len, min_len, dis_range_down, dis_range_up,  gap_len):
                    data_out[fasta_file][head].append(record)
    return data_out


def output_data_to_file(file_data_filtered, f_out):
    for fasta_file in file_data_filtered:
        print('>' + fasta_file, file=f_out)
        for head in file_data_filtered[fasta_file]:
            print('^' + head, file=f_out)
            for record in file_data_filtered[fasta_file][head]:
                string = '\t'.join(['\t'.join([str(ele[0]), str(ele[1]), ele[2]]) for ele in record[: -1]]) + '\t' + str(record[-1])
                print(string, file=f_out)
    return 0


def main(called_name='name', args=None):
    my_name = 'name'
    if my_name == called_name:
        f_name, max_len, min_len, dis_range_down, dis_range_up, gap_len = get_args(args)
        f_out = open('repeat_filtered.fasta', 'w')
        file_data = structure_file(f_name)
        file_data_filtered = filter_data(file_data, max_len, min_len, dis_range_down, dis_range_up, gap_len)
        output_data_to_file(file_data_filtered, f_out)
        f_out.close()
    return 0


if __name__ == '__main__':
    main()


