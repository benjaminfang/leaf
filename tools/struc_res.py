#! /usr/bin/env python3
# --*-- utf-8 code --*--


def structure_file(file_name):
    def structure_line(line):
        data_out = []
        tmp = line.split('\t')
        data_out = tmp
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


