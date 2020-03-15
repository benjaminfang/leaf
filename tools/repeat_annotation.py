#! /usr/bin/env python3

import argparse
from biolib.bioparser import Gff_parser
import os


def get_args(arg_list):
    parser = argparse.ArgumentParser()
    parser.add_argument('repeat_file', help='repeat file name.')
    parser.add_argument('gff_dir', help='gff directory path.')
    if arg_list:
        args = parser.parse_args(arg_list)
    else:
        args = parser.parse_args()
    repeat_file, gff_dir = args.repeat_file, args.gff_dir
    return repeat_file, gff_dir


def parse_gff_file(gff_file):
    def adept_dt_struc(gff_dt):
        dt_out = {}
        for contig in gff_dt:
            dt_out[contig] = {}
            for source in gff_dt[contig]:
                for seq_type in gff_dt[contig][source]:
                    if seq_type == 'region':
                        continue
                    for position in gff_dt[contig][source][seq_type]:
                        if position not in dt_out[contig]:
                            dt_out[contig][position] = {}
                        dt_out[contig][position][seq_type] = gff_dt[contig][source][seq_type][position]
        return dt_out


    dt = Gff_parser(gff_file)
    dt = dt.data['information']
    dt_adepted = adept_dt_struc(dt)
    return dt_adepted


def structure_repeat_data(file_name):
    def structure_line(line):
        data_out = line.split(';')
        tmp = []
        for ele in data_out[:-3]:
            ele = ele.split(',')
            ele[0] = int(ele[0])
            ele[1] = int(ele[1])
            tmp.append(ele)
        data_out = tmp + data_out[-3:]
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


def annotate_range(range_piece, gff_head):
    def extract_info(info_dict):
        item_priority = ['CDS', 'tRNA', 'rRNA', 'SRP_RNA', 'exon', 'gene', 'pseudogene', 'riboswitch', 'sequence_feature']
        tmp = []
        for info_type in info_dict:
            tmp.append([info_type, item_priority.index(info_type)])
        tmp.sort(key=lambda x:x[1])
        info_choose = tmp[0][0]
        info = info_dict[info_choose]
        strand = info['strand']
        attr = info['attributes']
        string = []
        for key in attr:
            string.append(key + '=' + attr[key])
        string = ';'.join(string)
        return {'info_type':info_choose, 'strand':strand, 'attributes':string}


    dt_out = {}
    range_start = range_piece[0]
    range_end = range_piece[1]
    for ele in gff_head:
        if ele[0] >= range_start and ele[1] <= range_end:
            info = extract_info(gff_head[ele])
            dt_out[ele] = info
    return dt_out


def output_to_file(data, f_out):
    def make_line(dt_in):
        string = []
        for range_piece in dt_in:
            start = str(range_piece[0])
            end = str(range_piece[1])
            info_type = dt_in[range_piece]['info_type']
            strand = dt_in[range_piece]['strand']
            attr = dt_in[range_piece]['attributes']
            string.append('\t'.join([start, end, strand, info_type, attr]))
        dt_out = '|'.join(string)
        return dt_out


    for fasta_file in data:
        print('>' + fasta_file, file=f_out)
        for head in data[fasta_file]:
            print('^' + head, file=f_out)
            for record in data[fasta_file][head]:
                p1 = record[:-4]
                p1 = ';'.join([','.join([str(ele[0]), str(ele[1]), ele[2]]) for ele in p1])
                p2 = ';'.join(record[-4:-1])
                print(p1, p2, sep=';', file=f_out)
                p3 = record[-1]
                p3 = make_line(p3)
                print(p3, file=f_out)
    return 0


def main(name='name', args=None):
    myname = 'name'
    if name == myname:
        repeat_file, gff_dir = get_args(args)
        f_out = open('repeat_annotated.fasta', 'w')
        repeat_data = structure_repeat_data(repeat_file)
        gff_file_dic = {ff:os.path.join(gff_dir, ff) for ff in os.listdir(gff_dir)}
        for fasta_file in repeat_data:
            gff_data = parse_gff_file(gff_file_dic[fasta_file[:-3] + 'gff'])
            for head in repeat_data[fasta_file]:
                gff_head = gff_data[head]
                for record in repeat_data[fasta_file][head]:
                    range_piece = record[0]
                    info = annotate_range(range_piece, gff_head)
                    record.append(info)
        output_to_file(repeat_data, f_out)
        f_out.close()
    return 0


if __name__ == '__main__':
    main()
