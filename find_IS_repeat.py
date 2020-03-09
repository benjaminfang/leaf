#! /usr/bin/env python3
# --*-- utf-8 code --*--

import os
import argparse
import statistics
import math
from biolib.bioparser import Fasta_parser
from biolib.biocodon import base_complement


def get_args(arg_list):
    parser = argparse.ArgumentParser()
    parser.add_argument('repeat_file', help='repeat file name.')
    parser.add_argument('fasta_dir', help='directory contain fasta file.')
    parser.add_argument('IS_file', help='IS file name.')
    if arg_list:
        args = parser.parse_args(arg_list)
    else:
        args = parser.parse_args()
    repeat_file, fasta_dir, IS_file = args.repeat_file, args.fasta_dir, args.IS_file
    return repeat_file, fasta_dir, IS_file


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
            fasta_file = os.path.basename(line[1:])
            data_out[fasta_file] = {}
        elif line[0] == '^':
            head_name = line[1:]
            data_out[fasta_file][head_name] = []
        else:
            line_struc = structure_line(line)
            data_out[fasta_file][head_name].append(line_struc)
    return data_out


def save_seq_as_fasta(seq, file_name, directory):
    f_path = os.path.join(directory, file_name)
    f_out = open(f_path, 'w')
    f_out.write('>seq_in_fasta\n')
    f_out.write(seq + '\n')
    f_out.close()
    return f_path


def makeblastdb(file_name, dbname, directory):
    cmdname = 'makeblastdb'
    opt_dbtype = '-dbtype nucl'
    opt_in = '-in ' + file_name
    out_path = os.path.join(directory, dbname)
    opt_out = '-out ' + out_path
    cmd = ' '.join([cmdname, opt_dbtype, opt_in, opt_out, '&>/dev/null'])
    os.system(cmd)
    return out_path


def runblast(item_file, task, db_path, directory):
    cmdname = 'blastn'
    opt_task = '-task ' + task
    opt_query = '-query ' + item_file
    opt_db = '-db ' + db_path
    blast_res = os.path.join(directory, 'blastres')
    opt_out = '-out ' + blast_res
    opt_outfmt = '-outfmt ' + '7'
    cmd = ' '.join([cmdname, opt_query, opt_db, opt_out, opt_outfmt, opt_task, '&>/dev/null'])
    os.system(cmd)
    return blast_res


def structure_blast_res(query_seq_length, blast_result_file):
    # this function is used to structure blast result in format 7.
    dt_out = []
    a_lines = [line.rstrip().split() for line in open(blast_result_file) if line[0] != '#']
    for line in a_lines:
        identity = round(float(line[2])/100, 5)
        q_s = int(line[6])
        q_e = int(line[7])
        s_s = int(line[8])
        s_e = int(line[9])
        if s_s <= s_e:
            origntation = '+'
        else:
            origntation = '-'
            s_s, s_e = s_e, s_s
        up_margine = q_s - 1
        down_margine = query_seq_length - q_e
        coverage = round((q_e - q_s + 1)/query_seq_length, 5)
        dt_out.append([s_s, s_e, origntation, up_margine, down_margine, coverage, identity])
    return dt_out



def append_sequence(repeat_data, fasta_dir):
    fasta_dir_data = {ff:os.path.join(fasta_dir, ff) for ff in os.listdir(fasta_dir)}
    for fasta_file in repeat_data:
        fasta_file_path = fasta_dir_data[fasta_file]
        fasta_data = Fasta_parser(fasta_file_path)
        fasta_data.join_lines()
        fasta_data = fasta_data.data
        for head in repeat_data[fasta_file]:
            for record in repeat_data[fasta_file][head]:
                piece = record[0]
                seq = fasta_data[head][piece[0] - 1 : piece[1]]
                if piece[2] == '-':
                    seq = list(seq)
                    seq.reverse()
                    seq = ''.join([base_complement[base] for base in seq])
                record.append(seq)
    return repeat_data


def judge_IS(repeat_data_appended, IS_blast_db):
    for fasta_file in repeat_data_appended:
        for head in repeat_data_appended[fasta_file]:
            for record in repeat_data_appended[fasta_file][head]:
                query_length = record[0][1] - record[0][0] + 1
                seq = record[-1]
                seq_fasta = save_seq_as_fasta(seq, 'seq_tmp.fasta', '.')
                blast_res = runblast(seq_fasta, 'megablast', IS_blast_db, '.' )
                blast_res_struced = structure_blast_res(query_length, blast_res)
                blast_res_struced = [ele for ele in blast_res_struced if ele[-1] >= 0.6 and ele[-2] >= 0.6]
                if len(blast_res_struced) > 0:
                    record.append('is_IS')
                else:
                    record.append('not_IS')
    return repeat_data_appended


def output_data_to_file(data, f_out):
    for fasta_file in data:
        print('>' + fasta_file, file=f_out)
        for head in data[fasta_file]:
            print('^' + head, file=f_out)
            for record in data[fasta_file][head]:
                ava_len = math.floor(statistics.mean([ele[1] - ele[0] + 1 for ele in record[:-2]]))
                string = ';'.join([','.join([str(ele[0]), str(ele[1]), ele[2]]) for ele in record[:-2]]) + ';' + str(ava_len) \
                    + ';' + record[-2] + ';' + record[-1]
                print(string, file=f_out)
    return 0

def main(name='name', args=None):
    myname = 'name'
    if name == myname:
        repeat_file, fasta_dir, IS_file = get_args(args)
        repeat_data = structure_file(repeat_file)
        repeat_data_appended = append_sequence(repeat_data, fasta_dir)
        IS_blast_db = makeblastdb(IS_file, 'IS_blast_db', '.')
        repeat_data_IS = judge_IS(repeat_data_appended, IS_blast_db)
        f_out = open('repeat_IS_judeged.fasta', 'w')
        output_data_to_file(repeat_data_IS, f_out)
        f_out.close()
    return 0


if __name__ == '__main__':
    main()


