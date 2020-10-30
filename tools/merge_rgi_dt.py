#! /usr/bin/env python3

import sys

def struc_rgi_head(dt_in):
    dt_out = [line.rstrip().split('\t') for line in open(dt_in)]
    tmp = {}
    for line in dt_out:
        id_num = line[4]
        tmp[id_num]= line[:4]
    dt_out = tmp
    return dt_out

 
def struc_rgi_anno(dt_in):
    dt_out = {}
    dt = [line.rstrip().split('\t') for line in open(dt_in)]
    dt = dt[1:]
    for line in dt:
        id_num = '_'.join(line[1].split('_')[:2])
        if id_num not in dt_out:
            dt_out[id_num] = []
        dt_out[id_num].append(line[2:5] + line[8:9] + line[14:17])
    return dt_out


def struc_pre_anno(dt_in):
    def struc_inner_lines(lines):
        dt_out = {}
        for line in lines:
            if line[0] != '\t':
                tmp = line.split('\t')
                key = tuple(tmp[:2])
                dt_out[key] = []
                dt_out[key].append(tmp)
                cont = []
                dt_out[key].append(cont)
            else:
                cont.append(line.strip())
        return dt_out


    dt_out = {}
    for line in open(dt_in):
        line = line.rstrip()
        if line[0] == '>':
            strain = line[1:]
            dt_out[strain] = {}
        elif line[0] == '^':
            contig = line[1:]
            dt_out[strain][contig] = []
        else:
            dt_out[strain][contig].append(line)
    for strain in dt_out:
        for contig in dt_out[strain]:
            lines = dt_out[strain][contig]
            lines_dt = struc_inner_lines(lines)
            dt_out[strain][contig] = lines_dt
            
    return dt_out


def merge_dt(rgi_head_dt, rgi_anno_dt, pre_anno):
    for id_num in rgi_anno_dt:
        pos_info = rgi_head_dt[id_num]
        pos = tuple(pos_info[2:4])
        for line in rgi_anno_dt[id_num]:
            line[0] = str(int(line[0]) + int(pos[0]) - 1)
            line[1] = str(int(line[1]) + int(pos[0]) - 1)
        strain = pos_info[0]
        contig = pos_info[1]
        try:
            pre_anno[strain][contig][pos].append(rgi_anno_dt[id_num])
        except Exception as ex:
            print(ex)
    return 0


def print_res(dt_in):
    f_out = open("repeat_final_anno.fasta", "w")
    for strain in dt_in:
        print('>' + strain, file=f_out)
        for contig in dt_in[strain]:
            print('^' + contig, file=f_out)
            for pos in dt_in[strain][contig]:
                dt = dt_in[strain][contig][pos]
                if len(dt) == 2:
                    pos_head = dt[0]
                    pos_head = pos_head[0:-2] + [pos_head[-1]] + ['No_AMR'] + [pos_head[-2]]
                    print('$' + '\t'.join(pos_head), file=f_out)
                    for line in dt[1]:
                        print('GFF' + '\t' + line, file=f_out)
                elif len(dt) == 3:
                    pos_head = dt[0]
                    pos_head = pos_head[0:-2] + [pos_head[-1]] + ['Yes_AMR'] + [pos_head[-2]]
                    print('$'+'\t'.join(pos_head), file=f_out)
                    for line in dt[1]:
                        print('GFF' + '\t' + line, file=f_out)
                    for line in dt[2]:
                        print('\t'.join(['RGI'] + line), file=f_out)
                else:
                    print('Error2')
                
    return 0


def main(myname='myname', args=None):
    rgi_head_dt = struc_rgi_head(sys.argv[1])
    rgi_anno_dt = struc_rgi_anno(sys.argv[2])
    pre_anno_dt = struc_pre_anno(sys.argv[3])
    merge_dt(rgi_head_dt, rgi_anno_dt, pre_anno_dt)
    #print(pre_anno_dt)
    print_res(pre_anno_dt)
    return 0


if __name__ == '__main__':
    main()
