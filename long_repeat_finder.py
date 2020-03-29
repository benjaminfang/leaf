#! /usr/bin/env python3

import os
import shutil
import time
import argparse
import multiprocessing as mp
import math
import statistics
from biolib.bioparser import Fasta_parser


def get_args(args_list):
    args_parser = argparse.ArgumentParser(prog='long_repeat_finder', description='This \
        utility used to find long repeat sequence among genome.')
    args_parser.add_argument('-D', '--directory_fasta_file', default=None,  help='directory path of fasta file.')
    args_parser.add_argument('-L', '--file_list', default=None, help='file list which want to process.')
    args_parser.add_argument('-S', '--seed_length', type=int, default=200, help='seed_length.')
    args_parser.add_argument('-I', '--identity_cutoff', type=float, default=0.9, help='identity cutoff used to filter seed sequence.')
    args_parser.add_argument('-C', '--coverage_cutoff', type=float, default=0.9, help='coverage cutoff used to filter seed sequence.')
    args_parser.add_argument('-E', '--expand_length', type=int, default=100, help='expanding length when expanding seed.')
    args_parser.add_argument('-T', '--threads', type=int, default=None, help='threads used.')
    if args_list:
        args = args_parser.parse_args(args_list)
    else:
        args = args_parser.parse_args()
    directory, file_list, seed_length, identity_cutoff, coverage_cutoff, expand_len, threads = args.directory_fasta_file, args.file_list, \
        args.seed_length, args.identity_cutoff, args.coverage_cutoff, args.expand_length, args.threads
    if (not (directory or file_list)) or (directory and file_list):
        raise Exception('the directory of file list must offer one.')
    if threads == None:
        threads = round(os.cpu_count() * 0.6)
    return directory, file_list, seed_length, identity_cutoff, coverage_cutoff, expand_len, threads


class Seed_provider:
    def __init__(self, data_source, seq_name=None, source_type='File'):
        if source_type == 'File':
            data = Fasta_parser(data_source)
            if len(data.data) > 1:
                raise Exception('too many item sequence in fasta file.')
            data.join_lines()
            for head in data.data:
                self.seq_name = head
                self.sequence = data.data[head]
            self.seq_length = len(self.sequence)
            self.seed_pointer = 0
            self.seed_island = [(1, self.seq_length)]
        elif source_type == 'Seq':
            self.seq_name = seq_name
            self.sequence = data_source
            self.seq_length = len(self.sequence)
            self.seed_pointer = 0
            self.seed_island = [(1, self.seq_length)]


    def trim_seed_island(self, start, end):
        new_seed_island = []
        for island in self.seed_island:
            island_start = island[0]
            island_end = island[1]
            if island_end < start or island_start > end:
                new_seed_island.append(island)
            elif island_start >= start and island_end <= end:
                pass
            elif island_start <= start and island_end >= end:
                first_segment = (island_start, start - 1) if island_start <= start-1 else None
                second_segment = (end + 1, island_end) if island_end >= end + 1 else None
                if first_segment:
                    new_seed_island.append(first_segment)
                if second_segment:
                    new_seed_island.append(second_segment)
            elif island_start < start and island_end <= end:
                new_segment = (island_start, start - 1)
                new_seed_island.append(new_segment)
            elif island_start >= start and island_end > end:
                new_segment = (end + 1, island_end)
                new_seed_island.append(new_segment)
            else:
                raise Exception('opps', island)
        new_seed_island.sort()
        self.seed_island = new_seed_island


    def judge_seed_on_island(self, seed_start, seed_end):
        for island in self.seed_island:
            if island[0] <= seed_start and island[1] >= seed_end:
                return True
        return False


    def next_seed(self, seed_length):
        while True:
            if self.seed_pointer > self.seq_length - seed_length:
                return None
            seed_start = self.seed_pointer + 1
            seed_end = self.seed_pointer + seed_length
            if self.judge_seed_on_island(seed_start, seed_end):
                seed_sequence = self.sequence[seed_start - 1: seed_end]
                self.trim_seed_island(seed_start, seed_end)
                self.seed_pointer = seed_end
                return (seed_start, seed_end, seed_sequence)
            else:
                self.seed_pointer += seed_length


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


def runblast(item_file, task, db_path, directory, out_f):
    cmdname = 'blastn'
    opt_task = '-task ' + task
    opt_query = '-query ' + item_file
    opt_db = '-db ' + db_path
    blast_res = os.path.join(directory, out_f)
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


def filter_seed_blast_res(blast_res_structed, coverage_cutoff, identity_cutoff):
    dt_out = []
    for entry in blast_res_structed:
        if entry[5] >= coverage_cutoff and entry[6] >= identity_cutoff:
            dt_out.append(entry)
    return dt_out


def check_blast_res_overlap(blast_res_structed):
    ranges = []
    i = 1
    for ele in blast_res_structed:
        ranges.append([i, ele[0]])
        ranges.append([i, ele[1]])
        i += 1
    ranges.sort(key=lambda x: x[1])
    for i in list(range(len(ranges)))[::2]:
        if ranges[i][0] != ranges[i+1][0]:
            return False
    return True


def judge_seed(seed_item, blast_res_structed, coverage_cutoff, identity_cutoff):
    dt_out = []
    judge_marker = True
    for entry in blast_res_structed:
        if entry[5] >= coverage_cutoff and entry[6] >= identity_cutoff:
            dt_out.append(entry)
    dt_out.sort(lambda x: x[0])
    if dt_out[0][0] != seed_item[0] or dt_out[0][1] != seed_item[1]:
        judge_marker = False
    front_ele = dt_out[0]
    tmp = []
    tmp.append(front_ele)
    for ele in dt_out[1:]:
        if ele[0] > front_ele[1]:
            tmp.append(ele)
            front_ele = ele
    if len(tmp) <= 1:
        judge_marker = False
    dt_out = tmp
    return dt_out, judge_marker


def offer_meanningful_seed(seed_instance, seed_length, blastdb, coverage_cutoff, identity_cutoff, directory):
    while True:
        seed_sequence = seed_instance.next_seed(seed_length)
        left_boundary = seed_sequence[0]
        if seed_sequence:
            seed_fasta_file = save_seq_as_fasta(seed_sequence[2], 'seed_seq.fasta', directory)
            blast_res = runblast(seed_fasta_file, 'blastn-short', blastdb, directory, 'blastres_seed')
            blast_res_structed = structure_blast_res(seed_length, blast_res)
            if len(blast_res_structed) == 0:
                print('Low complexity seed:', seed_sequence)
                continue
            blast_res_structed_filtered, judge_marker = judge_seed(seed_sequence, blast_res_structed, coverage_cutoff, identity_cutoff)
            if judge_marker:
                yield blast_res_structed_filtered, left_boundary
        else:
            break


def decide_expand_max_len(range_list, lock_up, lock_down, expand_length, whole_seq_length):
    tmp = []
    i = 1
    for ele in range_list:
        if ele[2] == '+':
            tmp.append([i, ele[0], 'l', 'u'])
            tmp.append([i, ele[1], 'r', 'd'])
        else:
            tmp.append([i, ele[0], 'l', 'd'])
            tmp.append([i, ele[1], 'r', 'u'])
        i += 1
    tmp.sort(key=lambda x: x[1])
    for i in list(range(len(tmp)))[::2]:
        if tmp[i][0] != tmp[i + 1][0]:
            print('WARNING: overlap fonund.')
            return {'up': 0, 'down': 0}

    expand_lens = {'up': [], 'down': []}
    first_point = tmp[0]
    last_point = tmp[-1]
    if lock_up and lock_down:
        return {'up': 0, 'down': 0}
    elif lock_up and not lock_down:
        expand_lens['up'].append(0)
        if first_point[3] == 'd':
            expand_lens['down'].append(first_point[1] - 1)
        if last_point[3] == 'd':
            expand_lens['down'].append(whole_seq_length - last_point[1])
        for i in list(range(1, len(tmp) - 1))[::2]:
            gap = tmp[i + 1][1] - tmp[i][1] - 1
            if gap < 0:
                gap = 0
            gap_half = math.floor(gap/2)
            gap_type = (tmp[i][3], tmp[i + 1][3])
            if gap_type == ('u', 'd') or gap_type == ('d', 'u'):
                expand_lens['down'].append(gap)
            elif gap_type == ('d', 'd'):
                expand_lens['down'].append(gap_half)
    elif not lock_up and lock_down:
        expand_lens['down'].append(0)
        if first_point[3] == 'u':
            expand_lens['up'].append(first_point[1] - 1)
        if last_point[3] == 'u':
            expand_lens['up'].append(whole_seq_length - last_point[1])
        for i in list(range(1, len(tmp) - 1))[::2]:
            gap = tmp[i + 1][1] - tmp[i][1] - 1
            if gap < 0:
                gap = 0
            gap_half = math.floor(gap/2)
            gap_type = (tmp[i][3], tmp[i + 1][3])
            if gap_type == ('u', 'd') or gap_type == ('d', 'u'):
                expand_lens['up'].append(gap)
            elif gap_type == ('u', 'u'):
                expand_lens['up'].append(gap_half)
    else:
        if first_point[3] == 'u':
            expand_lens['up'].append(first_point[1] - 1)
        else:
            expand_lens['down'].append(first_point[1] - 1)
        if last_point[3] == 'u':
            expand_lens['up'].append(whole_seq_length - last_point[1])
        else:
            expand_lens['down'].append(whole_seq_length - last_point[1])
        for i in list(range(1, len(tmp) - 1))[::2]:
            gap = tmp[i + 1][1] - tmp[i][1] - 1
            if gap < 0:
                gap = 0
            gap_half = math.floor(gap/2)
            expand_lens['up'].append(gap_half)
            expand_lens['down'].append(gap_half)
    up_max_expand = min(expand_lens['up'])
    down_max_expand = min(expand_lens['down'])
    if up_max_expand < 0 or down_max_expand < 0:
        print('opps')
    if up_max_expand < 0:
        up_max_expand = 0
    if down_max_expand < 0:
        down_max_expand = 0
    if up_max_expand > expand_length:
        up_max_expand = expand_length
    if down_max_expand > expand_length:
        down_max_expand = expand_length
    return {'up': up_max_expand, 'down': down_max_expand}


def expand_range(range_list, max_expand_length):
    dt_out = []
    up_expand = max_expand_length['up']
    down_expand = max_expand_length['down']
    for ele in range_list:
        if ele[2] == '+':
            dt_out.append([ele[0] - up_expand, ele[1] + down_expand, ele[2]])
        else:
            dt_out.append([ele[0] - down_expand, ele[1] + up_expand, ele[2]])
    # dt_out.sort(key=lambda x: x[0])
    return dt_out


def blast_expanded_seq(range_list, item_sequence, blastdb, directory):
    def filter_blast_res_in_range(range_list, blast_res_structed):
        dt_out = {}
        range_order = []
        for ele in range_list:
            dt_out[(ele[0], ele[1])] = []
            range_order.append((ele[0], ele[1]))
        for ele in blast_res_structed:
            for ele2 in dt_out:
                if ele[0] >= (ele2[0] - 1) and ele[1] <= (ele2[1] + 1):
                    dt_out[ele2].append([ele, ele[1] - ele[0] + 1])
        tmp = []
        for ele in range_order:
            if len(dt_out[ele]) >= 1:
                contain_ranged = dt_out[ele]
                contain_ranged.sort(key=lambda x: x[1])
                tmp.append(contain_ranged[-1][0])
            else:
                tmp.append(None)
        dt_out = tmp
        if None in dt_out:
            return None
        else:
            return dt_out


    query_seq_range = range_list[0]
    query_seq_length = query_seq_range[1] - query_seq_range[0] + 1
    query_seq = item_sequence[query_seq_range[0] - 1: query_seq_range[1]]
    query_seq_file = save_seq_as_fasta(query_seq, 'query_expande_seq.fasta', directory)
    blastres_file = runblast(query_seq_file, 'blastn', blastdb, directory, 'blastres_expanded')
    blast_res_structed = structure_blast_res(query_seq_length, blastres_file)
    blast_res_structed_filtered = filter_blast_res_in_range(range_list, blast_res_structed)
    if blast_res_structed_filtered == None:
        return None
    else:
        return blast_res_structed_filtered


def detect_and_lock_boundary(blast_res_structed, lock_up, lock_down):
    dt_out = []
    up_down_margine = {'up': [], 'down': []}
    margine_cutoff = 25
    for ele in blast_res_structed:
        up_down_margine['up'].append(ele[3])
        up_down_margine['down'].append(ele[4])
    up_down_margine['up'] = max(up_down_margine['up'])
    up_down_margine['down'] = max(up_down_margine['down'])

    if lock_up and lock_down:
        dt_out = [ele[:3] for ele in blast_res_structed]
    if (not lock_up) and up_down_margine['up'] > margine_cutoff:
        lock_up = True
        for ele in blast_res_structed:
            if ele[2] == '+':
                ele[0] += up_down_margine['up'] - ele[3]
            else:
                ele[1] -= up_down_margine['up'] - ele[3]
    if (not lock_down) and up_down_margine['down'] > margine_cutoff:
        lock_down = True
        for ele in blast_res_structed:
            if ele[2] == '+':
                ele[1] -= up_down_margine['down'] - ele[4]
            else:
                ele[0] += up_down_margine['down'] - ele[4]
    for ele in blast_res_structed:
        dt_out.append(ele[:3])
    return lock_up, lock_down, dt_out


def modify_blast_caused_overlap(range_list, lock_up, lock_down):
    dt_out = []
    tmp = []
    i = 1
    for ele in range_list:
        if ele[2] == '+':
            tmp.append([i, ele[0], 'l', 'u'])
            tmp.append([i, ele[1], 'r', 'd'])
        else:
            tmp.append([i, ele[0], 'l', 'd'])
            tmp.append([i, ele[1], 'r', 'u'])
        i += 1
    tmp.sort(key=lambda x: x[1])
    for i in list(range(1, len(tmp) - 1))[::2]:
        meet_type = (tmp[i][3], tmp[i+1][3])
        if (tmp[i][0] < tmp[i + 1][0]) and (tmp[i][1] == tmp[i + 1][1]):
            if meet_type == ('d', 'u'):
                if not lock_down:
                    tmp[i][1] -= 1
                elif not lock_up:
                    tmp[i + 1][1] += 1
                else:
                    tmp[i][1] -= 1
                lock_up = True
                lock_down = True
            elif meet_type == ('u', 'd'):
                if not lock_up:
                    tmp[i][1] -= 1
                elif not lock_down:
                    tmp[i + 1][1] += 1
                else:
                    tmp[i][1] -= 1
                lock_up = True
                lock_down = True
            elif meet_type == ('u', 'u'):
                tmp[i][1] -= 1
                lock_up = True
            elif meet_type == ('d', 'd'):
                tmp[i][1] -= 1
                lock_down = True
        elif (tmp[i][0] > tmp[i + 1][0]) and (tmp[i][1] == tmp[i + 1][1]):
            if meet_type == ('d', 'u'):
                if not lock_down:
                    tmp[i][1] += 1
                elif not lock_up:
                    tmp[i + 1][1] -= 1
                else:
                    tmp[i][1] += 1
                lock_up = True
                lock_down = True
            elif meet_type == ('u', 'd'):
                if not lock_up:
                    tmp[i][1] += 1
                elif not lock_down:
                    tmp[i + 1][1] -= 1
                else:
                    tmp[i][1] += 1
                lock_up = True
                lock_down = True
            elif meet_type == ('u', 'u'):
                tmp[i][1] += 1
                lock_up = True
            elif meet_type == ('d', 'd'):
                tmp[i][1] += 1
                lock_down = True
        elif tmp[i][0] > tmp[i + 1][0]:
            len_need_move = tmp[i + 1][1] - tmp[i][1] + 1
            half_move = len_need_move/2
            small_move = math.floor(half_move)
            big_move = len_need_move - small_move
            if meet_type == ('u', 'd'):
                if (not lock_up) and (not lock_down):
                    tmp[i][1] += big_move
                    tmp[i + 1][1] -= small_move
                elif not lock_up:
                    tmp[i][1] += len_need_move
                elif not lock_down:
                    tmp[i + 1][1] -= len_need_move
                else:
                    tmp[i][1] += big_move
                    tmp[i + 1][1] -= small_move
                lock_up = True
                lock_down = True
            elif meet_type == ('d', 'u'):
                if (not lock_up) and (not lock_down):
                    tmp[i][1] += small_move
                    tmp[i + 1][1] -= big_move
                elif not lock_up:
                    tmp[i + 1][1] -= len_need_move
                elif not lock_down:
                    tmp[i][1] += len_need_move
                else:
                    tmp[i][1] += small_move
                    tmp[i + 1][1] -= big_move
                lock_up = True
                lock_down = True
            elif meet_type == ('u', 'u'):
                tmp[i][1] += big_move
                tmp[i + 1][1] -= small_move
                lock_up = True
            elif meet_type == ('d', 'd'):
                tmp[i][1] += small_move
                tmp[i + 1][1] -= big_move
                lock_down = True

    tmp.sort(key=lambda x: x[1])
    for i in range(len(tmp))[::2]:
        oritation = (tmp[i][3], tmp[i+1][3])
        if oritation == ('u', 'd'):
            dt_out.append([tmp[i][1], tmp[i + 1][1], '+'])
        else:
            dt_out.append([tmp[i][1], tmp[i + 1][1], '-'])
    return lock_up, lock_down, dt_out


def print_data(data, f_out_name):
    f_out = open(f_out_name, 'a')
    print('>' + data['file_name'], file=f_out)
    for head in data['head']:
        print('^' + head, file=f_out)
        for record in data['head'][head]:
            string = []
            lengths = []
            for piece in record:
                string += [str(piece[0]), str(piece[1]), piece[2]]
                lengths.append(piece[1] - piece[0] + 1)
            string.append(str(round(statistics.mean(lengths))))
            string = '\t'.join(string)
            print(string, file=f_out)
    f_out.close()


def worker_func(args_in):
    file_name, seed_length, coverage_cutoff, identity_cutoff, expand_length, sub_tmp_dir, f_all_repeat, lock = args_in
    os.mkdir(sub_tmp_dir)
    max_expand_time = 12500
    data_print = {}
    file_basename = os.path.basename(file_name)
    data_print['file_name'] = file_basename
    data_print['head'] = {}
    fasta_dt = Fasta_parser(file_name)
    fasta_dt.join_lines()
    for head in fasta_dt.data:
        head_name = head.split()[0]
        data_print['head'][head_name] = []
        item_sequence = fasta_dt.data[head]
        whole_seq_length = len(item_sequence)
        sequence_file = save_seq_as_fasta(item_sequence, 'item_sequence.fasta', sub_tmp_dir)
        blastdb = makeblastdb(sequence_file, 'blastdb', sub_tmp_dir)
        seed_ins = Seed_provider(item_sequence, head, source_type='Seq')
        seeds_ok_blast_res = offer_meanningful_seed(seed_ins, seed_length, blastdb, coverage_cutoff, identity_cutoff, sub_tmp_dir)
        for seed in seeds_ok_blast_res:
            try:
                left_boundary = seed[1]
                range_list = [ele[:3] for ele in seed[0]]
                lock_up, lock_down = False, False
                expand_count = 0
                while True:
                    expand_count += 1
                    if expand_count > max_expand_time:
                        break
                    if lock_up and lock_down:
                        break
                    max_expand_length = decide_expand_max_len(range_list, lock_up, lock_down, expand_length, whole_seq_length)
                    if max_expand_length['up'] == 0:
                        lock_up = True
                    if max_expand_length['down'] == 0:
                        lock_down = True
                    range_list_expanded = expand_range(range_list, max_expand_length)
                    if not(lock_up and lock_down):
                        blast_res = blast_expanded_seq(range_list_expanded, item_sequence, blastdb, sub_tmp_dir)
                        if blast_res == None:
                            break
                        lock_up, lock_down, new_range = detect_and_lock_boundary(blast_res, lock_up, lock_down)
                        lock_up, lock_down, new_range = modify_blast_caused_overlap(new_range, lock_up, lock_down)
                        range_list = new_range
                    if range_list[0][0] < left_boundary - seed_length:
                        lock_up = True
                for piece in range_list:
                    seed_ins.trim_seed_island(piece[0], piece[1])
                data_print['head'][head_name].append(range_list)
            except Exception:
                lock.acquire()
                print(file_basename, head_name, seed)
                lock.release()
                for ele in seed:
                    seed_ins.trim_seed_island(ele[0], ele[1])
    lock.acquire()
    print_data(data_print, f_all_repeat)
    lock.release()
    shutil.rmtree(sub_tmp_dir)


def main(name='long_repeat_finder', args=None):
    myname = 'long_repeat_finder'
    if name == myname:
        directory, file_list, seed_length, coverage_cutoff, identity_cutoff, expand_length, threads = get_args(args)
        working_dir = 'working_dir_' + time.strftime('%Y%m%d%H%M%S')
        tmp_dir = os.path.join(working_dir, 'tmp')
        f_all_repeat = os.path.join(working_dir, 'all_repeat.fasta')
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
            os.mkdir(tmp_dir)
        if directory:
            item_file = [os.path.join(directory, file_name) for file_name in os.listdir(directory)]
        elif file_list:
            item_file = [line.rstrip() for line in open(file_list)]

        lock = mp.Manager().Lock()
        args_list = []
        i = 0
        for file_name in item_file:
            args_list.append([file_name, seed_length, coverage_cutoff, identity_cutoff, expand_length, \
                os.path.join(tmp_dir, 'tmp_' + str(i)), f_all_repeat, lock])
            i += 1
        pool = mp.Pool(threads)
        pool.map(worker_func, args_list)
    return 0


if __name__ == '__main__':
    main()
