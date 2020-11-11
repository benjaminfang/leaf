#! /usr/bin/env python3

import os
import time
import argparse
import multiprocessing as mp
from biolib import Fasta_parser, biocodon


def get_args(args_list):
    args_parser = argparse.ArgumentParser(prog='long_repeat_finder', description='This \
        utility is used to find long repeating sequence among genomes.')
    args_parser.add_argument('-D', '--directory_fasta_file', default=None, help='directory path of fasta file.')
    args_parser.add_argument('-L', '--file_list', default=None, help='file list which want to process.')
    args_parser.add_argument('-S', '--seed_length', type=int, default=200, help='initial seed length.')
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
    if threads is None:
        threads = round(os.cpu_count() * 0.6)
    return directory, file_list, seed_length, identity_cutoff, coverage_cutoff, expand_len, threads


class Seed:
    def __init__(self, sequence):
        self.sequence = sequence
        self.sequence_len = len(self.sequence)
        self.seed_pointer = 1
        self.seed_island = [[1, self.sequence_len]]
        self.current_seed = None
        self.seed_sequence = None
        self.living_island = None

    def trim_seed_island(self, black_block):
        # black_block is a list[start, end], which will be trimed form seed isoland.
        block_start = black_block[0]
        block_end = black_block[1]
        updated_seed_island = []
        for island in self.seed_island:
            island_start = island[0]
            island_end = island[1]
            if block_end < island_start or block_start > island_end:
                updated_seed_island.append(island)
            elif block_start <= island_start and block_end >= island_end:
                pass
            elif block_start <= island_start and block_end >= island_start:
                updated_seed_island.append([block_end + 1, island_end])
            elif block_start <= island_end and block_end >= island_end:
                updated_seed_island.append([island_start, block_start - 1])
            elif block_start > island_start and block_end < island_end:
                updated_seed_island.append([island_start, block_start - 1])
                updated_seed_island.append([block_end + 1, island_end])
            else:
                raise Exception(black_block, island)
        self.seed_island = updated_seed_island

    def whether_seed_on_island(self, seed_pos):
        # seed_pos: [seed_start, seed_end]
        seed_start = seed_pos[0]
        seed_end = seed_pos[1]
        for island in self.seed_island:
            if island[0] <= seed_start and island[1] >= seed_end:
                return True
        return False

    def give_living_island(self, seed_pos):
        seed_start = seed_pos[0]
        seed_end = seed_pos[1]
        for island in self.seed_island:
            if island[0] <= seed_start and island[1] >= seed_end:
                return island
        return None

    def next_seed(self, seed_length):
        while True:
            seed_start = self.seed_pointer
            seed_end = seed_start + seed_length - 1
            if seed_end > self.sequence_len:
                return None
            if self.whether_seed_on_island([seed_start, seed_end]):
                self.current_seed = [seed_start, seed_end]
                self.living_island = self.give_living_island([seed_start, seed_end])
                self.seed_pointer = seed_end + 1
                self.seed_sequence = self.sequence[seed_start - 1: seed_end]
                break
            else:
                self.seed_pointer = seed_end + 1
        return self


def makeblastdb(file_name, dbname, directory):
    cmdname = 'makeblastdb'
    opt_dbtype = '-dbtype nucl'
    opt_in = '-in ' + file_name
    out_path = os.path.join(directory, dbname)
    opt_out = '-out ' + out_path
    cmd = ' '.join([cmdname, opt_dbtype, opt_in, opt_out, '&>/dev/null'])
    os.system(cmd)
    return out_path


def save_seq_as_fasta(seq, file_name, head_name, directory):
    f_path = os.path.join(directory, file_name)
    f_out = open(f_path, 'w')
    f_out.write('>' + head_name + '\n')
    f_out.write(seq + '\n')
    f_out.close()
    return f_path


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
        identity = round(float(line[2]) / 100, 5)
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
        coverage = round((q_e - q_s + 1) / query_seq_length, 5)
        dt_out.append([s_s, s_e, origntation, up_margine, down_margine, coverage, identity])
    return dt_out


def modify_blast_res_s(seed_blast_res, seed):
    seed_blast_res.sort(key=lambda x: x[0])
    seed_index = None
    margine = {'up': [], 'down': []}
    for ele in seed_blast_res:
        margine['up'].append(ele[3])
        margine['down'].append(ele[4])
    up_margine_max = max(margine['up'])
    down_margine_max = max(margine['down'])
    i = 0
    for ele in seed_blast_res:
        if ele[0] == seed.current_seed[0] and ele[1] == seed.current_seed[1]:
            seed_index = i
        i += 1
        if ele[2] == '+':
            ele[0] += up_margine_max - ele[3]
            ele[1] -= down_margine_max - ele[4]
        else:
            ele[0] += down_margine_max - ele[4]
            ele[1] -= up_margine_max - ele[3]
    return seed_blast_res, seed_index


class Fragment_manager:
    def __init__(self, fragment_list, whole_seq_length, seed_index, seed):
        self.fragment = fragment_list
        self.whole_seq_length = whole_seq_length
        self.seed_index = seed_index
        self.seed_living_island = seed.living_island
        self.up_lock = False
        self.down_lock = False

    def lock_up(self):
        self.up_lock = True

    def lock_down(self):
        self.down_lock = True

    def check_overlap(self):
        tmp = []
        i = 0
        for ele in self.fragment:
            tmp.append([ele[0], i])
            tmp.append([ele[1], i])
            i += 1
        tmp.sort()
        for i in range(len(tmp))[::2]:
            if tmp[i][1] != tmp[i + 1][1]:
                return True
        return False

    def calcu_expand_len(self):
        seed_fragment = self.fragment[self.seed_index]
        living_island = self.seed_living_island
        tmp = []
        expand_len_coll = {'up': [], 'down': []}
        i = 1
        for ele in self.fragment:
            if ele[2] == '+':
                tmp.append([i, ele[0], 'l', 'u'])
                tmp.append([i, ele[1], 'r', 'd'])
            else:
                tmp.append([i, ele[0], 'l', 'd'])
                tmp.append([i, ele[1], 'r', 'u'])
        tmp.sort(key=lambda x: x[1])
        if tmp[0][3] == 'u':
            expand_len_coll['up'].append(tmp[0][1] - 1)
        else:
            expand_len_coll['down'].append(tmp[0][1] - 1)
        if tmp[-1][3] == 'u':
            expand_len_coll['up'].append(self.whole_seq_length - tmp[-1][1])
        else:
            expand_len_coll['down'].append(self.whole_seq_length - tmp[-1][1])
        tmp = tmp[1:-1]
        for i in range(len(tmp))[::2]:
            dis = tmp[i + 1][1] - tmp[i][1]
            if tmp[i][3] != tmp[i + 1][3]:
                expand_len_coll['up'].append(dis)
                expand_len_coll['down'].append(dis)
            else:
                if tmp[i][3] == 'u':
                    expand_len_coll['up'].append(dis // 2)
                else:
                    expand_len_coll['down'].append(dis // 2)
        expand_len_coll['up'].append(seed_fragment[0] - living_island[0])
        expand_len_coll['down'].append(living_island[1] - seed_fragment[1])
        return {'up': min(expand_len_coll['up']), 'down': min(expand_len_coll['down'])}

    def expand_fragment(self, side, length):
        expand_len = self.calcu_expand_len()
        up_expansable_len = expand_len['up']
        down_expansable_len = expand_len['down']
        up_expand = None
        down_expand = None
        if side == 'up':
            up_expand = min(up_expansable_len, length)
            for ele in self.fragment:
                if ele[2] == '+':
                    ele[0] -= up_expand
                else:
                    ele[1] += up_expand
        else:
            down_expand = min(down_expansable_len, length)
            for ele in self.fragment:
                if ele[2] == '+':
                    ele[1] += down_expand
                else:
                    ele[0] -= down_expand
        return self, up_expand or down_expand


def provide_init_seed_match_list(seed, blastdb_path, seed_length, seed_identity_cutoff, seed_coverage_cutoff, directory, whole_seq_length):
    while True:
        seed = seed.next_seed(seed_length)
        # seed: Seed class instance
        if not seed:
            break
        seed_sequence_fasta_file = save_seq_as_fasta(seed.seed_sequence, 'seed_sequence.fasta', 'seed_sequence', directory)
        seed_blast_res = runblast(seed_sequence_fasta_file, 'blastn-short', blastdb_path, directory, 'seed_blast_res')
        blast_res = structure_blast_res(seed_length, seed_blast_res)
        init_seed_match_list = []
        for ele in blast_res:
            if ele[5] >= seed_coverage_cutoff and ele[6] >= seed_identity_cutoff:
                init_seed_match_list.append(ele)
        if len(init_seed_match_list) > 2:
            blast_res, seed_index = modify_blast_res_s(init_seed_match_list, seed)
            framgment_manger_instance = Fragment_manager([ele[:3] for ele in init_seed_match_list], whole_seq_length, seed_index, seed)
            if not framgment_manger_instance.check_overlap():
                yield framgment_manger_instance


def get_seq(piece, sequence):
    if piece[2] == '+':
        return sequence[piece[0] - 1: piece[1]]
    else:
        seq = [biocodon.base_complement[ele] for ele in sequence[piece[0] - 1: piece[1]]]
        seq.reverse()
        return ''.join(seq)


def run_water(first_file, second_file, directory):
    # using local alignment method.
    cmd_name = 'water'
    opt_gapopen = '10'
    opt_gapextend = '0.5'
    opt_outfile = os.path.join(directory, 'water.res')
    cmd = ' '.join([cmd_name, '-gapopen', opt_gapopen, '-gapextend', opt_gapextend, '-outfile', opt_outfile, first_file, second_file, '&>/dev/null'])
    os.system(cmd)
    return opt_outfile


def struc_water_res(water_file):
    a_lines = []
    f_in = open(water_file)
    for line in f_in:
        line = line.rstrip('\n')
        if len(line) > 0 and line[0] != '#':
            a_lines.append(line)
    f_in.close()
    tmp = []
    for i in range(len(a_lines))[::3]:
        tmp.append([a_lines[i], a_lines[i + 1], a_lines[i + 2]])
    dt = {'first_pos': [], 'second_pos': [], 'first_seq': '', 'second_seq': '', 'mid_seq': ''}
    for ele in tmp:
        first = ele[0].split()
        dt['first_pos'].append([int(first[1]), int(first[-1])])
        dt['first_seq'] += first[2]
        # NOTE: note not same mid sequence length actually.
        mid = ele[1][21: 71]
        dt['mid_seq'] += mid
        second = ele[2].split()
        dt['second_pos'].append([int(second[1]), int(second[-1])])
        dt['second_seq'] += second[2]
    dt['first_pos'] = [dt['first_pos'][0][0], dt['first_pos'][-1][1]]
    dt['second_pos'] = [dt['second_pos'][0][0], dt['second_pos'][-1][1]]
    return dt


def detect_boundary(frag, water_res, side, expand_len_act, score=15):
    def locat_expand_start(first_seq, second_seq, side, indent_length):
        indent_str_len = 0
        first_base_len = indent_length
        second_base_len = 0
        if side == 'up':
            for ele in first_seq:
                if ele != '-':
                    indent_length -= 1
                indent_str_len += 1
                if indent_length == 0:
                    break
            for ele in second_seq[: indent_str_len]:
                if ele != '-':
                    second_base_len += 1
        else:
            for ele in first_seq[::-1]:
                if ele != '-':
                    indent_length -= 1
                indent_str_len += 1
                if indent_length == 0:
                    break
            for ele in second_seq[-indent_str_len:]:
                if ele != '-':
                    second_base_len += 1
        return indent_str_len, first_base_len, second_base_len

    def detect_boundary(water_expand_seq, side, score=15):
        first_seq = water_expand_seq['first_seq']
        first_pos = water_expand_seq['first_pos']
        second_seq = water_expand_seq['second_seq']
        second_pos = water_expand_seq['second_pos']
        mid_seq = water_expand_seq['mid_seq']
        if side == 'up':
            i = 0
            mis_match = 0
            gap = 0
            score_lose = 0
            for ele in mid_seq[::-1]:
                i += 1
                if ele == ' ':
                    mis_match = 0
                    gap += 1
                    score_lose += gap * 1
                elif ele == '.':
                    mis_match += 1
                    gap = 0
                    score_lose += mis_match * 1
                elif ele == '|':
                    gap = 0
                    mis_match = 0
                if score_lose > score:
                    break
            while True:
                if mid_seq[-i] == '.' or mid_seq[-i] == ' ':
                    i -= 1
                else:
                    break
            for ele in first_seq[: -i]:
                if ele != '-':
                    first_pos[0] += 1
            for ele in second_seq[: -i]:
                if ele != '-':
                    second_pos[0] += 1
            return {'first_pos': first_pos, 'second_pos': second_pos}
        else:
            i = 0
            mis_match = 0
            gap = 0
            score_lose = 0
            for ele in mid_seq:
                i += 1
                if ele == ' ':
                    mis_match = 0
                    gap += 1
                    score_lose += gap * 1
                elif ele == '.':
                    gap = 0
                    mis_match += 1
                    score_lose += mis_match * 1
                elif ele == '|':
                    gap = 0
                    mis_match = 0
                if score_lose > score:
                    break
            while True:
                if mid_seq[i - 1] == '.' or mid_seq[i - 1] == ' ':
                    i -= 1
                else:
                    break
            for ele in first_seq[i:]:
                if ele != '-':
                    first_pos[1] -= 1
            for ele in second_seq[i:]:
                if ele != '-':
                    second_pos[1] -= 1
            return {'first_pos': first_pos, 'second_pos': second_pos}

    data_out = []
    first_frag_len = frag.fragment[frag.seed_index][1] - frag.fragment[frag.seed_index][0] + 1
    for water_ele in water_res:
        first_seq = water_ele['first_seq']
        first_pos = water_ele['first_pos']
        second_seq = water_ele['second_seq']
        second_pos = water_ele['second_pos']
        mid_seq = water_ele['mid_seq']
        if side == 'up':
            indent_length = expand_len_act - first_pos[0] + 1
            indent_str_len, first_base_len, second_base_len = locat_expand_start(first_seq, second_seq, 'up', indent_length)
            water_res_expand = {'first_seq': first_seq[: indent_str_len], 'first_pos': [first_pos[0], first_pos[0] + first_base_len - 1],
                'second_seq': second_seq[: indent_str_len], 'second_pos': [second_pos[0], second_pos[0] + second_base_len - 1],
                'mid_seq': mid_seq[: indent_str_len]}
            detected_boundary = detect_boundary(water_res_expand, 'up', score)
            first_boundary = [detected_boundary['first_pos'][0], first_pos[1]]
            second_boundary = [detected_boundary['second_pos'][0], second_pos[1]]
            data_out.append([first_boundary, second_boundary])
        else:
            indent_length = expand_len_act - (first_frag_len - first_pos[1])
            indent_str_len, first_base_len, second_base_len = locat_expand_start(first_seq, second_seq, 'down', indent_length)
            water_res_expand = {'first_seq': first_seq[-indent_str_len:], 'first_pos': [first_pos[1] - first_base_len + 1, first_pos[1]],
                'second_seq': second_seq[-indent_str_len:], 'second_pos': [second_pos[1] - second_base_len + 1, second_pos[1]],
                'mid_seq': mid_seq[-indent_str_len:]}
            detected_boundary = detect_boundary(water_res_expand, 'down', score)
            first_boundary = [first_pos[0], detected_boundary['first_pos'][1]]
            second_boundary = [second_pos[0], detected_boundary['second_pos'][1]]
            data_out.append([first_boundary, second_boundary])
    return data_out


def whether_lock_down(frag, boundarys, lock_base_cutoff):
    # only judge down side lock_down bool.
    tmp = []
    i = 0
    seed_frag_len = frag.fragment[frag.seed_index][1] - frag.fragment[frag.seed_index][0] + 1
    for ele in boundarys:
        first_gap = seed_frag_len - ele[0][1]
        second_gap = frag.fragment[i][1] - frag.fragment[i][0] + 1 - ele[1][1]
        tmp.append(min(first_gap, second_gap))
        i += 1
    dis_max = max(tmp)
    if dis_max > lock_base_cutoff:
        return True
    else:
        return False


def transite_location(frag, boundarys, side):
    first_tmp = []
    second_tmp = []
    if side == 'up':
        for ele in boundarys:
            first_tmp.append(ele[0][0] - 1)
            second_tmp.append(ele[1][0] - 1)
        first_tmp = first_tmp[: frag.seed_index] + first_tmp[frag.seed_index + 1:]
        seed_subtract = min(first_tmp)
        second_tmp = second_tmp[: frag.seed_index] + [seed_subtract] + second_tmp[frag.seed_index + 1:]
        i = 0
        for ele in frag.fragment:
            if ele[2] == '+':
                ele[0] += second_tmp[i]
            else:
                ele[1] -= second_tmp[i]
            i += 1
    else:
        i = 0
        seed_frag_len = frag.fragment[frag.seed_index][1] - frag.fragment[frag.seed_index][0] + 1
        for ele in boundarys:
            first_tmp.append(seed_frag_len - ele[0][1])
            second_tmp.append(frag.fragment[i][1] - frag.fragment[i][0] + 1 - ele[1][1])
            i += 1
        first_tmp = first_tmp[: frag.seed_index] + first_tmp[frag.seed_index + 1:]
        seed_subtract = min(first_tmp)
        second_tmp = second_tmp[: frag.seed_index] + [seed_subtract] + second_tmp[frag.seed_index + 1:]
        i = 0
        for ele in frag.fragment:
            if ele[2] == '+':
                ele[1] -= second_tmp[i]
            else:
                ele[0] += second_tmp[i]
            i += 1
    return 0


def align_fragment_water(frag, directory, item_sequence):
    dt_out = []
    seed_frag = frag.fragment[frag.seed_index]
    seed_frag_seq = get_seq(seed_frag, item_sequence)
    seed_frag_file = save_seq_as_fasta(seed_frag_seq, 'seed_frag_seq_water.fasta', 'seed_frag_seq', directory)
    for piece in frag.fragment:
        piece_seq = get_seq(piece, item_sequence)
        second_file = save_seq_as_fasta(piece_seq, 'second_file_water.fasta', 'second_seq', directory)
        water_res = run_water(seed_frag_file, second_file, directory)
        water_res = struc_water_res(water_res)
        dt_out.append(water_res)
    return dt_out


def thread_worker(args_in):
    file_name, seed_length, coverage_cutoff, identity_cutoff, expand_length, sub_tmp_dir, repeat_result_file, lock = args_in
    sub_tmp_dir = os.path.join(sub_tmp_dir, os.path.basename(file_name))
    if not os.path.exists(sub_tmp_dir):
        os.mkdir(sub_tmp_dir)
    data = {}
    lock_base_cutoff = 10
    file_basename = os.path.basename(file_name)
    data['file_name'] = file_basename
    data['head'] = {}
    fasta_dt = Fasta_parser(file_name)
    fasta_dt.join_lines()
    for head in fasta_dt.data:
        print('===============================================================')
        print(head)
        head_name = head.split()[0]
        data['head'][head_name] = []
        # NOTE: save item sequence to a file and make blast datebase.
        item_sequence = fasta_dt.data[head]
        whole_seq_length = len(item_sequence)
        sequence_file = save_seq_as_fasta(item_sequence, 'item_sequence.fasta', 'whole_sequecne', sub_tmp_dir)
        blastdb = makeblastdb(sequence_file, 'blastdb', sub_tmp_dir)
        # NOTE: construct Seed class instance
        seed = Seed(item_sequence)
        # NOTE: fragment_res is iterable, and yeild by provide_init_seed_match_list function.
        fragment_res = provide_init_seed_match_list(seed, blastdb, seed_length, identity_cutoff, coverage_cutoff, sub_tmp_dir, whole_seq_length)
        for frag in fragment_res:
            print('1 initial fragment:', frag.fragment)
            # NOTE: First exapnd up side by length of seed and lock up side.
            fragment_ins, expand_len_act = frag.expand_fragment('up', seed_length)
            print('2 fragment up elonged actually elonged length:', fragment_ins.fragment, expand_len_act)
            fragment_ins.lock_up()
            water_compare_res = align_fragment_water(fragment_ins, sub_tmp_dir, item_sequence)
            print('3 water align res:', water_compare_res)
            boundarys = detect_boundary(fragment_ins, water_compare_res, 'up', expand_len_act, 15)
            print('4 up boundarys fixed:', boundarys)
            transite_location(fragment_ins, boundarys, 'up')
            print('5 transite to abs positon:', fragment_ins.fragment)
            print('--------------------')
            while True:
                if fragment_ins.down_lock:
                    break
                fragment_ins, expand_len_act = fragment_ins.expand_fragment('down', expand_length)
                print('l1 fragment down elonged', fragment_ins.fragment, expand_len_act)
                if expand_len_act < expand_length:
                    fragment_ins.lock_down()
                water_compare_res = align_fragment_water(fragment_ins, sub_tmp_dir, item_sequence)
                print('l2 water align res:', water_compare_res)
                boundarys = detect_boundary(fragment_ins, water_compare_res, 'down', expand_len_act, 15)
                print('l3 boundarys', boundarys)
                if whether_lock_down(fragment_ins, boundarys, lock_base_cutoff):
                    fragment_ins.lock_down()
            transite_location(fragment_ins, boundarys, 'down')
            print('6 transieted positon:', fragment_ins.fragment)
            for ele in fragment_ins.fragment:
                seed.trim_seed_island([ele[0], ele[1]])
            data['head'][head_name].append(fragment_ins)
    lock.acquire()
    f_out = open(repeat_result_file, 'a')
    print('>' + data['file_name'], file=f_out)
    for cont in data['head']:
        print('>>' + cont, file=f_out)
        for ele in data['head'][cont]:
            print(ele.fragment, file=f_out)
    f_out.close()
    lock.release()
    return 0


def main(name='long_repeat_finder', args=None):
    myname = 'long_repeat_finder'
    if name == myname:
        directory, file_list, seed_length, coverage_cutoff, identity_cutoff, expand_length, threads = get_args(args)
        working_dir = 'working_dir_' + time.strftime('%Y%m%d%H%M%S')
        tmp_dir = os.path.join(working_dir, 'tmp')
        repeat_result_file = os.path.join(working_dir, 'repeats_result.fasta')
        if not os.path.exists(working_dir):
            os.mkdir(working_dir)
            os.mkdir(tmp_dir)
        if directory:
            item_file = [os.path.join(directory, file_name) for file_name in os.listdir(directory)]
        elif file_list:
            item_file = [line.rstrip() for line in open(file_list)]

        lock = mp.Manager().Lock()
        args_list = []
        for file_name in item_file:
            args_list.append([file_name, seed_length, coverage_cutoff, identity_cutoff, expand_length,
                tmp_dir, repeat_result_file, lock])
        pool = mp.Pool(threads)
        pool.map(thread_worker, args_list)
    return 0


def test_thread_worker(args):
    thread_worker(args)
    return 0


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 3 and sys.argv[1] == 'test':
        if not os.path.exists('tmp_test'):
            os.mkdir('tmp_test')
        lock = mp.Manager().Lock()
        test_thread_worker([sys.argv[2], 200, 0.9, 0.9, 100, 'tmp_test', 'test_res', lock])
    else:
        main()
