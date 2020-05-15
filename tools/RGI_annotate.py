#!
import sys


def structure_repeat_data(file_name):
    def structure_line(line):
        data_out = []
        tmp = line.split('\t')
        for i in list(range(len(tmp) - 3))[::3]:
            data_out.append([int(tmp[i]), int(tmp[i + 1]), tmp[i + 2]])
        data_out.append(int(tmp[-3]))
        data_out += tmp[-2:]
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


def output_data(data):
    for strain in data:
        for cont in data[strain]:
            for rep in data[strain][cont]:
                print('>' + strain + '\t' + cont + '\t' + '\t'.join(rep[0][:2]))
                print(rep[-2])
    return 0


def annotate(rgi_fasta_file):
    pass


if __name__ == '__main__':
    file_name = sys.argv[1]
    data = structure_repeat_data(file_name)
    output_data(data)
    pass



