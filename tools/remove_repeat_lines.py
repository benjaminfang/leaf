import sys
a_line = [line.rstrip() for line in open(sys.argv[1])]

dt = {}
for line in a_line:
    if line[0] == '>':
        ff = line
        dt[ff] = {}
    elif line[0] == '^':
        head = line
        dt[ff][head] = []
    else:
        dt[ff][head].append(line)

dt_new = {}
for ff in dt:
    dt_new[ff] = {}
    for head in dt[ff]:
        dt_new[ff][head] = []
        block = set(dt[ff][head])
        dt_new[ff][head] = block


for ff in dt_new:
    print(ff)
    for head in dt_new[ff]:
        print(head)
        for line in dt_new[ff][head]:
            print(line)
                
