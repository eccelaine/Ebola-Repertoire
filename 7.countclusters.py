count = {}

with open('10Xout_formatted.csv') as fin, open('10Xout_counted.csv', 'w') as fout:
    for line in fin:
        ls = line.strip().split(',')
        if ls[0] not in count:
            count[ls[0]] = 0
        count[ls[0]] += 1

    fin.seek(0,0)
    FIRST = True
    for line in fin:
        if FIRST:
            fout.write('cluster_id,function,mongo_id,heavy_id,heavy_v,heavy_j,heavy_cdr3,light_id,light_v,light_j,light_cdr3,heavy_raw,light_raw,count\n')
#             fout.write(line.strip() + ',count\n')
            FIRST = False
        else:
            ls = line.strip().split(',')
            fout.write(line.strip() + ',' + str(count[ls[0]]) + '\n')
