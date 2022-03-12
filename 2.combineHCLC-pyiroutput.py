import argparse
import gzip
import json

INPUT_HEAVY = 'DavidHO-COV2-HC.json.gz'
INPUT_LIGHT = 'DavidHO-COV2-LC.json.gz'
OUT = 'DAVIDHO-paired.csv'
RUN = 'x'
SAMPLE = 'x'

#parser = argparse.ArgumentParser()
#parser.add_argument('heavy')
#parser.add_argument('light')
#parser.add_argument('outfile')
#parser.add_argument('--run')
#parser.add_argument('--sample')
#args = parse.parse_args()

d = {}


with gzip.open(INPUT_HEAVY, 'rt') as fin:
    for line in fin:
        j = json.loads(line.strip())
        if j['Sequence ID'] == 'ADI-18903':
            print('Hello heavy!!!!')
        d[j['Sequence ID']] = {
                'heavy': j
        }

num_lines = 0
with gzip.open(INPUT_LIGHT, 'rt') as fin:
    for line in fin:
        j = json.loads(line.strip())
        if j['Sequence ID'] == 'ADI-18903':
            print('Hello!!!!')
        d[j['Sequence ID']]['light'] = j
        num_lines += 1
        
print(num_lines)


with open(OUT, 'w') as fout:
    fout.write('run_id,sample_id,replicate_id,clonotype_id,heavy_label,heavy_v_gene,heavy_j_gene,heavy_isotype,heavy_cdr3_(aa),' +
                'heavy_cdr3_aa_length,heavy_percent_id,heavy_occurrences,heavy_mutations,heavy_d_gene,' +
                'heavy_cdr3_(dna),heavy_nt_trimmed,heavy_nt,heavy_aa,light_label,light_v_gene,light_j_gene,' +
                'light_isotype,light_cdr3_(aa),light_cdr3_aa_length,light_percent_id,light_occurrences,' +
                'light_mutations,light_d_gene,light_cdr3_(dna),light_nt_trimmed,light_nt,light_aa\n')
    for seqid in d:
        print(seqid)
        heavy = d[seqid]['heavy']
        light = d[seqid]['light']
        fout.write(RUN + ',' + SAMPLE + ',N/A,' + seqid + ',' + seqid + '_heavy,')
        fout.write(heavy['V family'] + ',')
        fout.write(heavy['J family'] + ',')
        fout.write('N/A,')
        fout.write(heavy.get('CDR3', {}).get('AA', "") + ',')
        fout.write(str(len(heavy.get('CDR3', {}).get('AA', ""))) + ',')
        fout.write(str(heavy.get('Total', {}).get('percent identity', "")) + ',')
        fout.write('1,')
        fout.write(str(heavy.get('Total', {}).get('mismatches', "")) + ',')
        fout.write(heavy.get('D family', 'None') + ',')
        fout.write(heavy.get('CDR3', {}).get('NT', "") + ',')
        fout.write(heavy['NT-Trimmed'] + ',')
        fout.write(heavy['Raw Sequence'] + ',')
        fout.write(heavy.get('AA', 'None') + ',')
        fout.write(seqid + '_light,')
        fout.write(light['V family'] + ',')
        fout.write(light['J family'] + ',')
        fout.write('N/A,')
        fout.write(light.get('CDR3', {}).get('AA', "") + ',')
        fout.write(str(len(light.get('CDR3', {}).get('AA', ""))) + ',')
        fout.write(str(light.get('Total', {}).get('percent identity', "")) + ',')
        fout.write('1,')
        fout.write(str(light.get('Total', {}).get('mismatches', "")) + ',')
        fout.write(light.get('D family', 'None') + ',')
        fout.write(light.get('CDR3', {}).get('NT', "") + ',')
        fout.write(light['NT-Trimmed'] + ',')
        fout.write(light['Raw Sequence'] + ',')
        fout.write(light.get('AA', 'None') + '\n')

