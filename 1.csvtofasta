import os
import re
import subprocess
import csv


for f in files:
    with open(f,'r') as fin, open(f.split('.')[0] + '-HC.fasta', 'w') as fout:
        reader = csv.DictReader(fin)
        for row in reader:
            fout.write('>' + row['ID'] + '\n')
            fout.write(row['HC'] + '\n')

    with open(f,'r') as fin, open(f.split('.')[0] + '-LC.fasta', 'w') as fout:
        reader = csv.DictReader(fin)
        for row in reader:
            fout.write('>' + row['ID'] + '\n')
            fout.write(row['LC'] + '\n')
