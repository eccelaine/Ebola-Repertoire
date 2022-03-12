# Ebola Repertoire
This repo contains the clustering code used for the manuscript: Systematic analysis of human antibody response to ebolavirus glycoprotein reveals high prevalence of neutralizing public clonotypes.

# Clustering to identify public clonotypes
1. Start with an input file of all sequences in a csv with headers "ID, HC, LC"

2. First run your input csv through "1.csvtofasta.py"
   This will convert your csv into a fasta file. 
   
3. Taking the fasta files of the heavy and light chain, run it through PyIR (https://github.com/crowelab/PyIR)
   Once you have PyIR working in your shell, use the command as follows
   > pyir-plus --legacy XXXX.fasta
   
4. Then using the two fasta files you just created, run them through "2.combineHCLC-pyiroutput.py". 
   This script should give you an ouput of a csv
   
5. Taking the output of script 2, now run the csv from step 4 through "3.makeclusteringfiles.py". 
   This script should give you two XXX.dat files as output
   
6. Run these two XXX.dat files through "4.clustering.py" using the below command: 
   > python2 Clusteringv2.1b_mod7.py --full-formatted-file=XXX.dat --keys-for-clustering 0 3 4 6 --clustered-formatted-outfile=clustered-XXX.dat.gz --sequence-identity=0.80 --clustering-type=complete

7. Take the two XXX.dat.gz files, unzip them, and run them through "5.combineHCLCclustering.py". 
   This should give you an output file of 10Xout.csv
   If just clustering on the heavy chain sequences (ex: bulk sequencing) run them through "5.combineHConlyclustering.py"
   
8. Take your 10Xout.csv and run it through "6.sortclusters"
   This should give you an output file of 10Xout_formatted.csv

9. Take your 10Xout_formateed.csv and run it through "7.countclusters.py"
   This should give you an output file of 10Xout_counted.csv as your final output.
