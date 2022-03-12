#!/usr/bin/env python
import time
import sys
import cPickle
import gzip
import csv
import os
import argparse
import numpy as np
from scipy.spatial import distance
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fclusterdata
from itertools import dropwhile
import timeit
import multiprocessing as mp, sys
from tempfile import NamedTemporaryFile
try:
        from Bio import SeqIO
        from Bio.Seq  import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.SeqIO.QualityIO import PairedFastaQualIterator
except ImportError:
        print("Biopython is not installed or is not in the path of python")
        sys.exit(1)

''' -------------------------------------------- '''
''' Function comes from Jordan Callicoat's post '''
''' http://code.activestate.com/recipes/363051/ '''
''' -------------------------------------------- '''
def flatten(l, ltypes=(list, tuple)):
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

class ClusterCdr3():

    def __init__(self, Commands):
         self.commands=Commands
         self.Cdr3Sequences=[]
         self.Cdr3Labels=[]
         self.distancematrix=[]
         self.distancetable=[]
         self.linkagematrix=[]
         self.SequenceClusters=[]
         if isinstance(self.commands['FastaSequences'],dict):
             for k,v in self.commands['FastaSequences'].items():
                 self.Cdr3Sequences.append(v)
                 self.Cdr3Labels.append(k)
             # self.distancematrix=np.array(self._distancematrix())
             #print('Before constructing distance matrix')
             #t0 = time.time()
             self.distancematrix = self._distancematrix_rss()
             #t1 = time.time()
             #print('After constructing distance matrix')
             #print('Distance matrix time = ', t1-t0)

             #print('Before constructing linkage matrix')
             #t0 = time.time()
             self.linkagematrix = linkage(self.distancematrix,method=self.commands['method'])
             #t1 = time.time()
             #print('After constructing linkage matrix')
             #print('Linkage matrix time = ', t1-t0)

             self._dumpLinkageTable()
             self._distancetable()
         elif isinstance(self.commands['FastaSequences'],list):
            for item in self.commands['FastaSequences']:
                 self.Cdr3Sequences.append(str(item.seq))
                 self.Cdr3Labels.append(item.name)
            self.distancematrix=np.array(self._distancematrix())
            self.linkagematrix = linkage(self.distancematrix,method=self.commands['method'])
            self._dumpLinkageTable()
            self._distancetable()
         else:
            pass

    def hamming(self,u, v):
        string=''
        if len(u) == len(v):
            string=u
        else:
            string=min(u,v, key=len)
        Sum=0
        try:
            for i,value in enumerate(string):
                if u[i]!=v[i]:
                    Sum=Sum+1
                else:
                    Sum=Sum+0
        except IndexError:
            print("CDR3s are not the same length")
        return float(Sum)
        
    def hamming_normalized(self,u, v):
        string=''
        if len(u) == len(v):
            string=u
        else:
            string=min(u,v, key=len)
        Sum=0
        Length=len(u)
        try:
            for i,value in enumerate(string):
                if u[i]!=v[i]:
                    Sum=Sum+1
                else:
                    Sum=Sum+0
        except IndexError:
            print("CDR3s are not the same length")
        return float(Sum)/float(Length)
        
    def _distancematrix(self):
        Temp=[]
        for i,cdr3a in enumerate(self.Cdr3Sequences):
            for j,cdr3b in enumerate(self.Cdr3Sequences):
                if j>i:
                    Temp.append(self.hamming_normalized(cdr3a,cdr3b))
        return Temp

    def _distancematrix_rss(self):
	    nseq = len(self.Cdr3Sequences)
	    nelem = (nseq * (nseq - 1)) / 2
	    Temp = np.zeros(int(nelem))
	    count = 0
	    for i,cdr3a in enumerate(self.Cdr3Sequences):
		    for j,cdr3b in enumerate(self.Cdr3Sequences):
			    if j>i:
				    Temp[count] = self.hamming_normalized(cdr3a,cdr3b)
				    count += 1
	    return Temp

    def _dumpLinkageTable(self):
        otu_table = open(self.commands['otu-table'],'w')
        N=len(self.Cdr3Sequences)
        M=N-1
        for i,entry in enumerate(self.linkagematrix):
            if self.linkagematrix[i][0] >= N:
                A=M-self.linkagematrix[i][0]
            else:
                A=self.linkagematrix[i][0]

            if self.linkagematrix[i][1]>=N:
                B=M-self.linkagematrix[i][1]
            else:
                B=self.linkagematrix[i][1]
            otu_table.write('%d,%d,%d,%.4f,%d\n'%(-1*i-1,A,B,self.linkagematrix[i][2],self.linkagematrix[i][3]))
        otu_table.close()	     

    def _distancetable(self):
        Temp=[]
        Nodes=[]
        Distances=[]
        Clusters=[]
        ClustersFlattened=[]
        Size=len(self.linkagematrix)
        N=Size
        '''---------------------------------------------------'''
        ''' Generate three arrays for indexing clustered data '''
        '''---------------------------------------------------'''
        for entry in self.linkagematrix:
            Size=Size+1
            if entry[0] > N and entry[1] > N:
                #print('entry[1] and entry[2] are larger than N')
                Index0=Nodes.index(entry[0])
                Index1=Nodes.index(entry[1])
                Clusters.append( (Clusters[Index0],Clusters[Index1]) )        
            elif entry[0] <=N and entry[1] > N:
                #print('entry[1] is larger than N')
                Index1=Nodes.index(entry[1])
                Clusters.append( (entry[0],Clusters[Index1]) )
            elif entry[0] > N and entry[1] <=N:
                #print('entry[0] is larger than N')
                Index0=Nodes.index(entry[0])
                Clusters.append( (Clusters[Index0],entry[1]) )
            else:
                Clusters.append( (entry[0],entry[1]) )
            Nodes.append(Size)
            Distances.append(entry[2])
        
        for entry in Clusters:
            ClustersFlattened.append( flatten(entry) )

        for i,j in zip(ClustersFlattened,Distances):
            self.distancetable.append((i,j))

    def _clustering(self):
        Temp=[]
        Clusters=[]
        Exclude=[]
        SequenceClusters=[]
        '''---------------------------------------------'''
        '''Determine all disances below some threshold  '''
        '''---------------------------------------------'''
        for i,value in enumerate(self.distancetable):
            if value[1] <=self.commands['cutoff']:
                Temp.append(value[1])
                Clusters.append(value[0])
        """
        Create an exlusion list to remove redundant clustering
        """
        Exclude=[]
        for i,v1 in enumerate(Clusters):
            if self.distancetable[i][1] <=self.commands['cutoff']:
                for j,v2 in enumerate(Clusters):
                    if j>i:
                        if len(set(v1).intersection(v2)):
                            if len(v1)>len(v2):
                                #print('Exclude redundant cluster',v2,j)
                                Exclude.append(j)
                            else:
                                #print('Exclude redundant cluster',v1,i)
                                Exclude.append(i)
                                
        '''-------------------------------'''
        ''' This is the end of our search '''
        '''-------------------------------'''
        ClusterNumber=0        
        for i,val in enumerate(Clusters):
            if i not in Exclude:
                try:
                    ''' I incremented the self.distancetable from [i] [i+1]'''
                    #print('Cluster',ClusterNumber,len(val), self.distancetable[i+1][1])
                    Sequences={}
                    for j in val:
                        idx=int(j)
                        #print(">%s\n%s"%(self.Cdr3Labels[idx],self.Cdr3Sequences[idx]))
                        Sequences[self.Cdr3Labels[idx]]=self.Cdr3Sequences[idx]
                    self.SequenceClusters.append(Sequences)
                    ClusterNumber=ClusterNumber+1
                except IndexError:
                    print('No clusters below threshold')
          
    def GetClustering(self):
        if self.SequenceClusters:
            return self.SequenceClusters
        else:
            self._clustering()
            return self.SequenceClusters
            
    def GetNonClustering(self):
        Unclustered={}
        if self.SequenceClusters:
            pass
        else:
            self._clustering
            
        for i,item in enumerate(self.Cdr3Labels):
            if any(item in k for k in self.SequenceClusters):
                pass
            else:
                Unclustered[self.Cdr3Labels[i]]=self.Cdr3Sequences[i]
        return Unclustered
        
    def getMatrix(self):
        return self.distancematrix
        
    def getDistanceTable(self):
        return self.distancetable
        
    def getLabels(self):
        return self.Cdr3Labels
        
    def getSequences(self):
        return self.Cdr3Sequences
        
    def getLinkageMatrix(self):
        return self.linkagematrix

#               0           1    2      3       4        5                6
#5cc8b743513b80c04c05206b year1 day0 IGHV3-23 IGHJ1 TNLKNGAYPH ACTAACCTGAAGAATGGGGCCTACCCCCAC
def cluster2(filename,chunk,sequence_identity,clustering_method):
    '''---------------------------------------'''
    ''' Function should operate on large file '''
    '''---------------------------------------'''
    if  os.path.splitext(filename)[-1] in [ '.gz', '.gzip']:
      file_handle = gzip.open(filename, "r")
    else:
      file_handle = open(filename, "r")
    otu_table=[]
    outfile=""
    RunningDict={}
    file_handle.seek(chunk[0])
    while(file_handle.tell() < chunk[1]):
        linestring=file_handle.readline().rstrip()
        data=linestring.split(' ')
        cdr3_length=len(data[3])
        RunningDict[data[0]]=data[3]
    ''' I do not check to ensure that each CDR3 is the same length '''
    if len(RunningDict)>1:
       normalized_hamming=(cdr3_length-int(round(sequence_identity*(cdr3_length))))/(float (cdr3_length))
       cdr3=ClusterCdr3({'FastaSequences':RunningDict,
                         'method':clustering_method,
                         'cutoff':normalized_hamming,
                         'otu-table':'erase.dat'})
       ''' ---------------------------------'''
       ''' Dump out table showing  clusters '''
       ''' ---------------------------------'''
       for cluster in cdr3.GetClustering():
          otu_table.append(cluster.keys())
       #otu_table.append(cdr3.GetNonClustering().keys())
    return  otu_table

def cluster(filename,chunk,sequence_identity,clustering_method):
    '''---------------------------------------'''
    ''' Function should operate on large file '''
    '''---------------------------------------'''
    if  os.path.splitext(filename)[-1] in [ '.gz', '.gzip']:
      file_handle = gzip.open(filename, "r")
    else:
      file_handle = open(filename, "r")

    '''-----------------------------------------'''
    ''' Store sequences in lists containing UMS '''
    '''-----------------------------------------'''
    otu_table=[]
    outfile=""
    RunningDict={}
    file_handle.seek(chunk[0])
    cdr3_length=0
    while(file_handle.tell() < chunk[1]):
        linestring=file_handle.readline().rstrip()
        data=linestring.split(' ')
        ''' -------------------------------------------------------------------'''
        ''' Store the sequences to avoid making routine dependent on Biopython '''
        '''--------------------------------------------------------------------'''
        cdr3_length=len(data[4])
        RunningDict[data[0]]=data[4]
    normalized_hamming=(cdr3_length-int(round(sequence_identity*(cdr3_length))))/(float (cdr3_length))
    cdr3=ClusterCdr3({'FastaSequences':RunningDict,
                      'method':clustering_method,
                      'cutoff':normalized_hamming,
                      'otu-table':'erase.dat'})
    ''' ---------------------------------'''
    ''' Dump out table showing  clusters '''
    ''' ---------------------------------'''
    for cluster in cdr3.GetClustering():
       otu_table.append(cluster.keys())
    #otu_table.append(cdr3.GetNonClustering().keys())
    file_handle.close()
    return otu_table

'''====================================================='''
''' Functions specific for handling the cluster files   ''' 
'''====================================================='''

def readin_lineage_file(table,filename):
   '''
    This function will read in a probability file. I have put in a numerical 
    precision check to ensure that the probabilities sum to 1.0 for each label.
    I did not include a header file, but this is something I can add later on.
   '''
   
   #5cc8b747513b80c04c052e54 year1 day0 IGHV3-15 IGHJ1 ATSLIQGHTLATD
   counter=0
   temp=[]
   if  os.path.splitext(filename)[-1] in [ '.gz', '.gzip']:
        with gzip.open(filename,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     if len(temp)>1:
                        table['cluster-%10.10d'%(counter)]=[]
                        table['cluster-%10.10d'%(counter)]=temp
                     else:
                        table['singleton-%10.10d'%(counter)]=[]
                        table['singleton-%10.10d'%(counter)]=temp
                     counter=counter+1
                     del temp
                     temp=[]
                else:
                     temp.append(' '.join(map(str,row)))
            csv_file.close()
   else:
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     if len(temp)>1:
                        table['cluster-%10.10d'%(counter)]=[]
                        table['cluster-%10.10d'%(counter)]=temp
                     else:
                        table['singleton-%10.10d'%(counter)]=[]
                        table['singleton-%10.10d'%(counter)]=temp
                     counter=counter+1
                     del temp
                     temp=[]
                else:
                     temp.append(' '.join(map(str,row)))
            csv_file.close()
   return table

def populate_table(table,filename):
   if  os.path.splitext(filename)[-1] in [ '.gz', '.gzip']:
        with gzip.open(filename,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     pass
                else:
                     label='%s %s %s'%(row[3],row[4],row[6])
                     if label in table:
                         table[label].append(row[0])
                     else:
                         table[label]=[]
                         table[label].append(row[0])
            csv_file.close()
   else:
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     pass
                else:
                     label='%s %s %s'%(row[3],row[4],row[6])
                     if label in table:
                         table[label].append(row[0])
                     else:
                         table[label]=[]
                         table[label].append(row[0])
            csv_file.close()
   return table

def table_of_mongoids(table,filename):
   if  os.path.splitext(filename)[-1] in [ '.gz', '.gzip']:
        with gzip.open(filename,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     pass
                else:
                     label='%s'%(row[0])
                     if label in table:
                         #table[label].append(' '.join(map(str,row)))
                         print('You must have unique MongoID')
                         pass
                     else:
                         #table[label]=[]
                         #table[label].append(' '.join(map(str,row)))
                         table[label]=' '.join(map(str,row))
            csv_file.close()
   else:
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     pass
                else:
                     label='%s'%(row[0])
                     if label in table:
                         #table[label].append(' '.join(map(str,row)))
                         pass
                     else:
                         #table[label]=[]
                         #table[label].append(' '.join(map(str,row)))
                         table[label]=' '.join(map(str,row))
            csv_file.close()
   return table


def read_in_formatted_file(filename):
   '''
    This function will read in a probability file. I have put in a numerical 
    precision check to ensure that the probabilities sum to 1.0 for each label.
    I did not include a header file, but this is something I can add later on.
   '''
   table={}
   if  os.path.splitext(filename)[-1] in [ '.gz', '.gzip']:
        with gzip.open(filename,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] != 'END':
                    if row[0] in table:
                        pass
                    else:
                        table[row[0]]=row[1:]
            csv_file.close()
   else:
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ') 
            for row in csv_reader:
                if row[0] != 'END':
                    if row[0] in table:
                        pass
                    else:
                        table[row[0]]=row[1:]
            csv_file.close()
   #print table
   #sys.exit(1)
   return table

def create_mongo_id_lookup_table(table):
   '''
    This function stores all the MongoID for each V3J and is tied to the starting table
   '''
   lookup_table={}
   for key in table.keys():
        label='%s %s %s'%(table[key][2],table[key][3],table[key][5])
        if label in lookup_table:
            lookup_table[label].append(key)
        else:
            lookup_table[label]=[]
            lookup_table[label].append(key)
   for key in lookup_table.keys():
       lookup_table[key]=list(set(lookup_table[key]))
   return lookup_table

def combine_formatted_files(table,filename):
   '''
    This function will read in a probability file. I have put in a numerical 
    precision check to ensure that the probabilities sum to 1.0 for each label.
    I did not include a header file, but this is something I can add later on.
   '''
   #5cc8b747513b80c04c052e54 year1 day0 IGHV3-15 IGHJ1 ATSLIQGHTLATD
   if  os.path.splitext(filename)[-1] in [ '.gz', '.gzip']:
        with gzip.open(filename,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     pass
                else:
                     label='%s_%s_%d'%(row[3],row[4],len(row[5]))
                     if label in table:
                         table[label].append(row)
                     else:
                         table[label]=[]
                         table[label].append(row)
            csv_file.close()
   else:
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     pass
                else:
                     label='%s_%s_%d'%(row[3],row[4],len(row[5]))
                     if label in table:
                         table[label].append(row)
                     else:
                         table[label]=[]
                         table[label].append(row)
            csv_file.close()
   return table


def deduplicate_table(table,integer_keys):
    reduced=NamedTemporaryFile(suffix='.cluster',dir='./',delete=False)
    deduplicated={}
    reduce_keys=integer_keys[1:]
    '''------------------------------------'''
    ''' You now have a table wich V-J-CDR3 '''
    '''------------------------------------'''
    for key in table.keys():
        '''--------------------------------'''
        ''' You now have a list of entries '''
        '''--------------------------------'''
        deduplicated[key]=[]
        for entry in table[key]:
            label='%s'%(' '.join(map(str,[entry[num] for num in reduce_keys])))
            if label not in deduplicated[key]:
                string=' '.join(map(str,entry))
                deduplicated[key].append(string)
            else:
                pass
    for key in deduplicated.keys():
        for line in deduplicated[key]:
            reduced.write('%s\n'%line)
        reduced.write('END\n')
    reduced.close()
    return reduced.name 

def dump_out_cluster_file(filename,integer_keys):
   '''
    This function will read in a probability file. I have put in a numerical 
    precision check to ensure that the probabilities sum to 1.0 for each label.
    I did not include a header file, but this is something I can add later on.
   '''
   clustered_outfile=NamedTemporaryFile(suffix='.cluster',dir='./',delete=False)
   #5cc8b747513b80c04c052e54 year1 day0 IGHV3-15 IGHJ1 ATSLIQGHTLATD
   if  os.path.splitext(filename)[-1] in [ '.gz', '.gzip']:
        with gzip.open(filename,'rb') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     clustered_outfile.write('END\n')
                else:
                     clustered_outfile.write('%s\n'%(' '.join(map(str,[row[num] for num in integer_keys]))))
            csv_file.close()
   else:
        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ')
            for row in csv_reader:
                if row[0] == 'END':
                     clustered_outfile.write('END\n')
                else:
                     clustered_outfile.write('%s\n'%(' '.join(map(str,[row[num] for num in integer_keys]))))
            csv_file.close()
   clustered_outfile.close()
   return clustered_outfile.name

def main():
    t0 = time.time()
    cluster={}
    usage = "usage: %prog [options] arg\n./cluster-cdr3.py --fasta-sequences=<file> --cdr3-length=<int> --sequence-identity-threshold=<float>"
    cluster = argparse.ArgumentParser(description="Will cluster cdr3 sequences of the same length")
    cluster.add_argument("--full-formatted-file",
                        action="store",
                        dest="full-formatted-file",
                        type=str,
                        help='FASTA formatted file containing either DNA or AA sequences')
    cluster.add_argument("--fasta-file",
                        action="store",
                        dest="fasta-file",
                        type=str,
                        help='FASTA formatted file containing either DNA or AA sequences')
    cluster.add_argument("--user-fasta-file",
                        action='store_true',
                        dest="user-fasta-file",
                        help='user defined their own formatted file',
                        default=False)
    cluster.add_argument("--user-formatted-file",
                        action='store_true',
                        dest="user-formatted-file",
                        help='user defined their own formatted file',
                        default=False)
    cluster.add_argument("--formatted-file",
                        action="store",
                        dest="formatted-file",
                        type=str,
                        help='FASTA formatted file containing either DNA or AA sequences')
    cluster.add_argument("--clustered-formatted-outfile",
                        action="store",
                        dest="clustered-formatted-outfile",
                        type=str,
                        help='FASTA formatted file containing either DNA or AA sequences')
    cluster.add_argument("--cdr3-length",
                        action="store",
                        type=int,
                        dest="cdr3-length",
                        help='length of sequence (all sequences must be of the same length)')
    cluster.add_argument("--sequence-identity",
                        action="store",
                        dest="sequence-identity",
                        type=float, 
                        help='sequence identity threshold',
                        default=0.90)
    cluster.add_argument("--linkage-file",
                        action="store",
                        dest="linkage-file",
                        type=str, 
                        help='linkage table',
                        default='linkage.table')
    cluster.add_argument("--number-processors",
                        action="store",
                        dest="number-processors",
                        type=int,
                        help='total number of processors to use',
                        default=40)
    cluster.add_argument('--keys-for-clustering',
                         action='store',
                         dest='keys-for-clustering',
                         nargs="*",
                         type=int,
                         default=[0,3,4,6])
    cluster.add_argument("--use-otu-file",
                  action='store_true',
                  help='do not dump out clustered files in FASTA file',
                  default=False)
    cluster.add_argument("--clustering-type",
                        action="store",
                        dest="clustering-type",
                        type=str,
                        help='type of clustering (single,complete or average)',
                        choices=['single', 'complete', 'average'],
                        default='complete')

    options = vars(cluster.parse_args())
    '''-----------------------------------'''
    ''' Read in sequences from FASTA file '''
    '''-----------------------------------'''
    #FastaSequences = SeqIO.index(options['fasta-file'], "fasta")

    #if len(FastaSequences)<2:
    #    sys.exit(1)

    '''--------------------------'''
    '''Read in the length of CDR '''
    '''--------------------------'''
    #Cdr3Length=options['cdr3-length']

    '''--------------------------------------------------'''
    ''' Identity will be converted to a Hamming Distance '''
    '''--------------------------------------------------'''
    SequenceIdentity=options['sequence-identity']

    '''----------------------------'''
    '''Normalized Hamming Distance '''
    '''----------------------------'''
    #NormalizedHamming=(Cdr3Length-int(round(SequenceIdentity*(Cdr3Length))))/(float (Cdr3Length))

    if options['user-formatted-file']==True:
        '''deuplicate the file '''
        cluster_file=options['formatted-file']
        sys.exit(1)
    elif options['user-fasta-file']==True:
        print('Currently not supported')
        sys.exit(1)
        FastaSequences = SeqIO.index(options['fasta-file'], "fasta")
        '''Convert these into the proper format'''
    else:
        table={}
        combine_formatted_files(table,options['full-formatted-file'])
        reduced_table=deduplicate_table(table,options['keys-for-clustering'])
        cluster_file=dump_out_cluster_file(reduced_table, options['keys-for-clustering'])
    '''-------------------------------------'''
    ''' Read in a correctly formatted table '''
    '''-------------------------------------'''
    global_position=0
    #fh=open(options['formatted-file'])
    fh=open(cluster_file)
    locations=[]
    anchor=0
    #filesize=os.path.getsize(options['formatted-file'])
    filesize=os.path.getsize(cluster_file)
    while( global_position <= filesize+1 ):
        row = fh.readline().rstrip()
        if row == 'END':
           locations.append((anchor,fh.tell()-4))
           anchor=fh.tell()
        else:
           global_position = global_position + len(row) + 1
    fh.close()

    '''-----------------------------------'''
    ''' Multiprocess and break up the job '''
    '''-----------------------------------'''
    jobs=[]  
    pool = mp.Pool(options['number-processors'])
    for chunk in locations:
        job = pool.apply_async(cluster2,(cluster_file,chunk,options['sequence-identity'],options['clustering-type']))
        jobs.append(job)

    '''------------------'''
    ''' Finish things up '''
    '''------------------'''
    output = []
    for job in jobs:
        try:
           output.append( job.get() )
        except ValueError:
           pass
    pool.close()
    '''------------------------------------------'''
    ''' Store clusters in a list of dictionaries '''
    '''------------------------------------------'''
    otu_table=[]
    for line in output:
        for entry in line:
            otu_table.append(entry)

    '''-------------------------------------------------------------------------------'''
    ''' Dump out DNA data into a flat file that is similar to the one used originally '''
    '''             The DNA sequences should be from each cluster (buggy if no clusters)                     '''
    '''-------------------------------------------------------------------------------'''
    if len(otu_table) > 0:
       clustered_outfile=gzip.open(options['clustered-formatted-outfile'],'wb')
       table=read_in_formatted_file(options['full-formatted-file'])
       lookup_table=create_mongo_id_lookup_table(table)
       clusteredlabels=[]
       for cluster in otu_table:
           for mongoid in cluster:
                clustered_outfile.write('%s %s\n'%(mongoid,' '.join(table[mongoid])))
                clusteredlabels.append(mongoid)
           clustered_outfile.write('END\n')
       clustered_outfile.close()
       print('====> Done finding all clusters\n')
       '''================================================'''
       ''' This next section is very messy and needs help '''
       '''================================================'''
       table={}
       #{'IGHV4-39 IGHJ3 GCGAGGCGCGGGAACTTCATGCCCCTTGATGCTTTTGATTTT': ['STAU-469']
       populate_table(table,options['full-formatted-file'])

       table_mongoids={}
       #'STAU-466': 'STAU-466 new new IGHV4-39 IGHJ4 ARRGNFMPLDAFDF
       table_of_mongoids(table_mongoids,options['full-formatted-file'])
       
       lineage={}
       readin_lineage_file(lineage,options['clustered-formatted-outfile'])
       
       '''=========================================================='''
       ''' Add identical v3j clonotypes to clusters (undeduplicate) '''
       '''=========================================================='''
       appear=[]
       outfile = gzip.open(options['clustered-formatted-outfile'], 'wb')
       for key in lineage.keys():
           duplicate_entry=[]
           for entry in lineage[key]:
               temp=entry.split(' ')
               v3j='%s'%' '.join(map(str,[temp[num] for num in options['keys-for-clustering'][1:]]))
               del temp
               for mongoid in table[v3j]:
                   appear.append(mongoid)
                   duplicate_entry.append(table_mongoids[mongoid])
               #remove v3j from table
               #del my_dict['key']
           '''--------------------------------------'''
           ''' Deduplciate to avoid double counting '''
           '''--------------------------------------'''
           for unique_record in list(set(duplicate_entry)):
                 outfile.write('%s\n'%(unique_record))
           outfile.write('%s\n'%'END')
           del duplicate_entry
       outfile.close()
       appear=list(set(appear))
       print('====> Done Adding identical v3j clonotypes to  all clusters\n')

       '''====================================================================='''
       ''' Now determine all labels that don't appear in the clustered file and '''
       ''' and create a list of lists for all labels not appearing in clusters '''
       '''====================================================================='''
       list_of_lists=[]
       for key in list(set(table_mongoids.keys()).difference(set(appear))):
           temp=table_mongoids[key].split(' ')
           temp2=[]
           v3j='%s'%' '.join(map(str,[temp[num] for num in options['keys-for-clustering'][1:]]))
           if key in table[v3j]:
               for entry in table[v3j]:
                   #print table['%s %s %s'%(temp[3],temp[4],temp[6])]  #provides list of identical v3js
                   temp2.append(table_mongoids[entry])
               list_of_lists.append(temp2)
               del temp
               del temp2
        
       '''===================================='''
       ''' Remove redundant copies of entries '''
       '''===================================='''
       new_k = []
       for elem in list_of_lists:
           if elem not in new_k:
               new_k.append(elem)

       '''============================='''
       ''' Dump out singleton clusters '''
       '''============================='''
       outfile = gzip.open(options['clustered-formatted-outfile'], 'ab')
       k = new_k
       for elements in k:
           for element in elements:
               outfile.write('%s\n'%element)
           outfile.write('%s\n'%'END')
       outfile.close()
       print('====> Done Adding identical v3j clonotypes to  all singleton clusters\n')
       sys.exit(1)

       ''' Try and place these into a lists of lists based on their identity with V3J clonotype'''
       #print(table)
       #print('Total number of lineages: %d'%(len(otu_table)))
       #print('Total number of sequences: %d'%(len(table)))
       #print('Total number of sequences in clusters: %d'%(len(clusteredlabels)))
       #print('Total number of singletons: %d'%(len(singletons)))
    else:
       print('No clusters!') 
       pass
    if os.path.exists(cluster_file):
       pass
       #os.remove(cluster_file)

if __name__ == "__main__":
     main()

