
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 00:17:13 2022

@author: jgrandinetti
"""
import os
from Bio import SeqIO
from Bio.Seq import Seq
from glob import glob

maindir = '/home/jgrandi7/Desktop/MiniProject/' #CHANGE THIS TO YOUR DIRECTORY

sratoolkit_path = maindir + 'sratoolkit.2.11.2-ubuntu64/bin/' 
sratoolkit_userrepo_path = maindir + 'my_repository'
SPAdes_path = maindir + 'SPAdes-3.15.4-Linux/bin/spades.py'
prokka_path = maindir + 'prokka/bin/'
#WARNING - these paths may change depending on which system and version of the dependencies you are running
#modify the above paths as needed

results = maindir + 'results' #do not touch this
os.system('mkdir ' + results)

def mp_retrieve():
    sratoolkit_command1 = sratoolkit_path + 'prefetch SRR8185310' + ' -O ' + results
    sratoolkit_command2 = sratoolkit_path + 'fastq-dump SRR8185310' + ' -O ' + results
    
    os.system(sratoolkit_command1)
    os.system(sratoolkit_command2)
    

def spadesexe(output):

    spadescommand = SPAdes_path + ' -k 55,77,99,127 -s ' + results + '/SRR8185310.fastq -o ' + results
    output.write('\n'+spadescommand)
    os.system(spadescommand)
    
def numcontigs(output):
         
    handle = open(results + '/contigs.fasta')  # open problem file
    listofseq = []  # new empty list to store sequences
    seqdict = {}
    greaterthan1000 = {}
    numbp = 0
    
    def get_key(val): #function to get a key given a dictionary value
       for key, value in seqdict.items():
          if val == value:
             return key
         
    for record in SeqIO.parse(handle, 'fasta'):
        listofseq.append(str(record.seq))  # append all sequences as string to list
        seqdict[record.id] = str(record.seq)
        
            
    for seq in seqdict.values():
        if len(seq) > 1000:
            greaterthan1000[get_key(seq)] = seq
            
    for seq in greaterthan1000.values():
        numbp += len(seq)
            
    output.write('\nThere are ' + str(len(greaterthan1000)) + ' contigs > 1000 in the assembly.')
    
    output.write('\nThere are ' + str(numbp) + ' bp in the assembly.')
    
    with open(results+'/contigs1000.fasta', 'w') as outfile:
        for seq in greaterthan1000.values():
            outfile.write('>'+get_key(seq)+'\n')
            outfile.write(seq+'\n')
            

def prokka(output):
    prokkacommand = prokka_path +'prokka --outdir ' + results + '/prokka --genus Escherichia --locustag ECOL ' + results +'/contigs1000.fasta'
    output.write('\n' + prokkacommand +'\n')
    os.system(prokkacommand)
    filename = glob(results + ('/prokka/PROKKA_********.txt'))[0]
    prokkadata = []
    with open(filename,'r') as anno:
        for line in anno:
            output.write(line)
            prokkadata.append(line)
    return prokkadata
    
def prokkacomp(output, prokkadata):
    CDS = int(prokkadata[3].strip('CDS: '))
    tRNA = int(prokkadata[5].strip('tRNA: '))
    statement = 'Prokka found ' + str(CDS - 4140) + ' additional CDS and ' + str(89 - tRNA) + ' less tRNA than the RefSeq.'
    outfile.write(statement)
    
with open(results + '/miniproject.log', 'w') as outfile:
    mp_retrieve()
    spadesexe(outfile)
    numcontigs(outfile)
    prokkadata = prokka(outfile)
    prokkacomp(outfile,prokkadata)
    print('Success! Thanks for using COMP-383-Mini-Project by jgrandi7')
    outfile.close()
    
    
