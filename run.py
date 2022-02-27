#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 00:17:13 2022

@author: jgrandinetti
"""
import os
from Bio import SeqIO
from Bio.Seq import Seq
from python_on_whales import docker

maindir = '/Users/jgrandinetti/Desktop/MiniProject/' #CHANGE THIS TO YOUR DIRECTORY

sratoolkit_path = maindir + 'sratoolkit.3.0.0-mac64/bin/'
sratoolkit_userrepo_path = maindir + 'user-repository'
SPAdes_path = maindir + 'SPAdes-3.15.4-Darwin/bin/spades.py'
prokka_path = maindir + 'prokka/bin/'

alternate_prokka = docker.run('staphb/prokka:latest')

results = maindir + 'results'
#createresults = 'mkdir ' + results
#os.system(createresults)



def mp_retrieve():
    sratoolkit_command1 = sratoolkit_path + 'prefetch SRR8185310' + ' -O ' + results
    sratoolkit_command2 = sratoolkit_path + 'fastq-dump SRR8185310' + ' -O ' + results
    
    os.system(sratoolkit_command1)
    os.system(sratoolkit_command2)
    

def spadesexe(output):

    spadescommand = SPAdes_path + ' -k 55,77,99,127 -s ' + results + '/SRR8185310.fastq -o ' + results
    output.write(spadescommand)
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
            
    output.write('There are ' + str(len(greaterthan1000)) + ' contigs > 1000 in the assembly.')
    
    output.write('\nThere are ' + str(numbp) + ' bp in the assembly.')
    
    with open(results+'/contigs1000.fasta', 'w') as outfile:
        for seq in greaterthan1000.values():
            outfile.write('>'+get_key(seq)+'\n')
            outfile.write(seq+'\n')
            

def prokka(output):
    prokkainit = docker.pull('staphb/prokka:latest')
    prokkacommand = docker.run('staphb/prokka:latest prokka --outdir ' + results + ' --genus Escherichia --locustag ECOL ' + results +'/contigs1000.fasta')
    output.write('\n' + prokkacommand)
    os.system(prokkainit)
    os.system(prokkacommand)
    
    
    
with open('miniproject.log', 'w') as outfile:
    #mp_retrieve()
    #spadesexe(outfile)
    numcontigs(outfile)
    prokka(outfile)
    outfile.close()
    
    
