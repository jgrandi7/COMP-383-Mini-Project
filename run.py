
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 00:17:13 2022

@author: jgrandinetti
"""

#This script was successfully tested on Ubuntu 20.04

import os
from Bio import SeqIO
from Bio.Seq import Seq
from glob import glob

maindir = '/home/jgrandi7/Desktop/MiniProject/' #CHANGE THIS TO YOUR DIRECTORY

sratoolkit_path = maindir + 'sratoolkit.2.11.2-ubuntu64/bin/' 
SPAdes_path = maindir + 'SPAdes-3.15.4-Linux/bin/spades.py'
prokka_path = maindir + 'prokka/bin/'
#WARNING - these paths may change depending on which system and version of the dependencies you are running
#modify the above paths as needed

results = maindir + 'results' #do not touch this line - sets results folder that will be created in your installation directory
os.system('mkdir ' + results)

def mp_retrieve():
    sratoolkit_command1 = sratoolkit_path + 'prefetch SRR8185310' + ' -O ' + results #prefetch command 
    sratoolkit_command2 = sratoolkit_path + 'fasterq-dump ' + results +'/SRR8185310/SRR8185310.sra' + ' -O ' + results #fasterq-dump command
    
    os.system(sratoolkit_command1)
    os.system(sratoolkit_command2)
    

def spadesexe(output):

    spadescommand = SPAdes_path + ' -k 55,77,99,127 -s ' + results + '/SRR8185310.fastq --only-assembler -o ' + results #spades command
    output.write('\n'+spadescommand) #spades command written to log file
    os.system(spadescommand)
    
def numcontigs(output):
         
    handle = open(results + '/contigs.fasta')  # open fasta file produced by spades
    listofseq = []  # new empty list to store sequences
    seqdict = {}
    greaterthan1000 = {}
    numbp = 0 #counter to total bp
    
    def get_key(val): #function to get a key given a dictionary value
       for key, value in seqdict.items():
          if val == value:
             return key
         
    for record in SeqIO.parse(handle, 'fasta'):
        listofseq.append(str(record.seq))  # append all sequences as string to list
        seqdict[record.id] = str(record.seq)
        
            
    for seq in seqdict.values(): #for sequences in dictionary
        if len(seq) > 1000: #if sequence is greater than 1000 bp
            greaterthan1000[get_key(seq)] = seq #set key and value in empty dictionary
            
    for seq in greaterthan1000.values(): #for every contig/sequence that is greater than 1000 bp
        numbp += len(seq) #add bp of that sequence to numbp
            
    output.write('\nThere are ' + str(len(greaterthan1000)) + ' contigs > 1000 in the assembly.')
    
    output.write('\nThere are ' + str(numbp) + ' bp in the assembly.')
    
    with open(results+'/contigs1000.fasta', 'w') as outfile:
        for seq in greaterthan1000.values(): #rewriting new fasta with contigs greater than 1000 bp
            outfile.write('>'+get_key(seq)+'\n')
            outfile.write(seq+'\n')
            

def prokka(output):
    prokkacommand = prokka_path +'prokka --outdir ' + results + '/prokka --genus Escherichia --locustag ECOL ' + results +'/contigs1000.fasta'
    #prokka tags genus Escherichia because we are looking for E.Coli, will output results in prokka folder in results folder, will take data from fasta we made in last step
    output.write('\n' + prokkacommand +'\n')
    os.system(prokkacommand)
    filename = glob(results + ('/prokka/PROKKA_********.txt'))[0]
    #we want to look at the prokka output txt file, but the name is based on the date that prokka is run and will therefore be different every day you run it
    #to get around this naming convention, I am using the glob library to look for the date as wildcard characters
    prokkadata = []
    with open(filename,'r') as anno:
        for line in anno:
            output.write(line)
            prokkadata.append(line)
    return prokkadata #returning prokka text output as a list
    
def prokkacomp(output, prokkadata): #we are inputting prokka results as a list, prokkadata, produced in previous function
    CDS = int(prokkadata[3].strip('CDS: ')) #getting CDS value from prokkadata list
    tRNA = int(prokkadata[5].strip('tRNA: ')) #getting tRNA value from prokkadata list
    statement = 'Prokka found ' + str(CDS - 4140) + ' additional CDS and ' + str(89 - tRNA) + ' less tRNA than the RefSeq.'
    #Prokka finds more CDS than in RefSeq assembly so subtract RefSeq from Prokka
    #Prokka finds less tRNA than in RefSeq assembly so subtract Prokka from RefSeq
    outfile.write(statement)
    
with open(results + '/miniproject.log', 'w') as outfile:
    mp_retrieve()
    spadesexe(outfile)
    numcontigs(outfile)
    prokkadata = prokka(outfile)
    prokkacomp(outfile,prokkadata)
    print('Success! Thanks for using COMP-383-Mini-Project by jgrandi7')
    outfile.close()
    
    
