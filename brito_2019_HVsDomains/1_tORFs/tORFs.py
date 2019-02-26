#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python version: Python 3
#
#   tORFs.py -> This code retrieves peptides encoded by all open
#               reading frames (ORFs) longer than 40 amino acids found
#               in genomes deposited on the NBCI repository.
#
# Usage: python tORFs.py workingDirectory auxFile
#
# Release date: 24/12/2017
# Last update: 26/02/2019
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os
from sys import *
from Bio import SeqIO
from Bio import Entrez
try:
    import urllib.request as urllib2
except ImportError:
    import urllib2

Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are

dir = argv[1]
auxFile = open(dir + argv[2], "r").readlines()

repDir = dir + 'seqDB/'
if 'seqDB' not in os.listdir(dir):
    os.system("mkdir %s" % repDir)

for line in auxFile:
    accNo = line.split("\t")[0].strip()
    taxa = line.split("\t")[1].strip()
    print('Fetching peptide sequences:', accNo, '(' + taxa + ')')
    net_handle = Entrez.efetch(db="nucleotide", id=accNo, rettype="fasta", retmode="text")
    for seq_record in SeqIO.parse(net_handle, "fasta"):
        seq = seq_record.seq

    c = 0
    outFile = open(dir + 'seqDB/' + taxa + '-' + accNo + '_peptides.fasta', 'w')
    for strand, nuc in [('forward', seq), ('reverse', seq.reverse_complement())]: # forward and reverse DNA strands
        for frame in range(3):
            length = 3 * ((len(seq) - frame) // 3)  # Multiple of three
            for pro in nuc[frame:frame + length].translate(table='Standard').split("*"): # list of all translated ORFs
                if len(pro) >= 40: # iterates over all peptides longer than 40 aminoacids
                    c += 1
                    entry = '>' + accNo + '-' + 'p' + '0' * (4 - len(str(c))) + str(c) + '.' + strand + '.frame' + str(frame+1) + '\n' + str(pro) + "\n"
                    outFile.write(entry)
    outFile.close()



