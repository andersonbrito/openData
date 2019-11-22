#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python version: Python 3
#
#   tORFs.py -> This code retrieves peptides encoded by all open
#               reading frames (ORFs) longer than 50 amino acids found
#               in genomes deposited on NBCI. This script requires
#               a tab-delimited file containing 'NCBI accession number'
#               and 'taxa name' (eg. NC_001806 \t HHV1).
#
# Usage: tORFs.py auxFile
#
# Release date: 24/12/2017
# Last update: 16/11/2019
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from Bio import SeqIO
from Bio import Entrez
import sys

Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are

auxFile = open(sys.argv[1], "r").readlines()

for line in auxFile:
    accNo = line.split("\t")[0].strip()
    taxa = line.split("\t")[1].strip()
    print('Fetching peptide sequences:', accNo, '(' + taxa + ')')
    net_handle = Entrez.efetch(db="nucleotide", id=accNo, rettype="fasta", retmode="text")
    for seq_record in SeqIO.parse(net_handle, "fasta"):
        seq = seq_record.seq

    c = 0
    outFile = open(taxa + '-' + accNo + '-peptides.fasta', 'w')
    genomeSize = len(str(seq))
    for strand, nuc in [('forward', seq), ('reverse', seq.reverse_complement())]:  # forward and reverse DNA strands
        for frame in range(3):
            length = 3 * ((len(seq) - frame) // 3)  # Multiple of three
            transNuc = str(nuc[frame:frame + length].translate(table='Standard'))
            dicStartEnd = {}
            lstPep = [p + '*' for p in transNuc.split("*")]
            lstPep[-1] = lstPep[-1][:-1]

            posRead = 0
            if strand == 'forward':
                posRead = frame + 1
            if strand == 'reverse':
                posRead = genomeSize - frame

            for num, pep in enumerate(lstPep, 1):
                if strand == 'forward':
                    start = posRead
                    end = posRead - 1 + len(pep) * 3
                    startEnd = (start, end)
                    posRead = end + 1

                if strand == 'reverse':
                    start = posRead
                    end = posRead + 1 - len(pep) * 3
                    startEnd = (start, end)
                    posRead = end - 1

                protein = pep.replace('*', '')
                if len(protein) >= 50:  # iterates over all peptides longer than 50 aminoacids
                    c += 1
                    entry = '>' + accNo + '-' + 'p' + '0' * (4 - len(str(c))) + str(c) + '.' + strand + '.frame' + str(
                        frame + 1) + '.' + str(min(startEnd)) + '-' + str(max(startEnd)) + '\n' + str(protein) + "\n"
                    outFile.write(entry)
print('\nDone! Translated ORFs saved in your working directory.\n')
