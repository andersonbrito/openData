#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python version: Python 3
#
#   domAnnot.py -> Given a set of peptides annotated with CDS positions,
#                  this code creates GFF annotations for Pfam domains
#                  detected in HMMER searchers. As input it requires a
#                  fasta file with the annotated peptides, and the output
#                  from hmmscan (domtblout.txt)
#
#
#   Usage: domAnnot.py fastaFile domtblout.txt
#
# Release date: 16/11/2019
# Last update: 16/11/2019
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from Bio import SeqIO
import sys

peptides = sys.argv[1]
domains = open(sys.argv[2]).readlines()

dThreshold = 1e-03
pThreshold = 1e-03

dPos = {}
fasta_sequences = SeqIO.parse(open(peptides),'fasta')
for fasta in fasta_sequences: # iterate over all fasta sequences
    id, seq = fasta.description, fasta.seq
    orf = id.split(".")[0]
    startOrf = id.strip().split(".")[-1].split("-")[0]
    endOrf = id.strip().split(".")[-1].split("-")[1]
    sense = id.strip().split(".")[1]
    dPos[orf] = startOrf + "_" + endOrf + "_" + sense

outFile = open("-".join(peptides.split("-")[:-1]) + '-domAnnot.gff', 'w')
header = "##gff-version 3\n# seqid	source	type	start	end	score	strand	frame	attributes\n"
print(header)
outFile.write(header)

genome = ''
for line in domains:
    if not line.startswith('#'):
        cevalue = line.split()[11]
        pevalue = line.split()[6]
        orf = line.split()[3].split(".")[0]
        if genome == '':
            genome = orf.split("-")[0]
        if orf in dPos.keys():
            if float(cevalue) <= dThreshold and float(pevalue) <= pThreshold:
                dcode = line.split()[1].split(".")[0]  # domain code
                dname = line.split()[0] # domain name
                sDom = int(line.split()[17]) # domain start position
                eDom = int(line.split()[18]) # domain end position
                sense = dPos[orf].split("_")[2]
                sorf = int(dPos[orf].split("_")[0]) # orf start position
                eorf = int(dPos[orf].split("_")[1]) # orf end position

                if sense == "forward":
                    sense = "+"
                    sanot = sorf + ((sDom * 3) - 2) - 1
                    eanot = sorf + (eDom * 3) - 1
                else:
                    sense = "-"
                    sanot = eorf - (eDom * 3) + 1
                    eanot = eorf - (sDom * 3) + 3

                entry = "\t".join([genome, "-", "misc_feature", str(sanot), str(eanot), ".", sense, ".", "Name="+dname+" ["+cevalue+"]; Note="+dcode])
                outFile.write(entry + "\n")
                print(entry)
