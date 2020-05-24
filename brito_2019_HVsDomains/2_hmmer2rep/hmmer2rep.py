#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python version: Python 3
#
#     hmmer2rep.py -> This code searches for domain hits on HMMER (hmmscan)
#                     outputs (domtblout), and convert it into a custom format,
#					  while retrieving structural information about such
#                     domains on Pfam.
#
# Usage: python hmmer2rep.py 'fullDirectoryPath' 'hmmscanOutput' 'auxFile'
#
# Release date: 23/12/2017
# Last update: 26/02/2019
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os
import re
from sys import *
try:
    import urllib.request as urllib2
except ImportError:
    import urllib2
from bs4 import BeautifulSoup as BS
import argparse

dThreshold = 1e-03
pThreshold = 1e-03

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="This script parses hmmscan outputs and generate as output a list of domains detected per taxa.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-input", required=True, help="hmmscan output file (domtblout.txt)")
    parser.add_argument("-list", required=True, help="TSV file containing accession number and taxa names")
    parser.add_argument("-output", required=True, help="TSV files with all domains detected per taxa")
    args = parser.parse_args()

    input = args.input
    list = args.list
    output = args.output


    domFile = open(input).readlines()
    auxFile = open(list, "r").readlines()

    # create the output directory 'repertoire'
    repDir = 'repertoire/'
    if 'repertoire' not in os.listdir():
        os.system("mkdir %s" % repDir)

    # dictionary of accession numbers and species acronyms
    dicNameAcc = {}
    for line in auxFile:
        dicNameAcc[line.split('\t')[0]] = line.split('\t')[1].strip()

    domOcc = {}
    # searche for PDB entries of domains on Pfam website, and returns the number of entries.
    def searchPDB(dAccNo):
        domPage = urllib2.urlopen("http://pfam.xfam.org/family/" + str(dAccNo))
        domPageData = domPage.read()
        soup = BS(domPageData, "html.parser")
        pdbEntries = str(soup.find('div', id="icons").find('span', id="structIcon").find('em').string).strip()
        if pdbEntries == '0':
            pdbEntries = '-'
        return pdbEntries

    # generate a dict of species and their proteins
    for genome in auxFile:
        taxa = genome.split('\t')[0]
        try:
            outputFile = open(repDir + dicNameAcc[taxa] + '-domRepertoire.txt', 'w')
            header = '# Results for ' + taxa + ' (' + dicNameAcc[taxa] + ')\n'
            outputFile.write(header)
            print(header.strip())

            count = 0
            # parse Hmmer Search results
            for line in domFile:
                if not line.startswith('#'):
                    spp = line.split()[3].split('-')[0]
                    prot = line.split()[3]
                    dcode = line.split()[1].split('.')[0]
                    cevalue = line.split()[11]
                    pevalue = line.split()[6]

                    # check if domain and sequence evalues are below the inclusion threshold
                    if spp == taxa and float(cevalue) <= dThreshold and float(pevalue) <= pThreshold:
                        if dcode in domOcc.keys():
                            entry = '\t'.join([prot, dcode, cevalue+'_'+pevalue, domOcc[dcode]])
                            count += 1
                        else:
                            gotPDB = searchPDB(dcode)
                            entry = '\t'.join([prot, dcode, cevalue+'_'+pevalue, gotPDB])
                            domOcc[dcode] = gotPDB
                            count += 1
                        outputFile.write(entry + '\n')
                        print(entry)
            print('\n')
            outputFile.write("END\n")
            outputFile.close()
        except:
            pass
