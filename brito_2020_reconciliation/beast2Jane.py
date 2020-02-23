#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python version: Python 3
#
#   beast2Jane.py -> This script converts beast trees into Jane
#                    cophylogenetic format. A tab delimited
#                    auxiliary file containing virus-host pairs
#                    (one pair pair line) is required. If any pair
#                    is missing from such list, it will be pruned from
#                    the final converted tree. The time zones are
#                    predefined in 5 Millions of years. Time zones
#                    with no nodes assigned will be filled by taxa in an
#                    artifical ladderized outgroup added for this purpose.
#
# Usage: python beast2Jane.py workingDirectory hostTree viralTree auxFile
#
# Release date: 22/08/2017
# Last update: 25/02/2019
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from sys import *
from Bio import Phylo
from io import StringIO

dir = argv[1]
hTree = open(dir + argv[2], 'r').readlines()
vTree = open(dir + argv[3], 'r').readlines()
auxFile = open(dir + argv[4], "r").readlines()
windowLength = int(argv[5])

# create a dictionary of viruses and hosts
vhPairs = {}
for pair in auxFile:
    vhPairs[pair.split("\t")[0]] = pair.split("\t")[1].strip()

allTaxa = list(vhPairs.keys()) + list(vhPairs.values())


def beast2nwk(treeFile, output):
    # to get numeric names, and create a dict with their clade names
    dicNames = {}
    for line in treeFile:
        if line.startswith('\t\t  '):
            num, spp = line.strip().replace(',', '').split()
            dicNames[num] = spp

    # to get the line containing the tree data
    treedata = ''
    for line in treeFile:
        if 'tree TREE1'.lower() in line.lower():
            line = line.replace("=", "£", 1)
            treedata = line.split("£")[1].strip()


    handle = StringIO(treedata)
    tree = Phylo.read(handle, "newick")

    # to rename clade names
    for clade in tree.find_clades():
        if str(clade.name) in dicNames.keys():
            clade.name = dicNames[clade.name]


    # get terminals and create a list
    listTerminal = []
    c = 0
    for clade in tree.find_clades():
        if str(clade.name) == 'None' and c == 0:
            for term in clade.get_terminals():
                listTerminal.append(term.name)
        c += 1


    # to remove taxa from original tree
    for tax in listTerminal:
        if tax not in allTaxa:
            tree.prune(target=tax.strip())

    # to output a pruned tree and fix potential problems of formatting
    Phylo.write([tree], dir + output + "_pruned.tree", 'nexus')
    treeLines = []
    for line in open(dir + output + "_pruned.tree").readlines():
        treeLines.append(line)

    outFile = open(dir + output + "_pruned.tree", "w")
    for line in treeLines:
        if line.startswith(" Tree tree1"):
            line = line.replace("):0.00000[", ")[")
            line = line.replace(";", ":0.00000;")
            outFile.write(line)
        else:
            outFile.write(line)
    outFile.close()

    return tree

# get viral and host tree lengths
allZones = []
def getTreeDistance(tree):
    treeDist = 0
    for clade in tree.find_clades():
        listComm = str(clade.comment).split(",")
        for num, com in enumerate(listComm):
            if 'height_95%_HPD' in com:
                smallNum = float(com[16:])
                bigNum = float(listComm[num + 1][:-1])
                if bigNum > treeDist:
                    treeDist = bigNum
                    break
    allZones.append(treeDist)
    return treeDist

getTreeDistance(beast2nwk(hTree, 'hTree'))
getTreeDistance(beast2nwk(vTree, 'vTree'))


# define the REVERSE time zones (Myr) boundaries based on the longest tree, and the given window length
limits = []
def zoneRanges(tree):
    limit = 0
    treeDist = float(getTreeDistance(tree))
    if treeDist >= max(allZones):
        while limit > treeDist*(-1):
            limits.append(limit)
            limit -= windowLength
zoneRanges(beast2nwk(hTree, 'hTree'))
zoneRanges(beast2nwk(vTree, 'vTree'))

limits = sorted(limits)
olderZone = int((limits[0]*(-1)) / windowLength)

#print(olderZone)

zCoverage = []
# set common numeric zones for both viral and host trees, where the last zone correspond to their tips
def setZones(tree):
    dicZones = {}
    for clade in tree.find_clades():
        terminals = "-".join(sorted([term.name for term in clade.get_terminals()]))

        listComm = str(clade.comment).split(",")

        # correct length difference between trees based on their delay
        zones = []
        for num, com in enumerate(listComm):
            if 'height_95%_HPD' in com:

                oldUB = float(listComm[num + 1][:-1])
                oldLB = float(com[16:])

                upperBoundary = (oldUB - float(oldUB-oldLB) * 1/4)*-1
                lowerBoundary = (oldLB + float(oldUB-oldLB) * 1/4)*-1

                for num, l in enumerate(limits):
                    try:
                        if l <= lowerBoundary < limits[num + 1]:
                            zones.append(num+1)
                        if l <= upperBoundary < limits[num + 1]:
                            zones.append(num + 1)
                        if upperBoundary < limits[0]: # in case one of the trees' root is older than the oldest limit
                            if 1 not in zones:
                                zones.append(1)
                    except:
                        pass

                if clade.is_terminal():
                    zones = [max(zones)]
                if min(zones) == max(zones):
                    zones = [max(zones)]

        # find gaps in host time zone
        for z in list(range(min(zones), max((zones)))):
            # z = int(z)
            if z not in zCoverage:
                zCoverage.append(z)

        zones = [str(z) for z in zones]
        dicZones[terminals] = ','.join(zones)

    # generate newick format jane trees
    for clade in tree.find_clades():
        terminals = "-".join(sorted([term.name for term in clade.get_terminals()]))
        for i in dicZones.keys():
            if terminals == i:
                clade.branch_length = None
                clade.confidence = None
                clade.comment = None
                clade.comment = dicZones[terminals]
            if len(terminals.split('-')) > 1:
                clade.name = None

    treeString = tree.format("newick").replace('0.00000', '')#.replace(':',':['))
    return(treeString)

setZones(beast2nwk(hTree, 'hTree')) # this must be active to run the lines below


gapZones = [zone for zone in list(range(1, olderZone)) if zone not in sorted(zCoverage)]
if 1 not in gapZones:
    gapZones.insert(0, 1)


# find gaps in host time zones to assign artificial ancestral outgroup to fill the gap
gNodes = []
solved = []
for n, z in enumerate(gapZones):
    if z == gapZones[-1]:
        if str(z) not in solved:
            gNodes.append((str(z),))
            solved.append(str(z))
        else:
            continue
    elif z == gapZones[0]:
        gNodes.append((str(z),))
        solved.append(str(z))
        continue
    else:
        if gapZones[n+1] - gapZones[n] > 1:
            if str(z) not in solved: # artificial nodes will only span two zones
                gNodes.append((str(z),))
                solved.append(str(z))
                continue
        else:
            zPair = (str(gapZones[n]), str(gapZones[n + 1]))
            if zPair not in gNodes and zPair[0] not in solved:
                gNodes.append(zPair)
                solved.extend(zPair)

# create list of artifical nodes to be grafted into the original tree
gNodes = gNodes[::-1]
vGraft = ',' + '(' * (len(gNodes)-1) + 'virus' + str(1) + ':[' + str(olderZone) + '],'
hGraft = ',' + '(' * (len(gNodes)-1)+ 'host' + str(1) + ':[' + str(olderZone) + '],'
vhPairs['virus' + str(1)] = 'host' + str(1)

for number, node in enumerate(gNodes):
    if node == gNodes[-1]:
        vGraft += '):[' + ','.join(node) + '];'
        hGraft += '):[' + ','.join(node) + '];'
    elif node == gNodes[-2]:
        vGraft += 'virus' + str(number + 2) + ':[' + str(olderZone) + ']):[' + ','.join(node) + ']'
        hGraft += 'host' + str(number + 2) + ':[' + str(olderZone) + ']):[' + ','.join(node) + ']'
        vhPairs['virus' + str(number + 2)] = 'host' + str(number + 2)
    else:
        vGraft += 'virus' + str(number + 2) + ':[' + str(olderZone) + ']):[' + ','.join(node) + '],'
        hGraft += 'host' + str(number + 2) + ':[' + str(olderZone) + ']):[' + ','.join(node) + '],'
        vhPairs['virus' + str(number + 2)] = 'host' + str(number + 2)


# output Jane format timed trees
janeOutput = open(dir + dir.split('/')[-2] + '_GraftTimed.nex', "w")
janeOutput.write('#NEXUS\nbegin host;\ntree host=(' + setZones(beast2nwk(hTree, 'hTree'))[:-2] + hGraft + '\nendblock;\n')
janeOutput.write('\nbegin parasite;\ntree parasite=(' + setZones(beast2nwk(vTree, 'vTree'))[:-2] + vGraft + '\nendblock;\n')

# output the list of virus-host pairs
janeOutput.write('\nbegin distribution;\nRange\n')
for virus, host in vhPairs.items():
    if host == list(vhPairs.values())[-1]:
        janeOutput.write('\t' + virus + ' : ' + host + ';\nendblock;')
    else:
        janeOutput.write('\t' + virus + ' : ' + host + ',\n')
janeOutput.close()


