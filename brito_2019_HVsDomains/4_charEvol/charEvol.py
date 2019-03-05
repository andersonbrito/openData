#!/usr/bin/python

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
# Python version: Python 3
#
#   charEvol.py -> This code processes matrices from Mesquite's ancestral
#                   state reconstruction analyses, and outputs files for
#                   visualizing the results on iTOL (itol.embl.de). It allows
#                   a closer look at the dynamics of domain gains, losses and
#                   duplications along the branches of a given tree.
#
# Usage: python charEvol.py workingDirectory beastTree originalMatrix mesquiteMatrix
#
# Release date: 29/11/2017
# Last update: 05/03/2019
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from Bio import Phylo
import networkx as nx
from collections import Counter
import random
from datetime import datetime
import os
from io import StringIO
from sys import *


dir = argv[1] # working directory
treeFile = argv[2] # phylogenetic tree
tree = open(dir + treeFile, 'r').readlines()
domMtxFile = open(dir + argv[3], "r").readlines() # matrix of domain counts
merMtxFile = open(dir + argv[4], "r").readlines() # 'Trace all characters' Mesquite output


# generate a unique folder name each time the code is run
day, time = str(datetime.today()).split()
uniqueDir = day + '_' + ''.join(time.split('.')[0].split(':'))

# convert Beast tree into Newick format
dicNames = {}
for line in tree:
    if line.startswith('\t\t  '):
        num, spp = line.strip().replace(',', '').split()
        dicNames[num] = spp

treedata = ''
for line in tree:
    sep = 'treetree1='
    if sep.lower() in line.lower().replace(' ',''): # to get the line containing the tree data
        line = line.replace(' ','').replace('Tree', 'tree').replace('TREE', 'tree')
        treedata = line.replace(' ','').split(sep)[1].strip()

handle = StringIO(treedata)
tree = Phylo.read(handle, "newick")

# to rename clade names
for clade in tree.find_clades():
    if str(clade.name) in dicNames.keys():
        clade.name = dicNames[clade.name]

    if str(clade.name) == 'None':
        listComm = str(clade.comment).split(",")
        for c in listComm:
            if 'posterior' in c:
                posterior = c.split("=")[-1]
                clade.confidence = float(posterior) * 100
        clade.comment = ""
    if str(clade.name) != 'None':
        clade.comment = ""

# name internal nodes
c = 2
for clade in tree.find_clades():
    if str(clade.name) == 'None':
        clade.comment = '&&NHX:name=n'+str(c)
        clade.name = 'n'+str(c)
        clade.confidence = None
    c += 1

# convert first tree into a network
vNet = Phylo.to_networkx(tree)

# create list of tree tips
lstTips = [term.name for term in tree.get_terminals()]

# change node names creating a new network
vGraph = nx.Graph()
for edge in vNet.edges():
    root = ('root', 'n2')
    if root not in vGraph:
        vGraph.add_edge('root', 'n2')
    vGraph.add_edge(edge[0].name, edge[1].name)

# check the edges
lstEdges = []
for edge in vGraph.edges():
    lstEdges.append(edge)


# print all node names and paths
evolPaths = []
for node in vGraph.nodes():
    # print(node)
    for path in nx.all_simple_paths(vGraph, source='root', target=node):
        if path[-1] in lstTips:
            evolPaths.append(','.join(path))

# export an iTol format tree with comments, and no 'node.name'
for clade in tree.find_clades():
    if str(clade.name) not in lstTips:
        clade.name = None
Phylo.write([tree], dir + treeFile.split(".")[0] + "_itol.tree", 'newick')


# get header of characters matrix
lstDomains = [dom.strip() for dom in domMtxFile[0].split('\t')[1:]]

# process and transpose original meristic matrix from Mesquite 'Trace all characters' analysis
merLst = []
start = ''
for line in merMtxFile:
    line = line.strip()
    if 'Char.\\Node' in line:
        start = 'Found'
        # print(line)
    if start == 'Found':
        if '(min.)' in line:
            line = line.replace(' (min.):', '', 10000).replace('(max.):', '', 10000).replace(';  ', '_', 10000).replace('; ', '', 10000)
            if line.endswith(';'):
                line = line[:-1]
            fixed = []
            for col in line.split('\t'):
                if 'character' in col:
                    domNum = int(col.split()[-1])-1
                    fixed.append(lstDomains[domNum])

                if len(col.split('_')) == 2:
                    minC, maxC = col.split('_')
                    if minC == maxC:
                        fixed.append(maxC)
                    else:
                        fixed.append(col)
            merLst.append(fixed)
        else:
            nNodeLst = []
            for nodeInfo in line.split('\t'):
                if nodeInfo not in lstTips:
                    nNodeLst.append('n' + nodeInfo) # add an 'n' to each internal node name
                else:
                    nNodeLst.append(nodeInfo)
            merLst.append(nNodeLst)


# transpose the fixed meristic matrix and output it
transpMerlist = []
# outTranspose = open(dir + 'traceAll_matrix_converted.txt', 'w')
for pos in list(map(list, zip(*merLst))):
    tLine = '\t'.join(pos)
    transpMerlist.append(tLine)
#     outTranspose.write('\t'.join(pos) + '\n')
# outTranspose.close()


# create a dictionary with parsimony estimates per node/domain
dicNodes = {}
for line in transpMerlist[1:]:
    dicNodes[line.split('\t')[0]] = [count.strip() for count in line.split('\t')[1:]]

for edge in lstEdges[1:]: # exclude ('root', 'n1') branch
    bfNode = edge[0]
    afNode = edge[1]

    if afNode not in lstTips:
        for num, profile in enumerate(zip(dicNodes[bfNode], dicNodes[afNode])):
            countBe, countAf = profile

            # choose randomly one of multiple parsimonious solutions
            if '_' in str(countBe): # pick a count for before node
                lstNum = list(range(int(countBe.split('_')[0]), int(countBe.split('_')[1]) + 1))
                dicNodes[bfNode][num] = random.choice(lstNum)

            if '_' in str(countAf):  # pick a count for after node
                lstNum = list(range(int(countAf.split('_')[0]), int(countAf.split('_')[1]) + 1))
                if dicNodes[bfNode][num] in lstNum:  # if ancestral node bfNode was uncertain, assign the same value here
                    dicNodes[afNode][num] = dicNodes[bfNode][num]
                else:  # if ancestral node bfNode was pre-defined, assign random choice from lstNum
                    dicNodes[afNode][num] = random.choice(lstNum)


# count the character changes (gain, loss, and duplication)
dicChanges = {}
for edge in lstEdges:
    bfNode = edge[0]
    afNode = edge[1]

    if bfNode == 'root':
        dicNodes['root'] = dicNodes['n2']

    dicEvents = {'gain': [], 'loss': [], 'dup': []}
    for pfam, countBe, countAf in zip(lstDomains, dicNodes[bfNode], dicNodes[afNode]):
        difference = int(countAf) - int(countBe)

        if int(countAf) > 1 and int(countBe) > 0 and difference > 0: # count multiple domain duplication events
            dicEvents['dup'].extend([pfam] * abs(difference))

        if int(countBe) == 0 and difference > 1:  # count domain gain followed by duplication events
            dicEvents['gain'].append(pfam)

            dicEvents['dup'].extend([pfam] * abs(difference - 1))

        if int(countBe) == 0 and difference == 1: # count single domain gain events
            dicEvents['gain'].append(pfam)

        if int(countBe) > 0 and difference == -1: # count single domain loss events
            dicEvents['loss'].append(pfam)

        if int(countBe) > 0 and difference < -1: # count multiple domain loss events
            dicEvents['loss'].extend([pfam] * abs(difference))

        if bfNode == 'root':
            if int(countAf) == 1:
                dicEvents['gain'].append(pfam)
            if int(countAf) > 1:
                dicEvents['gain'].append(pfam)
                dicEvents['dup'].append(pfam)
    dicChanges[edge] = dicEvents


# it prints a table with a list of event occurrences for the current run
if 'reconstructions' not in os.listdir(dir):
    os.system("mkdir %s%s" % (dir, 'reconstructions'))


subDir = 'reconstructions/' + uniqueDir + '/'
if uniqueDir not in os.listdir(dir + 'reconstructions/'):
    os.system("mkdir %s%s" % (dir, subDir))
    os.system("mkdir %s%s" % (subDir, 'raw_data'))
    os.system("mkdir %s%s" % (subDir, 'iTOL_annotations'))
    os.system("mv %s %s" % (dir + treeFile.split(".")[0] + "_itol.tree", subDir + 'iTOL_annotations'))


desFile = open(dir + subDir + 'raw_data/' + 'list_eventsPerBranch.txt', 'w')
desFile.write("branch\ttotal gains\tacquired domains\ttotal losses\tlost domains\ttotal duplications\tduplicated domains\n")
for edge, events in dicChanges.items():
    newDic = {}
    for type, doms in events.items():
        newDic[type] = []
        for pfam in doms:
            newDic[type] += [pfam]
    result = '\t'.join([(str(len(occ))) + '\t' + '; '.join(occ) for eve, occ in newDic.items()])
    desFile.write(edge[0] + ' → ' + edge[1] + '\t' + result + '\n')
print('\n')

# generate itol codes for labeling domains in branches
c = 1
newDic = {}
dicRes = {}
dicRes['domain'] = []
recEvents = {}
recEvents['gain'] = []
recEvents['loss'] = []
recEvents['dup'] = []
for edge, events in dicChanges.items():
    node = edge[1]
    if edge[1].isdigit():
        node = 'n' + edge[1]

    # define the spacing between elements mapped on a branch
    allDoms = [item for sublist in list(events.values()) for item in list(set(sublist))]
    step = ''
    if len(allDoms) > 1:
        step = 1 / (len(allDoms) + 1)
    elif len(allDoms) == 1:
        step = 0.5
    else:
        pass

    pos = step
    o = 1
    for type, doms in events.items():
        if node != 'n2':
            dicCount = dict(Counter(doms))
            duplicates = [item for item, count in dicCount.items() if count > 1]

            for d in list(set(doms)):
                # Using generic IDs as labels
                if d not in newDic.keys():
                    dId = 'd' + (3 - len(str(c))) * '0' + str(c)
                    newDic[d] = dId
                    dName = dId + '\t' + d
                    c += 1

                dLabel = newDic[d]
                if d in duplicates: # to show/identify only one domain representation when event happens multiple times
                    dLabel = dLabel + "*"

                if type == 'gain':
                    dicRes['domain'] += [','.join([node, '1', '1', '#54A800', '1', str(pos)[:4], d])]

                    if 'label' + str(o) not in dicRes.keys():
                        dicRes['label' + str(o)] = []
                    dicRes['label' + str(o)] += [','.join([node, dLabel, str(pos)[:4], '#54A800', 'normal', '1'])]

                    recEvents['gain'] += [dName]  # to add to the recurrent events list

                if type == 'loss':
                    dicRes['domain'] += [','.join([node, '1', '1', '#000000', '0', str(pos)[:4], d])]

                    if 'label' + str(o) not in dicRes.keys():
                        dicRes['label' + str(o)] = []
                    dicRes['label' + str(o)] += [','.join([node, dLabel, str(pos)[:4], '#000000', 'normal', '1'])]

                    recEvents['loss'] += [dName]  # to add to the recurrent events list

                if type == 'dup':
                    dicRes['domain'] += [','.join([node, '1', '1', '#7849A8', '1', str(pos)[:4], d])]

                    if 'label' + str(o) not in dicRes.keys():
                        dicRes['label' + str(o)] = []
                    dicRes['label' + str(o)] += [','.join([node, dLabel, str(pos)[:4], '#7849A8', 'normal', '1'])]

                    recEvents['dup'] += [dName]  # to add to the recurrent events list

                pos += step
                o += 1
        else:  # to add domains present at the root → n2 branch
            for d in doms:
                dName = 'dxx' + '\t' + d
                if type == 'gain':
                    recEvents['gain'] += [dName]  # to add to the recurrent events list
                if type == 'loss':
                    recEvents['loss'] += [dName]  # to add to the recurrent events list
                if type == 'dup':
                    recEvents['dup'] += [dName]  # to add to the recurrent events list


# outputs a list of recurrent events
recFile = open(dir + 'reconstructions/' + uniqueDir  + '/raw_data/' + 'reconstructedEvents.txt', 'w')
newHeaders = {'gain':'DOMAIN GAINS', 'loss':'DOMAIN LOSSES', 'dup':'DOMAIN DUPLICATIONS'}
for ev, domNames in recEvents.items():
    recFile.write('\n##### ' + newHeaders[ev] + ' #####\n')
    for d in domNames:
        recFile.write(d + '\n')


# outputs the legend for the current run
legFile = open(dir + 'reconstructions/' + uniqueDir  + '/raw_data/' + 'legend_domains.txt', 'w')
for pfam, id in newDic.items():
    legFile.write(id + ': ' + pfam + '\n')


# output annotation files for iTOL
for out, info in dicRes.items():
    outputFile = open(dir + 'reconstructions/' + uniqueDir  + '/iTOL_annotations/' +  out + '.txt', 'w')
    if 'domain' in out:
        outputFile.write('DATASET_SYMBOL\n\nSEPARATOR COMMA\n\nDATASET_LABEL,domain symbols\n\nCOLOR,#ffff00\n\nMAXIMUM_SIZE,3\n\nDATA\n\n')
    else:
        outputFile.write('DATASET_TEXT\n\nSEPARATOR COMMA\n\nDATASET_LABEL,'+out+'\n\nCOLOR,#ffff00\n\nLABEL_ROTATION,-45\n\nSIZE_FACTOR,1\n\nDATA\n\n')

    for line in info:
        outputFile.write(line + '\n')

print('Results saved at ' + dir + 'reconstructions/' + uniqueDir + '\n')
