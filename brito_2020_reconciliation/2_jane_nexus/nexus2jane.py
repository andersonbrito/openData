#!/usr/bin/python

# Created by: Anderson Brito
# Email: andersonfbrito@gmail.com
#
#
# Release date: 2020-07-19
# Last update: 2020-07-19

import argparse
from sys import *
from Bio import Phylo
from io import StringIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Conversion of virus and host time-calibrated trees into Jane nexus format",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--virus-tree", required=True, help="Virus phylogeny in nexus format")
    parser.add_argument("--host-tree", required=True, help="Host phylogeny in nexus format")
    parser.add_argument("--pairs", required=True, help="Virus-host pairs in TSV format")
    parser.add_argument("--window", required=True, type=int,  help="Size of time zones in unit of time")
    args = parser.parse_args()

    virus_tree = args.virus_tree
    host_tree = args.host_tree
    pairs = args.pairs
    windowLength = args.window

    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/past&future/PhD/works/reconcile/jane/spp-spp/reconcile_virevol1/run10_alpha_blunt/"
    # virus_tree = path + "virus_beast_mcc.tree"
    # host_tree = path + "host_baltic_ren.tree"
    # pairs = path + "vhPairs_alpha.tsv"
    # windowLength = 5

    # beast trees
    vTree = open(virus_tree, 'r').readlines()
    hTree = open(host_tree, 'r').readlines()
    auxFile = open(pairs, "r").readlines() # taxa to remove


    # create a dictionary of viruses and hosts
    vhPairs = {}
    for pair in auxFile:
        vhPairs[pair.split("\t")[0]] = pair.split("\t")[1].strip()

    viral_taxa = [taxon for taxon in vhPairs.keys()]
    host_taxa = [taxon for taxon in vhPairs.values()]
    allTaxa = viral_taxa + host_taxa
    
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


    # set common numeric zones for both viral and host trees, where the last zone correspond to their tips
    zCoverage = []
    def setZones(tree, output):
        # print(tree)
        dicZones = {}
        for clade in tree.find_clades():
            terminals = "-".join(sorted([term.name for term in clade.get_terminals()]))
            print('* Internal node leading to: ' + terminals)

            listComm = str(clade.comment).split(",")
            # correct length difference between trees based on their delay
            zones = []
            for num, com in enumerate(listComm):
                if 'height_95%_HPD' in com:
    
                    oldUB = float(listComm[num + 1][:-1])
                    oldLB = float(com[16:])
                    # print(oldUB, oldLB)

                    print('height_95%_HPD = ' + str(oldUB) + ' / ' + str(oldLB) + '\n')

                    leftBoundary = (oldUB)*-1
                    rightBoundary = (oldLB)*-1

                    # print(limits, len(limits))
                    for num, l in enumerate(limits):
                        # print(num, l)
                        try:
                            left = l
                            right = limits[num + 1]

                            if left <= rightBoundary < right and left <= leftBoundary < right:
                                zones.append(num + 1)
                                print('    >|< ' + str(leftBoundary) + ' is in zone ' + str(num + 1))

                            if left <= leftBoundary < right:
                                # print('    >>> ' + str(leftBoundary) + ' is zone ' + str(num + 1))
                                if leftBoundary < left + windowLength * 0.5:
                                    # print(leftBoundary, '<', left + windowLength * 0.5)

                                    # print(leftBoundary, left + windowLength * 0.5)
                                    # print('zone', num + 1, '- ', left, 'a', right)
                                    print('    <<< ' + str(leftBoundary) + ' is in zone ' + str(num+1))
                                    zones.append(num + 1)
                                else:
                                    if len(zones) == 0:
                                        if num + 2 > len(limits)-1:
                                            zones.append(num + 1)
                                            print('    ||< ' + str(leftBoundary) + ' is in zone ' + str(num + 1))
                                        else:
                                            zones.append(num + 2)
                                            # print(zones)
                                            print('    ||< ' + str(leftBoundary) + ' is in zone ' + str(num + 2))

                            if left <= rightBoundary < right:
                                # print('' + str(rightBoundary) + ' is zone ' + str(num + 1))
                                if rightBoundary > right - windowLength * 0.5:
                                    # print(rightBoundary, '>', right - windowLength * 0.5)

                                    # print(rightBoundary, right - windowLength * 0.5)
                                    # print('zone', num + 1, '- ', left, 'a', right)
                                    print('    >>> ' + str(rightBoundary) + ' is in zone ' + str(num+1))
                                    zones.append(num+1)
                                else:
                                    if len(zones) == 1:
                                        zones.append(num)
                                        # print(zones)
                                        print('    >|| ' + str(rightBoundary) + ' is in zone ' + str(num))

                            if leftBoundary < limits[0]: # in case one of the trees' root is older than the oldest limit
                                if 1 not in zones:
                                    zones.append(1)
                                    print('zone', num + 1, '-', left, 'a', right)
                                    print('   ||| ' + str(leftBoundary) + ' is in zone ' + str(num + 1))

                            # if clade.is_terminal():
                            #     zones = [max(zones)]


                        except:
                            pass
                    # print('{', leftBoundary, rightBoundary, '}')

                    if clade.is_terminal():
                        zones = [max(zones)]
                    if min(zones) == max(zones):
                        zones = [max(zones)]

            if clade.is_terminal() and len(zones) < 1:
                zones = [olderZone]
            # print(zones)
            print('')
            # find gaps in host time zone
            for z in list(range(min(zones), max(zones))):
                if z not in zCoverage:
                    zCoverage.append(z)

            zones = [str(z) for z in sorted(zones)]
            dicZones[terminals] = ','.join(zones)


        # export annotated trees in FigTree format
        new_tree = tree
        for clade in new_tree.find_clades():
            terminals = "-".join(sorted([term.name for term in clade.get_terminals()]))
            if terminals in dicZones.keys():
                clade.confidence = None
                clade.comment = None
                if str(clade.name) == 'None':
                    zones = dicZones[terminals].split(',')
                    if len(zones) > 1:
                        minimum = min(zones)
                        maximum = max(zones)
                        # print(minimum, maximum)
                        clade.name = '[' + str(minimum) + ',' + str(maximum) + ']'
                    else:
                        clade.name = '[' + str(zones[0]) + ']'
            else:
                clade.comment = None

        Phylo.write([new_tree], output + "_zones.tree", 'nexus')


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
        # print(treeString)
        return(treeString)


    # setZones(beast2nwk(hTree, 'hTree'), 'hTree')
    # setZones(beast2nwk(vTree, 'vTree'), 'vTree')

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
    janeOutput = open('virus-host_JaneTree.nex', "w")
    janeOutput.write('#NEXUS\nbegin host;\ntree host=(' + setZones(beast2nwk(hTree, 'hTree'), 'hTree')[:-2] + hGraft + '\nendblock;\n')
    janeOutput.write('\nbegin parasite;\ntree parasite=(' + setZones(beast2nwk(vTree, 'vTree'), 'vTree')[:-2] + vGraft + '\nendblock;\n')

    # output the list of virus-host pairs
    janeOutput.write('\nbegin distribution;\nRange\n')
    for virus, host in vhPairs.items():
        if host == list(vhPairs.values())[-1]:
            janeOutput.write('\t' + virus + ' : ' + host + ';\nendblock;')
        else:
            janeOutput.write('\t' + virus + ' : ' + host + ',\n')
    janeOutput.close()

    print('\n### Files successfully converted!\n')
    # print('...', gapZones)
    # print(',,,', zCoverage)