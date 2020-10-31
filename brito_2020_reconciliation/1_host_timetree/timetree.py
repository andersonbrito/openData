from Bio import Phylo
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--timetree", required=True, help="Tree downloaded from TimeTree")
    parser.add_argument("--list", required=True, help="Sublist of taxa to the extracted from tree")
    parser.add_argument("--clades", required=False, help="TSV files with list of clades in the new tree")
    parser.add_argument("--output", required=True, help="Pruned tree including only taxa in list")
    args = parser.parse_args()

    timetree = args.timetree
    list = args.list
    groups = args.clades
    output = args.output

    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/github/openData/brito_2020_reconciliation/host_timetree/"
    # # path = "/Users/anderson/GLab Dropbox/Anderson Brito/past&future/PhD/works/phylog/species_trees/viral_sppTrees/rhv04_virevol1/trees/ba/run0_host/"
    # timetree = path + "host_timetree.tree"
    # list = path + 'list_hosts.txt'
    # output = path + 'output2.tree'

    if groups == None:
        groups = 'clades.tsv'

    taxa_list = [target.strip() for target in open(list, "r").readlines() if target[0] not in ['\n', '#']]

    print('\n### Starting tree file processing...\n')
    tree = Phylo.read(timetree, 'newick')
    clades = [str(clade.name).replace('_', ' ') for clade in tree.find_clades() if clade.name != None]
    partial = []
    keep = []
    for taxon in taxa_list:
        if taxon in clades:
            # print(taxon)
            keep.append(taxon)
        else:
            genus_target = taxon.split(' ')
            found = ''
            for species in clades:
                species = species.split(' ')
                if len(genus_target) == 1:
                    if species[0] == genus_target[0] and found == '':
                        found = ' '.join(species)
                        # print('\t' + taxon + ' (' + found + ')')
                        keep.append(found)
                        partial.append(taxon)
                        continue

    not_found = [spp for spp in taxa_list if spp not in keep+partial]

    remove = [t for t in clades if t not in keep]
    c = 1
    for taxon in remove:
        taxon = taxon.replace(' ', '_')
        tree.prune(target=taxon)
        print(str(c) + ' - ' + taxon + ' was filtered out')
        c += 1


    if len(not_found) > 0:
        print('\n### These taxa were not found in TimeTree')
        for entry in not_found:
            print('\t- ' + entry)

    print('')
    c = 1
    clades = [str(clade.name).replace('_', ' ') for clade in tree.find_clades() if clade.name != None]
    outfile2 = open(output, 'w')
    for clade in tree.find_clades():
        if str(clade.name) == 'None':
            # print('None' + "\t" + str([term.name for term in clade.get_terminals()]))
            terminals = sorted([term.name for term in clade.get_terminals()])
            line = str(c) + '\t' + str(len(terminals)) + '\t' + ", ".join(terminals)
            outfile2.write(line + '\n')
            # print(str(c) + '\t' + str(len(terminals)) + '\t' + ", ".join(terminals))
            c += 1

    # save tree
    Phylo.write([tree], output, 'newick')
    print('\nTree file successfully filtered: \'' + output + '\n')

    # print(keep)
    # print(taxa_list)