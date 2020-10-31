import pandas as pd
import baltic as bt
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter nextstrain metadata files re-formmating and exporting only selected lines",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--timetree", required=True, help="Tree downloaded from TimeTree")
    parser.add_argument("--tmrcas", required=True, help="TSV file with clades and divergence times")
    parser.add_argument("--output", required=True, help="Pruned tree including only taxa in list")
    args = parser.parse_args()

    timetree = args.timetree
    tmrcas = args.tmrcas
    output = args.output

    # path = "/Users/anderson/GLab Dropbox/Anderson Brito/past&future/PhD/works/phylog/species_trees/viral_sppTrees/rhv04_virevol1/trees/ba/run0_host/"
    # timetree = path + "host_tree.tree"
    # tmrcas = path + 'tmrcas.txt'
    # output = path + 'new_host_tree.nexus'

    all_traits = ['node_name', 'height_95%_HPD']

    # load tree
    tree = bt.loadNewick(timetree)#, tip_regex='_([0-9\-]+)$')
    # print(tree)

    # tmrca dataframe
    df = pd.read_csv(tmrcas, encoding='utf-8', sep='\t', dtype='str')
    # print(df)

    # df['members'] = df['members'].apply(lambda x: ', '.join(sorted(x.split(','))))
    # print(df['members'].to_list())

    print('Starting tree file processing...')
    # transfer supporting value from a newick tree
    for k in sorted(tree.Objects, key=lambda q: q.height):  ## iterate over branches from most recent to oldest
        if k.branchType == 'node':  ## can only sort nodes
            terminals = ", ".join(sorted([leaf for leaf in k.leaves]))
            # print(terminals)
            if terminals in df['members'].to_list():
                traits = {}
                range = df.loc[df['members'] == terminals, 'range'].values[0].replace(' ', '')
                upper = range.split('-')[0]
                lower = range.split('-')[1]
                node_name = df.loc[df['members'] == terminals, 'group'].values[0]
                k.traits['node_name'] = node_name # add as a trait
                k.traits['height_95%_HPD'] = [float(upper), float(lower)]


    print('Number of objects in subtree: %d, annotations to include: %s, use number encoding for tips: %s'%(len(tree.Objects),'None','False'))
    # print('%s\n\n'%(tree.toString(traits=all_traits, numName=True, nexus=True,verbose=False))) ## can output to a minimal NEXUS format for viewing in figtree with correct trait parsing

    new_tree = '%s\n\n' % (tree.toString(traits=all_traits, numName=True, nexus=True,verbose=False)) ## can output to a minimal NEXUS format for viewing in figtree with correct trait parsing
    # print(new_tree)
    # print(tree.make_tree)

    outfile = open(output, 'w')
    outfile.write(new_tree)
    print('\nTree file successfully reformated: \'' + output + '\'\n')
