import os
import pandas as pd
import itertools
from io import StringIO
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Parse Jane files generated using CLI software",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--directory", required=True, help="Folder containing Jane CLI output file")
    args = parser.parse_args()

    path = './' + args.directory + '/'


    # path = '/Users/anderson/GLab Dropbox/Anderson Brito/github/openData/brito_2020_reconciliation/jane_cli/cli_output/'
    filenames = [path + f for f in os.listdir(path) if f.startswith("results_CR")]

    output1 = 'event_frequency.tsv'
    output2 = 'eventCounts.tsv'


    data1 = []
    data2 = {}
    for file in filenames:
        print(file)
        diff = {}
        num_regime = file.split('_CR')[1].split('.')[0]
        cost_regime = str('0' * (3 - len(num_regime))) + str(num_regime)
        start_tree = 0
        artificial_nodes = []

        if 'regime' not in data2:
            data2['regime'] = [cost_regime]
        else:
            data2['regime'] += [cost_regime]

        # identify and list artificial taxa to be ignored
        all_lines = open(file, 'r').readlines()
        for num, line in enumerate(all_lines):
            # print(line)
            if "Parasite Tree:" in line:
                # print('\t### found1')
                start_tree += 1
            if start_tree < 1:
                # print('not quite there yet')
                pass
            else:
                node_name = line.split(')')[0][1:] + '`'
                if "Name: virus" in line: # find artificial taxon and store node name
                    if node_name not in artificial_nodes:
                        artificial_nodes.append(node_name)
                else:
                    # print(artificial_nodes)
                    for node in artificial_nodes:
                        # print(node)
                        id = 'Name: ' + node
                        if id in line: # store extra artificial node names
                            if node_name not in artificial_nodes:
                                artificial_nodes.append(node_name)

            if "Best Solution:" in line:
                # print('\t### found2')
                start_tree += 1
                if 'CO' not in diff:
                    diff['Cospeciation'], diff['Duplication'], diff['Host Switch'], diff['Loss'] = [0, 0, 0, 0]

            if start_tree < 2:
                # print('not quite there yet again')
                pass
            else:
                for node in artificial_nodes:
                    if 'Parasite Node: ' + node in line:
                        # print(all_lines[num].strip())
                        event = all_lines[num + 1].split(': ')[1].strip()
                        diff[event] += 1
                        # print(node, event)

        # print(artificial_nodes)


        last_lines = [int(line.split(':')[1].strip()) for line in open(file, 'r').readlines()[-6:-2]]
        cospeciations, duplication, host_switch, loss = last_lines
        # print(last_lines)
        new_data2 = {}
        if 'CO' not in data2:
            data2['CO'], data2['IS'], data2['HS'], data2['LO'] = [], [], [], []
        new_data2['CO'], new_data2['IS'], new_data2['HS'] = 0, 0, 0

        # find co-phylogenetic events
        def isa_group_separator(line):
            return line == '--------------------------\n'

        # prepare input file
        s = open(file).read()
        f = s.replace('==================================', '--------------------------')
        f = StringIO(f)

        for key, group in itertools.groupby(f, isa_group_separator):
            # print(key,list(group))  # uncomment to see what itertools.groupby does.
            fields = list(group)
            if key == False:  # and len(list(group)) < 10:
                if len(fields) < 10: # ignore top content in input file
                    rec_data = {}
                    if 'Parasite' in ', '.join(fields):
                        for item in fields:
                            k = item.strip().split(':')[0].strip()
                            v = item.strip().split(':')[1].strip()
                            rec_data[k] = v

                        if "Node: virus" in ', '.join(fields) or rec_data['Parasite Node'] in artificial_nodes:
                            # print('IGNORE')
                            pass
                        else:
                            if rec_data['Association type'] == 'Host Switch':
                                node = rec_data['Parasite Node'].replace('`','')
                                host_parent = rec_data['Host'].split(',')[0].replace('`','').replace('(','')
                                host_child = rec_data['Host'].split(',')[1].replace('`','').replace(')','').strip()
                                virus_parent = rec_data['Switch Target'].split(',')[0].replace('`','').replace('(','')
                                virus_child = rec_data['Switch Target'].split(',')[1].replace('`','').replace(')','').strip()
                                event_id = node + ':' + host_parent + '-' + host_child + '..' + virus_parent + '-' + virus_child
                                entry = [cost_regime, event_id, rec_data['Association type'], rec_data['Event Time']]
                                print('\t'.join(entry))
                                data1 += [entry]
                                new_data2['HS'] += 1

                            elif rec_data['Association type'] == 'Duplication':
                                node = rec_data['Parasite Node'].replace('`','')
                                host_parent = rec_data['Host'].split(',')[0].replace('`','').replace('(','')
                                host_child = rec_data['Host'].split(',')[1].replace('`','').replace(')','').strip()
                                event_id = node + ':' + host_parent + '-' + host_child
                                entry = [cost_regime, event_id, rec_data['Association type'], rec_data['Event Time']]
                                print('\t'.join(entry))
                                data1 += [entry]
                                new_data2['IS'] += 1


                            elif rec_data['Association type'] == 'Cospeciation':
                                node = rec_data['Parasite Node'].replace('`','')
                                event_id = node + ':' + rec_data['Host'].replace('`','')
                                entry = [cost_regime, event_id, rec_data['Association type'], rec_data['Event Time']]
                                data1 += [entry]
                                print('\t'.join(entry))
                                new_data2['CO'] += 1


        if cospeciations - diff['Cospeciation'] != new_data2['CO']:
            data2['CO'] += [new_data2['CO']]
        else:
            data2['CO'] += [cospeciations - diff['Cospeciation']]

        if duplication - diff['Duplication'] != new_data2['IS']:
            data2['IS'] += [new_data2['IS']]
        else:
            data2['IS'] += [duplication - diff['Duplication']]

        if host_switch - diff['Host Switch'] != new_data2['HS']:
            data2['HS'] += [new_data2['HS']]
        else:
            data2['HS'] += [host_switch - diff['Host Switch']]

        if diff['Duplication'] > 1 or diff['Host Switch'] > 0:
            data2['LO'] += [loss - (diff['Duplication'] * 2) - diff['Host Switch']]
        else:
            data2['LO'] += [loss - diff['Loss']]

        print('\n############################################\n')
        # print(data2)

    # print(data)

    # write new metadata files
    df2 = pd.DataFrame(data1, columns=list(['cost_regime', 'event_id', 'event_type', 'time_zone']))
    df2 = df2.sort_values(by='cost_regime')
    df2.to_csv(output1, sep='\t', index=False)

    df3 = pd.DataFrame.from_dict(data2)
    df3['regime'] = df3['regime'].astype(int)
    df3 = df3.sort_values(by='regime')
    df3.to_csv(output2, sep='\t', index=False)

    print('\nJane reconciliation data successfully exported!\n')
