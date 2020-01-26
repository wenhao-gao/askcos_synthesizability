"""
The script to test the synthesizability of a batch of molecules using Tree Builder
section of ASKCOS (http://askcos.mit.edu). It will do retrosynthesis analysis for
the molecules in SMILES string under column "SMILES"(default). To see more guidance:
    $ python batch_TB.py --help

Written by Wenhao Gao in Jun. 2019
"""
import pandas as pd
import argparse
import requests
import json
import time

if __name__ == '__main__':

    # Parse the arguments from command
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help='The csv file containing smiles represented molecules')
    parser.add_argument("-n", "--n_start", default=0,
                        help='The csv file containing smiles represented molecules')
    parser.add_argument("--col", default='SMILES',
                        help='The column we want calculate in csv file')
    parser.add_argument("-o", "--output",  default='./output.txt',
                        help='The name for the output file that summarize the information of the task')
    args = parser.parse_args()

    # Replace this with the IP address of your ASKCOS server
    HOST = 'https://xx.xxx.xx.xxx'

    df = pd.read_csv(args.i)
    mols = df[args.col]

    t1 = time.time()
    pre = int(args.n_start)
    num_good = 0
    num_find = 0
    num_total = 0

    with open(args.output, 'wt') as f:
        for i in range(pre, len(mols)):
            print('Starting calculating No. %i' % (i))
            print(mols[i])
            
            # Check if is buyable
            params = {
                'smiles': mols[i] # required
            }
            resp = requests.get(HOST+'/api/price/', params=params, verify=False)
            
            if resp.json()['price'] == 0:
                # Parameters for Tree Builder
                params = {
                    'smiles': mols[i],  # required

                    # optional
                    'max_depth': 9,
                    'max_branching': 25,
                    'expansion_time': 60,
                    'max_ppg': 100,
                    'template_count': 1000,
                    'max_cum_prob': 0.999,
                    'chemical_property_logic': 'none',
                    'max_chemprop_c': 0,
                    'max_chemprop_n': 0,
                    'max_chemprop_o': 0,
                    'max_chemprop_h': 0,
                    'chemical_popularity_logic': 'none',
                    'min_chempop_reactants': 5,
                    'min_chempop_products': 5,
                    'filter_threshold': 0.1,
                    'return_first': 'true'
                }

                # For each entry, repeat to test up to 5 times if got error message
                for _ in range(5):
                    print('Trying to send the request, for the %i times now' % (_ + 1))
                    resp = requests.get(HOST + '/api/treebuilder/', params=params, verify=False)
                    if 'error' not in resp.json().keys():
                        print('A result is found!')
                        num_good += 1
                        if len(resp.json()['trees']) != 0:
                            num_find += 1
                            print(f'{len(resp.json()["trees"])} was found')
                        break
                    # pprint(resp.json())
            else:
                num_good += 1
                num_find += 1

            print()
            num_total += 1
            name = args.i.split('.')[0]
            with open(name + '_' + str(i) + '.json', 'w') as f_data:
                json.dump(resp.json(), f_data)

            # Print out intermediate results
            if i % 10 == 0:
                print('Finish calculating No. %i' % (i))
                print('Among them, %.3f are good' % (num_good / float(num_total)))
                print('Among them, %.3f finds a path' % (num_find / float(num_total)))
                print('Up to now, %.2f seconds in average were used' %((time.time() - t1)/num_total))
                print()

        # Print out summary of task to the output file
        t = time.time() - t1
        print('Finish tree building for %s, total time used: %.3f, average time used: %.3f'
              %(args.i, t, t/num_total), file=f)
        print('Among them, %.3f are good' % (num_good / float(num_total)), file=f)
        print('Among them, %.3f finds a path' % (num_find / float(num_total)), file=f)
