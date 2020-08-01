import pickle
import time
import json
import requests
import numpy as np
from hyperopt import fmin, tpe, hp, STATUS_OK, Trials
import graph_ga.sa_goal_directed_generation as ga


def evaluate(hps):
    ga.bo_main(hps['mu'], hps['sigma'])
    with open('./graph_ga/sa_goal_directed_results.json', 'r') as f:
        results = json.load(f)

    NUM = 30
    HOST = 'https://35.192.98.95'

    mols = []
    props = []

    for i in range(len(results['results'])):
        temp = []
        for j in range(NUM):
            temp.append(results['results'][i]['optimized_molecules'][j][1])
        props.append(np.mean(np.array(temp)))

    for i in range(len(results['results'])):
        for j in range(NUM):
            mols.append(results['results'][i]['optimized_molecules'][j][0])

    synth = []

    for _ in range(len(results['results'])):

        num_good = 0
        num_total = 0
        num_find = 0

        for i in range(len(mols)):
            # Check if is buyable
            params = {
                'smiles': mols[i]  # required
            }
            resp = requests.get(HOST + '/api/price/', params=params, verify=False)

            if resp.json()['price'] == 0:
                # Call tree builder
                params = {
                    'smiles': mols[i],  # required

                    # optional with defaults shown
                    'max_depth': 8,
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

                    'return_first': 'true'  # default is false
                }

                for _ in range(15):
                    print('Trying to send the request, for the %i times now' % (_ + 1))
                    resp = requests.get(HOST + '/api/treebuilder/', params=params, verify=False)
                    if 'error' not in resp.json().keys():
                        print(f'A result is found!, Total path: {len(resp.json()["trees"])}')
                        num_good += 1
                        if len(resp.json()['trees']) != 0:
                            num_find += 1
                        break
                    # pprint(resp.json())
            else:
                num_good += 1
                num_find += 1

            num_total += 1

        synth.append(num_find / num_total)

    loss = np.array(props) * np.array(synth)
    return 1 - loss.mean()


def main():

    SPACE = {
        'mu': hp.normal('mu', 3.5, 2),
        'sigma': hp.lognormal('sigma', 2, 4)
    }

    def objective(hps):
        return {
            'loss': evaluate(hps),
            'status': STATUS_OK,
            # -- store other results like this
            'eval_time': time.time(),
            'other_stuff': {'type': None, 'value': [0, 1, 2]},
            # -- attachments are handled differently
            'attachments':
                {'time_module': pickle.dumps(time.time)}
        }

    trials = Trials()

    best = fmin(objective,
                space=SPACE,
                algo=tpe.suggest,
                max_evals=30,
                trials=trials)

    print(best)


if __name__ == '__main__':
    main()
