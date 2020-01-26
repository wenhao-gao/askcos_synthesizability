import pandas as pd
import json
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', default=None,
                        help='Input JSON file.')
    args = parser.parse_args()
    
    print(f'Reading data ......')
    with open(args.input) as f:
        data = json.load(f)

    smiles = []
    source = []
    prop = []

    METHOD = args.input.split('.')[0]

    print(f'Extracting data ......')
    for result in data['results']:
        for smi in result['optimized_molecules']:
            smiles.append(smi[0])
            prop.append(result['benchmark_name'])
            source.append(METHOD)

    assert len(smiles) == len(source) == len(prop)
    print(f'Finish extratcing, pass the check with length!')

    df = pd.DataFrame({
        'SMILES': smiles,
        'property': prop,
        'source': source
    })

    print(f'Saving to csv file ......')
    df.to_csv(METHOD + '.csv', index=False)
    print(f'All finished!')
