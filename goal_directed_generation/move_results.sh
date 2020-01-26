#!/bin/bash

export TASK=chembl_wc

mv best_from_chembl/goal_directed_results.json hard/$TASK/best_from_chembl.json
mv smiles_ga/goal_directed_results.json hard/$TASK/smiles_ga.json
mv smiles_lstm_hc/goal_directed_results.json hard/$TASK/smiles_lstm.json
mv graph_ga/goal_directed_results.json hard/$TASK/graph_ga.json

cd hard

source activate rdkit
python json2csv.py -i $TASK/best_from_chembl.json
python json2csv.py -i $TASK/smiles_ga.json
python json2csv.py -i $TASK/smiles_lstm.json
python json2csv.py -i $TASK/graph_ga.json
