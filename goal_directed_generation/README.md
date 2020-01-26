# Synthesizability of generative methods

This repository contains codes and results to analyze synthesizability using 
[`ASKCOS`](http://askcos.mit.edu) for molecular generation and optimization algorithms.  
A more detailed explanation of the project
can be found in our [paper](https://arxiv.org/abs/1811.09621).

## Usage with ASKCOS API
If you want to test a set of molecules in SMILES (e.g. test.csv), please replace the HOST 
address in ```tb_analysis/batch_TB.py``` to your ASKCOS server IP, then
```bash
python batch_TB.py -i test.csv 
```
The results will be stored at the same directory with python file in json format.
One can also define the index of molecule to start with (default is 0) 
and the column name for SMILES strings (default is ""SMILES"). To see help message:
```bash
python batch_TB.py -h
```

## Usage with local ASKCOS

Under development

## Dependencies

### If you just want to test synthesizability
- json
- numpy
- pandas
- requests


### If you want to use the generative models and benchmarks to reproduce the result
Please refer to 
[`MOSES`](https://github.com/molecularsets/moses) for distribution learning implementations and
[`Guacamol Baselines`](https://github.com/BenevolentAI/guacamol_baselines) for goal-directed learning implementations

