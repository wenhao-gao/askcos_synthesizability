# Results

The whole results of analyses. Raw data of synthetic pathways in json format 
can be downloaded from [`figshare`](https://figshare.com/articles/Raw_data_of_ASKCOS_synthesizability_tests/11737092).

## Synthesizability of common databases

Summarized results of common databases can be found in ```dataset.csv```. 

The columns, form left to right, refers to:

- SMILES: the SMILES string of the enquiry molecule.
- dataset: the dataset it comes from.
- len_smiles: length of smiles string.
- mean_complexity: meanComplexity, an averaged expert scoring 
(only applicable for ones form Sheridan et al's).
- sa_score: SA_Score from Ertl, P., et al's.
- sc_score: SCScore from Coley, C. W., et al's.
- tb_depth: number of synthetic steps required (0 means compound is buyable, 11 means no pathway was found)
- tb_plausibility:plausibility of the shortest pathway (product of plausibility score of each reaction
predicted by a neural network)
- tb_price: total price of starting materials.
- tb_synthesizability: synthesizable or not.

Statistics of synthetic steps required to produce the set can be found in ```count_dataset.csv```.

## Synthesizability of unoptimized molecules

Summarized results of common databases can be found in ```distribution_learning.csv```. 

The column naming is similar as above, except:

- dataset: the dataset to be trained.
- method: the model used.

Statistics of synthetic steps required to produce the set can be found in ```count_dl_chembl.csv``` 
for models trained on ChEMBL, and ```count_dl_moses.csv``` for models trained on MOSES.

## Synthesizability of optimized molecules

Summarized results of common databases can be found in ```goal_*.csv```. 

```trivial``` or ```hard``` represents which set of benchmark it is about.

The last three letters, ```c/m+w+o/a/c```, refer to the type of task. 
```c/m``` refers to the datset to start with (ChEMBL or MOSES).
```o/a/c``` refers to heuristic biasing (without biasing, biased by SA_Score, biased by SCScore).

The column naming is similar as above, except:

- property: the name of the objective function.
- rank: the rank of this molecules among all candidates according to modified objective.
- objective: the value of objective functions without modifying.
- type: the type of the task (without biasing, biased by SA_Score, biased by SCScore)

The cases with improved objective of top-1 synthesizable candidate can be found in ```change_*.csv```
