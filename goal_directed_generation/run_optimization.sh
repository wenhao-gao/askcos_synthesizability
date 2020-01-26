#!/bash/bin

## Sample script of goal-directed generation

source activate rdkit

## Analyze data for mcts
# python -m graph_mcts_moses.analyze_dataset

## Pretrain model for LSTM
# CUDA_VISIBLE_DEVICES=7 python -m smiles_lstm_hc_moses.train_smiles_lstm_model

export SUITE=v5
export DATA=guacamol_v1_all
# export DATA=moses

python -m best_from_chembl.goal_directed_generation --smiles_file data/$DATA.smiles --suite $SUITE
python -m smiles_ga.goal_directed_generation --smiles_file data/$DATA.smiles --suite $SUITE
CUDA_VISIBLE_DEVICES= python -m smiles_lstm_hc.goal_directed_generation --smiles_file data/$DATA.smiles --suite $SUITE
# CUDA_VISIBLE_DEVICES= python -m smiles_lstm_hc_moses.goal_directed_generation --smiles_file data/$DATA.smiles --suite $SUITE
python -m graph_ga.goal_directed_generation --smiles_file data/$DATA.smiles --suite $SUITE
# python -m graph_mcts_guacamol.goal_directed_generation --suite $SUITE
