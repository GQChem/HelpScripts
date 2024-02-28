import argparse

parser = argparse.ArgumentParser(description='Copy best pdb from one cycle to the next one in PMPNN_AF2_Loop')
parser.add_argument('old_rank_output_csv_file', type=str)
parser.add_argument('pdb_file', type=str)

# Parse the arguments
args = parser.parse_args()

import shutil

#retrieve path of best one
with open(args.old_rank_output_csv_file,"r") as ranked_csv:
    data_keys = ranked_csv.readline().strip().split(',')
    seq_index = data_keys.index("sequence")
    path_index = data_keys.index("path")
    best_data = ranked_csv.readline().strip().split(',')
    best_path = best_data[path_index]
    best_seq = best_data[seq_index]
    print(f"Best PDB: {best_path}")
    print(f"Sequence: {best_seq}")
    shutil.copy(best_path,args.pdb_file)  