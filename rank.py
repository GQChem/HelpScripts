import argparse

parser = argparse.ArgumentParser(description='Makes csv files from ProteinMPNN.fa output ')
parser.add_argument('JOB_FOLDER', type=str, help = "folder containing af2_NN.log files")
parser.add_argument('af2_out_folder', type=str, help = "path to af2 generated pdbs")
parser.add_argument('pdb_file', type=str, help = "RMSD will be included as a metrics. Write - otherwise, don't leave empty!")
parser.add_argument('rank_output_csv_file', type=str, help = "where to save")
parser.add_argument('metric', type=str, help = "pLDDT or pTM or RMSD")
parser.add_argument('pymol_pse_file', type=str, help = "Path to pymol session to be created")
parser.add_argument('pymol_best_pse', type=int, help = "Create a pymol session contaning the N best models aligned with the original protein and colored by pLDDT")
parser.add_argument('only_first', type=bool, help = "Only compare the best folding of each sequence generated")

# Parse the arguments
args = parser.parse_args()

#Check prody is part of your environment, otherwise install it beforehand
import os

import pymol
from pymol import cmd

if args.pdb_file.endswith(".pdb"):
    cmd.load(args.pdb_file, "original")
def calculate_rmsd(folded_path):
    # Initialize PyMOL in headless mode (no GUI)
    pymol.pymol_argv = ['pymol', '-qc']  # -q for quiet, -c for no GUI
    pymol.finish_launching()
    # Load the two protein structures
    cmd.load(folded_path, "folded")
    # Align the proteins and calculate RMSD
    rmsd = cmd.align("original", "folded")[0]  # cmd.align returns a tuple, RMSD is the first element
    cmd.delete("folded")
    return rmsd

data = []
ranked_data = []
#Find log files in JOB folder
logfiles = [lf for lf in os.listdir(args.JOB_FOLDER) if lf.startswith("af2") and lf.endswith("log")]
for logfile in logfiles:
    full_path_log = os.path.join(args.JOB_FOLDER,logfile)
    with open(full_path_log,"r") as af2log:
        scores = dict()
        query_found = False
        reranking_found = False
        for logline in af2log:
            line = logline.strip()
            if line.startswith("Exception"):
                query_found = False
                reranking_found = False
                continue
            if not query_found:
                if "Query" in line:
                    comps = line.split(' ') #date time Query n/N NAME (length LENGTH)
                    scores["name"] = comps[4]
                    scores["length"] = comps[-1][:-1]
                    query_found = True
            else:
                if not reranking_found:
                    if "reranking" in line:
                        reranking_found = True
                else:
                    comps = line.split(' ') #date time rank_00N_..._model_M_seed_x pLDDT=... pTM=..
                    folded_pdb = comps[2]
                    rank = folded_pdb.split("rank_")[1][2]
                    scores["model"] = folded_pdb.split("model_")[1][0]
                    scores["pLDDT"] = comps[3].split("=")[1]
                    scores["pTM"] = comps[4].split("=")[1]    
                    folded_pdb_file = os.path.join(args.af2_out_folder,f"{name}_unrelaxed_{folded_pdb}"+".pdb")
                    if args.pdb_file.endswith(".pdb"):
                        name = scores["name"]
                        if os.path.exists(folded_pdb_file):
                            try:
                                scores["RMSD"] = "{:.4f}".format(calculate_rmsd(folded_pdb_file))
                            except Exception:
                                scores["RMSD"] = "-"
                        else:
                            scores["RMSD"] = "-"
                    if os.path.exists(folded_pdb_file):
                        scores["path"] = folded_pdb_file
                    else:
                        scores["path"] = "-"
                    data.append(scores)
                    scores = dict()
                    if args.only_first or rank == "5":
                        query_found = False
                        reranking_found = False

if len(data) == 0: 
    print("Ranking script couldn't parse output correctly")
else:
    if args.metric in list(data[0].keys()):
        small_to_big = args.metric == "RMSD"
        ranked_data = sorted(data,key=lambda x: float(x[args.metric]),reverse=not small_to_big)
    else:
        print(f"Error while ranking, key {args.metric} not found")
        ranked_data = data
    with open(args.rank_output_csv_file,"w") as csv_file:
        keys = list(ranked_data[0].keys())
        for i,key in enumerate(keys):
            if i > 0: csv_file.write(",")
            csv_file.write(key)
        for ranked_scores in ranked_data:
            csv_file.write("\n")
            for i,key in enumerate(keys):
                if i > 0: csv_file.write(",")
                csv_file.write(ranked_scores[key])

"""
PYMOL PSE CREATION
"""
def plddt(selection="all"):    
    blue_rgb = [0,76,202]
    blue = []
    for c in blue_rgb:
        blue.append(c/255.0)
    lightblue_rgb = [73, 196, 238]
    lightblue = []
    for c in lightblue_rgb:
        lightblue.append(c/255.0)
    yellow_rgb = [255, 213, 57]
    yellow = []
    for c in yellow_rgb:
        yellow.append(c/255.0)
    orange_rgb = [255, 113, 67]
    orange = []
    for c in orange_rgb:
        orange.append(c/255.0)      
    #select and color blue
    blue_upper = 100.0
    blue_lower = 90.0
    blue_sel_str = selection + " & ! b < " + str(blue_lower) + " & ! b > " + str(blue_upper)
    cmd.select('very_high', blue_sel_str)
    cmd.set_color('blue_plddt', blue)
    cmd.color('blue_plddt', 'very_high')
    #select and color lightblue
    lightblue_upper = 90.0
    lightblue_lower = 70.0
    lightblue_sel_str = selection + " & ! b < " + str(lightblue_lower) + " & ! b > " + str(lightblue_upper)
    cmd.select('confident', lightblue_sel_str)
    cmd.set_color('lightblue_plddt', lightblue)
    cmd.color('lightblue_plddt', 'confident')
    #select and color yellow
    yellow_upper = 70.0
    yellow_lower = 50.0
    yellow_sel_str = selection + " & ! b < " + str(yellow_lower) + " & ! b > " + str(yellow_upper)
    cmd.select('low', yellow_sel_str)
    cmd.set_color('yellow_plddt', yellow)
    cmd.color('yellow_plddt', 'low')
    #select and color orange
    orange_upper = 50.0
    orange_lower = 0.0
    orange_sel_str = selection + " & ! b < " + str(orange_lower) + " & ! b > " + str(orange_upper)
    cmd.select('very_low', orange_sel_str)
    cmd.set_color('orange_plddt', orange)
    cmd.color('orange_plddt', 'very_low')
    cmd.deselect()
cmd.extend('rank_plddt', plddt)

if len(ranked_data) > 0 and args.pymol_best_pse > 0:
    N_best = args.pymol_best_pse if args.pymol_best_pse < len(ranked_data) else len(ranked_data)
    for i in range(N_best):
        scores = ranked_data[i]
        base_name = scores["name"].split["_design_"][0][:-4]
        design,sequence = scores["name"].split["_design_"][1].split('_')
        model = scores["model"]
        short_name = f"{base_name}_d{design}s{sequence}m{model}"
        cmd.load(scores[i], short_name)
        cmd.do(f"rank_plddt {short_name}")
        # Align the proteins and calculate RMSD
        rmsd = cmd.align(short_name, "original")[0]  # cmd.align returns a tuple, RMSD is the first element
    cmd.save(args.pymol_pse_file)

cmd.quit()