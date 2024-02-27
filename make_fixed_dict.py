import argparse

parser = argparse.ArgumentParser(description='Reads RFD log to create fixed positions')
parser.add_argument('JOB_FOLDER', type=str, help="")
parser.add_argument('rfd_log_file', type=str, help = "Path to rfd.log file")
parser.add_argument('pLDDT_thr', type=float, help="Will consider residues with pLDDT > threshold as fixed ")
parser.add_argument('FIXED', type=str, help="rfd or pymol selection")
parser.add_argument('FIXED_CHAIN', type=str, help="A or B or whatever")
parser.add_argument('fixed_jsonl_file', type=str, help="output file")
parser.add_argument('sele_csv_file', type=str, help="output file with selections of fixed and mobile parts")

# Parse the arguments
args = parser.parse_args()

#Converts a pymol selection into an array
def sele_to_list(s):
    a = []
    if s == "": return a
    elif '+' in s:
        plus_parts = s.split('+')
        for pp in plus_parts:
            if '-' in pp:
                min,max = pp.split('-')
                for ri in range(int(min),int(max)+1):
                    a.append(ri)
            else:
                a.append(int(pp))
    else:
        if '-' in s:
            min,max = s.split('-')
            for ri in range(int(min),int(max)+1):
                a.append(ri)
        else:
            a.append(int(s))     
    return a
def list_to_sele(a):
    s = ""
    i = 0
    while i < len(a):   
        if i > 0: s += "+"
        s += f"{a[i]}"   
        #represent consecutive indeces with a dash
        if i < len(a) - 1:
            if int(a[i])+1 == int(a[i+1]):
                s += "-"
                j = i + 2
                while j < len(a):
                    if int(a[j]) != int(a[j-1])+1: break
                    j += 1
                i = j - 1
                s += f"{a[i]}" 
        i += 1        
    return s

import os
import pymol
from pymol import cmd

design_names = [d[:-4] for d in os.listdir(args.JOB_FOLDER) if d.endswith(".pdb")]
fixed_dict = dict()
mobile_dict = dict()
if args.FIXED == "rfd": #derive from logfile
    with open(args.rfd_log_file,"r") as rfdlog:
        all_log = [line.strip() for line in rfdlog.readlines()]
        for name in design_names:
            fixed_dict[name] = dict()
            mobile_dict[name] = dict()
            seq = ""
            #The sequence input is found a few lines after RFD says it is designing something
            log_found = False
            for line in all_log:
                if not log_found:
                    if "Making design" in line and line.endswith(name):
                        log_found = True
                        continue
                if log_found:
                    if "Sequence init" in line:
                        seq = line.split(" ")[-1]
                        break
            fixed_dict[name][args.FIXED_CHAIN] = [i+1 for i, x in enumerate(seq) if x != "-"]
            mobile_dict[name][args.FIXED_CHAIN] = [i+1 for i, x in enumerate(seq) if x == "-"]
else:
    #not rfd
    if args.pLDDT_thr < 100:
        # Initialize PyMOL in headless mode (no GUI)
        pymol.pymol_argv = ['pymol', '-qc']  # -q for quiet, -c for no GUI
        pymol.finish_launching()
    for name in design_names:
        fixed_dict[name] = dict()
        mobile_dict[name] = dict()
        if args.pLDDT_thr < 100:
            pdb_file = os.path.join(args.JOB_FOLDER,name+".pdb")
            cmd.load(pdb_file,"prot")
            fixed_residues = []
            mobile_residues = []
            atom_iterator = cmd.get_model("prot and name CA")
            for atom in atom_iterator.atom:
                if atom.b < args.pLDDT_thr:
                    if int(atom.resi) not in mobile_residues:
                        mobile_residues.append(int(atom.resi))
                else:
                    if int(atom.resi) not in fixed_residues:
                        fixed_residues.append(int(atom.resi))
            cmd.delete("prot")
            fixed_dict[name][args.FIXED_CHAIN] = fixed_residues[:]
            mobile_dict[name][args.FIXED_CHAIN] = mobile_residues[:]
        else:
            fixed_dict[name][args.FIXED_CHAIN] = sele_to_list(args.FIXED)
            mobile_dict[name][args.FIXED_CHAIN] = []

with open(args.fixed_jsonl_file,"w") as jsonl_file:
    #Python converts dictionaries to string having keys inside '', json only recognises ""
    jsonl_file.write(str(fixed_dict).replace("\'","\""))
with open(args.sele_csv_file,"w") as csv_file:
    csv_file.write("id,fixed,mobile")
    for id in fixed_dict.keys():
        fixed = list_to_sele(fixed_dict[id][args.FIXED_CHAIN])
        mobile = list_to_sele(mobile_dict[id][args.FIXED_CHAIN])
        csv_file.write("\n")
        csv_file.write(f"{id},{fixed},{mobile}")

#Dictionary of fixed positions looks like this
#{"5TTA": {"A": [1, 2, 3, 7, 8, 9, 22, 25, 33], "B": []}, "3LIS": {"A": [], "B": []}}