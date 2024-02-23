import argparse

parser = argparse.ArgumentParser(description='Reads RFD log to create fixed positions')
parser.add_argument('JOB_FOLDER', type=str, help="")
parser.add_argument('rfd_log_file', type=str, help = "Path to rfd.log file")
parser.add_argument('FIXED', type=str, help="rfd or pymol selection")
parser.add_argument('FIXED_CHAIN', type=str, help="A or B or whatever")
parser.add_argument('fixed_jsonl_file', type=str, help="output file")

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

import os
design_names = [d[:-4] for d in os.listdir(args.JOB_FOLDER) if d.endswith(".pdb")]
fixed_dict = dict()
if args.FIXED == "rfd": #derive from logfile
    with open(args.rfd_log_file,"r") as rfdlog:
        all_log = [line.strip() for line in rfdlog.readlines()]
        for name in design_names:
            fixed_dict[name] = dict()
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
else:
    for name in design_names:
        fixed_dict[name] = dict()
        fixed_dict[name][args.FIXED_CHAIN] = sele_to_list(args.FIXED)
with open(args.fixed_jsonl_file,"w") as jsonl_file:
    #Python converts dictionaries to string having keys inside '', json only recognises ""
    jsonl_file.write(str(fixed_dict).replace("\'","\""))

#Dictionary of fixed positions looks like this
#{"5TTA": {"A": [1, 2, 3, 7, 8, 9, 22, 25, 33], "B": []}, "3LIS": {"A": [], "B": []}}