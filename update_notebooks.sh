echo HelpScripts
rm -rf HelpScripts
mkdir HelpScripts
cd HelpScripts
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/af2_make_pse.py
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/fa_to_csv.py
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/make_fixed_dict.py
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/rank.py
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/pmpnn_af2_loop_find_best.py
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/pmpnn_af2_loop_fa_to_csv.py
cd ..

echo
echo Notebooks
cd ..
rm AlphaFold.ipynb
rm ProteinMPNN.ipynb
rm RFdiffusion.ipynb
rm PMPNN_AF2_Loop.ipynb
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/AlphaFold.ipynb
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/ProteinMPNN.ipynb
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/RFdiffusion.ipynb
wget https://raw.githubusercontent.com/GQChem/HelpScripts/main/PMPNN_AF2_Loop.ipynb
echo
echo UPDATE COMPLETE
