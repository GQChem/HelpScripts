#!/usr/bin/bash
echo A. ALPHAFOLD
echo
echo A.i Cloning repository...
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
echo
echo A.ii Installing...
bash install_colabbatch_linux.sh
rm install_colabbatch_linux.sh

echo
echo B. RFDIFFUSION
echo
echo B.i Cloning repository...
git clone https://github.com/RosettaCommons/RFdiffusion.git
echo
echo B.ii Downloading models...
cd RFdiffusion
mkdir models && cd models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt
echo
echo B.iii Setting up environment ProteinEnv
cd ../env
wget https://gitlab.uzh.ch/locbp/proteinmodels/-/blob/443bc562c91112ec2606435c48b0514ce6a5b249/ProteinEnv.yml
module load mamba
mamba env create -f ProteinEnv.yml 
mamba activate ProteinEnv
echo
echo B.iv Installing...
cd SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install
cd ../..
pip install -e . 
echo
echo B.v Installing Pymol...
conda install -c conda-forge -c schrodinger pymol-bundle

echo
echo C. PROTEINMPNN
echo
echo C.i Cloning repository...
cd ..
git clone https://github.com/dauparas/ProteinMPNN

echo
echo D. HELP SCRIPTS
echo
mkdir HelpScripts
cd HelpScripts
wget https://gitlab.uzh.ch/locbp/proteinmodels/-/blob/e5f8ae8c4c4e6ed64280cfb57877335a5d0cc87a/HelpScripts/af2_make_pse.py
wget https://gitlab.uzh.ch/locbp/proteinmodels/-/blob/e5f8ae8c4c4e6ed64280cfb57877335a5d0cc87a/HelpScripts/fa_to_csv.py
wget https://gitlab.uzh.ch/locbp/proteinmodels/-/blob/e5f8ae8c4c4e6ed64280cfb57877335a5d0cc87a/HelpScripts/make_fixed_dict.py
wget https://gitlab.uzh.ch/locbp/proteinmodels/-/blob/e5f8ae8c4c4e6ed64280cfb57877335a5d0cc87a/HelpScripts/rank.py
cd ..

echo
echo E. COPYING JUPYTER NOTEBOOKS in home
cd ..
wget https://gitlab.uzh.ch/locbp/proteinmodels/-/blob/e5f8ae8c4c4e6ed64280cfb57877335a5d0cc87a/AlphaFold.ipynb
wget https://gitlab.uzh.ch/locbp/proteinmodels/-/blob/e5f8ae8c4c4e6ed64280cfb57877335a5d0cc87a/ProteinMPNN.ipynb
wget https://gitlab.uzh.ch/locbp/proteinmodels/-/blob/e5f8ae8c4c4e6ed64280cfb57877335a5d0cc87a/RFdiffusion.ipynb
echo
echo DONE!
