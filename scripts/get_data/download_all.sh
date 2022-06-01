curr_dir=$(pwd)

# Download latest AlphaFold Protein Structure Database
cd ../../data/AFDB/
#1 sequences
wget http://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta
#2 structures
mkdir structures
cd structures
wget ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest/*
cd $curr_dir


# Download latest DisProt
# TO DO MANUALLY : https://www.disprot.org/download (consensus, fasta)
# save fasta file in : ../../data/DisProt/downloads/

# Download latest SCOPe (ASTRAL)
# CHOSE LAST VERSION MANUALLY : https://scop.berkeley.edu/astral/ (chose sequences with less than 95% identity to each other)
# To download latest to this date (31 May 2022):
cd ../../data/SCOPe/downloads/
wget https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-sel-gs-bib-95-2.08.fa
mv  astral-scopedom-seqres-gd-sel-gs-bib-95-2.08.fa astral_95.fasta
cd $curr_dir

# Download latest OPM
cd ../../data/OPM/downloads/
wget -c https://lomize-group-opm.herokuapp.com//primary_structures?fileFormat=csv -O proteins.csv
wget -c https://lomize-group-opm.herokuapp.com//structure_subunits?fileFormat=csv -O subunits.csv
cd structures
wget https://storage.googleapis.com/opm-assets/pdb/tar_files/Alpha-helical_polytopic.tar.gz
wget https://storage.googleapis.com/opm-assets/pdb/tar_files/Beta-barrel_transmembrane.tar.gz
wget https://storage.googleapis.com/opm-assets/pdb/tar_files/Bitopic_proteins.tar.gz
tar -zxvf Alpha-helical_polytopic.tar.gz -C polyotpic
tar -zxvf Beta-barrel_transmembrane.tar.gz -C beta
tar -zxvf Bitopic_proteins.tar.gz -C bitopic

cd $curr_dir
