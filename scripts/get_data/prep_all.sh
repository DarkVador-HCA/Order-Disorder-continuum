
# Prep DisProt (remove empty lines)
grep "\S" ../../data/DisProt/downloads/disprot_r2021_06.fasta > ../../data/DisProt/downloads/sequences.fasta

# Prep SCOPe (remove 'X', put sequences in capital letters, create id_cats.csv file, and keep only sequences from a,b,c,d classes)
python3 prep_SCOPe.py

#Prep OPM (get sequences from PDBs, remove long soluble segments (>30 aa), create id_cats.csv)
python3 prep_OPM.py