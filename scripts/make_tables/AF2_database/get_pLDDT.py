from Bio.PDB import *
import sys
import pickle as pk

pdb_file=sys.argv[1]
prot_name=pdb_file.split(".")[0].split('-model')[0].split("/")[-1]
genome=sys.argv[2]
data_dir=sys.argv[3]

parser=PDBParser()
structure=parser.get_structure(prot_name,pdb_file)
bfactors=[]
parent=""
for a in structure.get_atoms():
	par=a.get_parent()
	if par!=parent:
		bfactors.append(a.get_bfactor())
		parent=par

with open(data_dir+"/pLDDT/"+genome+"/"+prot_name+".pkl",'wb') as f:
	pk.dump(bfactors,f)
	
