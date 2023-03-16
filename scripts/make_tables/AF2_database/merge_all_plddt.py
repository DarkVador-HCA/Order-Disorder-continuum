import pickle as pk
import os,sys

data_dir=sys.argv[1]

all_p={}
for d in os.listdir(data_dir+"/pLDDT/"):
    for f in os.listdir(data_dir+"/pLDDT/"+d):
        prot=f.split(".")[0]            
        with open(data_dir+"/pLDDT/"+d+"/"+f,'rb') as f:
            plddt=pk.load(f)
        all_p[prot]=plddt

with open(data_dir+"/all_pLDDT.pkl","wb") as f:
    pk.dump(all_p,f)