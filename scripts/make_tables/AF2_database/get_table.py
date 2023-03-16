import pandas as pd
import pickle as pk
import os,sys

"""
Construit table table_AF2DB (article proteins, dataset S4) en utilisant infos de 2 dictionnaires générés par le script ?? 
"""

data_dir=sys.argv[1]

#df=pd.DataFrame(to_w,columns=["proteome","sequence id","annotation","contiguous region start","contiguous regions length (aa)","whole segment start","whole segment length (aa)","pLDDT category","structural state","HCA score","% of GATS"])


def get_plddt_cat(plddt_score):
    """
    Catégorie de confiance à partir du pLDDT
    """
    threshold_cat={50:"vl",70:"l",90:"c",100:"vh"}
    for threshold in [50,70,90,100]:
        if plddt_score<=threshold:
            return threshold_cat[threshold]

def get_ogns():
    """
    Organisme à partir de l'identifiant de séquence
    """
    ogn_prots={}
    for ogn in os.listdir("../ogn_prots"):
        with open("../ogn_prots/"+ogn,"r") as f:
            ogn_prots[ogn]=[l.strip() for l in f]
    ogn_prots["ALL_GENOMES"]=[]
    return ogn_prots

def get_info_from_fasta(fasta):
    """
    Infos tirées des fastas : annotations uniprot et nombre d'aa chargés
    """
    charged_aa=["D","E","K","R","H"]
    bin_charges,annots = {},{}
    with open(fasta,"r") as f:
        for l in f:
            if l[0]==">":
                l=l.strip().split()
                n=l[0].split(":")[1]
                annots[n]=" ".join(l[1:])
                bin_charges[n]=""
            else:
                seq=l.strip()
                for a in seq:
                    if a in charged_aa:
                        bin_charges[n]+="1"
                    else:
                        bin_charges[n]+="0"
    return bin_charges,annots


def define_contiguous_segments():
    """
    Définition des segments contigues de même pLDDT en segment foldable ou non foldable (à partir des dictionnaires hca_sum et all_pLDDT)
    """
    bornes=[]
    plddt=[]
    hca=[]
    sf_starts=[]
    prots=[]
    sf_sizes=[]
    for prot in plddt_all.keys():
        previousp=""
        previoush=20
        sf_start=0
        dSF_sizes={}
        tmp_sf_starts=[]
        if prot in dHCA.keys():
            for i in range(len(plddt_all[prot])):
                p=get_plddt_cat(plddt_all[prot][i])
                h=dHCA[prot][i]
                if h=='-':
                    h=-10.0
                h=float(h)
                if p==previousp and h==previoush:
                    bornes[-1][1]=i
                    dSF_sizes[sf_start]+=1
                else:
                    if h!=previoush:
                        dSF_sizes[i]=1
                        sf_start=i
                    else:
                        dSF_sizes[sf_start]+=1
                    prots.append(prot)
                    bornes.append([i,i])
                    plddt.append(p)
                    hca.append(h)
                    previousp=p
                    previoush=h
                    tmp_sf_starts.append(sf_start)  
            sf_sizes+=[dSF_sizes[v] for v in tmp_sf_starts]
            sf_starts+=tmp_sf_starts
    return prots,bornes,plddt,hca,sf_sizes,sf_starts


# input data#
with open(data_dir+"/hca_sum.pkl",'rb') as f:
    dHCA=pk.load(f)

with open(data_dir+"/all_pLDDT.pkl","rb") as f:
    plddt_all=pk.load(f)

if __name__=="__main__":    
    #
    prots,bornes,plddt,hca,sf_sizes,sf_starts=define_contiguous_segments()
    #
    fa=data_dir+"/downloads/sequences.fasta"
    charges,annot=get_info_from_fasta(fa)
    #
    ogn_prots=get_ogns()
    prots_ogns={}
    for o in ogn_prots.keys():
        for p in ogn_prots[o]:
            prots_ogns[p]=o
    #
    to_w=[]
    columns=["proteome","sequence id","annotation","contiguous region start","contiguous regions length (aa)","whole segment start","whole segment length (aa)","pLDDT category","HCA score","coverage by charged aa"]
    for i in range(len(prots)): # going through contiguous segments
        p=prots[i]
        length=bornes[i][1]-bornes[i][0]+1
        charge_cov=100*charges[p][bornes[i][0]:bornes[i][1]+1].count("1")/length
        to_w.append([ 
            prots_ogns[p],
            p,
            annot[p],
            bornes[i][0]+1,
            length,
            sf_starts[i]+1,
            sf_sizes[i],
            plddt[i],
            hca[i],
            100*charges[p][bornes[i][0]:bornes[i][1]+1].count("1")/length
            ])
    print(len(columns))
    print(len(to_w))
    print(len(to_w[0]))
    df=pd.DataFrame(to_w,columns=columns)
    df.style.set_precision(1)
    df.to_csv("../../../tables/datasetS4.csv",sep=";",index=False)
