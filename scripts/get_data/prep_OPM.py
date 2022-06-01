import os,sys,math
import pandas as pds
from pyHCA import HydroCluster,_transformSequence,_getAmas

"""Extraction of membrane domains using bounds and PDB files provided by OPM 
(sequences contain all residues in between the 1st the last that are in the membrane)"""
# os.chdir("../data/")
dAA = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

opmdir="../../data/OPM/"

##############################
### GET SEQUENCES FROM PDB ###
##############################

def parse_pdb(inpdb,bounds,subunits):
    """
    Get sequence from PDB file
    Arguments :
        inpdb = PDB file name (str)
        bounds = TM bounds (list[list[start,end]])
        subunits = subunit countaining TM domain (list[chain name])
    Return :
        seq = dict{"sequence name" : "amino acid sequence"}
    """
    seq={s:"" for s in subunits if s in bounds.keys()}
    go=False
    s0=""
    n0=-1
    rep_ch=False
    with open(inpdb,"r") as f:
        for l in f:
            if l[:6].strip()=='ATOM':
                n=int(l[22:26].strip())#resnum
                s=l[21]#chaine
                if s!=s0:
                    if s in seq.keys() and len(seq[s])!=0:
                        rep_ch=True
                    if go: # si on n'a pas pris en compte le dernier résidu de la chaine précédente => on ne prend pas en compte la chaine (séquence compromise)
                        del seq[s0]
                    go=False
                if s in seq.keys() and n!=n0:
                    if n==bounds[s][0]:
                        if len(seq[s])>0: # Si il y a plusieurs chaines du même nom => on ne prend aucune chaine en compte
                            del seq[s]
                        else: #sinon : c'est bon on peut commencer à prendre les résidus qui suivent en compte
                            go=True
                    if go: #si on a rencontré le premier résidu tm de la chaine, mais pas encore le dernier ==> on prend ajoute les résidus aux séquences
                        r=dAA[l[17:20].strip()] #r = code de l'aa à 1 lettre
                        
                        seq[s]+=r #ajout du résidu à la chaine
                    if n==bounds[s][1]: #si on en est au dernier résidu
                        go=False #on ne prend pas les résidus qui suivent en compte
                n0=n
                s0=s
    return seq

def get_seq_bounds():
    pdb_dirs=opmdir+"downloads/structures/"
    pdb_sub=get_pdb_bounds("../data/PDBbounds_per_sub.tsv",pdb_sub)
    #
    bounds={} #sequence bounds
    pbL=[]
    cpt=0
    for dir in pdb_dirs:
        for pdb_file in os.listdir(dir):
            pdb=pdb_file.split(".")[0]
            if pdb in pdb_sub.keys():
                sub_bounds=read_pdb(dir+pdb_file,pdb_sub[pdb])
                for s in sub_bounds.keys():

                    if len(sub_bounds[s])!=len(pdb_sub[pdb][s]):
                        cpt+=1
                        if pdb not in pbL:
                            pbL.append(pdb)
                bounds[pdb]=sub_bounds
    return bounds


def write_all_sequences():
    """
    From OPM tables (proteins and subunits) and PDB files, 
    extract sequences of transmembrane domains
    """
    #get list of pdbs to keep (only type 1 structures : transmembrane)
    df=pd.read_csv(opmdir+"downloads/proteins.csv")
    df["pdbid"]=df["pdbid"].map(lambda x: x.lstrip("=").strip("\"")).values
    dfkeep=df[["pdbid","classtype_id"]][df["type_id"]==1]
    tokeep=df["pdbid"][df["type_id"]==1]
    #get bounds on the selected structures
    df=pd.read_csv(opmdir+"downloads/subunits.csv")
    df["pdbid"]=df["pdbid"].map(lambda x: x.lstrip("=").strip("\"")).values
    segs=df[["pdbid","protein_letter","segment"]][df["pdbid"].isin(tokeep)]
    classes={1:'polytopic',2:'beta',11:'bitopic'}

    cpt_nosub=0#count of prots without subunits
    cpt_neg=0#count of structures with negative residues
    # cpt_unk=0#count of structures with unknown residues 
    cpt_seq=0

    all_seqs={} # aa sequences dict{"sequence name (pdb_chain)":"aa sequence"}
    all_bounds={}
    for pdb in tokeep: 
        c=dfkeep["classtype_id"][dfkeep["pdbid"]==pdb].values[0] # class (polytopic, bitopic or beta)
        dfpdb=segs[["protein_letter","segment"]][segs["pdbid"]==pdb] 
        subunits=dfpdb["protein_letter"].values
        if len(subunits)==0:
            cpt_nosub+=1
        #same as :
        # if dfpdb.shape==(0,2):
        #     cpt_nosub+=1
        bounds={}
        for s in subunits: # Get bounds for each chain
            l=dfpdb["segment"][dfpdb["protein_letter"]==s].values[0].split(",")
            if len(l[0].split("-"))>2:
                cpt_neg+=1
            else:
                bounds[s]=[[int(b.split("-")[0].split("(")[1].strip()),int(b.split("-")[1].split(")")[0].strip().replace('\r\n',''))] for b in l]
        max_bounds={s:[bounds[s][0][0],bounds[s][-1][1]] for s in subunits if s in bounds.keys()}
        if len(bounds)>0:
            all_bounds[pdb]=bounds
            seqs=parse_pdb("../data/"+classes[c]+"/"+pdb+".pdb",max_bounds,subunits)
            if len(seqs)>0:
                cpt_seq+=1
                all_seqs[pdb]=seqs

    with open(opmdir+"id_cats.csv","w") as f:
        for c in [1,2,11]:
            pdbs=dfkeep["pdbid"][dfkeep["classtype_id"]==c].values
            print(len(pdbs))
            for p in pdbs:
                if p in all_seqs.keys():
                    f.write("%s,%s\n"%(c,p))

    with open(opmdir+"downloads/PDB_bounds_per_sub.tsv","w") as f:
        f.write("#pdbid\tsubunit\tbounds\n")
        for p in all_bounds.keys():
            for s in all_bounds[p].keys():
                f.write("%s\t%s\t"%(p,s))
                sep=""
                for b in all_bounds[p][s]:
                    f.write("%s%i:%i"%(sep,b[0],b[1]))
                    if sep=="":
                        sep=","
                f.write("\n")
    
    #
    bounds=get_seq_bounds()
    cpt=0
    with open("../data/SEQ_bounds_per_sub.tsv","w") as f:
        f.write("#pdbid\tsubunit\tbounds\n")
        for pdb in bounds.keys():
            for sub in bounds[pdb].keys():
                if len(bounds[pdb][sub])==len(pdb_sub[pdb][sub]):
                    f.write("%s\t%s"%(pdb,sub))
                    s="\t"
                    for i in range(0,len(bounds[pdb][sub]),2):
                        f.write("%s%i:%i"%(s,bounds[pdb][sub][i],bounds[pdb][sub][i+1]))
                        if s=="\t":
                            s=","
                    f.write("\n")
                else:
                    cpt+=1
    # print("\n")
    # print("Prots absent from subunit file :",cpt_nosub)
    # print("Prots with res <0 :",cpt_neg)
    # print("Prots with unknown residues :",cpt_unk)
    # print("structures of type 1 :",len(tokeep), "structures dont on récupère la séquence :",len(tokeep)-cpt_nosub-cpt_neg-cpt_unk)
    # print("structures for which we get a sequence :",cpt_seq)

    with open(opmdir+'downloads/all_COMPLETEsequences.fasta','w') as f: #Save full sequences (before removing the soluble domains)
        for s in all_seqs.keys():
            for c in all_seqs[s].keys():
                if len(all_seqs[s][c])>0:
                    f.write(">%s_%s\n%s\n"%(s,c,all_seqs[s][c]))



##############################
### REMOVE SOLUBLE DOMAINS ###
##############################


def check_interTM_length(all_seq,bounds,length=30):
    seq=""
    new_bounds=[]
    diff=0
    marge=math.floor(length/2)
    amas=_getAmas(_transformSequence(all_seq))[1]
    lengths=[]
    for i in range(len(bounds)-1):
        tm1_start,tm1_end=bounds[i][0],bounds[i][1]
        tm2_start=bounds[i+1][0]
        seq+=all_seq[tm1_start:tm1_end+1] # segment TM
        new_bounds.append([bounds[i][0]-diff,bounds[i][1]-diff]) # bornes des semgments tm quand on soustrait parts solubles éliminées
        l=(tm2_start-1)-(tm1_end+1)+1#(bounds[i+1][0]-1)-(bounds[i][1]+1)+1
        lengths.append(l)
        if l<=length: # si le segment soluble fait moins que la taille limite (=> on considère que c'est une boucle et non un domaine)
            seq+=all_seq[tm1_end+1:tm2_start] # interTM
        else:
            j=-1
            if amas[tm1_end]==1 and amas[tm1_end+1]==1: # si la fin du sgment TM est à cheval sur un amas hydrophobe:
                j=1

                while amas[tm1_end+j]==1 and (tm1_end+j)<tm2_start:
                    j+=1
            kept_aa1=max(j+4,marge) # nb d'aa gardés après le premier passage TM
            j=-1
            if amas[tm2_start]==1 and amas[tm2_start-1]==1: # si la fin du segment TM est à cheval sur un amas hydrophobe:
                j=1
                while amas[tm2_start-j]==1 and (tm2_start-j)>tm1_end:
                    j+=1
            kept_aa2=max(j+4,marge)
            kept=kept_aa1+kept_aa2 # nb d'aa gardés avant le premier passage TM
            diff+=l-kept
            seq+=all_seq[tm1_end+1:tm1_end+kept_aa1+1] # aa gardés après le premier passage TM
            seq+=all_seq[tm2_start-kept_aa2:tm2_start] # aa gardés avant le deuxième passage TM
            # diff+=l-8 # taille du segment qu'on enlève - les aa qu'on laisse pour ne pas se retrouver avec 2 passages TM collés
            # seq+=all_seq[bounds[i][1]:bounds[i][1]+4] # aa gardés après le premier passage TM
            # seq+=all_seq[bounds[i+1][0]-4:bounds[i+1][0]-1] # aa gardés avant le deuxième passage TM

    seq+=all_seq[bounds[-1][0]:]
    new_bounds.append([bounds[-1][0]-diff,bounds[-1][1]-diff])
    return seq,new_bounds,lengths,diff


def read_fasta(fasta):
    seqs={}
    with open(fasta,'r') as f:
        for l in f:
            if l[0]==">":
                n=l[1:].strip()
                seqs[n]=""
            else:
                seqs[n]+=l.strip()
    return seqs

def read_bounds(bounds_file):
    bounds={}
    with open(bounds_file,'r') as f:
        for l in f:
            if l[0]!="#":
                l=l.strip().split("\t")
                n=l[0]+"_"+l[1]
                bounds[n]=[]
                for b in l[2].split(','):
                    b=b.split(':')
                    bounds[n].append([int(b[0]),int(b[1])])
    return bounds

def remove_long_soluble_parts():
    seqs=read_fasta(opmdir+"downloads/COMPLETEsequences.fasta")
    TMbounds=read_bounds(opmdir+"downloads/SEQ_bounds_per_sub.tsv")
    bounds={}
    n_seqs={}
    f=open("interTM_L"+str(max_interTM_length)+".csv","w")
    f.close()

    for n in seqs.keys():
        if n in TMbounds.keys():
            n_seqs[n],bounds[n],interTM_lengths,nb_rm_aa=check_interTM_length(seqs[n],TMbounds[n],length=max_interTM_length)
            with open("interTM_L"+str(max_interTM_length)+".csv","a") as f:
                for i in interTM_lengths:
                    f.write("%s,%i\n"%(n,i))
            with open("nb_rm_aa"+str(max_interTM_length)+".csv","a") as f:
                f.write("%s,%i\n"%(n,nb_rm_aa))

    with open('../sequences.fasta','w') as f:
        for s in n_seqs.keys():
            if len(n_seqs[s])>0:
                f.write(">%s\n%s\n"%(s,n_seqs[s]))


    if __name__=='__main__':
        write_all_sequences()
        remove_long_soluble_parts():