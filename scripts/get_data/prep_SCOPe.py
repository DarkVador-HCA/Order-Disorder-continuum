Sdir="../../data/SCOPe/"

with open(Sdir+"downloads/astral_95.fasta","r") as f:
    cats_id={"a":[],"b":[],"c":[],"d":[]} #categories we want to keep
    seqs={}
    want=True
    for l in f:
        if l[0]==">":
            want=False
            l=l[1:].split()
            nom=l[0]
            cat=l[1].split(".")[0]
            if cat in ["a","b","c","d"]:
                want=True
                cats_id[cat].append(nom)
                seqs[nom]=""
        elif want:
            for a in l.strip():
                if a!='X':
                    seqs[nom]+=a.upper()

#Saving each sequence category
with open(Sdir+"id_cats.csv","w") as f:
    for cat in cats_id.keys():
        for p in cats_id[cat]:
            f.write("%s,%s\n"%(cat,p))

#Saving sequences in right format
with open(Sdir+"downloads/sequences.fasta","w") as f:
    for p in seqs.keys():
        f.write(">SCOPe|%s\n%s\n"%(p,seqs[p]))
