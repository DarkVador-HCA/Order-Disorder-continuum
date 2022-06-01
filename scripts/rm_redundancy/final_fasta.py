#!/usr/bin/env python3
"""
Get fasta of non redundant sequences after MMSeqs2 clustering
"""

cpt=0
to_rm={db:[] for db in snakemake.config["DB"]}
with open(snakemake.input.clus,"r") as f:
    for l in f:
        l=l.strip().split()
        db1=l[0].split("|")[0]
        db2=l[1].split("|")[0]
        if db1!=db2:
            # print(l)
            p1="|".join(l[0].split("|")[1:])
            p2="|".join(l[1].split("|")[1:])
            if p1 not in to_rm[db1]:
                to_rm[db1].append(p1)
            if p2 not in to_rm[db2]:
                to_rm[db2].append(p2)
            cpt+=1

print("\n",cpt)
print("\n",len(to_rm))

for db in to_rm.keys():
    with open("temp/"+db+"_sequences.fasta","r") as f:
        want=False
        to_w=[]
        for l in f:
            if l[0]==">":
                want=False
                if l[1:].strip().split()[0] not in to_rm[db]:
                    want=True
                    to_w.append(l)
            elif want:
                to_w.append(l)

    with open(snakemake.config["DBdir"]+db+"/nr_sequences.fasta","w") as f:
        f.write("".join(to_w))
    with open(snakemake.config["DBdir"]+db+"/downloads/in_other_DB.txt","w") as f:
        f.write("\n".join(to_rm[db]))
