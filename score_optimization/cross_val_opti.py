#!/usr/bin/env python
from Bio import SeqIO
from pyHCA import HCA
import numpy as np
import os, sys, argparse, gzip
from sklearn.model_selection import KFold

NBINS=60#int(sys.argv[1])
arr=int(sys.argv[1])

import builtins
sys.stdout = open("stdouts/crossval_%i.txt"%arr, "w", buffering=1)
# sys.stdout = open("stdouts/compute_overlap.txt", "w", buffering=1)

def print(text):
    builtins.print(text)
    os.fsync(sys.stdout)

AA1 = ['A', 'C', 'D', 'E', 'F','G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

set_AA1 = set(AA1)

def compare_distrib(Z1, Z2,nbins=NBINS):
    """
    compare two distribution to find the minimum overlap
    """
    intervals = np.linspace(-10, 10, nbins)
    valZ1, _ = np.histogram(Z1, bins=intervals)
    valZ2, _ = np.histogram(Z2, bins=intervals)
    valZ1=np.true_divide(valZ1,np.sum(valZ1))
    valZ2=np.true_divide(valZ2,np.sum(valZ2))
    overlaps=np.sum(np.minimum(valZ1,valZ2))
    return overlaps

def prepare_sequence(seq):
    """
    compute score for a given sequence
    """
    size = len(seq)
    hcaprot = HCA(seq=seq)
    hcaclusters = hcaprot.get_clusters()
    clusters = list()
    residues_a, residues_b, residues_c = len(seq), 0, 0 # default, everything outside cluster
    for aa in seq:
        if aa not in set_AA1:
            aa = "X"
        
    for clust in hcaclusters: #Parcours des clusters hydrophobes
        if len(clust.hydro_cluster) > 2:
            for i, aa in enumerate(seq[clust.start: clust.stop]):
                if aa not in set_AA1:
                    aa = "X"
                    # print("Pourquoi?")
                if clust.hydro_cluster[i] == 1:
                    # in a cluster and is hydrophobe, increase b, decrease a
                    residues_b += 1
                    residues_a -= 1
                else:
                    # in a cluster but not hydrophobe, increase c, decrease a
                    residues_c +=  1
                    residues_a -= 1

    return (residues_a, residues_b, residues_c), size #

def compute_score(clusters, size, a, b, c): # a=poids aa hors HC, b=poids aa hydrophobes en HC, c=poids aa hydrophiles en HC
    """
    dcompute score for a given sequence
    """
    residues_a, residues_b, residues_c = clusters
    score = 0
    score += residues_a * a
    score += residues_b * b
    score += residues_c * c

    return score / size

def manual_opti(train1,train2):
    trial=0
    best_over=np.inf
    best_trial=0
    for a in range(-10,11):
        for b in range(-10,11):
            for c in range(-10,11):
                scores1 = [compute_score(clust, size, a, b, c) for (clust, size) in train1]
                # compute target
                scores2 = [compute_score(clust, size, a, b, c) for (clust, size) in train2]
                # evaluate
                overlap = compare_distrib(scores1, scores2)
                trial+=1
                if overlap<best_over:
                    best_over=overlap
                    best_trial=trial
                    best=[a,b,c]
                # print("Trial %s finished with value: %f and parameters: {'a': %i, 'b': %i, 'c': %i}. Best is trial %i with value: %f"%(trial,overlap,a,b,c,best_trial,best_over))
    print("%f overlap on train for : {'a': %i, 'b': %i, 'c': %i}"%(best_over,best[0],best[1],best[2]))
    return best


def kcross_val(data1,data2):
    k=10
    #data=data1+data2 Concaténer 2 df en ajoutant info de classe (desordre ou ordre)
    best_over=np.inf
    kfold = KFold(k,shuffle=True)
    index1=[]
    for d in data1:
        index1.append([i for i in kfold.split(d)])
    index2=[i for i in kfold.split(data2)]

    for i in range(k):
        train1=[data1[n][index1[n][i][0]] for n in range(len(data1))]
        train1=np.concatenate(train1)
        train2=data2[index2[i][0]]
        #
        test1=[data1[n][index1[n][i][1]] for n in range(len(data1))]
        test1=np.concatenate(test1)
        test2=data2[index2[i][1]]
        print("# KFold %i train1, %i train2, %i test1 %i test2"%(len(train1),len(train2),len(test1),len(test2)))
        a,b,c=manual_opti(train1,train2)
        
        scores1=[compute_score(clust, size, a, b, c) for (clust, size) in test1]
        scores2=[compute_score(clust, size, a, b, c) for (clust, size) in test2]
        overlap =  compare_distrib(scores1, scores2)
        print("# %f overlap on test for : {'a': %i, 'b': %i, 'c': %i}"%(overlap,a,b,c))
        if overlap<best_over:
            best_over=overlap
            best=[a,b,c]
    print("# Best value: {'a': %i, 'b': %i, 'c': %i} for overlap on testing data"%(best[0],best[1],best[2]))



##############
# Input Data #
##############


indir="opti_data/"
## read query
# print("Reading pdb")
path_to_scope = indir+"nrSCOPe.fasta"
with open(path_to_scope, "r") as inf:
    scopeseq = [str(record.seq) for record in SeqIO.parse(inf, "fasta")]
hcasco = [prepare_sequence(seq) for seq in scopeseq]

path_to_opm1 = indir+"nrOPM_1.fasta"
with open(path_to_opm1, "r") as inf:
    opm1seq = [str(record.seq) for record in SeqIO.parse(inf, "fasta")]
hcaopm1 = [prepare_sequence(seq) for seq in opm1seq]

path_to_opm11 = indir+"nrOPM_11.fasta"
with open(path_to_opm1, "r") as inf:
    opm11seq = [str(record.seq) for record in SeqIO.parse(inf, "fasta")]
hcaopm11 = [prepare_sequence(seq) for seq in opm11seq]

path_to_opm2 = indir+"nrOPM_2.fasta"
with open(path_to_opm1, "r") as inf:
    opm2seq = [str(record.seq) for record in SeqIO.parse(inf, "fasta")]
hcaopm2 = [prepare_sequence(seq) for seq in opm2seq]

## read targets
path_to_dis = indir+"nrDisProt_SF.fasta"
with open(path_to_dis, "r") as inf:
    disseq = [str(record.seq) for record in SeqIO.parse(inf, "fasta")]

hcadis = list()
for seq in disseq:
    clust, size = prepare_sequence(seq)
    if len(clust) > 0:
        hcadis.append((clust, size))

if __name__ == "__main__":
    kcross_val([np.array(hcasco),np.array(hcaopm1),np.array(hcaopm11),np.array(hcaopm2)],np.array(hcadis))

