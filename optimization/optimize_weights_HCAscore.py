#!/usr/bin/env python3
from Bio import SeqIO
from pyHCA import HCA
import numpy as np
import os, sys, argparse, gzip
import optuna

b=60

import builtins
sys.stdout = open("optimization.txt", "w", buffering=1)
def print(text):
    builtins.print(text)
    os.fsync(sys.stdout)

AA1 = ['A', 'C', 'D', 'E', 'F','G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
set_AA1 = set(AA1)

def compare_distrib(Z1, Z2):
    """
    compare two distribution to find the minimum overlap
    """
    intervals = np.linspace(-10, 10, b)
    valZ1, _ = np.histogram(Z1, bins=intervals)
    valZ2, _ = np.histogram(Z2, bins=intervals)
    valZ1=np.true_divide(valZ1,np.sum(valZ1))
    valZ2=np.true_divide(valZ2,np.sum(valZ2))
    overlaps=np.sum(np.minimum(valZ1,valZ2))
    return overlaps

def prepare_sequence(seq):
    """
    identify hydrophobic clusters in the sequence
    """
    size = len(seq)
    hcaprot = HCA(seq=seq)
    hcaclusters = hcaprot.get_clusters()
    clusters = list()
    residues_a, residues_b, residues_c = 0, 0, 0
    for aa in seq:
        if aa not in set_AA1:
            aa = "X"
        residues_a += 1 # outside of HC
    for clust in hcaclusters: #Parcours des clusters hydrophobes
        if len(clust.hydro_cluster) > 2:
            for i, aa in enumerate(seq[clust.start: clust.stop]):
                if aa not in set_AA1:
                    aa = "X"
                if clust.hydro_cluster[i] == 1:
                    # in a cluster and is hydrophobe, increase b, decrease a
                    residues_b += 1
                    residues_a -= 1
                else:
                    # in a cluster and but not hydrophobe, increase c, decrease a
                    residues_c +=  1
                    residues_a -= 1

    return (residues_a, residues_b, residues_c), size #

def compute_score(clusters, size, a, b, c): # a=weight aa outside HC, b=weight hydrophobic aa inside HC, c=weight hydrophilic aa inside HC
    """
    compute score for a given sequence
    """
    residues_a, residues_b, residues_c = clusters
    score = 0
    score += residues_a * a
    score += residues_b * b
    score += residues_c * c

    return score / size


def objective(trial):
    residues_a  = trial.suggest_discrete_uniform("a", -10, 10, 1) # poids résidus hors HC
    residues_b = trial.suggest_discrete_uniform("b", -10, 10, 1) # poids résidus hydrophobes en HC
    residues_c = trial.suggest_discrete_uniform("c", -10, 10, 1) # poids résidus hydrophiles en HC

    pdb_scores = [compute_score(clust, size, residues_a, residues_b, residues_c) for (clust, size) in hcapdb]
    # compute target
    dis_scores = [compute_score(clust, size, residues_a, residues_b, residues_c) for (clust, size) in hcadis]
    # evaluate
    overlap = compare_distrib(pdb_scores, dis_scores)
    return overlap

def optimize_study(study_name, storage, objective, n_trials):
    study = optuna.create_study(study_name=study_name, storage=storage, load_if_exists=True)
    study.optimize(objective, n_trials=n_trials)


# Input Data
indir="data/"
## read query
print("Reading pdb")
path_to_pdb = indir+"PDB/pdb.fasta"
with open(path_to_pdb, "r") as inf:
    pdbseq = [str(record.seq) for record in SeqIO.parse(inf, "fasta")]
hcapdb = [prepare_sequence(seq) for seq in pdbseq]

## read targets
print("Reading disprot")
path_to_dis = indir+"DisProtv7/disprot.fasta"
with open(path_to_dis, "r") as inf:
    disseq = [str(record.seq) for record in SeqIO.parse(inf, "fasta")]

hcadis = list()
for seq in disseq:
    clust, size = prepare_sequence(seq)
    if len(clust) > 0:
        hcadis.append((clust, size))

def main():
    trial=0
    best_over=np.inf
    best_trial=0
    for a in range(-10,11):
        for b in range(-10,11):
            for c in range(-10,11):
                pdb_scores = [compute_score(clust, size, a, b, c) for (clust, size) in hcapdb]
                # compute target
                dis_scores = [compute_score(clust, size, a, b, c) for (clust, size) in hcadis]
                # evaluate
                overlap = compare_distrib(pdb_scores, dis_scores)
                trial+=1
                if overlap<best_over:
                    best_over=overlap
                    best_trial=trial
                    best=[a,b,c]
                print("Trial %s finished with value: %f and parameters: {'a': %i, 'b': %i, 'c': %i}. Best is trial %i with value: %f"%(trial,overlap,a,b,c,best_trial,best_over))
    print("Best value: {'a': %i, 'b': %i, 'c': %i}"%(best[0],best[1],best[2]))


if __name__ == "__main__":
    main()