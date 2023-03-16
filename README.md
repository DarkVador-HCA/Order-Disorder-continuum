# Code and data from the manuscript

#### **A sequence-based foldability score combined with AlphaFold2 predictions to disentangle the protein order/disorder continuum** by *Apolline Bruley, Tristan Bitard-Feildel ,  Isabelle Callebaut, Elodie Duprat*
**Manuscript available *[here](http://doi.org/10.1002/prot.26441)*.**

# HCA score implementation
Data and scripts used for the HCA score implamentation can be found in **optimization**.

# Classification of order and disorder from advanced sequence databases 
The dvanced sequence databases are SCOPe v 2.0.7, DisProt v 8.0.2, OPM (date : August 30, 2021).
The data used from these databases can be found in **data**, and under a folder named after each database (**SCOPe**, **DisProt** and **OPM**).
### Downloads
Within that folder, **downloads** contains data directly downloaded from the databases. To update with the last available version of the database : 
```
cd scripts/get_data/
./download_all.sh
```
DisProt can only be downloaded manually from here : https://www.disprot.org/download)
SCOPe version should be changed directly in the script.  
### Pre-treatment of sequences
In the **downloads** repository, the file **sequence.fasta** contains all sequences extracted from each database, if you have updated your data, you can obtain this file by running : 
````
cd scripts/get_data/
./prep_all.sh
````
At this stage, OPM sequences will be scanned to identify and remove long soluble segments in between mambrane passages (all is encoded in **scripts/get_data/prep_OPM.py**).
### Redundancy filters
The file **nr_sequences.fasta** contains non redundant sequences from each dataset. To apply these filters : 
````
cd scripts/rm_redundancy/
snakemake
````
### Foldable segments and HCA score computation
**hca.out** contains the delineation of foldable segments in each non-redundant dataset and the HCA score of the sequence, and of each foldable segment. To obtain the result we used the pyHCA package : https://github.com/DarkVador-HCA/pyHCA and runned the following:
```
cd data/$DB/ # DB can take the value SCOPe, OPM or DisProt
hcatk segment -m domain -i nr_sequences.fasta -o hca.out
```
The scripts allowing the generation of the datasets S2 and S3 can be found in scripts/make_tables/ref_databases/.

The HCA score distributions for the reference database (figure 1) are represented on the notebook plot_HCAscore.ipynb. The detail behind this representation can be found in scripts/make_plots/hcascore_distributions.py.

# Leveraging AlphaFold2 predictions with pyHCA package

The scripts allowing the generation of the dataset S4 from the AlphaFold Protein Database (from the sequences.fasta and PDB files) can be found in scripts/make_tables/AF2_databases/.

The barplots (figure 5) are represented in the notebook plot_HCAscore.ipynb. The detail behind this representation can be found in scripts/make_plots/AF2_database/mk_barplot.py .
