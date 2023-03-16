curr_dir=$(pwd)

# Version of AFDB used in the article : v1 (https://ftp.ebi.ac.uk/pub/databases/alphafold/v1/)

# Download directory organization : 
# ls 
# UP000000803_7227_DROME.tar    UP000002296_353153_TRYCC.tar  UP000007305_4577_MAIZE.tar
# UP000000805_243232_METJA.tar  UP000002311_559292_YEAST.tar  UP000008153_5671_LEIIN.tar
# UP000000437_7955_DANRE.tar    UP000001450_36329_PLAF7.tar   UP000002485_284812_SCHPO.tar  UP000008816_93061_STAA8.tar
# UP000000559_237561_CANAL.tar  UP000001584_83332_MYCTU.tar   UP000002494_10116_RAT.tar     UP000008827_3847_SOYBN.tar
# UP000000589_10090_MOUSE.tar   UP000001940_6239_CAEEL.tar    UP000005640_9606_HUMAN.tar    UP000059680_39947_ORYSJ.tar
# UP000000625_83333_ECOLI.tar   UP000002195_44689_DICDI.tar   UP000006548_3702_ARATH.tar sequences.fasta



####################
### Extract data ###
####################
data_dir=$(../../../data/AFDB)
# 1 - get pLDDT from PDB files, creates 1 object by protein (by PDB file)
archs=($(ls -t $(data_dir)/downloads/*.tar))

mkdir $(data_dir)/pLDDT
for n in {0..20}; do
    arch=${archs[$n]}
    dir=${arch/.tar/}
    genome=${dir##*/}

    mkdir $dir
    tar -xf $arch -C $dir
    mkdir $(data_dir)/pLDDT/$genome

    for l in $(ls $dir/*.pdb.gz);do
        gunzip $l
        un=${l/.gz/}
        python3 get_pLDDT.py $un $genome $data_dir
        gzip $un
    done

    rm -r $dir
done

# 2 - merge informations collected for each protein during step 1 in 1 object
python3 merge_all_plddt.py $data_dir

# 3 - launch hcatk on sequences.fasta
gunzip $(data_dir)/downloads/sequences.fasta.gz
hcatk segment -i $(data_dir)/sequences.fasta -o $(data_dir)/hca.out -m domain

# 4 - Combine all to create dataset S4
python3 get_table.py $data_dir

