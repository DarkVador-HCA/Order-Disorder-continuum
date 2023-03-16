curr_dir=$(pwd)

# Download latest AlphaFold Protein Structure Database
cd ../downloads/
#1 sequences
wget http://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta
#2 structures
mkdir structures
cd structures
wget ftp://ftp.ebi.ac.uk/pub/databases/alphafold/latest/*
cd $curr_dir

# Version of AFDB used in the article : v1 (https://ftp.ebi.ac.uk/pub/databases/alphafold/v1/)

# Download directory organization : 
# ls 
# UP000000803_7227_DROME.tar    UP000002296_353153_TRYCC.tar  UP000007305_4577_MAIZE.tar
# UP000000805_243232_METJA.tar  UP000002311_559292_YEAST.tar  UP000008153_5671_LEIIN.tar
# UP000000437_7955_DANRE.tar    UP000001450_36329_PLAF7.tar   UP000002485_284812_SCHPO.tar  UP000008816_93061_STAA8.tar
# UP000000559_237561_CANAL.tar  UP000001584_83332_MYCTU.tar   UP000002494_10116_RAT.tar     UP000008827_3847_SOYBN.tar
# UP000000589_10090_MOUSE.tar   UP000001940_6239_CAEEL.tar    UP000005640_9606_HUMAN.tar    UP000059680_39947_ORYSJ.tar
# UP000000625_83333_ECOLI.tar   UP000002195_44689_DICDI.tar   UP000006548_3702_ARATH.tar sequences.fasta

