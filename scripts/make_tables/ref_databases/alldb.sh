touch TMPallSeq_resume.csv
touch TMPfoldSeg_resume.csv

for db in SCOPe DisProt OPM;do
    echo $db
    mkdir -p ../tables/$db ../tables/$db/
    python make_datasets_S2_S3.py ../data/$db/nr_sequences.fasta ../data/$db/hca.out $db
    tail -n+2 ../tables/$db/tables_S2.csv >> TMP_S2.csv
    tail -n+2 ../tables/$db/tables_S3.csv >> TMP_S3.csv
done


head -n1 ../tables/$db/allSeq_resume.csv > ../tables/allSeq_resume.csv
cat TMP_S2.csv >> ../results/datasetS2.csv

head -n1 ../tables/$db/tables/foldSeg_resume.csv > ../tables/foldSeg_resume.csv
cat TMP_S3.csv >> ../results/datasetS4.csv

rm TMP*