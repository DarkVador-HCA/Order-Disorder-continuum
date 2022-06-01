
# Definition des aguments :
while getopts c:i:m:s:n:o:f:a:d: flag
do
  case "${flag}" in
    c) cov=${OPTARG};;
    i) id=${OPTARG};;
    m) clumode=${OPTARG};;
    s)sensi=${OPTARG};;
    n)nbseqs=${OPTARG};;
    o) output=${OPTARG};;
    f) infasta=${OPTARG};;
    a) clus=${OPTARG};;
    d) wd=${OPTARG};;
  esac
done

echo $clus
# wd=mmseqs
  mkdir $wd
  #
   mmseqs createdb $infasta $wd/seqDB
  #
  mkdir $wd/tmp
   mmseqs cluster --min-seq-id $id --cluster-mode $clumode -c $cov -s $sensi --max-seqs $nbseqs $wd/seqDB $wd/cluDB $wd/tmp
   mmseqs result2repseq $wd/seqDB $wd/cluDB $wd/cluDB_rep
   mmseqs result2flat $wd/seqDB $wd/seqDB $wd/cluDB_rep $output --use-fasta-header
   mmseqs createtsv $wd/seqDB $wd/seqDB $wd/cluDB $clus
  # # cd ..
