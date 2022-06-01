cov=$1
id=$2
clumode=$3
output=$4
inputs=$5

echo $5

# Definition des aguments :
while getopts c:i:m:s:n:o:d: flag
do
  case "${flag}" in
    c) cov=${OPTARG};;
    i) id=${OPTARG};;
    m) clumode=${OPTARG};;
    s)sensi=${OPTARG};;
    n)nbseqs=${OPTARG};;
    o) output=${OPTARG};;
    d) indbs=${OPTARG};;
    \?)echo ${OPTARG};tmp=$infastas${OPTARG};infastas=$tmp;;
    # "") tmp=$indir${OPTARG};;
  esac
done

wd=$(pwd)
# infiles=$(echo $indir | tr ";" "\n")
mkdir mmseqs
cd mmseqs
# touch tmpseqs.txt
for db in $indbs; do
  echo $wd"/"$db
  sed  "s/\(>\)/&$db|/" ../temp/${db}_sequences.fasta >> tmpseqs.fasta
done
# #

cd mmseqs
mmseqs createdb tmpseqs.fasta seqDB
mkdir tmp
mmseqs cluster --min-seq-id $id --cluster-mode $clumode -c $cov -s $sensi --max-seqs $nbseqs seqDB cluDB tmp
mmseqs createtsv seqDB seqDB cluDB $wd"/"$output
cd ..
rm -r mmseqs

