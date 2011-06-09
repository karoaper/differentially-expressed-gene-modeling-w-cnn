IFS=","

for folder in $1
do
	for ((i=$2;i<$3;i++)); do ./gene_expression /Users/khovsep/Projects/GeneTrackAllign/Data/PSU\ data/${folder}/gene_expression_cut$i.conf > /Users/khovsep/Projects/GeneTrackAllign/Data/PSU\ data/${folder}/gene_expression_cut$i.out & done;
done
