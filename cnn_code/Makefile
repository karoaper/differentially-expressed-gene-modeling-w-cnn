all: weights genes

weights:
	g++ -I/usr/include/eblearn -L/usr/lib extract_weights.cpp -lm -leblearntools -leblearn -lidx -o extract_weights

genes:
	
	g++ -I/usr/include/eblearn -L/usr/lib train.cpp -lm -leblearntools -leblearn -lidx -o gene_expression