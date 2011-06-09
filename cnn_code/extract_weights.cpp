/* ==========================
 Code for training a convolutional neural network for regression, used for the DREAM protein binding affinity prediction challenge
 */

#ifndef __WINDOWS__
#include <fenv.h>
#endif


#define Malloc(type,n) (type *) malloc(sizeof(type)*n)

#include "libeblearn.h"
#include "libeblearntools.h"

typedef double t_net;

using namespace std;
using namespace ebl; // all eblearn objects are under the ebl namespace

int main(int argc, char **argv) 
{
	// regular main without gui
	cout << "Extracting DREAM 5 weights for conf file " << argv[1];
	configuration conf(argv[1]); // configuration file
	char * weightsfn = argv[2];
	
	
	
	
	int tracklen;// = conf.get_uint("track_length");	
	int numtracks;// = conf.get_uint("num_tracks");
		
	
	char temp [20000];
	//	fprintf(stderr,"ran until here\n");
	sprintf(temp,"%s/%s",conf.get_cstring("root"),conf.get_cstring("train_name"));
	FILE *trainingdata_file = fopen(temp,"r");
	if (trainingdata_file==NULL) perror ("Error opening file");
	
	
	
	char firstchar;
	int paramval;
	char paramname [100];
	while ((firstchar = getc(trainingdata_file)) == '#')
	{
		fscanf(trainingdata_file," %s = ",paramname);
		
		if(strcmp(paramname,"numtracks") == 0)
		{
			fscanf(trainingdata_file,"%d\n",&paramval);
			if(numtracks == 0)
				numtracks = paramval;
		}
		else if(strcmp(paramname,"featuresize") == 0)
		{
			fscanf(trainingdata_file,"%d\n",&paramval);
			if(tracklen == 0)
				tracklen = paramval;
		}
		else
		{
			fscanf(trainingdata_file,"%*[^\n]\n");
		}
	}
	
	fprintf(stderr,"%d %d\n",numtracks,tracklen);
	
	ungetc(firstchar,trainingdata_file);
	
	init_drand(time(NULL)); // initialize random seed
	
	int ks1 = conf.get_uint("size_kernel_l1");
	int ks2 = conf.get_uint("size_kernel_l2");
	int ss1 = conf.get_uint("subsample_size_l1");
	int ss2 = conf.get_uint("subsample_size_l2");
	int numhid = conf.get_uint("num_hidden");
	
	
	sprintf(temp,"%s",conf.get_cstring("layer_connections_file"));	
	FILE *connections_file = fopen(temp,"r");
	if (connections_file==NULL) perror ("Error opening file");
	
	fprintf(stderr,"Reading the connections file\n");
	int nf1,nf2;
	
	fscanf(connections_file,"====\n");
	fscanf(connections_file,"%d\n",&nf1);
	int * table0_rows = new int[nf1]; 
	int totalrows = 0;
	for(int i=0;i<nf1;i++)
	{
		fscanf(connections_file,"%d %*[^\n]\n",&table0_rows[i]);
		totalrows += table0_rows[i];
	}
	//	for(int i=0;i<nf1;i++)
	//		fprintf(stderr,"%d ",table0_rows[i]);
	
	idx<intg> table0 = idx<intg>(totalrows,2);
	fscanf(connections_file,"====\n");
	fscanf(connections_file,"%d\n",&nf2);
	int * table1_rows = new int[nf2]; 
	totalrows = 0;
	for(int i=0;i<nf2;i++)
	{
		fscanf(connections_file,"%d %*[^\n]\n",&table1_rows[i]);
		totalrows += table1_rows[i];
	}
	//	for(int i=0;i<nf2;i++)
	//		fprintf(stderr,"%d ",table1_rows[i]);
	
	idx<intg> table1 = idx<intg>(totalrows,2);
	
	rewind(connections_file);
	
	int ** out_in_map_0 = Malloc(int *,nf1);
	int ** out_in_map_1 = Malloc(int *,nf2);
	
	
	fscanf(connections_file,"====\n");
	fscanf(connections_file,"%*d\n");
	int tablerowind = 0;
	for(int i=0;i<nf1;i++)
	{
		fscanf(connections_file,"%*d");
		out_in_map_0[i] = Malloc(int,table0_rows[i]);

		int temp;
		for(int j=0;j<table0_rows[i];j++)
		{
			fscanf(connections_file," %d",&temp);
			table0.set(temp,tablerowind,0);
			table0.set(i,tablerowind++,1);
			out_in_map_0[i][j] = temp;
		}
		fscanf(connections_file,"\n");
	}
	
	tablerowind = 0;
	fscanf(connections_file,"====\n");
	fscanf(connections_file,"%*d\n");
	for(int i=0;i<nf2;i++)
	{
		fscanf(connections_file,"%*d");
		out_in_map_1[i] = Malloc(int,table1_rows[i]);
		
		int temp;
		for(int j=0;j<table1_rows[i];j++)
		{
			fscanf(connections_file," %d",&temp);
			table1.set(temp,tablerowind,0);
			table1.set(i,tablerowind++,1);

			out_in_map_1[i][j] = temp;
		}
		fscanf(connections_file,"\n");
	}
	
	fclose(connections_file);
	
	
//	fprintf(stderr,"Constructing the gene_expression network\n");
	
//	parameter<t_net> theparam(60000); // create trainable parameter
	
	
//	idx<ubyte> trainingclasses(conf.get_uint("training_size"));
//	idx<ubyte> testingclasses(conf.get_uint("testing_size"));
	
	
//	sprintf(temp,"%s_weights.dat",argv[1]);
	
	
	//	cout << temp;
	idx<t_net> weights = load_matrix<t_net>(weightsfn);
	weights.pretty();
	double *weights_vals = weights.idx_ptr();

	
	sprintf(temp,"%s_layer1_weights.dat",weightsfn);
	fprintf(stderr,"%s\n",temp);
	FILE * outputfile = fopen(temp,"w");
	fprintf(stderr,"%s\n",temp);
	
	int ind = 0;
	//	printf("\n");
	for(int i=0;i<nf1;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		fprintf(stderr,"%d\n",ind);
	
		for(int j=0;j<table0_rows[i];j++)
		{
			fprintf(stderr,"%d\n",ind);
			fprintf(outputfile,"%d:",out_in_map_0[i][j]);
			for(int k=0;k<ks1;k++)
				fprintf(outputfile," %lf",weights_vals[ind++]);
			fprintf(outputfile,"\n");
		}
		//		ind++;   NO BIAS TERM IN THIS MIX
	}
	fprintf(stderr,"%d\n",ind);
	
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<nf1;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);

	fclose(outputfile);
	fprintf(stderr,"%d\n",ind);
	
	sprintf(temp,"%s_layer1_subsamplingweights.dat",weightsfn);
	outputfile = fopen(temp,"w");
//	printf("nf1 = %d\n",nf1);
//	printf("%lf %lf %lf\n",weights_vals[ind],weights_vals[ind-1],weights_vals[ind+1]);
	for(int i=0;i<nf1;i++)
	{
		fprintf(outputfile,"%d: %lf %lf\n",i,weights_vals[ind],weights_vals[nf1+ind]);
		ind++;
	}
	
	ind += nf1;
	
	fprintf(stderr,"%d\n",ind);
	fclose(outputfile);

	sprintf(temp,"%s_layer2_weights.dat",weightsfn);
	outputfile = fopen(temp,"w");

	for(int i=0;i<nf2;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		
		for(int j=0;j<table1_rows[i];j++)
		{
			fprintf(outputfile,"%d:",out_in_map_1[i][j]);
			for(int k=0;k<ks2;k++)
				fprintf(outputfile," %lf",weights_vals[ind++]);
			fprintf(outputfile,"\n");
		}
		//		ind++;   NO BIAS TERM IN THIS MIX
	}
	fprintf(stderr,"%d\n",ind);
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<nf2;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	
	fclose(outputfile);
	fprintf(stderr,"%d\n",ind);
	
	sprintf(temp,"%s_layer2_subsamplingweights.dat",weightsfn);
	outputfile = fopen(temp,"w");
	for(int i=0;i<nf2;i++)
	{
		fprintf(outputfile,"%d: %lf %lf\n",i,weights_vals[ind],weights_vals[nf2+ind]);
		ind++;
	}
	
	ind += nf2;
	fclose(outputfile);
	fprintf(stderr,"%d\n",ind);

	
	sprintf(temp,"%s_layer3_weights.dat",weightsfn);
	outputfile = fopen(temp,"w");
	
	int ks3 = ((tracklen-ks1+1)/ss1 - ks2 + 1)/ss2;
	
	for(int i=0;i<numhid;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		
		for(int j=0;j<nf2;j++)
		{
			fprintf(outputfile,"%d:",j);
			for(int k=0;k<ks3;k++)
				fprintf(outputfile," %lf",weights_vals[ind++]);
			fprintf(outputfile,"\n");
		}
		//		ind++;   NO BIAS TERM IN THIS MIX
	}
	fprintf(stderr,"%d\n",ind);
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<numhid;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	
	fclose(outputfile);
	fprintf(stderr,"%d\n",ind);
				
	sprintf(temp,"%s_layer4_weights.dat",weightsfn);
	outputfile = fopen(temp,"w");
	
	for(int i=0;i<2;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		
		for(int j=0;j<numhid;j++)
			fprintf(outputfile," %lf",weights_vals[ind++]);
		fprintf(outputfile,"\n");
		//		ind++;   NO BIAS TERM IN THIS MIX
	}
	fprintf(stderr,"%d\n",ind);
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<2;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	fprintf(stderr,"%d\n",ind);
	
	fclose(outputfile);

	
	//weights.pretty();
	
	return 0;
}
