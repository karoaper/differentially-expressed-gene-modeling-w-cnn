/* ==========================
 Code for training a convolutional neural network for regression, used for the DREAM protein binding affinity prediction challenge
 */


#ifndef __WINDOWS__
#include <fenv.h>
#endif

#include "libeblearn.h"
#include "libeblearntools.h"

#ifdef __GUI__
#include "libeblearngui.h"
#endif
//#include <iostream>
//#include <fstream>

#define UBYTE_INPUT		0
#define DOUBLE_INPUT	1
#define DEFAULT_MAX_NUMPARAMS	200000

#define Malloc(type,n) (type *) malloc(sizeof(type)*n)

#define convert_classes(n)	(n == 1 ? 0 : 1)

typedef double t_net;
void store_weights(char *fn,double *weights_vals,int numlayers, int nf1,int nf2, int nf3, int numout, int ks1, int ks2, int ks3, int **out_in_map_0, int **out_in_map_1, int * table0_rows, int * table1_rows);
void store_weights_use_one_conv_layer(char *fn,double *weights_vals,int numlayers, int nf1,int nf2, int numout, int ks1, int ks2, int **out_in_map_0, int * table0_rows);
void plot_results(float *correct,int num,int total,float maxtrain, int maxtrainind, float maxtest, int maxtestind,const char *fn);



using namespace std;
using namespace ebl; // all eblearn objects are under the ebl namespace

template <class T> void retrieve_network_states(supervised_trainer<t_net, T, ubyte> &thetrainer, labeled_datasource<t_net, T, ubyte> &ds, infer_param &infp, char * fn);
void store_states(bbstate_idx<t_net> &data,FILE * fin, const char *layername);

int main(int argc, char **argv)
{
	cout << "* GENE EXPRESSION: learning to classify gene expression direction (UP vs DOWN) using eblearn C++ library *" << endl;
	configuration conf(argv[1]); // configuration file

	int tracklen = (conf.exists("track_length") ? conf.get_uint("track_length"): 0);  //conf.get_uint("track_length");	
	int numtracks = (conf.exists("num_tracks") ? conf.get_uint("num_tracks"): 0);  //conf.get_uint("num_tracks");
	
	int typeofinput = (conf.exists("input_type") ? conf.get_uint("input_type"): UBYTE_INPUT);
	
	
	
	
	char temp [20000];
	char *results_plot_fn;
	
	
	bool storeweights = conf.exists_true("storeweights");
	
	if(storeweights)
		printf("Storing weights");
	
	bool plotresults = conf.exists_true("plotresults");
	bool store_states = conf.exists_true("store_states");
	if(plotresults)
	{
		if(conf.exists("results_plot_fn"))
			sprintf(results_plot_fn,"%s", conf.get_cstring("results_plot_fn"));
//			results_plot_fn = conf.get_cstring("results_plot_fn");
		else 
		{
			results_plot_fn = new char[20000];
			sprintf(results_plot_fn,"%s_resultsplot.png", argv[1]);
		}

		printf("Plotting incremental accuracy curves to %s\n", results_plot_fn);
	}
	
	
	
	sprintf(temp,"%s/%s",conf.get_cstring("root"),conf.get_cstring("train_name"));
	FILE *trainingdata_file = fopen(temp,"r");
	if (trainingdata_file==NULL) perror ("Error opening training file");
	
	
	char firstchar;
	int paramval;
	char paramname [200];
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
	
	
	bool use_one_conv_layer = conf.exists_true("use_one_conv_layer");
	int ks1 = conf.get_uint("size_kernel_l1");
	int ss1 = conf.get_uint("subsample_size_l1");

	int ks2,ss2;
	if (!use_one_conv_layer)
	{
		ks2 = conf.get_uint("size_kernel_l2");
		ss2 = conf.get_uint("subsample_size_l2");
	}
	
	int numhid = conf.get_uint("num_hidden");
	
	

			   
	sprintf(temp,"%s",conf.get_cstring("layer_connections_file"));	
	FILE *connections_file = fopen(temp,"r");
	if (connections_file==NULL) perror ("Error opening connections file");
	
	fprintf(stderr,"Reading the connections file\n");
	int nf1,nf2;
	
	fscanf(connections_file,"====\n");
	fscanf(connections_file,"%d\n",&nf1);
	int * table0_rows = new int[nf1]; 
	int * table1_rows;

	
	int totalrows = 0;
	for(int i=0;i<nf1;i++)
	{
		fscanf(connections_file,"%d %*[^\n]\n",&table0_rows[i]);
		totalrows += table0_rows[i];
	}
	
	idx<intg> table0 = idx<intg>(totalrows,2);
	
	totalrows = 0;
	if (!use_one_conv_layer)
	{
		fscanf(connections_file,"====\n");
		fscanf(connections_file,"%d\n",&nf2);
		table1_rows = new int[nf2]; 

		for(int i=0;i<nf2;i++)
		{
			fscanf(connections_file,"%d %*[^\n]\n",&table1_rows[i]);
			totalrows += table1_rows[i];
		}
	}
	idx<intg> table1 = idx<intg>(totalrows,2);
	
	rewind(connections_file);
	
	int ** out_in_map_0 = Malloc(int *,nf1);
	int ** out_in_map_1;
	
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
	
	if (!use_one_conv_layer)
	{
		out_in_map_1 = Malloc(int *,nf2);
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
	}
	
	
	fclose(connections_file);
	
	
	fprintf(stderr,"Constructing the gene_expression network\n");
	
	parameter<t_net> theparam(conf.exists("maxnumparams")? conf.get_uint("maxnumparams"): DEFAULT_MAX_NUMPARAMS); // create trainable parameter
	fprintf(stderr,"fixed numparams\n");
		
	idx<ubyte> trainingclasses(conf.get_uint("training_size"));
	idx<ubyte> testingclasses(conf.get_uint("testing_size"));
   
	if(typeofinput == UBYTE_INPUT)
	{
		idx<ubyte> trainingdata(conf.get_uint("training_size"),numtracks,tracklen,1);
		idx<ubyte> testingdata(conf.get_uint("testing_size"),numtracks,tracklen,1);
		fprintf(stderr,"before reading in data");
		for(int i=0;i<conf.get_uint("training_size");i++)
		{
			int temp2;
			fscanf(trainingdata_file,"%d",&temp2);
			trainingclasses.set(convert_classes(temp2),i);
			
			for(int j=0;j<numtracks;j++)
			{
				for(int k=0;k<tracklen;k++)
				{
					ubyte temp3;
					fscanf(trainingdata_file," %u",&temp3);
					trainingdata.set(temp3,i,j,k,0);
				}
			}
			fscanf(trainingdata_file,"\n");
		}
		fclose(trainingdata_file);
		
		sprintf(temp,"%s/%s",conf.get_cstring("root"),conf.get_cstring("test_name"));	
		FILE *testingdata_file = fopen(temp,"r");
		
		while ((firstchar = getc(testingdata_file)) == '#')
			fscanf(testingdata_file,"%*[^\n]\n");
		
		ungetc(firstchar,trainingdata_file);
		
		for(int i=0;i<conf.get_uint("testing_size");i++)
		{
			int temp2;
			fscanf(testingdata_file,"%d",&temp2);
			testingclasses.set(convert_classes(temp2),i);
			for(int j=0;j<numtracks;j++)
			{
				for(int k=0;k<tracklen;k++)
				{
					ubyte temp3;
					fscanf(testingdata_file," %u",&temp3);
					testingdata.set(temp3,i,j,k,0);
				}
			}
			fscanf(testingdata_file,"\n");
		}
		fclose(testingdata_file);
		fprintf(stderr,"read all the data no prob\n");
		
		
		
		labeled_datasource<t_net, ubyte, ubyte>
		train_ds(trainingdata,trainingclasses,conf.get_cstring("train_name")),
		test_ds(testingdata,testingclasses,conf.get_cstring("test_name"));
		
		train_ds.set_weigh_samples(conf.exists_true("wsamples"), true, 
								   conf.exists_true("wnorm"),
								   conf.exists("min_sample_weight") ?
								   conf.get_double("min_sample_weight") : 0.0);
		train_ds.set_shuffle_passes(conf.exists_true("shuffle_passes"));
		//! randomization
		if (conf.exists_true("fixed_randomization")) {
			cout << "Random seed is fixed (0)." << endl;
			init_drand(0); // fixed random seed
		} else {
			init_drand(time(NULL)); // initialize random seed
			cout << "Random seed is variable." << endl;
		}

		
		
		test_ds.set_test(); 
		
		idxdim dims(train_ds.sample_dims()); // get order and dimensions of sample		
		idx<t_net> targets = create_target_matrix<t_net>(train_ds.get_nclasses(), 1.0);

		// learning parameters
		gd_param gdp(/* double leta*/ conf.get_double("eta"),
					 /* double ln */ 	0.0,
					 /* double l1 */ 	0.0,
					 /* double l2 */ 	0.0,
					 /* int dtime */ 	0,
					 /* double iner */0.0, 
					 /* double a_v */ 0.0,
					 /* double a_t */ 0.0,
					 /* double g_t*/ 	0.0);
		infer_param infp;
		classifier_meter trainmeter, testmeter;
		forget_param_linear fgp(1, 0.5);
		
		float * correct;
		int max_train_acc_it = 0,max_test_acc_it=0;
		float max_train_acc = 0, max_test_acc = 0;
		
		
		if (!use_one_conv_layer)
		{

			idx<intg> tblmax = table1.select(1, 1);
			intg nf1 = 1 + idx_max(tblmax);
			idx<intg> table2 = full_table(nf1, numhid);
			
			int ki2 = ((tracklen-ks1+1)/ss1 - ks2 + 1)/ss2;


			net_cscscf<t_net> mygenes =  net_cscscf<t_net>(theparam,tracklen,1,
														   ks1,1,table0,
														   ss1,1,
														   ks2,1,table1,
														   ss2,1,
														   ki2,1,table2,
														   train_ds.get_nclasses());

			
			supervised_euclidean_machine<t_net, ubyte> thenet((module_1_1<t_net>&)mygenes, targets, dims);
			supervised_trainer<t_net, ubyte, ubyte> thetrainer(thenet, theparam);
			
			thenet.forget(fgp);
			
			// estimate second derivative on 100 iterations, using mu=0.02
			thetrainer.compute_diaghessian(train_ds, 100, 0.02);
			
			// first show classification results without training
			thetrainer.test(train_ds, trainmeter, infp);
			thetrainer.test(test_ds, testmeter, infp);
			
			// now do training iterations 
			cout << "Training network on REAL GENE EXPRESSION data set with " << train_ds.size();
			cout << " training samples and " << test_ds.size() <<" test samples:" << endl;
			
			
			
			
			
			if (plotresults)
			{
				correct = new float[2*(conf.get_uint("iterations")+1)];
				
				
				correct[0] = trainmeter.class_normalized_average_success();
				correct[conf.get_uint("iterations")+1] = testmeter.class_normalized_average_success();
			}
			for (uint i = 1; i <= conf.get_uint("iterations"); ++i) 
			{
				thetrainer.train(train_ds, trainmeter, gdp, 1,infp);	         // train
				thetrainer.test(train_ds, trainmeter, infp);         // test
				thetrainer.test(test_ds, testmeter, infp);	                 // test
				float cur_train_acc = trainmeter.class_normalized_average_success();
				float cur_test_acc = testmeter.class_normalized_average_success();
				
				if (cur_train_acc > max_train_acc)
				{
					max_train_acc = cur_train_acc;
					max_train_acc_it = i;
				}
				if (cur_test_acc > max_test_acc)
				{
					max_test_acc = cur_test_acc;
					max_test_acc_it = i;
				}
				
				if(storeweights)
				{
					sprintf(temp,"%s_weights_%d.txt",argv[1],i);
					int ks3 = ((tracklen-ks1+1)/ss1 - ks2 + 1)/ss2;
					
					store_weights(temp,theparam.x.idx_ptr(),4, nf1,nf2,numhid,2,ks1,ks2,ks3,out_in_map_0,out_in_map_1,table0_rows,table1_rows);
				}
				if (plotresults)
				{
					correct[i] = trainmeter.class_normalized_average_success();
					correct[i+conf.get_uint("iterations")+1] = testmeter.class_normalized_average_success();
					
					plot_results(correct,i,conf.get_uint("iterations"),max_train_acc,max_train_acc_it, max_test_acc,max_test_acc_it,results_plot_fn);
				}
				if (store_states)
				{
					sprintf(temp,"%s_states_%d.txt",argv[1],i);
					retrieve_network_states<ubyte>(thetrainer, test_ds, infp, temp);	
				}
				
				printf("MAX TRAIN CORRECT = %f (Iteration = %d)\nMAX TEST CORRECT = %f (Iteration = %d)\n",max_train_acc,max_train_acc_it, max_test_acc,max_test_acc_it);
				
				thetrainer.compute_diaghessian(train_ds, 100, 0.02); // recompute 2nd der
			}
			
			
		}
		else 
		{
			
			idx<intg> tblmax = table0.select(1, 1);
			intg nf1 = 1 + idx_max(tblmax);
			idx<intg> table1 = full_table(nf1, numhid);
			
			int ki1 = (tracklen-ks1+1)/ss1;
			

			net_cscf<t_net> mygenes =  net_cscf<t_net>(theparam,tracklen,1,ks1,1,
													   table0,ss1,1,ki1,1,
													   table1,train_ds.get_nclasses());

			
			supervised_euclidean_machine<t_net, ubyte> thenet((module_1_1<t_net>&)mygenes, targets, dims);
			supervised_trainer<t_net, ubyte, ubyte> thetrainer(thenet, theparam);
			
			thenet.forget(fgp);
			
			// estimate second derivative on 100 iterations, using mu=0.02
			thetrainer.compute_diaghessian(train_ds, 100, 0.02);
			
			// first show classification results without training
			thetrainer.test(train_ds, trainmeter, infp);
			thetrainer.test(test_ds, testmeter, infp);
			
			// now do training iterations 
			cout << "Training network on GENE EXPRESSION data set with " << train_ds.size();
			cout << " training samples and " << test_ds.size() <<" test samples:" << endl;
			
			
			if (plotresults)
			{
				correct = new float[2*(conf.get_uint("iterations")+1)];
								
				correct[0] = trainmeter.class_normalized_average_success();
				correct[conf.get_uint("iterations")+1] = testmeter.class_normalized_average_success();
			}
			for (uint i = 1; i <= conf.get_uint("iterations"); ++i) 
			{
				thetrainer.train(train_ds, trainmeter, gdp, 1,infp);	         // train
				thetrainer.test(train_ds, trainmeter, infp);         // test
				thetrainer.test(test_ds, testmeter, infp);	                 // test
				float cur_train_acc = trainmeter.class_normalized_average_success();
				float cur_test_acc = testmeter.class_normalized_average_success();
				
				if (cur_train_acc > max_train_acc)
				{
					max_train_acc = cur_train_acc;
					max_train_acc_it = i;
				}
				if (cur_test_acc > max_test_acc)
				{
					max_test_acc = cur_test_acc;
					max_test_acc_it = i;
				}
				
				if(storeweights)
				{
					sprintf(temp,"%s_weights_%d.txt",argv[1],i);
					int ks3 = (tracklen-ks1+1)/ss1;
					
					store_weights_use_one_conv_layer(temp,theparam.x.idx_ptr(),3,nf1,numhid,2,ks1,ks3,out_in_map_0,table0_rows);
				}
				if (plotresults)
				{
					correct[i] = trainmeter.class_normalized_average_success();
					correct[i+conf.get_uint("iterations")+1] = testmeter.class_normalized_average_success();
					
					plot_results(correct,i,conf.get_uint("iterations"),max_train_acc,max_train_acc_it, max_test_acc,max_test_acc_it,results_plot_fn);
				}
				if (store_states)
				{
					sprintf(temp,"%s_states_%d.txt",argv[1],i);
					retrieve_network_states<ubyte>(thetrainer, test_ds, infp, temp);	
				}
				printf("MAX TRAIN CORRECT = %f (Iteration = %d)\nMAX TEST CORRECT = %f (Iteration = %d)\n",max_train_acc,max_train_acc_it, max_test_acc,max_test_acc_it);
				
				thetrainer.compute_diaghessian(train_ds, 100, 0.02); // recompute 2nd der
			}
		}
	}
	else
	{
		idx<double> trainingdata(conf.get_uint("training_size"),numtracks,tracklen,1);
		idx<double> testingdata(conf.get_uint("testing_size"),numtracks,tracklen,1);

		fprintf(stderr,"before reading in data");
		for(int i=0;i<conf.get_uint("training_size");i++)
		{
			int temp2;
			fscanf(trainingdata_file,"%d",&temp2);
			trainingclasses.set(convert_classes(temp2),i);
			
			for(int j=0;j<numtracks;j++)
			{
				for(int k=0;k<tracklen;k++)
				{
					double temp3;
					fscanf(trainingdata_file," %lf",&temp3);
					trainingdata.set(temp3,i,j,k,0);
				}
			}
			fscanf(trainingdata_file,"\n");
		}
		fclose(trainingdata_file);
		
		sprintf(temp,"%s/%s",conf.get_cstring("root"),conf.get_cstring("test_name"));	
		FILE *testingdata_file = fopen(temp,"r");
		
		while ((firstchar = getc(testingdata_file)) == '#')
			fscanf(testingdata_file,"%*[^\n]\n");
		
		ungetc(firstchar,trainingdata_file);
		
		for(int i=0;i<conf.get_uint("testing_size");i++)
		{
			int temp2;
			fscanf(testingdata_file,"%d",&temp2);
			testingclasses.set(convert_classes(temp2),i);
			for(int j=0;j<numtracks;j++)
			{
				for(int k=0;k<tracklen;k++)
				{
					double temp3;
					fscanf(testingdata_file," %lf",&temp3);
					testingdata.set(temp3,i,j,k,0);
				}
			}
			fscanf(testingdata_file,"\n");
		}
		fclose(testingdata_file);
		fprintf(stderr,"read all the data no prob\n");
				
		labeled_datasource<t_net, double, ubyte>
		train_ds(trainingdata,trainingclasses,conf.get_cstring("train_name")),
		test_ds(testingdata,testingclasses,conf.get_cstring("test_name"));
		test_ds.set_test(); 
		train_ds.set_weigh_samples(conf.exists_true("wsamples"), true, 
								   conf.exists_true("wnorm"),
								   conf.exists("min_sample_weight") ?
								   conf.get_double("min_sample_weight") : 0.0);
		train_ds.set_shuffle_passes(conf.exists_true("shuffle_passes"));
		//! randomization
		if (conf.exists_true("fixed_randomization")) {
			cout << "Random seed is fixed (0)." << endl;
			init_drand(0); // fixed random seed
		} else {
			init_drand(time(NULL)); // initialize random seed
			cout << "Random seed is variable." << endl;
		}

		idxdim dims(train_ds.sample_dims()); // get order and dimensions of sample		
		idx<t_net> targets = create_target_matrix<t_net>(train_ds.get_nclasses(), 1.0);
		
		// learning parameters
		gd_param gdp(/* double leta*/ conf.get_double("eta"),
					 /* double ln */ 	0.0,
					 /* double l1 */ 	0.0,
					 /* double l2 */ 	0.0,
					 /* int dtime */ 	0,
					 /* double iner */0.0, 
					 /* double a_v */ 0.0,
					 /* double a_t */ 0.0,
					 /* double g_t*/ 	0.0);
		infer_param infp;
		classifier_meter trainmeter, testmeter;
		forget_param_linear fgp(1, 0.5);
		
		float * correct;
		int max_train_acc_it = 0,max_test_acc_it=0;
		float max_train_acc = 0, max_test_acc = 0;


		if (!use_one_conv_layer)
		{
			idx<intg> tblmax = table1.select(1, 1);
			intg nf1 = 1 + idx_max(tblmax);
			idx<intg> table2 = full_table(nf1, numhid);
			
			int ki2 = ((tracklen-ks1+1)/ss1 - ks2 + 1)/ss2;
			
			net_cscscf<t_net> mygenes =  net_cscscf<t_net>(theparam,tracklen,1,
														   ks1,1,table0,
														   ss1,1,
														   ks2,1,table1,
														   ss2,1,
														   ki2,1,table2,
														   train_ds.get_nclasses());
			
			
			supervised_euclidean_machine<t_net, ubyte> thenet((module_1_1<t_net>&)mygenes, targets, dims);
			supervised_trainer<t_net, double, ubyte> thetrainer(thenet, theparam);
			thenet.forget(fgp);
			
			// estimate second derivative on 100 iterations, using mu=0.02
			thetrainer.compute_diaghessian(train_ds, 100, 0.02);
			
			// first show classification results without training
			thetrainer.test(train_ds, trainmeter, infp);
			thetrainer.test(test_ds, testmeter, infp);
			
			// now do training iterations 
			cout << "Training network on REAL GENE EXPRESSION data set with " << train_ds.size();
			cout << " training samples and " << test_ds.size() <<" test samples:" << endl;
			
			
			if (plotresults)
			{
				correct = new float[2*(conf.get_uint("iterations")+1)];
				
				
				correct[0] = trainmeter.class_normalized_average_success();
				correct[conf.get_uint("iterations")+1] = testmeter.class_normalized_average_success();
			}
			for (uint i = 1; i <= conf.get_uint("iterations"); ++i) 
			{
				thetrainer.train(train_ds, trainmeter, gdp, 1,infp);	         // train
				thetrainer.test(train_ds, trainmeter, infp);         // test
				thetrainer.test(test_ds, testmeter, infp);	                 // test
				float cur_train_acc = trainmeter.class_normalized_average_success();
				float cur_test_acc = testmeter.class_normalized_average_success();
				
				if (cur_train_acc > max_train_acc)
				{
					max_train_acc = cur_train_acc;
					max_train_acc_it = i;
				}
				if (cur_test_acc > max_test_acc)
				{
					max_test_acc = cur_test_acc;
					max_test_acc_it = i;
				}
				
				if(storeweights)
				{
					sprintf(temp,"%s_weights_%d.txt",argv[1],i);
					int ks3 = ((tracklen-ks1+1)/ss1 - ks2 + 1)/ss2;
					
					store_weights(temp,theparam.x.idx_ptr(),4, nf1,nf2,numhid,2,ks1,ks2,ks3,out_in_map_0,out_in_map_1,table0_rows,table1_rows);
				}
				if (plotresults)
				{
					correct[i] = trainmeter.class_normalized_average_success();
					correct[i+conf.get_uint("iterations")+1] = testmeter.class_normalized_average_success();
					
					plot_results(correct,i,conf.get_uint("iterations"),max_train_acc,max_train_acc_it, max_test_acc,max_test_acc_it,results_plot_fn);
				}
				if (store_states)
				{
					sprintf(temp,"%s_states_%d.txt",argv[1],i);
					retrieve_network_states<double>(thetrainer, test_ds, infp, temp);	
				}
				printf("MAX TRAIN CORRECT = %f (Iteration = %d)\nMAX TEST CORRECT = %f (Iteration = %d)\n",max_train_acc,max_train_acc_it, max_test_acc,max_test_acc_it);
				
				thetrainer.compute_diaghessian(train_ds, 100, 0.02); // recompute 2nd der
			}
		}
		else 
		{
			idx<intg> tblmax = table0.select(1, 1);
			intg nf1 = 1 + idx_max(tblmax);
			idx<intg> table1 = full_table(nf1, numhid);
			
			int ki1 = (tracklen-ks1+1)/ss1;
			
			net_cscf<t_net> mygenes =  net_cscf<t_net>(theparam,tracklen,1,ks1,1,
													   table0,ss1,1,ki1,1,
													   table1,train_ds.get_nclasses());
			supervised_euclidean_machine<t_net, ubyte> thenet((module_1_1<t_net>&)mygenes, targets, dims);
			supervised_trainer<t_net, double, ubyte> thetrainer(thenet, theparam);
			
			thenet.forget(fgp);
			
			// estimate second derivative on 100 iterations, using mu=0.02
			thetrainer.compute_diaghessian(train_ds, 100, 0.02);
			
			// first show classification results without training
			thetrainer.test(train_ds, trainmeter, infp);
			thetrainer.test(test_ds, testmeter, infp);
			
			// now do training iterations 
			cout << "Training network on REAL GENE EXPRESSION data set with " << train_ds.size();
			cout << " training samples and " << test_ds.size() <<" test samples:" << endl;
			
			if (plotresults)
			{
				correct = new float[2*(conf.get_uint("iterations")+1)];
				
				
				correct[0] = trainmeter.class_normalized_average_success();
				correct[conf.get_uint("iterations")+1] = testmeter.class_normalized_average_success();
			}
			for (uint i = 1; i <= conf.get_uint("iterations"); ++i) 
			{
				thetrainer.train(train_ds, trainmeter, gdp, 1,infp);	         // train
				thetrainer.test(train_ds, trainmeter, infp);         // test
				thetrainer.test(test_ds, testmeter, infp);	                 // test
				float cur_train_acc = trainmeter.class_normalized_average_success();
				float cur_test_acc = testmeter.class_normalized_average_success();
				
				if (cur_train_acc > max_train_acc)
				{
					max_train_acc = cur_train_acc;
					max_train_acc_it = i;
				}
				if (cur_test_acc > max_test_acc)
				{
					max_test_acc = cur_test_acc;
					max_test_acc_it = i;
				}
				
				if(storeweights)
				{
					sprintf(temp,"%s_weights_%d.txt",argv[1],i);
					int ks3 = (tracklen-ks1+1)/ss1;
					
					store_weights_use_one_conv_layer(temp,theparam.x.idx_ptr(),3,nf1,numhid,2,ks1,ks3,out_in_map_0,table0_rows);
				}
				if (plotresults)
				{
					correct[i] = trainmeter.class_normalized_average_success();
					correct[i+conf.get_uint("iterations")+1] = testmeter.class_normalized_average_success();
					
					plot_results(correct,i,conf.get_uint("iterations"),max_train_acc,max_train_acc_it, max_test_acc,max_test_acc_it,results_plot_fn);
				}
				if (store_states)
				{
					sprintf(temp,"%s_states_%d.txt",argv[1],i);
					retrieve_network_states<double>(thetrainer, test_ds, infp, temp);	
				}
				
				printf("MAX TRAIN CORRECT = %f (Iteration = %d)\nMAX TEST CORRECT = %f (Iteration = %d)\n",max_train_acc,max_train_acc_it, max_test_acc,max_test_acc_it);
				
				thetrainer.compute_diaghessian(train_ds, 100, 0.02); // recompute 2nd der
			}
		}
	}
	
	
	
		

	//	int typenet = conf.get_uint("type_network");
	return 0;
}



template <class T> void retrieve_network_states(supervised_trainer<t_net, T, ubyte> &thetrainer, labeled_datasource<t_net, T, ubyte> &ds, infer_param &infp, char * fn)
{
	bbstate_idx<t_net> *input = new bbstate_idx<t_net>(ds.sample_dims());
	bbstate_idx<ubyte> label;
	bbstate_idx<t_net>* hi; //! temporary buffer pointer
	bbstate_idx<t_net>* ho; //! temporary buffer pointer
	layers <t_net, bbstate_idx<t_net> > &fmodloc = (layers <t_net, bbstate_idx<t_net> >&) thetrainer.machine.fmod;

	FILE *fin = fopen(fn,"w");
	const char * layersnames [4] = {"convolution","subsamples","full","output"};
	int outputlayerinds [4] = {2,3,6,7};
	ds.seek_begin();
	for (unsigned int i = 0; i < ds.size(); ++i) 
	{
		ds.fprop(*input, label);
		bool correct = thetrainer.test_sample(*input, label, infp);
		fprintf(fin,"============\nlabel: %d\npredicted label: %d\n",label.x.get(),(correct ? label.x.get() : (label.x.get() +1) % 2 ));
		store_states(*input,fin,"input");

		int outputind = 0;
		if(1)
		{
			hi = input;
			// loop over modules
			uint i1;
			for(i1 = 0; i1 < fmodloc.modules->size(), outputind<4; i1++)
			{
				// if last module, output into out
				if (i1 == fmodloc.modules->size() - 1)
					ho = new bbstate_idx<t_net>(ds.sample_dims());
				else 
				{ // not last module, use hidden buffers
					ho = (bbstate_idx<t_net>*)(*(fmodloc.hiddens))[i1];
					// allocate hidden buffer if necessary
					if (ho == NULL)
					{
						// create idxdim of same order but sizes 1
						idxdim d = hi->x.get_idxdim();
						for (int k = 0; k < d.order(); ++k)
							d.setdim(k, 1);
						// assign buffer
						(*(fmodloc.hiddens))[i1] = new bbstate_idx<t_net>(d);
						ho = (bbstate_idx<t_net>*)(*(fmodloc.hiddens))[i1];
					}
				}
				// run module
				(*(fmodloc.modules))[i1]->fprop(*hi,*ho);
				if(i1 == outputlayerinds[outputind])
				{
					store_states(*ho,fin,layersnames[outputind]);
					outputind ++;
				}
				hi = ho;
			}
			if (i1 == fmodloc.modules->size())
				delete ho;
		}
		ds.next();
	}
	delete input;
	fclose(fin);
}


void store_states(bbstate_idx<t_net> &data,FILE * fin, const char *layername)
{
	double * xvals = data.x.idx_ptr();
	
	int numfeatures = data.x.dim(0);
	int length = data.x.dim(1);

	int ind = 0;
	fprintf(fin,"%s:",layername);
	for(int i=0;i<numfeatures;i++)
	{
		if (i>0)
			fprintf(fin," |");
		
		for(int j=0;j<length;j++)
			fprintf(fin," %lf",xvals[ind++]);
	}
	fprintf(fin,"\n");
}

						
void plot_results(float *correct,int num,int total,float maxtrain, int maxtrainind, float maxtest, int maxtestind, const char *fn)
{
	FILE *p = popen("gnuplot -persist","w");
	if (p)
	{
		fprintf(p, "set xrange [0:%d]\n",num);
		fprintf(p, "set terminal png\n");
		fprintf(p, "set output \"%s\"\n",fn);
		fprintf(p, "set mytics 5\n");
		fprintf(p, "set mxtics 5\n");
		fprintf(p, "set grid xtics ytics mxtics mytics\n");
		fprintf(p, "set title \"Hit-rate graph: max train=%f(it=%d); max test=%f(it=%d)\"\n",maxtrain,maxtrainind,maxtest,maxtestind);
		fprintf(p, "plot \"-\" binary endian=little array=%d format=\"%%float\" title \"Training\" with lines, \"-\" binary endian=little array=%d format=\"%%float\" title \"Test\" with lines\n",total+1,total+1);
		
		fwrite(correct, sizeof(float), 2*(1+total), p );
		fflush(p);
		
		fflush(stderr);
		fprintf(p,"exit \n");
		pclose(p);
	}		
}



void store_weights(char *fn,double *weights_vals,int numlayers, int nf1,int nf2, int nf3, int numout, int ks1, int ks2, int ks3, int **out_in_map_0, int **out_in_map_1, int * table0_rows, int * table1_rows)
{
	FILE * outputfile = fopen(fn,"w");
	int ind = 0;

	fprintf(outputfile,"######LAYER 0######%\n");
	for(int i=0;i<nf1;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		for(int j=0;j<table0_rows[i];j++)
		{
			fprintf(outputfile,"%d:",out_in_map_0[i][j]);
			for(int k=0;k<ks1;k++)
				fprintf(outputfile," %lf",weights_vals[ind++]);
			fprintf(outputfile,"\n");
		}
	}
	
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<nf1;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	
	fprintf(outputfile,"######LAYER 0S######%\n");

	for(int i=0;i<nf1;i++)
	{
		fprintf(outputfile,"%d: %lf %lf\n",i,weights_vals[ind],weights_vals[nf1+ind]);
		ind++;
	}
	
	ind += nf1;
	
	fprintf(outputfile,"######LAYER 1######%\n");

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
	}

	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<nf2;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	

	fprintf(outputfile,"######LAYER 1S######%\n");

	for(int i=0;i<nf2;i++)
	{
		fprintf(outputfile,"%d: %lf %lf\n",i,weights_vals[ind],weights_vals[nf2+ind]);
		ind++;
	}
	
	ind += nf2;

	fprintf(outputfile,"######LAYER 2######%\n");

	for(int i=0;i<nf3;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		
		for(int j=0;j<nf2;j++)
		{
			fprintf(outputfile,"%d:",j);
			for(int k=0;k<ks3;k++)
				fprintf(outputfile," %lf",weights_vals[ind++]);
			fprintf(outputfile,"\n");
		}
	}

	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<nf3;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	
	fprintf(outputfile,"######LAYER 3######%\n");
	
	for(int i=0;i<numout;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		
		for(int j=0;j<nf3;j++)
			fprintf(outputfile," %lf",weights_vals[ind++]);
		fprintf(outputfile,"\n");
	}
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<numout;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	
	fclose(outputfile);
}


void store_weights_use_one_conv_layer(char *fn,double *weights_vals,int numlayers, int nf1,int nf2, int numout, int ks1, int ks2, int **out_in_map_0, int * table0_rows)
{
	FILE * outputfile = fopen(fn,"w");
	int ind = 0;
	
	fprintf(outputfile,"######LAYER 0######%\n");
	for(int i=0;i<nf1;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		for(int j=0;j<table0_rows[i];j++)
		{
			fprintf(outputfile,"%d:",out_in_map_0[i][j]);
			for(int k=0;k<ks1;k++)
				fprintf(outputfile," %lf",weights_vals[ind++]);
			fprintf(outputfile,"\n");
		}
	}
	
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<nf1;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	
	fprintf(outputfile,"######LAYER 0S######%\n");
	
	for(int i=0;i<nf1;i++)
	{
		fprintf(outputfile,"%d: %lf %lf\n",i,weights_vals[ind],weights_vals[nf1+ind]);
		ind++;
	}
	
	ind += nf1;
		
	fprintf(outputfile,"######LAYER 1######%\n");
	
	for(int i=0;i<nf2;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		
		for(int j=0;j<nf1;j++)
		{
			fprintf(outputfile,"%d:",j);
			for(int k=0;k<ks2;k++)
				fprintf(outputfile," %lf",weights_vals[ind++]);
			fprintf(outputfile,"\n");
		}
	}
	
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<nf2;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	
	fprintf(outputfile,"######LAYER 2######%\n");
	
	for(int i=0;i<numout;i++)
	{
		fprintf(outputfile,"=====Output feature %d=========\n",i);
		
		for(int j=0;j<nf2;j++)
			fprintf(outputfile," %lf",weights_vals[ind++]);
		fprintf(outputfile,"\n");
	}
	fprintf(outputfile,"=====Output feature BIASES=========\n");
	for(int i=0;i<numout;i++)
		fprintf(outputfile,"%d: %lf\n",i,weights_vals[ind++]);
	
	fclose(outputfile);
}
