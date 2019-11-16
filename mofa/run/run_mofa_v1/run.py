from mofapy.core.entry_point import entry_point
import pandas as pd
import numpy as np
import sys
import argparse

######################
## Define arguments ##
######################

p = argparse.ArgumentParser( description='' )

# I/O options
p.add_argument( '--input_folder',			   type=str,				required=True,		   help='Input data file (matrix format)' )
p.add_argument( '--outfile',			   type=str,				required=True,		   help='Output file to store the model (.hdf5)' )

# Model options
p.add_argument( '--factors',			type=int,			  default=25,			 help='Number of factors' )

# Training options
p.add_argument( '--seed',				  type=int,				default=0,			   help='Random seed' )
p.add_argument( '--verbose',  action="store_true",							  help='Do stochastic inference?' )
p.add_argument( '--start_elbo',			type=int,			  default=1,			 help='Start of ELBO computation' )
p.add_argument( '--elbo_freq',			type=int,			  default=1,			 help='Frequency of ELBO computation' )
p.add_argument( '--iterations',	   type=int,			  default=100,			  help='Number of iterations')
p.add_argument( '--convergence_mode',	   type=str,			  default="medium",			  help='Convergence mode')

# Stochastic inference options
p.add_argument( '--stochastic_inference',  action="store_true",							  help='Do stochastic inference?' )
p.add_argument( '--batch_size',			type=float,			  default=0.5,			 help='Batch size (fraction of samples)' )
p.add_argument( '--learning_rate',			type=float,			  default=0.75,			 help='Batch size (fraction of samples)' )
p.add_argument( '--forgetting_rate',	   type=float,			  default=0.,			  help='Forgetting rate for stochastic inference')

args = p.parse_args()

###############
## Load data ##
###############

views = [""]
sample_groups = [ "E6.5_1", "E6.5_5", "E7.0_10", "E7.0_14", "E7.25_26", "E7.25_27"]
data = [None]*len(views)
for m in range(len(views)):
	data[m] = [None]*len(sample_groups)
	for g in range(len(sample_groups)):
		datafile = "%s/%s.txt.gz" % (args.input_folder, sample_groups[g])
		data[m][g] = pd.read_csv(datafile, header=0, sep='\t')
		# Center features per group
		data[m][g] -= np.nanmean(data[m][g], axis=0)
	data[m] = np.concatenate(data[m], axis=0)
	data[m] = data[m][:,data[m].std(axis=0)>0]


#####################
## Train the model ##
#####################

# initialise biofam	
ent = entry_point()

# Set data
ent.set_data(data)

# Set model options
ent.set_model_options(factors=args.factors, likelihoods=["gaussian"], sparsity=True)

# Set data options
ent.set_data_options(center_features=True, scale_views=False)

# Parse the data
ent.parse_data()

# Set training options
if args.seed == 0: args.seed = None
ent.set_train_options(iter=args.iterations, verbose=args.verbose, elbofreq=args.elbo_freq, seed=args.seed, tolerance=0.5)

# Define prior distributions
ent.define_priors()

# Define initialisations of variational distributions
ent.define_init()

# Parse intercept factor
ent.parse_intercept()

# Train the model
ent.train_model()


# Save the model
ent.save_model(args.outfile)
