from biofam.run.entry_point import entry_point
from time import time
import pandas as pd
import numpy as np
import sys
import resource
import argparse

p = argparse.ArgumentParser( description='' )
p.add_argument( '--datafile',               type=str,                required=True,           help='Input data file (matrix format)' )
p.add_argument( '--outfile',               type=str,                required=True,           help='Output file to store the model (.hdf5)' )
p.add_argument( '--seed',                  type=int,                default=0,               help='Random seed' )
p.add_argument( '--stochastic_inference',  action="store_true",                              help='Do stochastic inference?' )
p.add_argument( '--tau',            type=float,              default=0.75,             help='Batch size (fraction of samples)' )
p.add_argument( '--verbose',  action="store_true",                              help='Do stochastic inference?' )
p.add_argument( '--nostop',  action="store_true",                              help='Do not stop after reaching convergence?' )
p.add_argument( '--batch_size',            type=float,              default=0.5,             help='Batch size (fraction of samples)' )
p.add_argument( '--elbo_freq',            type=int,              default=1,             help='Frequency of ELBO computation' )
p.add_argument( '--forgetting_rate',       type=float,              default=0.,              help='Forgetting rate for stochastic inference')
p.add_argument( '--iterations',       type=int,              default=100,              help='Number of iterations')
p.add_argument( '--factors',            type=int,              default=25,             help='Number of factors' )
p.add_argument( '--start_elbo',            type=int,              default=1,             help='Start of ELBO computation' )
p.add_argument( '--convergence_mode',       type=str,              default="medium",              help='Convergence mode')

args = p.parse_args()


# datafile = "/Users/ricard/data/gastrulation10x/mofa/data.txt.gz"
data = pd.read_csv(args.datafile, delimiter="\t", header=0)

# initialise biofam    
ent = entry_point()

# Set data options
lik="gaussian"
ent.set_data_options(likelihoods=lik, center_features_per_group=True, scale_features=False, scale_views=False)

# Set data
data = [[data]]
features_names_dict = {k: list(v[0].columns.values) for (k,v) in zip(["RNA"],data)}
samples_names_dict = {k: list(v.index) for (k,v) in zip(["E6.5_to_E7.5"],data[0])}
ent.set_data_matrix(data, samples_names_dict=samples_names_dict, features_names_dict=features_names_dict)


# Set model options
ent.set_model_options(factors=args.factors, likelihoods=lik, sl_z=True, sl_w=True, ard_z=False, ard_w=True)

# Set training options
if args.seed == 0: args.seed = None
ent.set_train_options(iter=args.iterations, convergence_mode=args.convergence_mode, startELBO=args.start_elbo, elbofreq=args.elbo_freq, gpu_mode=True, verbose=args.verbose, seed=args.seed, nostop=args.nostop, startSparsity=1)

# Set stochastic vb options
if args.stochastic_inference:
    ent.set_stochasticity_options(tau=args.tau, forgetting_rate=args.forgetting_rate, batch_size=args.batch_size)

# Build the model
ent.build()

# Train the model
ent.run()

# Save the model
ent.save(args.outfile)