#!/g/stegle/ricard/anaconda/envs/mofa2/bin/python
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        # number of cores
#SBATCH --mem 10G                   # memory pool for all cores
#SBATCH -t 0-10:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o slurm.%N.%j.out          # STDOUT
#SBATCH -e slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=ricard@ebi.ac.uk

# sbatch -p gpu --gres gpu:1 run.py

from mofapy2.run.entry_point import entry_point
import pandas as pd
import numpy as np
import argparse

################################
## Initialise argument parser ##
################################

p = argparse.ArgumentParser( description='' )

# I/O options
p.add_argument( '--input_folder',          type=str,              required=True,          help='Input data file (matrix format)' )
p.add_argument( '--outfile',               type=str,              required=True,          help='Output file to store the model (.hdf5)' )

# Model options
p.add_argument( '--factors',               type=int,              default=25,             help='Number of factors' )

# Training options
p.add_argument( '--seed',                  type=int,              default=42,             help='Random seed' )
p.add_argument( '--verbose',               action="store_true",                           help='Do stochastic inference?' )
p.add_argument( '--start_elbo',            type=int,              default=1,              help='Start of ELBO computation' )
p.add_argument( '--elbo_freq',             type=int,              default=1,              help='Frequency of ELBO computation' )
p.add_argument( '--iterations',            type=int,              default=100,            help='Number of iterations')
p.add_argument( '--convergence_mode',      type=str,              default="medium",       help='Convergence mode')

# Stochastic inference options
p.add_argument( '--stochastic_inference',  action="store_true",                           help='Do stochastic inference?' )
p.add_argument( '--batch_size',            type=float,              default=0.5,          help='Batch size (fraction of samples)' )
p.add_argument( '--learning_rate',         type=float,              default=0.75,         help='Learning rate' )
p.add_argument( '--forgetting_rate',       type=float,              default=0.,           help='Forgetting rate for stochastic inference')

args = p.parse_args()

###############
## Load data ##
###############

views = [""]
sample_groups = [ "E6.5_1", "E6.5_5", "E7.0_10", "E7.0_14", "E7.25_26", "E7.25_27"]
# sample_groups = [ "E6.5_1"]
data = [None]*len(views)
for m in range(len(views)):
    data[m] = [None]*len(sample_groups)
    for g in range(len(sample_groups)):
        datafile = "%s/%s.txt.gz" % (args.input_folder, sample_groups[g])
        # data[m][g] = np.load(datafile)
        data[m][g] = pd.read_csv(datafile, header=0, sep='\t')

# data = pd.read_csv(args.input_folder, delimiter="\t", header=0)

########################
## Create MOFA object ##
########################

# initialise entry point    
ent = entry_point()

# Set data
features_names = [ x[0].columns.values for x in data ]
samples_names = [ x.index for x in data[0] ]
ent.set_data_matrix(data, views_names=["RNA"], groups_names=sample_groups, samples_names=samples_names, features_names=features_names)

# Set model options
# ent.set_model_options(factors=args.factors, spikeslab_factors=False, spikeslab_weights=True) # original
ent.set_model_options(factors=args.factors, spikeslab_factors=False, spikeslab_weights=False)

# Set training options
ent.set_train_options(iter=args.iterations, convergence_mode=args.convergence_mode, startELBO=args.start_elbo, freqELBO=args.elbo_freq, verbose=args.verbose, seed=args.seed)

# Set stochastic inference options
if args.stochastic_inference:
    ent.set_stochastic_options(learning_rate=args.learning_rate, forgetting_rate=args.forgetting_rate, batch_size=args.batch_size)

###############################
## Build and train the model ##
###############################


# Build the model
ent.build()

# Train the model
ent.run()

####################
## Save the model ##
####################

ent.save(args.outfile, save_data=True)
