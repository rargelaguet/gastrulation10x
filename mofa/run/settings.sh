job() {
	command="bsub -M $(( $1 * 1024 )) -n $2 \
	-q $3 ${@:4}"
	echo $command
	eval $command
}

gpujob() {
	command="bsub -M $(( $1 * 1024 )) -n $2 \
	-q $3 -P gpu -gpu - ${@:4}"
	echo $command
	eval $command
}

#########
## I/O ##
#########

if echo $HOSTNAME | grep -q "ricard"; then
    input_folder="/Users/ricard/data/gastrulation10x_mofa/data"
    output_folder="/Users/ricard/data/gastrulation10x_mofa/hdf5"
elif echo $HOSTNAME | grep -q "ebi"; then
	input_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/data'
	output_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/hdf5'
elif echo $HOSTNAME | grep -q "embl"; then
	input_folder='/g/stegle/ricard/gastrulation10x_mofa/data'
	output_folder='/g/stegle/ricard/gastrulation10x_mofa/hdf5'
else
	echo "Computer not recognised"; exit
fi

###################
## Model options ##
###################

# Nomber of factors
num_factors=15

#######################
## Training settings ##
#######################

# Convergence mode ("fast", "medium", "slow")
convergence_mode="fast"

# Maximum number of iterations
maxiter=250

# frequency of ELBO computation
elbofreq=1

# iteration to start ELBO computations
start_elbo=1
