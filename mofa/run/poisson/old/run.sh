job() {
	command="bsub -M $(( $1 * 1024 )) -n $2 \
	-q $3 ${@:4}"
	echo $command
	eval $command
}

# gpujob() {
# 	command="bsub -M $(( $1 * 1024 )) -n $2 \
# 	-q $3 -P gpu -gpu - ${@:4}"
# 	echo $command
# 	eval $command
# }

#########
## I/O ##
#########

if echo $HOSTNAME | grep -q "ricard"; then
    input_folder="/Users/ricard/data/gastrulation10x_mofa/data"
    output_folder="/Users/ricard/data/gastrulation10x_mofa/hdf5"
elif echo $HOSTNAME | grep -q "ebi"; then
	input_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/data'
	output_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/hdf5_largeK'
else
	echo "Computer not recognised"; exit
fi

###################
## Model options ##
###################

# Nomber of factors
num_factors=10

#######################
## Training settings ##
#######################

# Convergence mode ("fast", "medium", "slow")
convergence_mode="fast"

# Maximum number of iterations
maxiter=1000

# iteration to start and frequency of ELBO computations
start_elbo=1
freq_elbo=1

# Number of trials
ntrials=1


#########
## Run ##
#########

for j in $(seq 1 $ntrials); do
	outfile="$output_folder/foo_$j.hdf5"
	cmd="python ../run.py --input_folder $input_folder --outfile $outfile --factors $num_factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $freq_elbo --seed $j"
	# job 15 2 research-rh7 $cmd
	# jobid=$(sbatch -p gpu --gres gpu:1 --exclude=gpu8,gpu9 $cmd | awk '{print $NF}')
	# eval `python $cmd`
	eval $cmd
done