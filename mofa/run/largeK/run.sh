# job() {
# 	command="bsub -M $(( $1 * 1024 )) -n $2 \
# 	-q $3 ${@:4}"
# 	echo $command
# 	eval $command
# }

# gpujob() {
# 	command="bsub -M $(( $1 * 1024 )) -n $2 \
# 	-q $3 -P gpu -gpu - ${@:4}"
# 	echo $command
# 	eval $command
# }

#########
## I/O ##
#########

if [ "$HOSTNAME" = ricard ]; then
    input_folder="/Users/ricard/data/gastrulation10x_mofa/data"
    output_folder="/Users/ricard/data/gastrulation10x_mofa/hdf5"
elif [ "$HOSTNAME" = login.cluster.embl.de ]; then
	input_folder='/g/stegle/ricard/gastrulation10x_mofa/data'
	output_folder='/g/stegle/ricard/gastrulation10x_mofa/hdf5'
elif [ "$HOSTNAME" = noah-login-03.ebi.ac.uk ]; then
	input_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/data'
	output_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/hdf5'
else
	echo "Computer not recognised"; exit
fi

###################
## Model options ##
###################

# Nomber of factors
num_factors=25

#######################
## Training settings ##
#######################

# Convergence mode ("fast", "medium", "slow")
convergence_mode="fast"

# Maximum number of iterations
iter=1000

# frequency of ELBO computation
elbofreq=1

# iteration to start ELBO computations
start_elbo=1

# Number of trials
ntrials=1


#########
## Run ##
#########

for j in $(seq 1 $ntrials); do
	outfile="$output_folder/model_$j.hdf5"
	cmd="../run.py --input_folder $input_folder --outfile $outfile --factors $factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq --seed $j"
	job 15 1 research-rh7 $cmd
	# jobid=$(sbatch -p gpu --gres gpu:1 --exclude=gpu8,gpu9 $cmd | awk '{print $NF}')
	# eval $cmd
done