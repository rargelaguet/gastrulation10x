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
    input_folder="/Users/ricard/data/gastrulation10x/mofa/data"
    output_folder="/Users/ricard/data/gastrulation10x/mofa/hdf5"
elif [ "$HOSTNAME" = login.cluster.embl.de ]; then
	input_folder='/g/stegle/ricard/gastrulation10x/mofa/data'
	output_folder='/g/stegle/ricard/gastrulation10x/mofa/hdf5'
elif [ "$HOSTNAME" = noah-login-03.ebi.ac.uk ]; then
	input_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/data'
	output_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/hdf5'
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
convergence_mode="medium"

# Maximum number of iterations
iter=2500

# frequency of ELBO computation
elbofreq=1

# iteration to start ELBO computations
start_elbo=1

# Number of trials
ntrials=1

###################################
## Stochastic inference settings ##
###################################

# Iteration to start stochastic inference
start_stochastic=1

# batch size (a multiple of 0.05)
# batch_size=( 0.10 0.25 0.50 )
batch_size=( 0.25 )

# Starting learning rate
# learning_rate=( 0.10 0.25 0.50 0.75 1.0 )
learning_rate=( 0.75 )

# Forgetting rate
# forgetting_rate=( 0.50 1.0 1.5 2.0 )
forgetting_rate=( 0.5 )

##############################
## Run stochastic inference ##
##############################

# for i in "${batch_size[@]}"; do
# 	for j in "${forgetting_rate[@]}"; do
# 		for k in "${tau[@]}"; do
# 			for trial in $(seq 1 $ntrials); do
# 				outfile="$output_folder/E6.5_k15_${i}_${j}_${k}_${trial}.hdf5"
# 				cmd="python run.py --datafile $datafile --outfile $outfile --seed $trial --stochastic_inference --batch_size $i --forgetting_rate $j --tau $k --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq"
# 				# gpujob 10 1 research-rh7 $cmd
# 				# job 10 1 research-rh7 $cmd
# 				# eval $cmd
# 			done
# 		done
# 	done
# done

############################
## Run standard inference ##
############################

# for j in $(seq 1 $ntrials); do
# 	# outfile="$output_folder/E6.5-E7.25_$j.hdf5"
# 	outfile="$output_folder/model_$j.hdf5"
# 	cmd="run.py --input_folder $input_folder --outfile $outfile --factors $factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq --seed $j"
# 	# gpujob 15 1 research-rh7 $cmd
# 	jobid=$(sbatch -p gpu --gres gpu:1 --exclude=gpu8,gpu9 $cmd | awk '{print $NF}')
# 	# eval $cmd
# done