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

# I/O
# input_folder="/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/data"
# output_folder="/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/hdf5"

input_folder="/Users/ricard/data/gastrulation10x/mofa/data"
output_folder="/Users/ricard/data/gastrulation10x/mofa/hdf5"

# Number of trials
ntrials=1

# Number of factors
factors=5

# ELBO Convergence settings
convergence_mode="fast"
maxiter=15
elbofreq=1
start_elbo=1


for j in $(seq 1 $ntrials); do
	outfile="$output_folder/model_$j.hdf5"
	cmd="python run.py --input_folder $input_folder --outfile $outfile --factors $factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq --seed $j --verbose"
	eval $cmd
done