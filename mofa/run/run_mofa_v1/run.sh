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
input_folder="/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/data"
output_folder="/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/hdf5/mofa_v1"

# input_folder="/Users/ricard/data/gastrulation10x/mofa/data"
# output_folder="/Users/ricard/data/gastrulation10x/mofa/hdf5/mofa_v1"

# Number of trials
ntrials=1

# Number of factors
factors=10

# ELBO Convergence settings
maxiter=11
elbofreq=5


for j in $(seq 1 $ntrials); do
	outfile="$output_folder/E6.5-E7.25_$j.hdf5"
	cmd="python run.py --input_folder $input_folder --outfile $outfile --factors $factors --iterations $maxiter --elbo_freq $elbofreq --seed $j"
	job 20 1 research-rh7 $cmd
	# eval $cmd
done