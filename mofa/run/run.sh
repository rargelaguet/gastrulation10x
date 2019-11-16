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
output_folder="/hps/nobackup2/research/stegle/users/ricard/gastrulation10x/mofa/hdf5"

# input_folder="/Users/ricard/data/gastrulation10x/mofa/data"
# output_folder="/Users/ricard/data/gastrulation10x/mofa/hdf5"

# Number of trials
ntrials=1

# Number of factors
factors=10

# ELBO Convergence settings
convergence_mode="medium"
maxiter=5000
elbofreq=5
start_elbo=500

# Stochastic settings
# batch_size=( 0.25 0.5 1.0 )
# tau=( 0.05 0.15 0.25 0.5 0.75 1.0 )
# forgetting_rate=( 0.5 0.25 0 )

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

for j in $(seq 1 $ntrials); do
	# outfile="$output_folder/E6.5-E7.25_$j.hdf5"
	outfile="$output_folder/model_$j.hdf5"
	cmd="python run.py --input_folder $input_folder --outfile $outfile --factors $factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq --seed $j --verbose"
	gpujob 15 1 research-rh7 $cmd
	# eval $cmd
done