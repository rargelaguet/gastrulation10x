
###################
## Load settings ##
###################

if echo $HOSTNAME | grep -q "ricard"; then
	source /Users/ricard/gastrulation10x/mofa/run/settings.sh
elif echo $HOSTNAME | grep -q "ebi"; then
	source /homes/ricard/gastrulation10x/mofa/run/settings.sh
	output_folder='/hps/nobackup2/research/stegle/users/ricard/gastrulation10x_mofa/hdf5/downsampling_v3'
else
	echo "Computer not recognised"; exit
fi

downsampling_fraction=($(seq 0.00 0.025 0.25))
maxiter=100
start_elbo=999
num_factors=10

#########
## Run ##
#########

for i in "${downsampling_fraction[@]}"; do
	outfile="$output_folder/model_downsample_${i}.hdf5"
	cmd="python run.py --input_folder $input_folder --outfile $outfile --downsampling_fraction $i --factors $num_factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq"
	job 7 1 research-rh7 $cmd
	# eval $cmd
done