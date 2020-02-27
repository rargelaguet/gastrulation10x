
###################
## Load settings ##
###################

if echo $HOSTNAME | grep -q "ricard"; then
	source /Users/ricard/gastrulation10x/mofa/run/settings.sh
elif echo $HOSTNAME | grep -q "ebi"; then
	source /homes/ricard/gastrulation10x/mofa/run/settings.sh
else
	echo "Computer not recognised"; exit
fi

downsampling_fraction=($(seq 0.00 0.025 0.25))

#########
## Run ##
#########

for i in $(seq 1 $downsampling_fraction); do
	outfile="$output_folder/model_downsample_$i.hdf5"
	cmd="python run.py --input_folder $input_folder --outfile $outfile --downsampling_fraction $downsampling_fraction --factors $num_factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq"
	job 15 2 research-rh7 $cmd
	eval $cmd
done