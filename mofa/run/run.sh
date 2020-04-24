
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

maxiter=50
start_elbo=999


#########
## Run ##
#########

outfile="$output_folder/test.hdf5"
cmd="python run.py --input_folder $input_folder --outfile $outfile --factors $num_factors --iterations $maxiter --convergence_mode $convergence_mode --start_elbo $start_elbo --elbo_freq $elbofreq"
# gpujob 15 1 research-rh7 $cmd
# jobid=$(sbatch -p gpu --gres gpu:1 --exclude=gpu8,gpu9 $cmd | awk '{print $NF}')
eval $cmd
