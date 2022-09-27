#! /bin/bash


######## login 
#SBATCH --job-name=1e4post_18
#SBATCH --output=./job-outs/rep_20dmgparal/xia_10000._xip_10000._PSI0_0.005_PSI1_0.5/mercury_post_18.out
#SBATCH --error=./job-outs/rep_20dmgparal/xia_10000._xip_10000._PSI0_0.005_PSI1_0.5/mercury_post_18.err


#SBATCH --account=pi-lhansen
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00

####### load modules
module load python/booth/3.8/3.8.5  gcc/9.2.0

echo "$SLURM_JOB_NAME"

echo "Program starts $(date)"

python3 /home/bcheng4/TwoCapital_Shrink/abatement/postdamage_spe_xi_psi_gammalist_name2.py --num_gamma 20 --xi_a 10000. --xi_g 10000.  --epsilonarr 0.1 0.1 0.1  --fractionarr 0.1 0.01 0.01   --maxiterarr 20000 20000 20000  --id 18 --psi_0 0.005 --psi_1 0.5 --name rep_20dmgparal --hK 0.1 --hY  	0.1 --hL  0.1

echo "Program ends $(date)"

