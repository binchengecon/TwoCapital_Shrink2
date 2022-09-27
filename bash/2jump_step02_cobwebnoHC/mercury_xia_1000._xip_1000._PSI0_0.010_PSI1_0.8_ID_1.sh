#! /bin/bash


######## login 
#SBATCH --job-name=post15_4
#SBATCH --output=./job-outs/2jump_step02_cobwebnoHC/xia_1000._xip_1000._PSI0_0.010_PSI1_0.8/mercury_post_1.out
#SBATCH --error=./job-outs/2jump_step02_cobwebnoHC/xia_1000._xip_1000._PSI0_0.010_PSI1_0.8/mercury_post_1.err


#SBATCH --account=pi-lhansen
#SBATCH --partition=standard
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=7-00:00:00

####### load modules
module load python/booth/3.8/3.8.5  gcc/9.2.0

echo "$SLURM_JOB_NAME"

echo "Program starts $(date)"

python3 /home/bcheng4/TwoCapital_Shrink/abatement/postdamage_2jump_cobweb_smartguess.py --num_gamma 3 --xi_a 1000. --xi_g 1000.  --epsilonarr 0.1 0.01  --fractionarr 0.1 0.01   --maxiterarr 60000 80000  --id 1 --psi_0 0.010 --psi_1 0.8 --name 2jump_step02_cobwebnoHC --hXarr 0.2 0.2 0.2 --Xminarr 4.00 0.0 -5.5 0.0 --Xmaxarr 9.00 4.0 0.0 3.0

echo "Program ends $(date)"

