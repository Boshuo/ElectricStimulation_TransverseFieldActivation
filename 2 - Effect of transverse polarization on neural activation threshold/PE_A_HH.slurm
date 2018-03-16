#!/bin/bash
#
#SBATCH --array=1-475%150
#SBATCH --job-name=PE_A_HH
#SBATCH --error=/PE_A_HH/Err_Out/%a.err
#SBATCH --output=/PE_A_HH/Err_Out/%a.out
#SBATCH --mem=4G

uname -n 1>&2

i=$(($SLURM_ARRAY_TASK_ID))

/opt/apps/matlabR2016a/bin/matlab -r "mod_prmtr.model_name = 'PE_A_HH'; mod_prmtr.id = ${i}; out_ctrl.if_save_data = 1; out_ctrl.if_write_log = 1; out_ctrl.if_plot = 1; main_ES_HH( mod_prmtr, out_ctrl ); quit;"