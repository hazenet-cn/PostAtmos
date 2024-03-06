#!/bin/bash
#SBATCH --job-name=draw_PNG
#SBATCH --partition=64c512g
#SBATCH --mail-type=end
#SBATCH --time=10:00:00
#SBATCH --mail-user=
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=/dssg/home/acct-esehazenet/hazenet-pg5/public_code/post_process.log  # Confirm the modification
#SBATCH --error=/dssg/home/acct-esehazenet/hazenet-pg5/public_code/post_process.log # Confirm the modification

/dssg/home/acct-esehazenet/share/.conda/envs/env_air/bin/python -c "
import sys;
sys.path.append('/Path/To/PostAtmos');
from post_process import Post_Process;
extracted_layers = [1] # # Confirm the modification of altitude layers. The first layer is mandatory, with the format like [1,5,10]
sites_includedprovinces = ['all'] # Confirm the modification, for example, 'Shanghai', 'Jiangsu Province', all need to be in lowercase
sites_includedcities =['all']
# Confirm the modification. The four parameters are processing type (WRF,CMAQ_regularsites,CMAQ_supersites,EMISSION,adjoint_sensitivity,adjoint_semi_normalized_sensitivity),
# start and end dates, data source paths. Multiple paths should be separated by semicolons. Multiple paths can have different grid settings.
objProcess = Post_Process('CMAQ_regularsites','2019-07-01','2019-07-31',\
'/Path/To/Forward_China_12km_2019-06-25_2019-07-31_BAU',\
 extracted_layers,sites_includedprovinces,sites_includedcities);
objProcess.run()
 