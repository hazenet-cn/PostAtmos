#!/bin/bash
#SBATCH --job-name=draw_PNG
#SBATCH --partition=64c512g
#SBATCH --mail-type=end
#SBATCH --time=10:00:00
#SBATCH --mail-user=
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=/Path/To/PostAtmos/samples/post_process.log  # Confirm the modification
#SBATCH --error=/Path/To/PostAtmos/samples/post_process.log # Confirm the modification
source /Path/To/.conda/envs/env_air/bin/activate # Confirm the modification
python -c "
import sys;
root_dir = '/Path/To/PostAtmos'# Confirm the modification
sys.path.append(root_dir);
from post_process import Post_Process;
extracted_layers = [1] # # Confirm the modification of altitude layers. The first layer is mandatory, with the format like [1,5,10]
sites_includedprovinces = ['all'] # Confirm the modification, for example, 'Shanghai', 'Jiangsu Province', all need to be in lowercase
sites_includedcities =['all']
# Confirm the modification. The eight parameters are 
# processing type (WRF,CMAQ_regularsites,CMAQ_supersites,EMISSION,adjoint_sensitivity,adjoint_semi_normalized_sensitivity),
# start and end dates, path of the PostAtmos,
# paths of the simulation data to be verified. Multiple paths should be separated by semicolons. Multiple paths can have different grid settings.
# data filtering configurations
objProcess = Post_Process('CMAQ_regularsites','2019-01-14','2019-01-14',\
 root_dir, root_dir+'/data',\
 extracted_layers,sites_includedprovinces,sites_includedcities);
objProcess.run();"
 