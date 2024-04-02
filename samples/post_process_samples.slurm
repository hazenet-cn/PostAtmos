#!/bin/bash
#SBATCH --job-name=draw_PNG
#SBATCH --partition=64c512g
#SBATCH --mail-type=end
#SBATCH --time=10:00:00
#SBATCH --mail-user=
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=/Path/To/PostAtmos/samples/post_process.log  # This is your log path
#SBATCH --error=/Path/To/PostAtmos/samples/post_process.log # This is your log path
source /Path/To/.conda/envs/env_air/bin/activate # This is your python path
python -c "
import sys;
root_dir = '/Path/To/PostAtmos'# This is your PostAtmos path
sys.path.append(root_dir);
from post_process import Post_Process;
extracted_layers = [1] # # This is destination altitude layers. The first layer is mandatory, with the format like [1,5,10]
sites_includedprovinces = ['all'] # Administrative division I, for example, 'all'(need to be in lowercase), 'Shanghai', 'Jiangsu', 
sites_includedcities =['all'] # Administrative division II, for example, 'all', 'Nanjing', 'Hangzhou'
# Confirm the modification. The eight parameters are 
# processing type (WRF,CMAQ),
# start and end dates, path of the PostAtmos,
# paths of the simulation data to be verified. Multiple paths should be separated by semicolons. Multiple paths can have different grid settings.
# data filtering configurations
objProcess = Post_Process('CMAQ','2019-01-14','2019-01-14',\
 root_dir, root_dir+'/data/simulation',\
 extracted_layers,sites_includedprovinces,sites_includedcities);
objProcess.run();"
 