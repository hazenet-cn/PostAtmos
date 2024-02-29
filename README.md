# PostAtmos

Post-processor PostAtmos from the paper “Rapid graphical verification of atmospheric modeling relies on high-performance computing” in Python is built for rapid, high-quality and flexible verification of atmospheric model simulation results using high-performance computing.

<img src="https://github.com/hazenet-cn/PostAtmos/blob/main/imgs/post_process.png"  width = "70%" height = "70%"/>

## Data
The data directory stores the background and font files needed to run the program. And the program currently can recognise CMAQ, WRF result output files.

## Requirements
+ MPI
+ Pymysql
+ Pandas  
+ Numpy
+ Matplotlib

## Usage
/PostAtmos/environment/env_air.yml may be used to create a conda environment.

$ conda env create --file env_air.yml --name target_env

Edit path in /PostAtmos/post_process_samples.slurm and run:  

$ sbatch /PostAtmos/post_process_samples.slurm

or run /PostAtmos/post_process_samples.ipynb in Jupyter Notebook.

platform: linux-64
