# PostAtmos

Post-processor PostAtmos from the paper “Rapid graphical validation of atmospheric modeling relies on high-performance computing” is built in Python. By harnessing the parallel computational capabilities of high-performance computing platforms, it can successfully accelerate the validation process by approximately 165 times compared to traditional serial techniques.

<img src="https://github.com/hazenet-cn/PostAtmos/blob/main/docs/imgs/post_process.png"  width = "70%" height = "70%"/>

The data directory stores a set of test data required to run the program,  
including simulation data (.nc files) and observation data (.csv files).  


 and then modify the database connection in post_process.py  
The program currently can recognise CMAQ, WRF result output files.

## Requirements
+ MPI
+ Pymysql
+ Pandas  
+ Numpy
+ Matplotlib

## Usage

### Step 1. Create a conda environment
/PostAtmos/config/env_air.yml may be used to create a conda environment.  
```shell
$ conda env create --file env_air.yml --name target_env
```

### Step 2. Import the csv data into your database
/PostAtmos/environment/env_air.yml may be used to create a conda environment.



Edit path in /PostAtmos/post_process_samples.slurm and run:  

$ sbatch /PostAtmos/post_process_samples.slurm

or run /PostAtmos/post_process_samples.ipynb in Jupyter Notebook.

platform: linux-64
