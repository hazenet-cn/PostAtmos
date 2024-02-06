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
/environment/requirements.txt may be used to create a conda environment.

## Usage
Edit path in 模拟结果后处理作图示例文件test.slurm and run:
$ sbatch /dssg/home/acct-esehazenet/hazenet-pg5/PostAtmos/模拟结果后处理作图示例文件test.slurm
platform: linux-64