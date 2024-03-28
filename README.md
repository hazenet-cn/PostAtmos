# PostAtmos

Post-processor "PostAtmos" from the paper “Rapid graphical validation of atmospheric modeling relies on high-performance computing” built in Python. By harnessing the parallel computational capabilities of high-performance computing platforms, it can successfully accelerate the validation process by approximately 165 times compared to traditional serial techniques.

<img src="https://github.com/hazenet-cn/PostAtmos/blob/main/docs/imgs/post_process.png"  width = "70%" height = "70%"/>


The config directory stores conda environment configuration file.  
The data directory stores a set of test data required to run the program,including simulation data (.nc files) and observation data 
(.csv files).  
The samples directory stores the slurm sample files for running PostAtmos, as well as a display of some of the results of running PostAtmos.
The resources dicrectory stores the background and font files needed in running.


 and then modify the database connection in post_process.py  
The program currently can recognise CMAQ, WRF result output files.

## Requirements
+ Linux based HPCs
+ MySQL Community Server
+ OpenMPI
+ Python 

## Usage

### Step 1. Create a conda environment
/PostAtmos/config/env_air.yml may be used to create a conda environment.  
```shell
$ conda env create --file env_air.yml --name target_env  
```

### Step 2. Install Open MPI  
Refer to https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html  

### Step 3. Import the csv data into your database
```sql
-- create database in MySQL Client
CREATE DATABASE observation;

-- use database
USE observation;

-- create table for observation
CREATE TABLE Supersite_Sites (
    station_code VARCHAR(255) PRIMARY KEY,
    station_name VARCHAR(255),
    station_province VARCHAR(255),
    station_city VARCHAR(255),
    station_lon FLOAT,
    station_lat FLOAT
);  

CREATE TABLE YourTableName (
    ID BIGINT(20) UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    station_code VARCHAR(20),
    SO2 DOUBLE,
    NO2 DOUBLE,
    CO DOUBLE,
    O3 DOUBLE,
    PM10 DOUBLE,
    PM25 DOUBLE,
    pubtime TIMESTAMP,
    ctime DATETIME
);  

-- import CSV files
-- 注意：这需要在MySQL客户端或命令行中执行，不是标准SQL语句
LOAD DATA INFILE '/path/to/your/file.csv'
INTO TABLE obs_PM25
FIELDS TERMINATED BY ',' -- 或你的分隔符
ENCLOSED BY '"' -- 如果字段被引号包围
LINES TERMINATED BY '\n'
IGNORE 1 ROWS; -- 如果CSV文件包含标题行
```

### Step 4. Download Gridded Global Relief Data
click https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/ETOPO2v2g/netCDF/ETOPO2v2g_f4_netCDF.zip
and move the file "ETOPO2v2g_f4.nc" into resources directory.

### Step 5. Job submission using Slurm
Edit the configurations in /PostAtmos/samples/post_process_samples.slurm that are labeled in the sample file and run:  
```shell
$ sbatch /PostAtmos/samples/post_process_samples.slurm
```
