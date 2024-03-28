# PostAtmos

Post-processor "PostAtmos" from the paper “Rapid graphical validation of atmospheric modeling relies on high-performance computing” built in Python. By harnessing the parallel computational capabilities of high-performance computing platforms, it can successfully accelerate the validation process by approximately 165 times compared to traditional serial techniques.

<img src="https://github.com/hazenet-cn/PostAtmos/blob/main/docs/imgs/post_process.png"  width = "70%" height = "70%"/>


The config directory stores conda environment configuration file.  
The data directory stores a set of test data required to run the program,including simulation data (.nc files) and observation data 
(.csv files).  
The samples directory stores the slurm sample files for running PostAtmos, as well as a display of some of the results of running PostAtmos.  
The resources dicrectory stores the background and font files needed in running.

The program currently can recognise output files of CMAQ and WRF.

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
CREATE TABLE t_na_stations (
    station_code VARCHAR(255) PRIMARY KEY,
    station_name VARCHAR(255),
    station_province VARCHAR(255),
    station_city VARCHAR(255),
    station_lon FLOAT,
    station_lat FLOAT
);  

CREATE TABLE t_na_station_realtime (
    ID BIGINT(20) UNSIGNED AUTO_INCREMENT PRIMARY KEY,
    station_code VARCHAR(20),
    NO2 DOUBLE,
    pubtime TIMESTAMP
);  
```
Import CSV files in MySQL Client: take Mysql workbench as an example:  

Click the destination table  
<img src="https://github.com/hazenet-cn/PostAtmos/blob/main/docs/imgs/step_img1.png"  width = "70%" height = "70%"/>

Select the table path you want to import  
<img src="https://github.com/hazenet-cn/PostAtmos/blob/main/docs/imgs/step_img2.png"  width = "70%" height = "70%"/>

Select destination  
<img src="https://github.com/hazenet-cn/PostAtmos/blob/main/docs/imgs/step_img3.png"  width = "70%" height = "70%"/>

Set import configuration  
<img src="https://github.com/hazenet-cn/PostAtmos/blob/main/docs/imgs/step_img4.png"  width = "70%" height = "70%"/>

The obs_data.csv in the data/observation folder corresponds to the t_na_station_realtime table in the database and the obs_station.csv corresponds to the t_na_stations table in the database.

### Step 4. Download Gridded Global Relief Data
click https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/ETOPO2v2g/netCDF/ETOPO2v2g_f4_netCDF.zip
and move the file "ETOPO2v2g_f4.nc" into resources directory.

### Step 5. Edit sql database information in post_process.py(row 41)
```python
airdb_engine = sqlalchemy.create_engine("dialect+driver://username:password@host:port/database")  #observation station and data
```

### Step 6. Job submission using Slurm
Edit the configurations in /PostAtmos/samples/post_process_samples.slurm that are labeled in the sample file and run:  
```shell
$ sbatch /PostAtmos/samples/post_process_samples.slurm
```
