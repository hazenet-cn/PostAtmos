# PostNC

Post-processor PostNC from the paper “Rapid graphical validation of NetCDF-based modeling relies on high-performance computing: case study of atmospheric models” built in Python. By harnessing the parallel computational capabilities of high-performance computing platforms, it can successfully accelerate the validation process by approximately 165 times compared to traditional serial techniques.

<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/post_process.png"  width = "70%" height = "70%"/>

Software Input:  
Output from WRF or CMAQ runs;  
Observations stored in a database.  

Software Output:  
CSV folder-matched simulated values and observations  
PNG folder - spatial distribution without background (can be displayed on an online map)  
Scatter folder - scatter plot of all observations and simulated values within the simulation range  
Timeseries folder - time series graphs for each site  

The config directory stores conda environment configuration file.  
The data directory stores a set of test data required to run the program,including simulation data (.nc files) and observation data 
(.csv files).  
The samples directory stores the slurm sample files for running PostNC, as well as a display of some of the results of running PostNC.  
The resources dicrectory stores the background and font files needed in running.

The program currently can recognise output files of CMAQ and WRF.

## Requirements
+ Linux based HPCs
+ MySQL Community Server
+ OpenMPI
+ Python 

## Usage

### Step 1. Create a conda environment
/PostNC/config/env_air.yml may be used to create a conda environment.  
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
<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/step_img1.png"  width = "70%" height = "70%"/>

Select the table path you want to import  
<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/step_img2.png"  width = "70%" height = "70%"/>

Select destination  
<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/step_img3.png"  width = "70%" height = "70%"/>

Set import configuration  
<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/step_img4.png"  width = "70%" height = "70%"/>

The obs_data.csv in the data/observation folder corresponds to the t_na_station_realtime table in the database and the obs_station.csv corresponds to the t_na_stations table in the database.

### Step 4. Download Gridded Global Relief Data
click https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/ETOPO2v2g/netCDF/ETOPO2v2g_f4_netCDF.zip
and move the file "ETOPO2v2g_f4.nc" into resources directory.

### Step 5. Edit sql database and python information in post_process.py(row 41,43)
```python
engine_str = "dialect+driver://username:password@host:port/database" #observation station and data  
python_path = "/Path/To/.conda/envs/env_air/bin"
```
dialect: the type of database, e.g. mysql, postgresql, sqlite, etc.  
driver: the database driver to use, e.g. pymysql, psycopg2, sqlite3. For some databases, the driver is optional.  
username: The user name to connect to the database.  
password: The user's password.  
host: Host name or IP address of the database server.  
port: The port number of the database server, the default MySQL port is 3306.  
database: Name of the database to connect to.  

### Step 6. Job submission using Slurm
Edit the configurations in /PostNC/samples/post_process_samples.slurm that are labeled in the sample file
```
# Confirm the modification. The eight parameters are 
# processing type (WRF or CMAQ),
# start and end dates, path of the PostNC,
# paths of the simulation data to be verified. Multiple paths should be separated by semicolons. Multiple paths can have different grid settings.
# data filtering configurations
objProcess = Post_Process('CMAQ','2019-01-14','2019-01-14',\
 root_dir, root_dir+'/data/simulation',\
 extracted_layers,sites_includedprovinces,sites_includedcities);
```

and run:  
```shell

$ sbatch /PostNC/samples/post_process_samples.slurm
```

If the following error occurs:  
sbatch: error: Batch script contains DOS line breaks (\r\n)  
sbatch: error: instead of expected UNIX line breaks (\n).  
run:  
```shell
$ dos2unix /PostNC/samples/post_process_samples.slurm  
```
## Customized Guidelines
Variables recommended for modification in post_process.py and their meaning:
| Number | Variable             | Default Value           | Explanation                                   |
|:------:|:--------------------:|:-----------------------:|:---------------------------------------------:|
|   1    | `self.factors`       | list in source codes    | list of target species                        |
|   2    | `self.output_dir`    | root_dir/PostProcess_xx | Validation result output path                 |
|   3    | `dem_file`           | "ETOPO2v2g_f4.nc"       | Spatial distribution figure background        |
|   4    | `self.lcc_proj`      | ccrs.LambertConformal   | Spatial distribution plot projection          |
|   5    | `df`                 | "asia_city.csv"         | Urban points labelled on the background       |
|   6    | `conc_color`         | list in source codes    | Colorbar color                                |
|   7    | `c_levels`           | [minvalue, maxvalue]    | Color mapping and range                       |
|   8    | `self.lon_tickrange` | numpy.arange(80,140,10) | Range of longitude to be ticked               |
|   9    | `self.lat_tickrange` | numpy.arange(10,60,10)  | Range of latitude to be ticked                |
|  10    | `dpi`                | 600                     | Image resolution                              |

## Contact and Support
- Please report any bugs through 'Issues'.  

- Contact to us:  
chengz88@sjtu.edu.cn  (Zhen Cheng)  
liujiayin@sjtu.edu.cn  (Jiayin Liu)  