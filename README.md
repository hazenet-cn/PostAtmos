# PostNC: Post-Processor of NetCDF-Based Model Validation

PostNC is a post-processor built using Python, designed for the validation of numerical models with output in NetCDF. It's from the paper "Rapid graphical validation of NetCDF-based modeling relies on high-performance computing: case study of atmospheric models". The following diagram illustrates the overall workflow of PostNC.  

<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/post_process.png"  width = "70%" height = "70%"/>

Currently, the program can recognize output files from CMAQ and WRF models. To validate other files output in netCDF format, manual modifications to the code's data processing procedures are necessary, such as changing variable names to adapt to different models.  

By leveraging the parallel computing capabilities of high-performance computing platforms, PostNC can accelerate the validation process for CMAQ case by approximately 165 times compared to traditional serial techniques.  

## Directory Structure
- config: Stores the conda environment configuration file.  
- data: Contains test data required to run the program, including simulation data (.nc files) and observation data (.csv files).  
- samples: Stores the Slurm sample files for running PostNC and some results of running PostNC.  
- resources: Stores the background and font files needed in running.  

## Inputs 
- Output from WRF or CMAQ runs;  
- Observations stored in a database.  

Sample input files for testing can be found in the /PostNC/data folder.  

## Outputs
- CSV folder: Tables storing matched simulated values and observations. 

| source_index | SiteCode | TimePoint           | factor_index | PM25        | PM25_obs |
|--------------|----------|---------------------|--------------|-------------|----------|
| 0            | 1001A    | 2019-01-06 08:00:00 | 4            | 57.8730850  | 38.0     |
| 0            | 1002A    | 2019-01-06 08:00:00 | 4            | 33.8960838  | 41.0     |
| 0            | 1003A    | 2019-01-06 08:00:00 | 4            | 57.8730850  | 28.0     |
| 0            | 1004A    | 2019-01-06 08:00:00 | 4            | 57.8730850  | 34.0     |
| 0            | 1005A    | 2019-01-06 08:00:00 | 4            | 72.5394439  | 37.0     |
| 0            | 1006A    | 2019-01-06 08:00:00 | 4            | 57.8730850  | 34.0     |
| 0            | 1007A    | 2019-01-06 08:00:00 | 4            | 57.8730850  | 25.0     |
| 0            | 1008A    | 2019-01-06 08:00:00 | 4            | 51.3097457  | 23.0     |
| 0            | 1009A    | 2019-01-06 08:00:00 | 4            | 43.6032104  | 16.0     |
| 0            | 1011A    | 2019-01-06 08:00:00 | 4            | 54.7686958  | 20.0     |
| 0            | 1012A    | 2019-01-06 08:00:00 | 4            | 51.1413612  | 47.0     |
| 0            | 1015A    | 2019-01-06 08:00:00 | 4            | 76.5963745  | 52.0     |
| 0            | 1017A    | 2019-01-06 08:00:00 | 4            | 79.1437759  | 43.0     |
| 0            | 1018A    | 2019-01-06 08:00:00 | 4            | 79.1437759  | 48.0     |

- MP4 folder: Videos depicting offline spatial distribution of data over time(with background).  

<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/CO.gif"  width = "70%" height = "70%"/>

- PNG folder: Pictures showing spatial distribution of data without background (can be displayed on an online map).  

<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/web_spatial_graph.png"  width = "70%" height = "70%"/>

- Scatter folder: scatter plots of all observations and simulated values of target species within the simulation domain.

<img src="https://github.com/hazenet-cn/PostNC/blob/main/samples/Sample_output/Scatter/Scatter_NO2_CMAQ_regularsites_CHINA_12km_Assimilation_2019-01-06_2019-01-31.png"  width = "70%" height = "70%"/>

- Timeseries folder: time series graphs for each site within the simulation domain.  

<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/timeseries.png"  width = "70%" height = "70%"/>

These outputs can be found in the /PostNC/samples/Sample_output folder.  

## Requirements
+ Linux based HPCs
+ MySQL Community Server
+ Intel oneAPI
+ Python3
+ Git

## Installation and Usage

### Step 1. Clone PostNC repository
```
git clone https://github.com/hazenet-cn/PostNC.git
```

### Step 2. Create a conda environment
/PostNC/config/env_air.yml may be used to create a conda environment.  
```shell
$ conda env create --file env_air.yml --name target_env  
```

### Step 3. Install Intel oneAPI and Intel MPI Library  
Visit the Intel oneAPI Download Page at https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html and choose the Intel® HPC Toolkit suitable for your operating system.  
For installation instructions, refer to the Intel oneAPI HPC Toolkit Getting Started Guide for Linux at https://www.intel.com/content/www/us/en/docs/oneapi-hpc-toolkit/get-started-guide-linux/2024-0/overview.html.  

### Step 4. Prepare observation data in database
You can either use your own database of observations or the sampsle observations provided in the data catalog for testing. The /PostNC/data/observation folder contains CSV tables of observation data that you can import into a MySQL database as follows:  
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
To import CSV files into a MySQL database, you can use a MySQL client such as MySQL Workbench. In the /PostNC/data/observation folder, the obs_data.csv file corresponds to the t_na_station_realtime table in the database, and the obs_station.csv file corresponds to the t_na_stations table in the database. Here's how you can import these files using MySQL Workbench:  

Click the destination table  
<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/step_img1.png"  width = "70%" height = "70%"/>

Select the table path you want to import  
<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/step_img2.png"  width = "70%" height = "70%"/>

Select destination  
<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/step_img3.png"  width = "70%" height = "70%"/>

Set import configuration  
<img src="https://github.com/hazenet-cn/PostNC/blob/main/docs/imgs/step_img4.png"  width = "70%" height = "70%"/>

### Note:
To minimize complications, the database table names are imported as t_na_stations for station information and t_na_station_realtime for station observations. If you use your own database tables, you will need to update the database names in the code accordingly(row 124, 1301 in post_process.py); otherwise, they will not be recognized.  

### Step 5. Download Gridded Global Relief Data
click https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/ETOPO2v2g/netCDF/ETOPO2v2g_f4_netCDF.zip
and move the file "ETOPO2v2g_f4.nc" into resources directory.

### Step 6. Edit sql database and python information in post_process.py(row 41,43)
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

### Step 7. Job submission using Slurm
Edit the configurations in /PostNC/samples/post_process_samples.slurm that are labeled in the sample file
```
# Confirm the modification. The eight parameters are 
# processing type(self.data_type in post_process.py) (WRF or CMAQ),
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
- You can modify variables in post_process.py to customize the validation process, such as target species, output directory, spatial distribution figure background, and image resolution.  

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

- To optimize the Post_Process class for supporting additional models beyond CMAQ and WRF, follow these steps:

1. Introduce a new data_type to represent the new model.
2. In the __init__ method, add a branch to identify the variables and attributes specific to the new model's output.
3. In the extract_2dsimu method, add a branch to handle the identification of the new model's output file and the unit conversion of its variables.
4. Ensure that the draw_offline_PNG, draw_online_PNG, and submit_sites methods also include the same branch to support the new model.  
For ease of implementation, it is recommended to replicate the existing branches in the aforementioned methods and modify the relevant fields accordingly for the new model.

## Contact and Support
- Please report any bugs through 'Issues'.  
- Contact to us:  
Zhen Cheng: chengz88@sjtu.edu.cn  
Jiayin Liu: liujiayin@sjtu.edu.cn  