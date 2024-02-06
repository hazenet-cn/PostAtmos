#!/bin/bash
#SBATCH --job-name=draw_PNG
#SBATCH --partition=64c512g
#SBATCH --mail-type=end
#SBATCH --time=10:00:00
#SBATCH --mail-user=
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=/dssg/home/acct-esehazenet/hazenet-pg5/public_code/post_process.log  # 确认修改
#SBATCH --error=/dssg/home/acct-esehazenet/hazenet-pg5/public_code/post_process.log # 确认修改

/dssg/home/acct-esehazenet/share/.conda/envs/env_air/bin/python -c "
import sys;
sys.path.append('/dssg/home/acct-esehazenet/hazenet-pg5/public_code/Class_files@Hazenet');
from post_process import Post_Process;
# 确认修改,四个参数分别为处理类型（ WRF,CMAQ_BAU,EMISSION,adjoint_sensitivity,adjoint_semi_normalized_sensitivity)
# 开始和结束日期，数据源（多个数据源用;隔开，多个数据源需要是同一套网格）
extracted_layers = [1] # # 确认修改，高度层，第一层必选，格式如 [1,5,10]
sites_includedprovinces = ['all'] # 确认修改，例如 '上海市','江苏省', all需要小写
sites_includedcities =['all']
objProcess = Post_Process('CMAQ_BAU','2019-01-06','2019-01-31',\
'/dssg/home/acct-esehazenet/hazenet-pg5/CHINA_12km_Assimilation_2019-01-06_2019-01-31',\
 extracted_layers,sites_includedprovinces,sites_includedcities);
objProcess.run();"

 