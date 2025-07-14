# this is the demo for how to run the code
import scipy.io as sio
import numpy as np
import os
from scipy.io import savemat
import scipy.stats as stats
from GLM import generate_covariates, glm

cwd = os.getcwd()
# make your own data dir if not exist
if not os.path.exists(cwd + '\PANSS_All'):
    os.makedirs(cwd + '\PANSS_All')
data_dir = cwd + '\PANSS_All'

name1 = '\PANSS_All\FG1_P'
if not os.path.exists(cwd + name1):
    	os.makedirs(cwd + name1)
glm_dir = cwd + name1

# set the number of the unrelated subjects
# set the number of brain regions in your atlas
unrelated_subj_num = 1118
region_num = 210

glm(unrelated_subj_num, fc_type)



