import scipy.io as sio
import numpy as np
import os
import scipy.stats as stats
from scipy.io import savemat, loadmat
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

cwd = os.getcwd()
# make your own data dir if not exist
if not os.path.exists(cwd + '\PANSS_All'):
    os.makedirs(cwd + '\PANSS_All')
data_dir = cwd + '\PANSS_All'

def glm(unrelated_subj_num, fc_type):
	name1 = '\PANSS_All\FG1_P1'
	if not os.path.exists(cwd + name1):
		os.makedirs(cwd + name1)
	glm_dir = cwd + name1 + '\\'

	# get the age, sex... info
	age = np.genfromtxt(data_dir + r'\age.txt')
	sex = np.genfromtxt(data_dir + r'\sex.txt')
	cog = np.genfromtxt(data_dir + r'\P_Scores.txt')
	motion = np.genfromtxt(data_dir + r'\FD.txt')

	# center all variables
	cen_age = age - np.mean(age)
	cen_cog = cog - np.mean(cog)
	cen_motion = motion - np.mean(motion)

	# compute the interactions
	sexcog = sex * cen_cog
	agecog = cen_age * cen_cog

	# concatenate all the covariates
	x = np.concatenate((cen_age.reshape(-1,1), sex.reshape(-1,1), cen_cog.reshape(-1,1), cen_motion.reshape(-1,1), agecog.reshape(-1,1), sexcog.reshape(-1,1)), axis=1)

	regionalcp = loadmat(data_dir + r'\data1_1.mat')['data1_1']
	
	roi_number = regionalcp.shape[1]

	# regionalcp_z = np.arctanh(regionalcp)

	# GLM
	exog = sm.add_constant(x)
	endog = regionalcp

	p_value = []
	params = []
	conf_int = []
	for i in range(roi_number):
		model = sm.GLM(endog[:,i], exog, sm.families.Gaussian())
		res = model.fit()
		p_value.append(res.pvalues)
		params.append(res.params)
		conf_int.append(res.conf_int())
	p_value = np.array(p_value)
	params = np.array(params)

	savemat(glm_dir + 'GLM_parameters.mat', {'parameters': params})

	for i in range(1,exog.shape[1]):
		p_corrected = multipletests(p_value[:,i], alpha = 0.001, method = 'fdr_bh')

		allroi_association = np.zeros(roi_number)
		sigroi_association = np.zeros(roi_number)
		for j in range(roi_number):
			if params[j,i] > 0:
				allroi_association[j] = (-np.log(p_corrected[1]))[j]
			else:
				allroi_association[j] = -(-np.log(p_corrected[1]))[j]

		savemat(glm_dir + 'GLMvar' + str(i) + fc_type + '_allroi.mat', {'allroi_association': allroi_association})

		for j in np.where(p_corrected[0] == True)[0]:
			if params[j,i] > 0:
				sigroi_association[j] = (-np.log(p_corrected[1]))[j]
			else:
				sigroi_association[j] = -(-np.log(p_corrected[1]))[j]

		savemat(glm_dir + 'GLMvar' + str(i) + fc_type + '_sigroi.mat', {'sigroi_association': sigroi_association})

		mdic = {'conf_int': conf_int}
		savemat(glm_dir + 'GLM_' + fc_type + '_confint.mat', mdic)
