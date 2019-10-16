clear all;close all;clc;

load('D:\MA_Manuscript_project\multimodality_adas_prediction_020319\data\multimodality_adas_prediction.mat')

all_mri_score = cmeas{1,2};
mri_derived_scores = cmeas{1,6};
mmse = cmeas{1,1};
avf_PET = cmeas{1,3};
amyloid_tau = cmeas{1,4};
fdg_PET = cmeas{1,5};
var_importance = [];

age;
apoe;
RID;
edu;

clear cmeas

X_complex_concat = [zscore(age) zscore(edu) zscore(apoe) zscore(mri_derived_scores) zscore(mmse) zscore(avf_PET) zscore(amyloid_tau) zscore(fdg_PET)];

regress_explorer(X_complex_concat, ADAS);
univariate_feeder;
bridge_results;


