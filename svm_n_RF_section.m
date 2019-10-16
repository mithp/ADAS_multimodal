%Code similar to PLSR works with samme folds to enable comparisions
%% Performing SVM regression and Random Forest

Mdl_SVM = fitrsvm(x_model,y_model,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','UseParallel',true));
Mdl_RF = fitrensemble(x_model,y_model,'OptimizeHyperparameters','auto',...
    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName',...
    'expected-improvement-plus','UseParallel',true));

yfit_model_svm = predict(Mdl_SVM,x_model);
yfit_val_svm = predict(Mdl_SVM,x_val);

yfit_model_rf = predict(Mdl_RF,x_model);
yfit_val_rf = predict(Mdl_RF,x_val);

[r2_model_folds_svm(fold) RMSECV_t_svm] = rsquare(y_model,yfit_model_svm);
[r2_val_folds_svm(fold) RMSEP_t_svm] = rsquare(y_val,yfit_val_svm);

RMSECV_fold_ssvm(fold) = RMSECV_t_svm/range(y_model)*100;
RMSEP_folds_svm(fold) = RMSEP_t_svm/range(y_val)*100;

[r2_model_folds_rf(fold) RMSECV_t_rf] = rsquare(y_model,yfit_model_rf);
[r2_val_folds_rf(fold) RMSEP_t_rf] = rsquare(y_val,yfit_val_rf);

RMSECV_folds_rf(fold) = RMSECV_t_rf/range(y_model)*100;
RMSEP_folds_rf(fold) = RMSEP_t_rf/range(y_val)*100;

close all;
clc;