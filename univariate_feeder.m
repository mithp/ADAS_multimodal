
%Code for running the codes on the fold from Multivariate example
ccc;
load('D:\MA_Manuscript_project\Public_code\Data\all_runs.mat')

model_complex = model_timeline_wise_data;
eval_complex = eval_timeline_wise_data;

clear model_timeline_wise_data eval_timeline_wise_data

%extracting the data sizes for univariate modelling
time_line = size(model_complex,1);
number_of_runs = size(model_complex{1,1},1);
K_folds = size(model_complex{1,1}{1,1},1);

for time_line_counter = 1:time_line
    
    
    for number_of_runs_counter = 1: number_of_runs
        
        %model_complex_per_run =
        
        for K_fold_counter = 1:K_folds
            
            model_predictors = model_complex{time_line_counter,1}{number_of_runs_counter,1}{K_fold_counter,1};
            model_response = model_complex{time_line_counter,1}{number_of_runs_counter,1}{K_fold_counter,2};
            eval_predictors = eval_complex{time_line_counter,1}{number_of_runs_counter,1}{K_fold_counter,1};
            eval_reponse = eval_complex{time_line_counter,1}{number_of_runs_counter,1}{K_fold_counter,2};
            
            %to select individual modalities lets hardcode the indices so
            %we can index and send only portion of predictors to regression
            %modelling, thus we have same folds but different modality all
            %at once. First coloum is age, then edu, apoe etc
            %modality =[0,1,2,3,12,17,20,25];
            start_index = [0,1,2,3,12,17,20];
            end_index =   [1,2,3,12,17,20,25];
            
            for modality_indexer = 5 : size(start_index,2)
                
                %%SEction for regression
                start_index_ind = start_index(modality_indexer)+1; % increment by 1 to avoid copying previous modality
                end_index_ind = end_index(modality_indexer);
                
                % The portion for modelling
                x_model = model_predictors(:,start_index_ind:end_index_ind);
                y_model  = model_response;
                
                %Let's keep this portion for testing
                x_val  = eval_predictors(:,start_index_ind:end_index_ind);
                y_val  = eval_reponse;
                
                % initialize maxinmum number of PLS components
                if size(x_model,2) <= 15
                    ncompmax = size(x_model,2); %if there are fewer variables than 15
                else
                    ncompmax = 15;
                end
                
                %4 key components to assess the performance of the model
                r2train = zeros(ncompmax,1);
                r2test = zeros(ncompmax,1);
                rmsecvC = zeros(ncompmax,1);
                rmsepC = zeros(ncompmax,1);
                Error_Percent = zeros(ncompmax,1); %not in use at the momemnt
                ncomp = 0;
                
                %To find the optimal number of components lets further divide
                %the model data so the validation data is truely untouched
                
                pls_divide_count = size(y_model,1);
                [trainInd,testInd,~] = dividerand(pls_divide_count,0.80,0.20,0);
                
                % Inital model is based on this data
                x_train = x_model(trainInd,:);
                y_train  = y_model(trainInd,:);
                
                %Used to select the number of parameters
                x_test  = x_model(testInd,:);
                y_test  = y_model(testInd,:);
                
                
                clear r2_train r2_test rmse_train rmse_test
                
                %% Simple Partial least squares regression model (with 1-10 components)
                for ncomp=1:ncompmax
                    
                    [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
                    yfit_PLS_train = [ones(size(x_train,1),1) x_train]*betaPLS;
                    yfit_PLS_test = [ones(size(x_test,1),1) x_test]*betaPLS;
                    [r2_train(ncomp) rmse_train(ncomp)] = rsquare(y_train,yfit_PLS_train);
                    [r2_test(ncomp) rmse_test(ncomp)] = rsquare(y_test,yfit_PLS_test);
                    
                end
                
                %Finding the optimal number of components for this split
                number_of_components = find(rmse_test == min(rmse_test));
                
                %lets evaluate the validation set to find out the truth
                [~,~,~,~,betaPLS_model,~,~, ~]= plsregress(x_model,y_model,number_of_components);%,'CV',number_of_components);
                
                yfit_PLS_model = [ones(size(x_model,1),1) x_model] * betaPLS_model; %model predition training data
                yfit_PLS_val = [ones(size(x_val,1),1) x_val]*betaPLS_model; %model predition testing data
                
                [r2_model_modality(modality_indexer) RMSECV_t] = rsquare(y_model,yfit_PLS_model);
                [r2_val_modality(modality_indexer) RMSEP_t] = rsquare(y_val,yfit_PLS_val);
                
                %figure;scatter(y_val,yfit_PLS_val);
                %pause(0.2)
                
                %need to add code to save the training, testing and prediction
                %data: this is still pending
                
                N_fold(modality_indexer) = number_of_components;
                RMSECV_modality(modality_indexer) = RMSECV_t/range(y_model)*100;
                RMSEP_modality(modality_indexer) = RMSEP_t/range(y_val)*100;
                %close all
                
                %The code for SVM and RF in a seerate file to keeps things
                %simple
                svm_n_RF__univariate_section;
                
                results_per_modality_model(:,modality_indexer) = {y_model,yfit_PLS_model,yfit_model_svm,yfit_model_rf} ;
                results_per_modality_eval(:,modality_indexer) = {y_val,yfit_PLS_val,yfit_val_svm,yfit_val_rf} ;
                
                
                %%REgression section ends
                
            end
            
            
            
            collate_R2_RMSECV = {r2_model_modality,r2_val_modality,RMSECV_modality,RMSEP_modality,...
                r2_model_modality_svm,r2_val_modality_svm,RMSECV_modality_ssvm,RMSEP_modality_svm,...
                r2_model_modality_rf,r2_val_modality_rf,RMSECV_modality_rf,RMSEP_modality_rf};
            
            results_fold_model (:, K_fold_counter) = {results_per_modality_model};
            results_fold_eval (:, K_fold_counter) = {results_per_modality_eval,collate_R2_RMSECV};
            
            clear results_per_modality_model results_per_modality_eval collate_R2_RMSECV
        end
        
        runs_model(:, number_of_runs_counter) = {results_fold_model};
        runs_eval(:, number_of_runs_counter) = {results_fold_eval};
        
        clear results_fold_model results_fold_eval
        
    end
    
    time_line_model(:, time_line_counter) = {runs_model};
    time_line_eval(:, time_line_counter) = {runs_eval};
    
    clear runs_model runs_eval
    
end

location = strcat('&\Results_matfiles\univariate_runs_',string(date),'.mat');
save(location,'time_line_model','time_line_eval');
