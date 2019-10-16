function [] = regress_explorer(X_complex_concat,ADAS)


for time_line = 1:4 % Controlls the timeline periods
    
    clear ADAS_score X_complex
    
    %Naive contcatenation
    X_complex = X_complex_concat;
    
    %load respectice timline ADAS score
    ADAS_score = ADAS(:,time_line);
    %Diagnosis_DX = Dx(:,time_line);
    
    y_data = ADAS_score;
    %remove patients with no score
    
    find_index_with_only_score = find(y_data ~= -1);
    X_complex= X_complex(find_index_with_only_score,:);
    y_data= y_data(find_index_with_only_score,:);
    
    %this calls Genetic algorithm code to find out the impotance of the
    %each paramter contributing to PLSR's RMSECV
    var_importance(time_line,:) = ga_plot_MRI (X_complex,y_data,time_line);
    
    %Y data normalization
    %ADAS_score_normalized = (y_data);
    
    %Runs controll the repetition for cross-validation
    for runs = 1 : 10 % This takes aprox 45 hours, so for testing purpose change this to 1 or 2
        
        
        % Get the size to feed to k-folds
        num_points = size(y_data,1);
        
        %specify the folds
        K = 5;
        
        %this function assigns the fold indices randlomly
        indices = crossvalind('Kfold',num_points,K);
        
        for fold = 1 : K
            
            
            % Logical operation to mark the indices
            valInd = (indices == fold);
            modelInd = ~valInd;
            
            % The portion for modelling
            x_model = X_complex(modelInd,:);
            y_model  = y_data(modelInd,:);
            
            %Let's keep this portion for testing
            x_val  = X_complex(valInd,:);
            y_val  = y_data(valInd,:);
            
            % initialize maxinmum number of PLS components
            if size(X_complex,2) <= 15
                ncompmax = size(X_complex,2); %if there are fewer variables than 15
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
            
            [r2_model_folds(fold) RMSECV_t] = rsquare(y_model,yfit_PLS_model);
            [r2_val_folds(fold) RMSEP_t] = rsquare(y_val,yfit_PLS_val);
            
            %figure;scatter(y_val,yfit_PLS_val);
            %pause(0.2)
            
            %need to add code to save the training, testing and prediction
            %data: this is still pending
            
            N_fold(fold) = number_of_components;
            RMSECV_folds(fold) = RMSECV_t/range(y_model)*100;
            RMSEP_folds(fold) = RMSEP_t/range(y_val)*100;
            %close all
            
            %The code for SVM and RF in a seerate file to keeps things
            %simple
            svm_n_RF_section;
            
            %univariate_section;
            
            %Save the modelling data
            model_data_fold (fold,:) = {x_model,y_model,yfit_PLS_model,yfit_model_svm,yfit_model_rf};
            eval_data_fold (fold,:)= {x_val,y_val,yfit_PLS_val,yfit_val_svm,yfit_val_rf};
            
            
            
        end
        close all;
        clc;
        
        %collect all the parameters and save
        R2_model(time_line,runs) = median(r2_model_folds) ;
        R2_validation (time_line,runs) = median(r2_val_folds);
        RMSECV(time_line,runs) = median(RMSECV_folds);
        RMSEP(time_line,runs) = median(RMSEP_folds);
        N (time_line,runs) = median(N_fold);
        
        R2_model_svm(time_line,runs) = median(r2_model_folds_svm) ;
        R2_validation_svm (time_line,runs) = median(r2_val_folds_svm);
        RMSECV_svm(time_line,runs) = median(RMSECV_fold_ssvm);
        RMSEP_svm(time_line,runs) = median(RMSEP_folds_svm);
        
        R2_model_rf(time_line,runs) = median(r2_model_folds_rf) ;
        R2_validation_rf (time_line,runs) = median(r2_val_folds_rf);
        RMSECV_rf(time_line,runs) = median(RMSECV_folds_rf);
        RMSEP_rf(time_line,runs) = median(RMSEP_folds_rf);
        
        %Collect_timeline_wise
        model_run_wise_data(runs,:)  = {model_data_fold,indices};
        eval_runs_wise_data(runs,:)  = {eval_data_fold};
        
    end
    
    % collect_data_per_runs
    model_timeline_wise_data(time_line,:)  = {model_run_wise_data};
    eval_timeline_wise_data(time_line,:)  = {eval_runs_wise_data};
    
    %end of time_line_for_loop
end
disp('Complete')

location = strcat('&multimodality_adas_prediction_020319\Results_matfiles\all_runs_',string(date),'.mat')
save(location,'model_timeline_wise_data','eval_timeline_wise_data');

end

