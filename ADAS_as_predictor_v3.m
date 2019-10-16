
% ADAS and ADAS + CMEAS on score change
% predicted new score = baseline score + predicted score change  and  actual new score = baseline score + actual score change
% + score or scores change

ccc;
addpath D:\MA_Manuscript_project\multimodality_adas_prediction_020319\genpls
load('D:\MA_Manuscript_project\multimodality_adas_prediction_020319\data\multimodality_adas_prediction_020319.mat')

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

for trails = 1 : 4
    
    for time_line = 2:4
        
        clear ADAS_score X_complex indices
        clear X_complex
        
        if trails == 1 || trails == 3
            X_complex = [ADAS(:,1)];
        elseif trails == 2 || trails == 4
            X_complex = [zscore(apoe) zscore(mri_derived_scores) zscore(mmse) zscore(avf_PET) zscore(amyloid_tau) zscore(fdg_PET) ADAS(:,1)];
        end
        
        %load respectice timline ADAS score
        ADAS_score = ADAS(:,time_line);
        %Diagnosis_DX = Dx(:,time_line);
        y_data = ADAS_score;
        y_data_baseline = ADAS(:,1);
        
        
        
        %X_complex = ADAS(:,1);
        find_index_with_only_score = find(y_data ~= -1);
        
        %remove patients with no score
        
        find_index_with_only_score = find(y_data ~= -1);
        X_complex= X_complex(find_index_with_only_score,:);
        y_data= y_data(find_index_with_only_score,:);
        
        y_data_baseline = y_data_baseline(find_index_with_only_score,:);
        y_change = y_data - y_data_baseline ;
        
        y_TL_score = y_data;
        
        if trails == 1 || trails == 2
            y_data = (y_change);
        end
        
        %Diagnosis_DX = Diagnosis_DX(find_index_with_only_score,:);
        
        %Runs controll the repetition for cross-validation
        for runs = 1 : 10
            
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
                
                y_BS_M_score = y_data_baseline(modelInd,:);
                y_TS_M_score = y_TL_score(modelInd,:);
                
                %Let's keep this portion for testing
                x_val  = X_complex(valInd,:);
                y_val  = y_data(valInd,:);
                
                y_BS_V_score = y_data_baseline(valInd,:);
                y_TS_V_score = y_TL_score(valInd,:);
                
                % initialize maxinmum number of PLS components
                if size(X_complex,2) <= 15
                    ncompmax = size(X_complex,2); %if there are fewer variables than 15
                else
                    ncompmax = 15;
                end
                
                ncomp = 0;
                
                %To find the optimal number of components lets further divide
                %the model data so the validation data is truely untouched
                
                pls_divide_count = size(y_model,1);
                [trainInd,testInd,~] = dividerand(pls_divide_count,0.80,0.20,0);
                
                % Inital model is based on this data
                x_train = x_model(trainInd,:);
                y_train  = y_model(trainInd,:);
                y_BS_Train_score = y_BS_M_score(trainInd,:);
                y_TS_Train_score = y_TS_M_score(trainInd,:);
                
                %Used to select the number of parameters
                x_test  = x_model(testInd,:);
                y_test  = y_model(testInd,:);
                y_BS_Test_score = y_BS_M_score(testInd,:);
                y_TS_Test_score = y_TS_M_score(testInd,:);
                
                %% Simple Partial least squares regression model (with 1-10 components)
                for ncomp=1:ncompmax
                    
                    [XL,~,Xs,~,betaPLS,PctVar,msep, stats]= plsregress(x_train,y_train,ncomp,'CV',10);
                    yfit_PLS_train = [ones(size(x_train,1),1) x_train]*betaPLS;
                    yfit_PLS_test = [ones(size(x_test,1),1) x_test]*betaPLS;
                    
                    if trails == 3 || trails == 4
                        [r2_train(ncomp) rmse_train(ncomp)] = rsquare(y_train,yfit_PLS_train);
                        [r2_test(ncomp) rmse_test(ncomp)] = rsquare(y_test,yfit_PLS_test);
                        
                        
                    elseif trails == 1 || trails == 2
                        [r2_train(ncomp) rmse_train(ncomp)] = rsquare(y_TS_Train_score,y_BS_Train_score + yfit_PLS_train );
                        [r2_test(ncomp) rmse_test(ncomp)] = rsquare(y_TS_Test_score,y_BS_Test_score + yfit_PLS_test );
                    end
                end
                
                %Finding the optimal number of components for this split
                number_of_components = find(rmse_test == min(rmse_test));
                
                %lets evaluate the validation set to find out the truth
                [~,~,~,~,betaPLS_model,~,~, ~]= plsregress(x_model,y_model,number_of_components);%,'CV',number_of_components);
                
                yfit_PLS_model = [ones(size(x_model,1),1) x_model] * betaPLS_model; %model predition training data
                yfit_PLS_val = [ones(size(x_val,1),1) x_val]*betaPLS_model; %model predition testing data
                
                %figure; scatter(y_val,yfit_PLS_val);
                %close all;
                
                if trails == 1 || trails == 2
                    eval_data_fold (fold,:)= {y_val,yfit_PLS_val, y_BS_V_score,y_TS_V_score};
                elseif trails == 3 || trails == 4
                    eval_data_fold (fold,:)= {y_val,yfit_PLS_val, [],[]};
                end
                
                
                clear rmse_test
                
            end
            
            %need to unshuffle the values
            
            clear fold_1_actual fold_2_actual fold_3_actual fold_4_actual fold_5_actual
            clear y_unshuffled plsr_unshuffled
            %collting data from folds
            fold_1_actual = [eval_data_fold{1,1}, eval_data_fold{1,2}, eval_data_fold{1,3}, eval_data_fold{1,4} ];
            fold_2_actual = [eval_data_fold{2,1}, eval_data_fold{2,2}, eval_data_fold{2,3}, eval_data_fold{2,4} ];
            fold_3_actual = [eval_data_fold{3,1}, eval_data_fold{3,2}, eval_data_fold{3,3}, eval_data_fold{3,4} ];
            fold_4_actual = [eval_data_fold{4,1}, eval_data_fold{4,2}, eval_data_fold{4,3}, eval_data_fold{4,4} ];
            fold_5_actual = [eval_data_fold{5,1}, eval_data_fold{5,2}, eval_data_fold{5,3}, eval_data_fold{5,4} ];
            
            %counters for folds
            f1 = 0;f2 = 0;f3 = 0;f4 = 0;f5 = 0;
            y_unshuffled = nan * ones(size(indices,1),1);
            plsr_unshuffled = nan * ones(size(indices,1),1);
            
            for f_id = 1 : size(indices,1)
                
                fold_value = indices(f_id);
                
                switch fold_value
                    case 1
                        f1 = f1+1;
                        original_pos_value = fold_1_actual(f1,:);
                    case 2
                        f2 = f2+1;
                        original_pos_value = fold_2_actual(f2,:);
                    case 3
                        f3 = f3+1;
                        original_pos_value = fold_3_actual(f3,:);
                    case 4
                        f4 = f4+1;
                        original_pos_value = fold_4_actual(f4,:);
                    case 5
                        f5 = f5+1;
                        original_pos_value = fold_5_actual(f5,:);
                end
                
                
                if trails == 1 || trails == 2
                    y_unshuffled(f_id,:) = original_pos_value(4);
                elseif trails == 3 || trails == 4
                    y_unshuffled(f_id,:) = original_pos_value(1);
                end
                
                %plsr_unshuffled(f_id,:) = original_pos_value(3) + original_pos_value(2);
                
                if trails == 1 || trails == 2
                    plsr_unshuffled(f_id,:) = original_pos_value(3) + original_pos_value(2);
                elseif trails == 3 || trails == 4
                    plsr_unshuffled(f_id,:) = original_pos_value(2);
                end
                
                clear original_pos_value
            end
            y_aggregate(:,runs) = y_unshuffled ;
            y_pred_aggregate(:,runs) = plsr_unshuffled ;
            
        end
        
        
        [maeci_P corrci_P] = nihpd_bs_ci_mae(y_pred_aggregate,y_aggregate(:,1),1000,0.05);
        
        %p = nihpd_permutation_test(yhat,yhat2,y,nboot)
        
        results(time_line,:) = [mean(maeci_P) std(maeci_P) mean(corrci_P) std(corrci_P)];
        
        clear y_aggregate y_pred_aggregate
        
    end
    
    results_change (trails) = {results};
    clear results
end


disp('Complete')

counter = 1;

for time_ = 2 : 4
    
    subplot(3,2, counter)
    y = [results_change{1,1}(time_,3) , results_change{1,2}(time_,3), results_change{1,3}(time_,3) , results_change{1,4}(time_,3)];         % random y values (3 groups of 4 parameters)
    errY = [ results_change{1,1}(time_,4) , results_change{1,2}(time_,4), results_change{1,3}(time_,4) , results_change{1,4}(time_,4)  ];          % 10% error
    h = barwitherr(errY, y);% Plot with errorbars
    set(gca,'XTickLabel',{'ADAS pred change','CMEAS+ADAS','ADAS pred score','CMEAS+ADAS'})
    ylabel('Mean Correlation')
    set(h(1),'FaceColor','k');
    grid on
    counter = counter +1;
    subplot(3,2, counter)
    clear y errY h
    
    y = [results_change{1,1}(time_,1) , results_change{1,2}(time_,1),results_change{1,3}(time_,1) , results_change{1,4}(time_,1) ];         % random y values (3 groups of 4 parameters)
    errY = [ results_change{1,1}(time_,2) , results_change{1,2}(time_,2) ,results_change{1,3}(time_,2) , results_change{1,4}(time_,2)  ];          % 10% error
    h = barwitherr(errY, y);% Plot with errorbars
    set(gca,'XTickLabel',{'ADAS pred change','CMEAS+ADAS','ADAS pred score','CMEAS+ADAS'})
    ylabel('Mean MAE')
    set(h(1),'FaceColor','k');
    switch time_
        case 2
            title ('12 months');
        case 3
            title ('24 months');
        case 4
            title ('36 months');
    end
    grid on
    counter = counter +1;
    clear y errY h
    
end
suptitle('PLSR prediction based on previous scores')
