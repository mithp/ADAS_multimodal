
%ccc;

function [results_modality_wise] = unshuffler_univariate()

%load('D:\MA_Manuscript_project\multimodality_adas_prediction_020319\saved_results\univariate_runs_30-Mar-2019.mat')
load('U:\MA_Manuscript_project\full_code\saved_results\univariate_runs_05-Apr-2019.mat')%  change this once the runs are comeplete
load('D:\MA_Manuscript_project\multimodality_adas_prediction_020319\saved_results\all_runs_28-Mar-2019.mat')
clear time_line_model
load('D:\MA_Manuscript_project\multimodality_adas_prediction_020319\data\multimodality_adas_prediction_020319.mat')
clear age apoe cmeas edu RID

%time_line_eval = {runs_eval};
clear runs_eval

%time_line_eval = full_data;

for modality_counter = 1 : size(time_line_eval{1,1}{1,1}{1,1},2)
    
    for time_line = 1 : size(time_line_eval,2)
        
        for runs = 1 : size(time_line_eval{1,time_line},2)
            
            %foldid = model_timeline_wise_data{1,1}{1,2};
            foldid = model_timeline_wise_data{time_line,1}{runs,2};
            
            y_unshuffled = nan * ones(size(foldid,1),1);
            plsr_unshuffled = nan * ones(size(foldid,1),1);
            svr_unshuffled = nan * ones(size(foldid,1),1);
            rf_unshuffled = nan * ones(size(foldid,1),1);
            
            fold_1_actual = [time_line_eval{1,time_line}{1,runs}{1,1}{1,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,1}{2,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,1}{3,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,1}{4,modality_counter}];
            fold_2_actual = [time_line_eval{1,time_line}{1,runs}{1,2}{1,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,2}{2,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,2}{3,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,2}{4,modality_counter}];
            fold_3_actual = [time_line_eval{1,time_line}{1,runs}{1,3}{1,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,3}{2,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,3}{3,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,3}{4,modality_counter}];
            fold_4_actual = [time_line_eval{1,time_line}{1,runs}{1,4}{1,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,4}{2,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,4}{3,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,4}{4,modality_counter}];
            fold_5_actual = [time_line_eval{1,time_line}{1,runs}{1,5}{1,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,5}{2,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,5}{3,modality_counter},...
                time_line_eval{1,time_line}{1,runs}{1,5}{4,modality_counter}];
            
            %counters for folds
            f1 = 0;f2 = 0;f3 = 0;f4 = 0;f5 = 0;
            
            for f_id = 1 : size(foldid,1)
                
                fold_value = foldid(f_id);
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
                
                y_unshuffled(f_id,:) = original_pos_value(1);
                plsr_unshuffled(f_id,:) = original_pos_value(2);
                svr_unshuffled(f_id,:) = original_pos_value(3);
                rf_unshuffled(f_id,:) = original_pos_value(4);
                clear original_pos_value
                
            end
            
            y_data = ADAS(:,time_line);
            
            find_index_with_only_score = find(y_data ~= -1);
            y_data= y_data(find_index_with_only_score,:);
            dignosis = Dx(find_index_with_only_score,time_line);
            %remove patients with no score
            
            assert(isequal(y_data,y_unshuffled(:,1)));
            clear f1 f2 f3 f4 f5
            
            unshuffled_y_runs (:,runs) = y_unshuffled;
            unshuffled_plsr_runs (:,runs) = plsr_unshuffled;
            unshuffled_svr_runs (:,runs) = svr_unshuffled;
            unshuffled_rf_runs (:,runs) = rf_unshuffled;
        end
        
        p_PS = nihpd_permutation_test(unshuffled_plsr_runs,unshuffled_svr_runs,unshuffled_y_runs(:,1),1000);
        p_PR = nihpd_permutation_test(unshuffled_plsr_runs,unshuffled_rf_runs,unshuffled_y_runs(:,1),1000);
        p_SR = nihpd_permutation_test(unshuffled_svr_runs,unshuffled_rf_runs,unshuffled_y_runs(:,1),1000);
        
        [maeci_P corrci_P] = nihpd_bs_ci_mae(unshuffled_plsr_runs,unshuffled_y_runs(:,1),1000,0.05);
        [maeci_S corrci_S] = nihpd_bs_ci_mae(unshuffled_svr_runs,unshuffled_y_runs(:,1),1000,0.05);
        [maeci_R corrci_R] = nihpd_bs_ci_mae(unshuffled_rf_runs,unshuffled_y_runs(:,1),1000,0.05);
        
        %if time_line == 1
            %Stratified by diagnosis
            [maeci_P_1 corrci_P_1] = nihpd_bs_ci_mae(unshuffled_plsr_runs(dignosis == 1),unshuffled_y_runs(dignosis == 1,10),1000,0.05);
            [maeci_P_2 corrci_P_2] = nihpd_bs_ci_mae(unshuffled_plsr_runs(dignosis == 2),unshuffled_y_runs(dignosis == 2,10),1000,0.05);
            [maeci_P_3 corrci_P_3] = nihpd_bs_ci_mae(unshuffled_plsr_runs(dignosis == 3),unshuffled_y_runs(dignosis == 3,10),1000,0.05);
            
            [maeci_S_1 corrci_S_1] = nihpd_bs_ci_mae(unshuffled_svr_runs(dignosis == 1),unshuffled_y_runs(dignosis == 1,10),1000,0.05);
            [maeci_S_2 corrci_S_2] = nihpd_bs_ci_mae(unshuffled_svr_runs(dignosis == 2),unshuffled_y_runs(dignosis == 2,10),1000,0.05);
            [maeci_S_3 corrci_S_3] = nihpd_bs_ci_mae(unshuffled_svr_runs(dignosis == 3),unshuffled_y_runs(dignosis == 3,10),1000,0.05);
            
            [maeci_R_1 corrci_R_1] = nihpd_bs_ci_mae(unshuffled_rf_runs(dignosis == 1),unshuffled_y_runs(dignosis == 1,10),1000,0.05);
            [maeci_R_2 corrci_R_2] = nihpd_bs_ci_mae(unshuffled_rf_runs(dignosis == 2),unshuffled_y_runs(dignosis == 2,10),1000,0.05);
            [maeci_R_3 corrci_R_3] = nihpd_bs_ci_mae(unshuffled_rf_runs(dignosis == 3),unshuffled_y_runs(dignosis == 3,10),1000,0.05);
            
            results_UNI(time_line,:) = [p_PS p_PR p_SR maeci_P corrci_P maeci_S corrci_S maeci_R corrci_R maeci_P_1 corrci_P_1 maeci_P_2 corrci_P_2 maeci_P_3 corrci_P_3 ...
                maeci_S_1 corrci_S_1 maeci_S_2 corrci_S_2 maeci_S_3 corrci_S_3...
                maeci_R_1 corrci_R_1 maeci_R_2 corrci_R_2 maeci_R_3 corrci_R_3];
            
        %else
            
         %   results_UNI(time_line,:) = [p_PS p_PR p_SR maeci_P corrci_P maeci_S corrci_S maeci_R corrci_R nan nan nan nan nan nan nan nan nan nan nan nan ...
          %      nan nan nan nan nan nan nan nan nan nan nan nan ...
           %     nan nan nan nan nan nan nan nan nan nan nan nan];
            
        %end
        
        clear unshuffled_y_runs unshuffled_plsr_runs unshuffled_svr_runs unshuffled_rf_runs
        
    end
    
    results_modality_wise(:,:,modality_counter) = results_UNI;
    
    
end

results_modality_wise;

end


