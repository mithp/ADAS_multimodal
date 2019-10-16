%ccc;

function [results_MV] = unshuffler()
load('D:\MA_Manuscript_project\multimodality_adas_prediction_020319\data\multimodality_adas_prediction_020319.mat')
clear age apoe cmeas edu RID
load('D:\MA_Manuscript_project\multimodality_adas_prediction_020319\saved_results\all_runs_28-Mar-2019.mat')

load('D:\MA_Manuscript_project\multimodality_adas_prediction_020319\saved_results\univariate_runs_30-Mar-2019.mat')
clear time_line_model

for time_line = 1 :  size(eval_timeline_wise_data,1)
    
    for runs = 1: size(eval_timeline_wise_data{1,1},1)
        
        % code to unshuffle the fold indies to normal
        foldid = model_timeline_wise_data{time_line,1}{runs,2};
        
        %collecting fold chucks from folds
        fold_1_actual = [eval_timeline_wise_data{time_line,1}{runs,1}{1,2},...
            eval_timeline_wise_data{time_line,1}{runs,1}{1,3},...
            eval_timeline_wise_data{time_line,1}{runs,1}{1,4},...
            eval_timeline_wise_data{time_line,1}{runs,1}{1,5}];
        fold_2_actual = [eval_timeline_wise_data{time_line,1}{runs,1}{2,2},...
            eval_timeline_wise_data{time_line,1}{runs,1}{2,3},...
            eval_timeline_wise_data{time_line,1}{runs,1}{2,4},...
            eval_timeline_wise_data{time_line,1}{runs,1}{2,5}];
        fold_3_actual = [eval_timeline_wise_data{time_line,1}{runs,1}{3,2},...
            eval_timeline_wise_data{time_line,1}{runs,1}{3,3},...
            eval_timeline_wise_data{time_line,1}{runs,1}{3,4},...
            eval_timeline_wise_data{time_line,1}{runs,1}{3,5}];
        fold_4_actual = [eval_timeline_wise_data{time_line,1}{runs,1}{4,2},...
            eval_timeline_wise_data{time_line,1}{runs,1}{4,3},...
            eval_timeline_wise_data{time_line,1}{runs,1}{4,4},...
            eval_timeline_wise_data{time_line,1}{runs,1}{4,5}];
        fold_5_actual = [eval_timeline_wise_data{time_line,1}{runs,1}{5,2},...
            eval_timeline_wise_data{time_line,1}{runs,1}{5,3},...
            eval_timeline_wise_data{time_line,1}{runs,1}{5,4},...
            eval_timeline_wise_data{time_line,1}{runs,1}{5,5}];
        
        %counters for folds
        f1 = 0;f2 = 0;f3 = 0;f4 = 0;f5 = 0;
        
        y_unshuffled = nan * ones(size(foldid,1),1);
        plsr_unshuffled = nan * ones(size(foldid,1),1);
        svr_unshuffled = nan * ones(size(foldid,1),1);
        rf_unshuffled = nan * ones(size(foldid,1),1);
        
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
        
        unshuffled_y_runs (:,runs) = y_unshuffled;
        unshuffled_plsr_runs (:,runs) = plsr_unshuffled;
        unshuffled_svr_runs (:,runs) = svr_unshuffled;
        unshuffled_rf_runs (:,runs) = rf_unshuffled;
        
        clear y_unshuffled plsr_unshuffled rf_unshuffled
        
        
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
        
        results_MV(time_line,:) = [p_PS p_PR p_SR maeci_P corrci_P maeci_S corrci_S maeci_R corrci_R maeci_P_1 corrci_P_1 maeci_P_2 corrci_P_2 maeci_P_3 corrci_P_3 ...
            maeci_S_1 corrci_S_1 maeci_S_2 corrci_S_2 maeci_S_3 corrci_S_3...
            maeci_R_1 corrci_R_1 maeci_R_2 corrci_R_2 maeci_R_3 corrci_R_3];
        
    %else
     %   results_MV(time_line,:) = [p_PS p_PR p_SR maeci_P corrci_P maeci_S corrci_S maeci_R corrci_R nan nan nan nan nan nan nan nan nan nan nan nan ...
      %      nan nan nan nan nan nan nan nan nan nan nan nan ...
       %     nan nan nan nan nan nan nan nan nan nan nan nan];
    %end
    
    clear unshuffled_y_runs unshuffled_plsr_runs unshuffled_svr_runs unshuffled_rf_runs
    
end
results_MV;

end


