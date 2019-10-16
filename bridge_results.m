
% Bringing univariate and multivariate on the same plaatform
MV_DATA = unshuffler;
UNI_DATA = unshuffler_univariate;

ALL_DATA = UNI_DATA;
ALL_DATA(:,:,9) = MV_DATA; % Change this after all runs

time_lines = size(ALL_DATA,1);
modalities = size(ALL_DATA,2);

timeline_data = transpose(squeeze(ALL_DATA(1,:,:))); % modalities as rows

for time_ = 1 : time_lines
    
    clear timeline_data
    timeline_data = transpose(squeeze(ALL_DATA(time_,:,:))); % modalities as rows
    
    for i = 1: size(timeline_data,1) %modality
        
        maeci_P(time_,i)  = mean(timeline_data(i,4:5));
        maeci_P_STD(time_,i)  = std(timeline_data(i,4:5));
        
        corrci_P(time_,i) = mean(timeline_data(i,6:7));
        corrci_P_STD(time_,i) = std(timeline_data(i,6:7));
        
        maeci_S(time_,i)  = mean(timeline_data(i,8:9));
        maeci_S_STD(time_,i)  = std(timeline_data(i,8:9));
        
        corrci_S(time_,i) = mean(timeline_data(i,10:11));
        corrci_S_STD(time_,i) = std(timeline_data(i,10:11));
        
        maeci_R(time_,i)  = mean(timeline_data(i,12:13));
        maeci_R_STD(time_,i)  = std(timeline_data(i,12:13));
        
        corrci_R(time_,i) = mean(timeline_data(i,14:15));
        corrci_R_STD(time_,i) = std(timeline_data(i,14:15));
        
        PS = timeline_data(i,1);
        PR = timeline_data(i,2);
        SR = timeline_data(i,3);
        
        if abs(PS) < 0.05
            stars = 1;
            %continue;
        elseif abs(PS) < 0.01
            stars = 2;
            %continue;
        elseif abs(PS) < 0.001
            stars = 3;
        else
            stars = 0;
        end
        
        p_PS(time_,i) = stars;
        clear stars
        if abs(PR) < 0.05
            stars = 1;
            %continue;
        elseif abs(PR) < 0.01
            stars = 2;
            %continue;
        elseif abs(PR) < 0.001
            stars = 3;
        else
            stars = 0;
        end
        
        p_PR(time_,i) = stars;
        clear stars
        
        if abs(SR) < 0.05
            stars = 1;
            %continue;
        elseif abs(SR) < 0.01
            stars = 2;
            %continue;
        elseif abs(SR) < 0.001
            stars = 3;
        else
            stars = 0;
        end
        
        p_SR(time_,i) = stars;
        
        %if time_ == 1
        
        maeci_P_1(time_,i)     = mean(timeline_data(i,16:17));
        maeci_P_1_STD(time_,i) = std(timeline_data(i,16:17));
        corrci_P_1(time_,i)     = mean(timeline_data(i,18:19));
        corrci_P_1_STD(time_,i) = std(timeline_data(i,18:19));
        maeci_P_2(time_,i)     = mean(timeline_data(i,20:21));
        maeci_P_2_STD(time_,i) = std(timeline_data(i,20:21));
        corrci_P_2(time_,i)    = mean(timeline_data(i,22:23));
        corrci_P_2_STD(time_,i) = std(timeline_data(i,22:23));
        maeci_P_3(time_,i)     = mean(timeline_data(i,24:25));
        maeci_P_3_STD(time_,i) = std(timeline_data(i,24:25));
        corrci_P_3(time_,i)     = mean(timeline_data(i,26:27));
        corrci_P_3_STD(time_,i) = std(timeline_data(i,26:27));
        
        maeci_S_1(time_,i)     = mean(timeline_data(i,28:29));
        maeci_S_1_STD(time_,i) = std(timeline_data(i,28:29));
        corrci_S_1(time_,i)     = mean(timeline_data(i,30:31));
        corrci_S_1_STD(time_,i) = std(timeline_data(i,30:31));
        maeci_S_2(time_,i)     = mean(timeline_data(i,32:33));
        maeci_S_2_STD(time_,i) = std(timeline_data(i,32:33));
        corrci_S_2(time_,i)     = mean(timeline_data(i,34:35));
        corrci_S_2_STD(time_,i) = std(timeline_data(i,34:35));
        maeci_S_3(time_,i)     = mean(timeline_data(i,36:37));
        maeci_S_3_STD(time_,i) = std(timeline_data(i,36:37));
        corrci_S_3(time_,i)     = mean(timeline_data(i,38:39));
        corrci_S_3_STD(time_,i) = std(timeline_data(i,38:39));
        
        maeci_R_1(time_,i)     = mean(timeline_data(i,40:41));
        maeci_R_1_STD(time_,i) = std(timeline_data(i,40:41));
        corrci_R_1(time_,i)     = mean(timeline_data(i,42:43));
        corrci_R_1_STD(time_,i) = std(timeline_data(i,42:43));
        maeci_R_2(time_,i)     = mean(timeline_data(i,44:45));
        maeci_R_2_STD(time_,i) = std(timeline_data(i,44:45));
        corrci_R_2(time_,i)     = mean(timeline_data(i,46:47));
        corrci_R_2_STD(time_,i) = std(timeline_data(i,46:47));
        maeci_R_3(time_,i)     = mean(timeline_data(i,48:49));
        maeci_R_3_STD(time_,i) = std(timeline_data(i,48:49));
        corrci_R_3(time_,i)     = mean(timeline_data(i,50:51));
        corrci_R_3_STD(time_,i) = std(timeline_data(i,50:51));
        
        %end
        clear PS PR SR stars
    end
end



%legend('AGE','EDU','APOE','MRI','MMSE','AVF','Tau','FDG','Multimodal')

clear time_ ;
counter = 1;
for time_ = 1 : time_lines
    if time_ == 1
        flag = 1;
    else
        flag = 0;
    end
    subplot(time_lines,2,counter)
    
    y = [corrci_P(time_,:);corrci_S(time_,:);corrci_R(time_,:)];         % random y values (3 groups of 9 parameters)
    errY = [ corrci_P_STD(time_,:); corrci_S_STD(time_,:); corrci_R_STD(time_,:) ];          % 10% error
    h = barwitherr(errY', y');% Plot with errorbars
    h(1).FaceColor = [.1 .1 .1];h(2).FaceColor = [.5 .5 .5];h(3).FaceColor = [.8 .8 .8];
    %set(gca,'XTickLabel',{'PLSR','SVR','RF'})
    set(gca,'XTickLabel',{'AGE','EDU','APOE','MRI','NPB','AVF45-PET','CSF','FDG','Multimodal'})
    set(gca, 'fontsize', 12, 'linewidth', 2);
    set(gca,'FontSize',12,'TickDir','in','fontWeight','bold', 'FontName', 'Times New Roman');
    
    if flag == 1
        %legend('AGE','EDU','APOE','MRI','NPB','AVF45-PET','CSF','FDG','Multimodal')
        legend('PLSR','SVR','RF')
    end
    ylabel('Mean Correlation')
    set(gca, 'fontsize', 12, 'linewidth', 2);
    set(gca,'FontSize',12,'TickDir','in','fontWeight','bold', 'FontName', 'Times New Roman');
    set(h(1),'FaceColor','k');
    clear y errY
    
    grid on;
    clear y errY
    counter = counter +1;
    
    switch time_
        case 1
            title ('Baseline');
        case 2
            title ('12 months');
        case 3
            title ('24 months');
        case 4
            title ('36 months');
            
    end
    
    
    subplot(time_lines,2,counter)
    
    y = [maeci_P(time_,:);maeci_S(time_,:);maeci_R(time_,:)];         % random y values (3 groups of 4 parameters)
    errY = [ maeci_P_STD(time_,:); maeci_S_STD(time_,:); maeci_R_STD(time_,:) ];          % 10% error
    h = barwitherr(errY', y');% Plot with errorbars
    h(1).FaceColor = [.1 .1 .1];h(2).FaceColor = [.5 .5 .5];h(3).FaceColor = [.8 .8 .8];
    set(gca,'XTickLabel',{'AGE','EDU','APOE','MRI','NPB','AVF45-PET','CSF','FDG','Multimodal'})
    set(gca, 'fontsize', 12, 'linewidth', 2);
    set(gca,'FontSize',12,'TickDir','in','fontWeight','bold', 'FontName', 'Times New Roman');
    %legend('AGE','EDU','APOE','MRI','MMSE','AVF','Tau','FDG','Multimodal')
    ylabel('Mean MAE (ADAS Score)')
    set(gca, 'fontsize', 12, 'linewidth', 2);
    set(gca,'FontSize',12,'TickDir','in','fontWeight','bold', 'FontName', 'Times New Roman');
    set(h(1),'FaceColor','k');
    
    counter = counter +1;
    grid on
    
    
end
suptitle('Regression performance at different timelines');
clear y errY h

figure;
clear time_ ;
counter = 1;
for time_ = 1 : time_lines
    if time_ == 1
        flag = 1;
    else
        flag = 0;
    end
    subplot(time_lines,2,counter)
    
    y = [maeci_P_1(time_,:);maeci_P_2(time_,:);maeci_P_3(time_,:)];         % random y values (3 groups of 4 parameters)
    errY = [ maeci_P_1_STD(time_,:);maeci_P_2_STD(time_,:);maeci_P_3_STD(time_,:) ];          % 10% error
    h = barwitherr(errY', y');% Plot with errorbars
    h(1).FaceColor = [.1 .1 .1];h(2).FaceColor = [.5 .5 .5];h(3).FaceColor = [.8 .8 .8];
    set(gca,'XTickLabel',{'AGE','EDU','APOE','MRI','NPB','AVF45-PET','CSF','FDG','Multimodal'})
    set(gca, 'fontsize', 12, 'linewidth', 2);
    set(gca,'FontSize',12,'TickDir','in','fontWeight','bold', 'FontName', 'Times New Roman');
    %legend('AGE','EDU','APOE','MRI','MMSE','AVF','Tau','FDG','Multimodal')
    ylabel('Mean MAE (ADAS Score)')
    set(gca, 'fontsize', 12, 'linewidth', 2);
    set(gca,'FontSize',12,'TickDir','in','fontWeight','bold', 'FontName', 'Times New Roman');
    
    switch time_
        case 1
            title ('Baseline');
        case 2
            title ('12 months');
        case 3
            title ('24 months');
        case 4
            title ('36 months');
            
    end
    set(h(1),'FaceColor','k');
    clear y errY h
    counter = counter +1;
    grid on
    
    subplot(time_lines,2,counter)
    
    y = [corrci_P_1(time_,:);corrci_P_2(time_,:);corrci_P_3(time_,:)];         % random y values (3 groups of 4 parameters)
    errY = [ corrci_P_1_STD(time_,:);corrci_P_2_STD(time_,:);corrci_P_3_STD(time_,:)];          % 10% error
    h = barwitherr(errY', y');% Plot with errorbars
    h(1).FaceColor = [.1 .1 .1];h(2).FaceColor = [.5 .5 .5];h(3).FaceColor = [.8 .8 .8];
    %set(gca,'XTickLabel',{'Normal','MCI','AD'})
    set(gca,'XTickLabel',{'AGE','EDU','APOE','MRI','NPB','AVF45-PET','CSF','FDG','Multimodal'})
    if flag == 1
        %legend('AGE','EDU','APOE','MRI','NPB','AVF45-PET','CSF','FDG','Multimodal')
        legend('Normal','MCI','AD')
    end
    ylabel('Mean correlation')
    set(gca, 'fontsize', 12, 'linewidth', 2);
    set(gca,'FontSize',12,'TickDir','in','fontWeight','bold', 'FontName', 'Times New Roman');
    set(h(1),'FaceColor','k');
    clear y errY h
    counter = counter +1;
    grid on
end

suptitle('PLSR Test Results');
