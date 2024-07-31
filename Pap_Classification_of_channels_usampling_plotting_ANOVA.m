%% BF Anova 4 factors
clc;
clear all;
close all;
addpath(genpath('bayesFactor-master'))
participants_info=tdfread(['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\participants.tsv']);
subjects_analysed=[1:25 27:40 42:45 47:58]; % already remove 41 becuase "data" does not exist

ictal_or_inter='ictal';
t=0;
if strcmp(ictal_or_inter,'ictal')
    load(['Connectiv_classif_oversamp_pooled_Seizure.mat']);
elseif strcmp(ictal_or_inter,'interictal')
    load(['Connectiv_classif_oversamp_pooled_Interictal.mat']);
end

t=t+1;
ff=1; % 1= all feaures included
iter_equalis=1;
iter_rand=1;
for p=1:size(Ground_truth,1) % patients
    for fld=1:size(Ground_truth,4) % folds
            Performance_within(p,1:7)=Evaluate(Ground_truth{p,iter_equalis,iter_rand,fld},Predictions{p,iter_equalis,iter_rand,fld});
            try
                [~,~,~,Performance_within(p,8)]=perfcurve(Ground_truth{p,iter_equalis,iter_rand,fld},Predictions{p,iter_equalis,iter_rand,fld},1);
            catch
                Performance_within(p,8)=nan;
            end
    end
end


Measures={'accuracy','sensitivity','specificity','precision','recall','f1-measure','gmean','auc'};
metric=8; %AUC
data_tmp=Performance_within(:,metric);

pt=0;
for P=subjects_analysed
    pt=pt+1;
    g1{pt,1}=participants_info.outcome(P,1);
    g2{pt,1}=participants_info.target(P,:);
    g3{pt,1}=participants_info.lesion_status(P,:);
    g4{pt,1}=participants_info.implant(P,:);

    % if contains(participants_info.outcome(P,1),'S')
    %     g_outcome(pt,1)=1;
    % else
    %     g_outcome(pt,1)=0;
    % end
    % if contains(participants_info.lesion_status(P,1),'N')
    %     g_outcome(pt,3)=1;
    % else
    %     g_outcome(pt,3)=0;
    % end
    % if contains(participants_info.implant(P,1),'S')
    %     g_outcome(pt,4)=1;
    % else
    %     g_outcome(pt,4)=0;
    % end
    
end
clearvars data_ready
nans=isnan(data_tmp);
nans([1])=1; % bad labeling of lesion (1)
nans([4 25 30 37 38])=1; % resection targets with fewer than 5 patients
data_ready=table(data_tmp,g1,g2,g3,g4);
data_ready(nans,:)=[];
data_ready.Properties.VariableNames = {'Effect' 'Outcome' 'ROI' 'Lesion' 'Recording'};
bfFull= bf.anova(data_ready,'Effect ~ Outcome+ROI+Lesion+Recording');
bfRestricted= bf.anova(data_ready,'Effect ~ ROI+Lesion+Recording');
effects(1) = bfFull/bfRestricted;
bfRestricted= bf.anova(data_ready,'Effect ~ Outcome+Lesion+Recording');
effects(2) = bfFull/bfRestricted;
bfRestricted= bf.anova(data_ready,'Effect ~ Outcome+ROI+Recording');
effects(3) = bfFull/bfRestricted;
bfRestricted= bf.anova(data_ready,'Effect ~ Lesion+Outcome+ROI');
effects(4) = bfFull/bfRestricted;


% effects(1) = bf.ttest2(data_tmp(logical(g_outcome(:,1)==1)),data_tmp(logical(g_outcome(:,1)==0)));
% effects(3) = bf.ttest2(data_tmp(logical(g_outcome(:,3)==1)),data_tmp(logical(g_outcome(:,3)==0)));
% effects(4) = bf.ttest2(data_tmp(logical(g_outcome(:,4)==1)),data_tmp(logical(g_outcome(:,4)==0)));

% Plotting accuracies across variables
for condition=1:4
    clearvars data data_random
    figure
    feats=1;
    data(feats,:)=data_tmp;
    cond=nan(size(data,2),1);
    
    if condition==1
        % %         strings={'1','2','3','4'};
        ss={'S','F'};
        strings={'Engel I','Engel II-IV'};
        col=cool;
        shiftee=0.1;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(ss)
                %                 if contains(participants_info.engel(s,1),strings{c})
                if contains(participants_info.outcome(s,1),ss{c})
                    cond(s)=(c);
                end
            end
        end
    elseif condition==2
        %         strings={'FRONTAL','TEMPORAL','FP','MTL','PARIETAL','MFL','INSULAR'};
        strings={'FRONTAL','TEMPORAL','MTL'};
        col=jet;
        shiftee=0.2;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(strings)
                if contains(participants_info.target(s,:),strings{c})
                    cond(s)=(c);
                end
            end
        end
        cond(isnan(cond))=[];
    elseif condition==3
        strings={'LESIONAL','NON-LESIONAL'};
        col=cool;
        shiftee=0.1;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(strings)
                if contains(participants_info.lesion_status(s,:),strings{c})
                    cond(s)=(c);
                end
            end
        end
        cond(1)=[]; % unused
        data(1)=[];
    elseif condition==4
        strings={'SEEG','ECOG'};
        col=cool;
        shiftee=0.1;
        s=0;
        for i=subjects_analysed
            s=s+1;
            for c=1:length(strings)
                if contains(participants_info.implant(s,:),strings{c})
                    cond(s)=(c);
                end
            end
        end
    end
    feats=1;
    bar_width=0.2;
    tmp=round(linspace(1,length(col),length(strings)));
    cols={col(tmp,:)};
    multip=floor(length(unique(cond))/2)*0.2;
    uniques=unique(cond);
    datas=nan(55,length(uniques));
    c=0;
    datas=nan(55,length(uniques));
    for i=uniques'
        c=c+1;
        cond(40)=0; % unused
        data_tmp_tmp=rmoutliers(data(feats,cond==i));
        datas(1:length(data_tmp_tmp),c)=data_tmp_tmp;
        if strcmp(ictal_or_inter,'Interictal')
            colours=[0 0 0];
        else
            colours=[0.3 0.3 0.3];
        end
        swarmchart([c*ones(length(data_tmp_tmp),1)]',data_tmp_tmp,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
        hold on;
    end
    [nanmean(datas) nanstd(datas)]
    boxplot(datas,'Whisker',inf,'color','k');
    
    box off
    grid on
    set(gca,'TickDir','out','Fontsize',16)
    xticks([1:size(datas,2)])
    xlim([0.5 size(datas,2)+0.5])
    plot([0.5 size(datas,2)+0.5],[0.5 0.5],'--k')
    title(['ANOVA BF = ',sprintf('%0.2f',effects(condition))])
    ylabel('\it Localisation (AUC)')
    ylim([0.35 1.0])
    yticks([0.4:0.1:1])

    if condition==2
        strings={'FRT','TPR','MTL'};
    end
    xticklabels(strings)
    if condition==1
        xlabel('Outcome')
    elseif condition==2
        xlabel('Region of Resection')
    elseif condition==3
        xlabel('Pathology')
    elseif condition==4
        xlabel('Recording')
    end
    
    
    clearvars datas
    % BFs across conditions
    pt=0;
    for P=subjects_analysed
        pt=pt+1;
        if condition==1
            g{pt,1}=participants_info.outcome(P,1);
        elseif condition==2
            g{pt,1}=participants_info.target(P,:);
        elseif condition==3
            g{pt,1}=participants_info.lesion_status(P,:);
        elseif condition==4
            g{pt,1}=participants_info.implant(P,:);
        end
    end
    
    data_bf=data;
    nans=isnan(data_bf);
    nans([1])=1; % bad labeling of lesion (1)
    nans([4 25 30 37 38])=1; % resection targets with fewer than 5 patients
    data_bf(nans)=[];
    g(nans)=[];
    uniq_conds_tmp=unique(g);
    
    combs_tmp=nchoosek(1:length(uniq_conds_tmp),2);
    for comb=1:size(combs_tmp,1)
        try
            conditions_BFtmp(comb,1)=bf.ttest2(data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,1)))),data_bf(strcmp(g,uniq_conds_tmp(combs_tmp(comb,2)))));
        catch
            conditions_BFtmp(comb,1)=nan;
        end
        combinations_tmp{comb,1}=uniq_conds_tmp(combs_tmp(comb,1));
        combinations_tmp{comb,2}=uniq_conds_tmp(combs_tmp(comb,2));
    end
    combinations_BF{condition}=conditions_BFtmp;
    combinations{condition}=combinations_tmp;
    clearvars data_bf conditions_BFtmp combinations_tmp
end