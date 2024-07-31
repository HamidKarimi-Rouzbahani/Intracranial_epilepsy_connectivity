% This code classifies resected and non-resected contacts within each patient
clc
clear all
close all
addpath(genpath('bayesFactor-master'))

%% Featrues
% all
% without corr, te, gpfit and ccm
connectivities={'cce','di',...
    'dtf','dcoh','pdcoh',...
    'sgc','gd','psi','lmfit',...
    'anm','cds','reci',...
    'igci'};

net_feats={'_InStrgth','_OutStrgth','_SrcPassTim','_ClustCoef','_Eccent','_NodBtw'};
net_feats_detailed={'in strength','out strength','source passage time','clustering coefficient','eccentricity','node betweenness'};
ppr_ordrd_connectivities={'ANM','IGCI','CDS','RECI','CCE','DI',...    
                          'GD','PSI','DTF','DCOH','PDCOH','SGC',...
                          'LMFIT'};

ppr_ordrds=[10 13 11 12 1 2 7 8 3 4 5 6 9];


ff=0;
for i=1:length(connectivities)
    for j=1:length(net_feats)
        ff=ff+1;
        Features_labels{ff}={[connectivities{i} net_feats{j}]};
    end
end

%% Evaluation of features
ictal_or_inter='ictal';
all_patients=[1:39 41:56]; % patient indices

X_all_norm=nan(30,104,2,length(all_patients),5,3);

p=0;
for Patient=all_patients % patient indices
    p=p+1;
    if strcmp(ictal_or_inter,'interictal')
        if Patient==2
            trls=2;
        else
            trls=1:2;
        end
    else
        if (Patient==2 || Patient==6 || Patient==27)
            trls=1;
        elseif (Patient==9 || Patient==10 || Patient==15 || Patient==28 || Patient==45 || Patient==50)
            trls=1:2;
        elseif (Patient==1 || Patient==3 || Patient==4 || Patient==5 || Patient==11 || Patient==13 || Patient==16 || Patient==21 || Patient==22 || Patient==30 || Patient==31 || Patient==33 || Patient==35 || Patient==41 || Patient==43 || Patient==44 || Patient==49 || Patient==56)
            trls=1:3;
        elseif (Patient==7 || Patient==12 || Patient==20 || Patient==23 || Patient==24 || Patient==46)
            trls=1:4;
        else
            trls=1:5;
        end
    end

    for trial=trls
        for epoch=1:3

            % Classificaiton is performed for each patient separately
            if strcmp(ictal_or_inter,'ictal')
                load([num2str(Patient),'_Project2_connect_features_SOZ_NON_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.mat']);
                load([num2str(Patient),'_Project2_data_for_PyConnectivity_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.mat'],'channels_for_connectiv_inds','channels_resctd_inds','channels_soz_inds');
            elseif strcmp(ictal_or_inter,'interictal')
                load([num2str(Patient),'_Project2_connect_features_SOZ_NON_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.mat']);
                load([num2str(Patient),'_Project2_data_for_PyConnectivity_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.mat'],'channels_for_connectiv_inds','channels_resctd_inds','channels_soz_inds');
            end

            % Labels
            labels=[1 0]; % SOZ followed by NON channels as set in previous code

            data_targ=Features_all(labels==1,:);
            data_non=Features_all(labels==0,:);

            %% Normalise between 0 and 1
            
            X=abs([data_targ;data_non]);
            y=[ones(size(data_targ,1),1);zeros(size(data_non,1),1)];
            X_norm=(X-min(X))./(max(X)-min(X));
            Xtmp=X_norm(y==1,:);
            X_all_norm(1:size(Xtmp,1),1:size(Xtmp,2),1,p,trial,epoch)=Xtmp;
            Xtmp=X_norm(y==0,:);
            X_all_norm(1:size(Xtmp,1),1:size(Xtmp,2),2,p,trial,epoch)=Xtmp;

            %% Unnormalised
            % X_all_norm(1:size(data_targ,1),1:size(data_targ,2),1,p,trial,epoch)=data_targ;
            % X_all_norm(1:size(data_non,1),1:size(data_non,2),2,p,trial,epoch)=data_non;

        end
    end
end

%% Over subjects separated by trials and epochs
% figure;
trials=1:5; % which trials to average
epoch=1:3;    % which epochs to average
net_feats={'_InStrgth','_OutStrgth','_SrcPassTim','_ClustCoef','_Eccent','_NodBtw'};
net_feats_detailed={'in strength','out strength','source passage time','clustering coefficient','eccentricity','node betweenness'};

for nf=1:1
    feat_to_plot=net_feats{nf};
    for i=1:length(Features_labels)
        if contains(Features_labels{i},feat_to_plot)
            targ_feat(i)=1;
        else
            targ_feat(i)=0;
        end
    end
    targ_feat=logical(targ_feat);

    % subplot(3,2,nf)
    % bias=0.1;

    X_all_norm_1=squeeze(nanmean(nanmean(nanmean(X_all_norm(:,targ_feat,1,:,trials,epoch),6),5),1));
    X_all_norm_1=X_all_norm_1(ppr_ordrds,:);
    
    X_all_norm_2=squeeze(nanmean(nanmean(nanmean(X_all_norm(:,targ_feat,2,:,trials,epoch),6),5),1));
    X_all_norm_2=X_all_norm_2(ppr_ordrds,:);
    
  
    figure;
    c=0;
    datas_tmp_1_outed=nan(55,13);
    datas_tmp_2_outed=nan(55,13);
    shifting=0.7;
    for feat=1:13
        c=c+4;
        colours=[0.7 0.3 0.3];
        data_tmp1=X_all_norm_1(feat,~isnan(X_all_norm_1(feat,:)));
        datas_tmp_1_outed(1:length(data_tmp1),feat)=data_tmp1';
        c1=[(c-shifting)*ones(sum(~isnan(data_tmp1)),1)]';
        cs1(feat)=c1(1);
        swarmchart(c1,data_tmp1,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
        hold on;
        colours=[0.3 0.3 0.7];
        data_tmp2=X_all_norm_2(feat,~isnan(X_all_norm_2(feat,:)));
        datas_tmp_2_outed(1:length(data_tmp2),feat)=data_tmp2';
        c2=[(c+shifting)*ones(sum(~isnan(data_tmp2)),1)]';
        cs2(feat)=c2(1);
        swarmchart(c2,data_tmp2,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
        BFs=bf.ttest2(data_tmp1',data_tmp2');
        y_position=0.02;
        if BFs<100
            t=text(c,y_position,sprintf('%0.2f',BFs));
        elseif BFs>100
            t=text(c,y_position,'>100');
        end
        t.Rotation=90;
        t.FontSize=14;
        if BFs>10
            t.FontWeight='Bold';
            t.Color=[0.3 0.7 0.3];
        end
        clearvars data_tmp1 data_tmp2
    end

    boxplot(datas_tmp_1_outed,cs1,'Positions',cs1,'Whisker',inf,'color',[0.7 0.3 0.3]);
    hold on;
    boxplot(datas_tmp_2_outed,cs2,'Positions',cs2,'Whisker',inf,'color',[0.3 0.3 0.7]);
    box off
    set(gca,'TickDir','out','Fontsize',16)
    ylabel(['\it ',net_feats_detailed{1,nf}])
    xticks([cs1+shifting])
    xticklabels([ppr_ordrd_connectivities])
    title([ictal_or_inter])    
    grid on
    ylim([0 1])
    set(gca,'fontsize',14)
    xlim([2 54])
end
ccc
%% BF Anova 4 factors
addpath(genpath('bayesFactor-master'))
participants_info=tdfread(['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\participants.tsv']);
subjects_analysed=[1:25 27:40 42:45 47:58]; % already remove 41 becuase "data" does not have that
trials=1:5; % which trials to average
epoch=1:3;    % which epochs to average
net_feats={'_InStrgth','_OutStrgth','_SrcPassTim','_ClustCoef','_Eccent','_NodBtw'};

for nf=1:1
    feat_to_plot=net_feats{nf};
    for i=1:length(Features_labels)
        if contains(Features_labels{i},feat_to_plot)
            targ_feat(i)=1;
        else
            targ_feat(i)=0;
        end
    end
    targ_feat=logical(targ_feat);

    X_all_norm_1=squeeze(nanmean(nanmean(nanmean(X_all_norm(:,targ_feat,1,:,trials,epoch),6),5),1));
    X_all_norm_2=squeeze(nanmean(nanmean(nanmean(X_all_norm(:,targ_feat,2,:,trials,epoch),6),5),1));

    for conect=1:size(X_all_norm_1,1)
        data_tmp=abs([X_all_norm_1(conect,:)-X_all_norm_2(conect,:)]');

        pt=0;
        for P=subjects_analysed
            pt=pt+1;
            g1{pt,1}=participants_info.outcome(P,1);
            g2{pt,1}=participants_info.target(P,:);
            g3{pt,1}=participants_info.lesion_status(P,:);
            g4{pt,1}=participants_info.implant(P,:);
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
        Effects_all(nf,conect,:)=effects;
        [nf conect]
    end
end
