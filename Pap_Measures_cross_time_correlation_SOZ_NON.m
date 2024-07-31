% This code classifies resected and non-resected contacts within each patient
clc
clear all
close all
addpath(genpath('bayesFactor-master'))

%% Featrues
% without corr, te, gpfit and ccm
connectivities={'cce','di',...
    'dtf','dcoh','pdcoh',...
    'sgc','gd','psi','lmfit',...
    'anm','cds','reci',...
    'igci'};

net_feats={'_InStrgth','_OutStrgth','_SrcPassTim','_ClustCoef','_Eccent','_NodBtw'};

ff=0;
for i=1:length(connectivities)
    for j=1:length(net_feats)
        ff=ff+1;
        Features_labels{ff}={[connectivities{i} net_feats{j}]};
    end
end

%% Evaluation of features

all_patients=[1:39 41:56]; % patient indices

X_all_ictal=nan(78,2,length(all_patients),5,3);
X_all_inter=nan(78,2,length(all_patients),5,3);

for wind=1:2
    if wind==1
        ictal_or_inter='interictal';
    elseif wind==2
        ictal_or_inter='ictal';
    end
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
            
            % X=abs([data_targ;data_non]);
            % y=[ones(size(data_targ,1),1);zeros(size(data_non,1),1)];
            % X_norm=(X-min(X))./(max(X)-min(X));
            % 
            % if wind==1
            %     Xtmp=X_norm(y==1,:);
            %     X_all_inter(1:78,1,p,trial,epoch)=Xtmp;
            %     Xtmp=X_norm(y==0,:);
            %     X_all_inter(1:78,2,p,trial,epoch)=Xtmp;
            % else
            %     Xtmp=X_norm(y==1,:);
            %     X_all_ictal(1:78,1,p,trial,epoch)=Xtmp;
            %     Xtmp=X_norm(y==0,:);
            %     X_all_ictal(1:78,2,p,trial,epoch)=Xtmp;
            % end

            %% Unnormalised
            if wind==1
                X_all_inter(1:78,1,p,trial,epoch)=data_targ;
                X_all_inter(1:78,2,p,trial,epoch)=data_non;
            else
                X_all_ictal(1:78,1,p,trial,epoch)=data_targ;
                X_all_ictal(1:78,2,p,trial,epoch)=data_non;
            end

            end
        end
    end
end

%% Over subjects separated by trials and epochs
trials=1:5; % which trials to average
epoch=1:3;    % which epochs to average
net_feats={'_InStrgth','_OutStrgth','_SrcPassTim','_ClustCoef','_Eccent','_NodBtw'};
for wind=1:2
    for nf=1
        feat_to_plot=net_feats{nf};
        for i=1:length(Features_labels)
            if contains(Features_labels{i},feat_to_plot)
                targ_feat(i)=1;
            else
                targ_feat(i)=0;
            end
        end
        targ_feat=logical(targ_feat);
        if wind==1
            X_all_SOZ_inter=squeeze(nanmean(nanmean(X_all_inter(targ_feat,1,:,trials,epoch),5),4));
            X_all_NON_inter=squeeze(nanmean(nanmean(X_all_inter(targ_feat,2,:,trials,epoch),5),4));
        else
            X_all_SOZ_ictal=squeeze(nanmean(nanmean(X_all_ictal(targ_feat,1,:,trials,epoch),5),4));
            X_all_NON_ictal=squeeze(nanmean(nanmean(X_all_ictal(targ_feat,2,:,trials,epoch),5),4));
        end
    end
end

%% Correlations
connect_to_plot=1:length(connectivities);
ppr_ordrd_connectivities={'ANM','IGCI','CDS','RECI','CCE','DI',...    
                          'GD','PSI','DTF','DCOH','PDCOH','SGC',...
                          'LMFIT'};

ppr_ordrds=[10 13 11 12 1 2 7 8 3 4 5 6 9];
net_feats_detailed={'in strength','out strength','source passage time','clustering coefficient','eccentricity','node betweenness'};

for i=1:length(connectivities)
    subplot(3,5,i)
    x=X_all_SOZ_inter(ppr_ordrds(i),:)'-X_all_NON_inter(ppr_ordrds(i),:)';
    y=X_all_SOZ_ictal(ppr_ordrds(i),:)'-X_all_NON_ictal(ppr_ordrds(i),:)';

    [rs,ps]=corr(x,y);

    scatter(x,y,100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
    xlabel('Interictal')
    ylabel('Ictal')
    hold on;

    [p,S] = polyfit(x,y,1);
    [y_fit,delta] = polyval(p,x,S);
    plot(x,y_fit,'k-')

    colour_str=['\color[rgb]{', num2str(0.3), ',', num2str(0.8), ',', num2str(0.3), '}'];
    if ps<0.01
        title({[ppr_ordrd_connectivities{i},':\it',colour_str,' r = ',sprintf('%.2f',rs), ',\it',colour_str,' P < 0.01']})
    else
        title({[ppr_ordrd_connectivities{i}, ':\it r = ',sprintf('%.2f',rs), ',\it P = ',sprintf('%.2f',ps)]})
    end
    set(gca,'fontsize',14)
end
    sgtitle(['\it ',net_feats_detailed{1}])
