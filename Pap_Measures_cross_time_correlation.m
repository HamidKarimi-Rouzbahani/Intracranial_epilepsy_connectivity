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

ff=0;
for i=1:length(connectivities)
    for j=1:length(net_feats)
        ff=ff+1;
        Features_labels{ff}={[connectivities{i} net_feats{j}]};
    end
end

%% Evaluation of features
all_patients=[1:39 41:56]; % patient indices
X_all_norm_int=nan(30,104,2,length(all_patients),5,3);
X_all_norm_ict=nan(30,104,2,length(all_patients),5,3);

for inter_or_ict=1:2
    if inter_or_ict==1
        ictal_or_inter='interictal';
    elseif inter_or_ict==2
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
                % classificaiton is performed for each patient separately
                if strcmp(ictal_or_inter,'ictal')
                    load([num2str(Patient),'_Project2_connect_features_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.mat']);
                    load([num2str(Patient),'_Project2_data_for_PyConnectivity_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.mat'],'channels_for_connectiv_inds','channels_resctd_inds','channels_soz_inds');
                elseif strcmp(ictal_or_inter,'interictal')
                    load([num2str(Patient),'_Project2_connect_features_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.mat']);
                    load([num2str(Patient),'_Project2_data_for_PyConnectivity_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.mat'],'channels_for_connectiv_inds','channels_resctd_inds','channels_soz_inds');
                end

                % labels
                labels=channels_soz_inds(logical(channels_for_connectiv_inds)); % SOZ channels

                data_targ=Features_all(labels==1,:);
                data_non=Features_all(labels==0,:);

                X=abs([data_targ;data_non]);
                y=[ones(size(data_targ,1),1);zeros(size(data_non,1),1)];

                X_norm=(X-min(X))./(max(X)-min(X)); % normalisation within patient
                if inter_or_ict==1
                    Xtmp=X_norm(y==1,:);
                    X_all_norm_int(1:size(Xtmp,1),1:size(Xtmp,2),1,p,trial,epoch)=Xtmp;
                    Xtmp=X_norm(y==0,:);
                    X_all_norm_int(1:size(Xtmp,1),1:size(Xtmp,2),2,p,trial,epoch)=Xtmp;
                else
                    Xtmp=X_norm(y==1,:);
                    X_all_norm_ict(1:size(Xtmp,1),1:size(Xtmp,2),1,p,trial,epoch)=Xtmp;
                    Xtmp=X_norm(y==0,:);
                    X_all_norm_ict(1:size(Xtmp,1),1:size(Xtmp,2),2,p,trial,epoch)=Xtmp;
                end
            end
        end
    end
end
%% Plotting within metric
clc
close all
trials=1:5; % which trials to average
epoch=1:3;    % which epochs to average

ppr_ordrd_connectivities={'ANM','IGCI','CDS','RECI','CCE','DI',...
    'GD','PSI','DTF','DCOH','PDCOH','SGC',...
    'LMFIT'};

ppr_ordrds=[10 13 11 12 1 2 7 8 3 4 5 6 9];
net_feats_detailed={'in strength','out strength','source passage time','clustering coefficient','eccentricity','node betweenness'};

f=0;
net_feats={'_InStrgth','_OutStrgth','_SrcPassTim','_ClustCoef','_Eccent','_NodBtw'};

for nf=1:6
    f=f+1;
    feat_to_plot=net_feats{nf};
    for i=1:length(Features_labels)
        if contains(Features_labels{i},feat_to_plot)
            targ_feat(i)=1;
        else
            targ_feat(i)=0;
        end
    end
    targ_feat=logical(targ_feat);

    X_all_norm_SOZ_int=squeeze(nanmean(nanmean(nanmean(X_all_norm_int(:,targ_feat,1,:,trials,epoch),6),5),1));
    X_all_norm_NON_int=squeeze(nanmean(nanmean(nanmean(X_all_norm_int(:,targ_feat,2,:,trials,epoch),6),5),1));

    X_all_norm_SOZ_ict=squeeze(nanmean(nanmean(nanmean(X_all_norm_ict(:,targ_feat,1,:,trials,epoch),6),5),1));
    X_all_norm_NON_ict=squeeze(nanmean(nanmean(nanmean(X_all_norm_ict(:,targ_feat,2,:,trials,epoch),6),5),1));
    figure;

    for conect=1:13
        subplot(3,5,conect)
        differ_int=(X_all_norm_SOZ_int(ppr_ordrds(conect),:)-X_all_norm_NON_int(ppr_ordrds(conect),:));
        differ_ict=(X_all_norm_SOZ_ict(ppr_ordrds(conect),:)-X_all_norm_NON_ict(ppr_ordrds(conect),:));

        subplot(3,5,conect)


        x=differ_int';
        y=differ_ict';
        [rs,ps]=corr(x,y,'rows','complete');
        scatter(x',y',100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
        xlabel('\it Interictal')
        ylabel('\it Ictal')
        hold on;


        both_valued=~isnan(x.*y);
        [p,S] = polyfit(x(both_valued),y(both_valued),1);
        [y_fit,delta] = polyval(p,x(both_valued),S);
        plot(x(both_valued),y_fit,'k-')

        colour_str=['\color[rgb]{', num2str(0.3), ',', num2str(0.8), ',', num2str(0.3), '}'];
        if ps<0.01
            title({[ppr_ordrd_connectivities{conect},':\it',colour_str,' r = ',sprintf('%.2f',rs), ',\it',colour_str,' P < 0.01']})
        else
            title({[ppr_ordrd_connectivities{conect}, ':\it r = ',sprintf('%.2f',rs), ',\it P = ',sprintf('%.2f',ps)]})
        end
        set(gca,'fontsize',14)
    end
    sgtitle(['\it \Delta ',net_feats_detailed{nf},' (SOZ - NON-SOZ)'])
end

%% Plotting between interictal in-strength and ictal outstrength
clc
close all
trials=1:5; % which trials to average
epoch=1:3;    % which epochs to average

ppr_ordrd_connectivities={'ANM','IGCI','CDS','RECI','CCE','DI',...
    'GD','PSI','DTF','DCOH','PDCOH','SGC',...
    'LMFIT'};

ppr_ordrds=[10 13 11 12 1 2 7 8 3 4 5 6 9];
net_feats_detailed={'in strength','out strength','source passage time','clustering coefficient','eccentricity','node betweenness'};

f=0;
net_feats={'_InStrgth','_OutStrgth','_SrcPassTim','_ClustCoef','_Eccent','_NodBtw'};

f=f+1;
nf=1;
feat_to_plot=net_feats{nf};
for i=1:length(Features_labels)
    if contains(Features_labels{i},feat_to_plot)
        targ_feat(i)=1;
    else
        targ_feat(i)=0;
    end
end
targ_feat=logical(targ_feat);

X_all_norm_SOZ_int=squeeze(nanmean(nanmean(nanmean(X_all_norm_int(:,targ_feat,1,:,trials,epoch),6),5),1));
X_all_norm_NON_int=squeeze(nanmean(nanmean(nanmean(X_all_norm_int(:,targ_feat,2,:,trials,epoch),6),5),1));


nf=2;
feat_to_plot=net_feats{nf};
for i=1:length(Features_labels)
    if contains(Features_labels{i},feat_to_plot)
        targ_feat(i)=1;
    else
        targ_feat(i)=0;
    end
end
targ_feat=logical(targ_feat);

X_all_norm_SOZ_ict=squeeze(nanmean(nanmean(nanmean(X_all_norm_ict(:,targ_feat,1,:,trials,epoch),6),5),1));
X_all_norm_NON_ict=squeeze(nanmean(nanmean(nanmean(X_all_norm_ict(:,targ_feat,2,:,trials,epoch),6),5),1));
figure;

for conect=1:13
    subplot(3,5,conect)
    differ_int=(X_all_norm_SOZ_int(ppr_ordrds(conect),:)-X_all_norm_NON_int(ppr_ordrds(conect),:));
    differ_ict=(X_all_norm_SOZ_ict(ppr_ordrds(conect),:)-X_all_norm_NON_ict(ppr_ordrds(conect),:));

    subplot(3,5,conect)


    x=differ_int';
    y=differ_ict';
    [rs,ps]=corr(x,y,'rows','complete');
    scatter(x',y',100,'MarkerEdgeColor',[0.3 0.3 0.95],'MarkerFaceColor',[0 0 0])
    xlabel('\it Interictal (\Delta in strength)')
    ylabel('\it Ictal (\Delta out strength)')
    hold on;


    both_valued=~isnan(x.*y);
    [p,S] = polyfit(x(both_valued),y(both_valued),1);
    [y_fit,delta] = polyval(p,x(both_valued),S);
    plot(x(both_valued),y_fit,'k-')

    colour_str=['\color[rgb]{', num2str(0.3), ',', num2str(0.8), ',', num2str(0.3), '}'];
    if ps<0.01
        title({[ppr_ordrd_connectivities{conect},':\it',colour_str,' r = ',sprintf('%.2f',rs), ',\it',colour_str,' P < 0.01']})
    else
        title({[ppr_ordrd_connectivities{conect}, ':\it r = ',sprintf('%.2f',rs), ',\it P = ',sprintf('%.2f',ps)]})
    end
    set(gca,'fontsize',14)
end
