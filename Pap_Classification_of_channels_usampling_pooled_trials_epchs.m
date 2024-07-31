clc
clear all;
close all;
%% Classification: Within ictal or interictal
ictal_or_inter='interictal';
all_patients=[1:39 41:56]; % patient indices

addpath(genpath('F:\RESEARCH\Hamid\Matlab\Tools\dkbsl-matlab_smote-41e4193'))

p=0;
for Patient=all_patients % patient indices
    if strcmp(ictal_or_inter,'interictal')
        if Patient==2
            trials=2;
        else
            trials=1:2;
        end
    else
        if (Patient==2 || Patient==6 || Patient==27)
            trials=1;
        elseif (Patient==9 || Patient==10 || Patient==15 || Patient==28 || Patient==45 || Patient==50)
            trials=1:2;
        elseif (Patient==1 || Patient==3 || Patient==4 || Patient==5 || Patient==11 || Patient==13 || Patient==16 || Patient==21 || Patient==22 || Patient==30 || Patient==31 || Patient==33 || Patient==35 || Patient==41 || Patient==43 || Patient==44 || Patient==49 || Patient==56)
            trials=1:3;
        elseif (Patient==7 || Patient==12 || Patient==20 || Patient==23 || Patient==24 || Patient==46)
            trials=1:4;
        else
            trials=1:5;
        end
    end
    p=p+1;
    X_all=[];
    y_all=[];
    for trial=trials
        for epoch=1:3
            % classificaiton is performed for each patient separately
            if strcmp(ictal_or_inter,'ictal')
                load([num2str(Patient),'_Project2_connect_features_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.mat']);
                load([num2str(Patient),'_Project2_data_for_PyConnectivity_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.mat'],'channels_for_connectiv_inds','channels_resctd_inds','channels_soz_inds');
            elseif strcmp(ictal_or_inter,'interictal')
                load([num2str(Patient),'_Project2_connect_features_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.mat']);
                load([num2str(Patient),'_Project2_data_for_PyConnectivity_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.mat'],'channels_for_connectiv_inds','channels_resctd_inds','channels_soz_inds');
            end

            % data and labels
            labels=channels_soz_inds(logical(channels_for_connectiv_inds)); % labels: SOZ channels

            data_targ=Features_all(labels==1,:);
            data_non=Features_all(labels==0,:);
            X=abs([data_targ;data_non]);
            y=[ones(size(data_targ,1),1);zeros(size(data_non,1),1)];
            X_all=vertcat(X_all,X);
            y_all=vertcat(y_all,y);
        end
    end
% ccc
    X=normalize(X_all);
    % X=(X_all);
    y=y_all;

    for iter_equalis=1:10
        %% Classification
        % randomising class labels (to generate null distribution for statistical
        % testing)?: No, we do it in another file; so iter_rand=1
        for iter_rand=1:1
            if iter_rand~=1
                y_r=randsample(y,length(y));
            else
                y_r=y;
            end
            % classify the contacts using the decision tree
            CVO = cvpartition(y_r,'k',5); % 5 fold cross-validation
            for fld = 1:CVO.NumTestSets
                trIdx = CVO.training(fld);
                teIdx = CVO.test(fld);

                X_train = X(trIdx, :);
                Y_train = y_r(trIdx);
                X_test = X(~trIdx, :);
                Y_test = y_r(~trIdx);

                if length(unique(Y_train))>1

                    % balancing classes through upsampling
                    classes=[0 1];
                    c=0;
                    for classLabel = classes
                        c=c+1;
                        classIndices(c) = sum(Y_train == classLabel);
                    end
                    % balancing classses through upsampling
                    if classIndices(1)~=classIndices(2)
                        [~,tmp]=min(classIndices);
                        min_class_label=classes(tmp);
                        classes(tmp)=[];
                        maj_class_label=classes;
                        maj_class=X_train(Y_train==maj_class_label,:);
                        min_class_upsampled=repmat(X_train(find(Y_train==min_class_label),:),[ceil(max(classIndices)./min(classIndices)) 1]);
                        min_class=min_class_upsampled(1:size(maj_class,1),:);
                        X_train_f=[min_class;maj_class];
                        Y_train_f=[repmat(min_class_label,[size(min_class,1) 1]);repmat(maj_class_label,[size(maj_class,1) 1])];
                    else
                        Y_train_f=Y_train;
                        X_train_f=X_train;
                    end

                    % normalising the data within the training ad testing
                    % sets separately to avoid circular analysis
                    Classifier_Model = TreeBagger(50,normalize(X_train_f),Y_train_f,...
                        Method="classification",...
                        OOBPrediction="on",OOBPredictorImportance="on");
                    impCART{p,iter_equalis,iter_rand,fld} = Classifier_Model.OOBPermutedPredictorDeltaError;
                    preds_tmp=predict(Classifier_Model,normalize(X_test));

                    for i=1:length(preds_tmp)
                        preds(i,1)=str2num(preds_tmp{i});
                    end
                    Predictions{p,iter_equalis,iter_rand,fld}=preds;
                    classif=nanmean(preds==Y_test);
                    Ground_truth{p,iter_equalis,iter_rand,fld}=Y_test;
                    clearvars preds
                else
                    Predictions{p,iter_equalis,iter_rand,fld}=nan;
                    Ground_truth{p,iter_equalis,iter_rand,fld}=nan;
                    impCART{p,iter_equalis,iter_rand,fld} = nan;
                end
            end
        end
        clearvars y_r
    end
    if strcmp(ictal_or_inter,'ictal')
        save(['Connectiv_classif_oversamp_pooled_Seizure.mat'],'Predictions','Ground_truth','impCART','-v7.3');
    elseif strcmp(ictal_or_inter,'interictal')
        save(['Connectiv_classif_oversamp_pooled_Interictal.mat'],'Predictions','Ground_truth','impCART','-v7.3');
    end
    [Patient]
end
ccc
%% Plotting classifications
clc;
clear all;
close all;
addpath(genpath('bayesFactor-master'))
ictal_or_inter='ictal';
t=0;
if strcmp(ictal_or_inter,'ictal')
    load(['Connectiv_classif_oversamp_pooled_Seizure.mat']);
elseif strcmp(ictal_or_inter,'interictal')
    load(['Connectiv_classif_oversamp_pooled_Interictal.mat']);
end
clearvars -except trls trial epoch t ictal_or_inter Predictions Ground_truth impCART
t=t+1;
ff=1; % 1= all feaures included
iter_rand=1; % 1= number of iterations
for p=1:size(Ground_truth,1) % patients
    for fld=1:size(Ground_truth,4) % folds
        for iter_equalis=1:size(Ground_truth,2) % iterations
            Performance_within(p,iter_equalis,iter_rand,fld,1:7)=Evaluate(Ground_truth{p,iter_equalis,iter_rand,fld},Predictions{p,iter_equalis,iter_rand,fld});
            try
                [~,~,~,Performance_within(p,iter_equalis,iter_rand,fld,8)]=perfcurve(Ground_truth{p,iter_equalis,iter_rand,fld},Predictions{p,iter_equalis,iter_rand,fld},1);
            catch
                Performance_within(p,iter_equalis,iter_rand,fld,8)=nan;
            end
        end
    end
end



Performance_within=squeeze(nanmean(nanmean(Performance_within(:,:,:,:,:),4),2));

Measures={'accuracy','sensitivity','specificity','precision','recall','f1-measure','gmean','auc'};
metric=8; %AUC

datas=Performance_within(:,metric);
colours={[0 0 0]};
datas_tmp_outed=nan(55,1);
tmpp=rmoutliers(datas(:,1));
datas_tmp_outed(1:length(tmpp),1)=tmpp;
swarmchart([ones(length(tmpp),1)],tmpp,'MarkerFaceColor',colours{1},'MarkerEdgeColor',colours{1},'MarkerFaceAlpha',0.3,'MarkerEdgeAlpha',0.3);
hold on;
boxplot(datas_tmp_outed,'Whisker',inf,'color','k');
xticks([1])
xticklabels([])
ylabel('\it Localisation (AUC)')
if strcmp(ictal_or_inter,'ictal')
    title([,'Ictal median = ',sprintf('%0.2f',nanmedian(datas_tmp_outed))])
else
    title([,'Interictal median = ',sprintf('%0.2f',nanmedian(datas_tmp_outed))])
end
box off
grid on
plot([0.5 2.5],[0.5 0.5],'--k')
clearvars datas datas_tmp_outed
ylim([0.4 1])
set(gca,'TickDir','out','Fontsize',18)

%% Features importance all
% clc
% clear all
% close all
% ictal_or_inter='ictal';
% 
% connectivities={'di',...
%     'dtf','dcoh','pdcoh',...
%     'sgc','gd','psi','lmfit',...
%     'gpfit','anm','cds','reci',...
%     'igci'};
% 
% 
% net_feats={'-InS','-OutS','-Tim','-ClsCof','-Eccnt','-NodBtw'};
% 
% ff=0;
% for i=1:length(connectivities)
%     for j=1:length(net_feats)
%         ff=ff+1;
%         Features_labels{ff}=[connectivities{i} net_feats{j}];
%     end
% end
% 
% t=0;
% if strcmp(ictal_or_inter,'ictal')
%     load(['Connectiv_classif_oversamp_pooled_Seizure.mat']);
% elseif strcmp(ictal_or_inter,'interictal')
%     load(['Connectiv_classif_oversamp_pooled_Interictal.mat']);
% end
% 
% feats=1:size(impCART{1,1,1,1},2);
% 
% ff=1;
% iter_rand=1;
% iter_equalis=1;
% feat_imp=nan(size(impCART,1),size(Features_labels,2),size(impCART,4));
% for p=1:size(impCART,1)
%     for fld=1:size(impCART,4)
%         for iter_equalis=1:size(impCART,2) % iterations
%             if size(impCART{p,iter_equalis,iter_rand,fld},1)~=0
%                 feat_imp(p,:,fld)=impCART{p,iter_equalis,iter_rand,fld};
%             else
%                 feat_imp(p,:,fld)=nan;
%             end
%         end
%     end
% end
% 
% datas=squeeze(nanmean(feat_imp,3));
% 
% datas_tmp_outed=nan(size(impCART,1),size(Features_labels,2));
% rng(6)
% a=jet;
% colours_tmp=a(randi(size(a,1),[1 13]),:);
% colours=[];
% j=0;
% for i=1:size(datas,2)
%     if mod(i,size(net_feats,2))==1
%         j=j+1;
%         colours=vertcat(colours,repmat(colours_tmp(j,:),[size(net_feats,2) 1]));
%     end
% end
% 
% c=0;
% for i=1:size(datas,2)
%     c=c+1;
%     tmpp=rmoutliers(datas(:,i));
%     datas_tmp_outed(1:length(tmpp),c)=tmpp;
%     swarmchart([c*ones(size(tmpp,1),1)]',tmpp,20,'MarkerFaceColor',colours(i,:),'MarkerEdgeColor',colours(i,:),'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
%     hold on;
%     clearvars tmpp
% end
% 
% boxplot(datas_tmp_outed,'Whisker',inf,'color','k');
% box off
% grid on
% ylabel('Feature Contribution (A.U.)')
% xticks([1:length(feats)])
% xticklabels([Features_labels])
% xlim([0 size(Features_labels,2)+1])
% title([ictal_or_inter])
% set(gca,'TickDir','out','Fontsize',11)

%% Features importance collapsed 
clc
clear all
close all
ictal_or_inter='ictal';

connectivities={'cce','di',...
    'dtf','dcoh','pdcoh',...
    'sgc','gd','psi','lmfit',...
    'anm','cds','reci',...
    'igci'};

ppr_ordrd_connectivities={'ANM','IGCI','CDS','RECI','CCE','DI',...    
                          'GD','PSI','DTF','DCOH','PDCOH','SGC',...
                          'LMFIT'};

ppr_ordrds=[10 13 11 12 1 2 7 8 3 4 5 6 9];

net_feats={'InS','OutS','Tim','ClsCof','Eccnt','NodBtw'};
net_feats_detailed={'in strength','out strength','source passage time','clustering coefficient','eccentricity','node betweenness'};

ff=0;
for i=1:length(connectivities)
    for j=1:length(net_feats)
        ff=ff+1;
        Features_labels{ff}=[connectivities{i} net_feats{j}];
    end
end

t=0;
if strcmp(ictal_or_inter,'ictal')
    load(['Connectiv_classif_oversamp_pooled_Seizure.mat']);
elseif strcmp(ictal_or_inter,'interictal')
    load(['Connectiv_classif_oversamp_pooled_Interictal.mat']);
end

feats=1:size(impCART{1,1,1,1},2);

ff=1;
iter_rand=1;
iter_equalis=1;
feat_imp=nan(size(impCART,1),size(Features_labels,2),size(impCART,4));
for p=1:size(impCART,1)
    for fld=1:size(impCART,4)
        for iter_equalis=1:size(impCART,2) % iterations
            if size(impCART{p,iter_equalis,iter_rand,fld},1)~=0
                feat_imp(p,:,fld)=impCART{p,iter_equalis,iter_rand,fld};
            else
                feat_imp(p,:,fld)=nan;
            end
        end
    end
end

datas=squeeze(nanmean(feat_imp,3));

subplot(121)
for cat=1:length(net_feats)
    tmpp=[];
    inds_inc=[];
    for i=1:size(datas,2)
        if mod(i,length(net_feats))==cat
            inds_inc=horzcat(inds_inc,i);
        elseif cat==6 && mod(i,length(net_feats))==0
            inds_inc=horzcat(inds_inc,i);
        end
    end
    datas_new(:,cat)=reshape(datas(:,inds_inc),[],1);
end
datas_tmp_outed=nan(size(impCART,1),size(datas_new,2));

colours=[0.5 0.5 0.5];
c=0;
for i=1:size(datas_new,2)
    c=c+1;
    tmpp=rmoutliers(datas_new(:,i));
    datas_tmp_outed(1:length(tmpp),c)=tmpp;
    swarmchart([c*ones(size(tmpp,1),1)]',tmpp,20,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
    hold on;
    clearvars tmpp
end

boxplot(datas_tmp_outed,'Whisker',inf,'color','k');
box off
grid on
ylabel('\it Feature contribution (a.u.)')
xticks([1:length(feats)])
xticklabels([net_feats_detailed])
xlim([0 size(datas_new,2)+1])
ylim([-0.2 1])
title([ictal_or_inter])
set(gca,'TickDir','out','Fontsize',18)

% BF matrices
BFs_metrics=nan(size(datas_tmp_outed,2));
for i=1:size(datas_tmp_outed,2)
    for j=i+1:size(datas_tmp_outed,2)
        BFs_metrics(i,j)=bf.ttest2(datas_tmp_outed(:,i),datas_tmp_outed(:,j));
    end
end


% Collapsing across graph metrics
datas_new=[];
datas=squeeze(nanmean(nanmean(nanmean(feat_imp,5),4),3));
subplot(122)
for cat=1:length(connectivities)
    inds_inc=length(net_feats)*(cat-1)+1:length(net_feats)*(cat);
    datas_new(:,cat)=reshape(datas(:,inds_inc),[],1);    
end

datas_tmp_outed=nan(size(impCART,1),size(datas_new,2));
colours=[0.5 0.5 0.5];
c=0;
for i=1:size(datas_new,2)
    c=c+1;
    tmpp=rmoutliers(datas_new(:,ppr_ordrds(i)));
    datas_tmp_outed(1:length(tmpp),c)=tmpp;
    swarmchart([c*ones(size(tmpp,1),1)]',tmpp,20,'MarkerFaceColor',colours,'MarkerEdgeColor',colours,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.3);
    hold on;
    clearvars tmpp
end

boxplot(datas_tmp_outed,'Whisker',inf,'color','k');

box off
grid on
ylabel('\it Feature contribution (a.u.)')
xticks([1:length(feats)])
xticklabels([ppr_ordrd_connectivities])
xlim([0 size(datas_new,2)+1])
ylim([-0.2 1])
title([ictal_or_inter])
set(gca,'TickDir','out','Fontsize',18)

% BF matrices
BFs_connects=nan(size(datas_tmp_outed,2));
for i=1:size(datas_tmp_outed,2)
    for j=i+1:size(datas_tmp_outed,2)
        BFs_connects(i,j)=bf.ttest2(datas_tmp_outed(:,i),datas_tmp_outed(:,j));
    end
end


