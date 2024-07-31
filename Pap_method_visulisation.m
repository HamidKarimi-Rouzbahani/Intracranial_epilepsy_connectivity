% This code extracts features from the interictal and ictal recordings
% and saves the extracted features in new files
% Some features are extracted using external functions which should be
% available in paths accissible to Matlab

% INPUTS: epoched and cleaned data from C1_ds4100_ictal_interictal_epoching
% OUTPUTS: features data to be used by C3_Separating_target_non_target_contacts_all_feats
%% Feature extraction
clc
clear all;
close all;
% adding eeglab toolbox
addpath(genpath('F:\Toolbox\eeglab2021.1'))
ictal_or_inter='interictal'; % 'ictal' or 'interictal'

Patient=23
if Patient==1
    Patient_initials='060';
    modality='seeg';
    ictal_trials=[1:3];
    inter_ictal_trials=[1:2];
elseif Patient==2
    Patient_initials='064';
    ictal_trials=[1];
    modality='ecog';
    inter_ictal_trials=[2];
elseif Patient==3
    Patient_initials='065';
    modality='ecog';
    ictal_trials=[1:3];
    inter_ictal_trials=[1:2];
elseif Patient==4
    Patient_initials='070';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==5
    Patient_initials='074';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==6
    Patient_initials='075';
    ictal_trials=[1];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==7
    Patient_initials='080';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==8
    Patient_initials='082';
    ictal_trials=[1:5];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==9
    Patient_initials='086';
    ictal_trials=[1:2];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==10
    Patient_initials='087';
    ictal_trials=[1:2];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==11
    Patient_initials='088';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==12
    Patient_initials='089';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==13
    Patient_initials='094';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==14
    Patient_initials='097';
    ictal_trials=[1:5];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==15
    Patient_initials='105';
    ictal_trials=[1:2];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==16
    Patient_initials='106';
    ictal_trials=[1:3];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==17
    Patient_initials='107';
    ictal_trials=[1:5];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==18
    Patient_initials='111';
    ictal_trials=[1:5];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==19
    Patient_initials='112';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==20
    Patient_initials='114';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==21
    Patient_initials='116';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==22
    Patient_initials='117';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==23
    Patient_initials='123';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==24
    Patient_initials='126';
    ictal_trials=[1:4];
    modality='ecog';
    inter_ictal_trials=[1:2];
elseif Patient==25
    Patient_initials='130';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==26
    Patient_initials='133';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==27
    Patient_initials='134';
    ictal_trials=[1];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==28
    Patient_initials='135';
    ictal_trials=[1:2];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==29
    Patient_initials='138';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==30
    Patient_initials='139';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==31
    Patient_initials='140';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==32
    Patient_initials='141';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==33
    Patient_initials='142';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==34
    Patient_initials='144';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==35
    Patient_initials='146';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==36
    Patient_initials='148';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==37
    Patient_initials='150';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==38
    Patient_initials='151';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==39
    Patient_initials='157';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==40
    Patient_initials='158';
    ictal_trials=[4];
    modality='seeg';
    inter_ictal_trials=[];
elseif Patient==41
    Patient_initials='160';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==42
    Patient_initials='162';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==43
    Patient_initials='163';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==44
    Patient_initials='164';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==45
    Patient_initials='166';
    ictal_trials=[1:2];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==46
    Patient_initials='171';
    ictal_trials=[1:4];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==47
    Patient_initials='172';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==48
    Patient_initials='173';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==49
    Patient_initials='177';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==50
    Patient_initials='179';
    ictal_trials=[1:2];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==51
    Patient_initials='180';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==52
    Patient_initials='181';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==53
    Patient_initials='185';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==54
    Patient_initials='187';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==55
    Patient_initials='188';
    ictal_trials=[1:5];
    modality='seeg';
    inter_ictal_trials=[1:2];
elseif Patient==56
    Patient_initials='190';
    ictal_trials=[1:3];
    modality='seeg';
    inter_ictal_trials=[1:2];
end

% Loading the data file from C0
if strcmp(ictal_or_inter,'ictal')
    trials=ictal_trials;
elseif strcmp(ictal_or_inter,'interictal')
    trials=inter_ictal_trials;
end

% loading each trial/epoch/recodings
trial=1;
if strcmp(ictal_or_inter,'ictal')
    load([Patient_initials,'_SEEG_Seizure',num2str(trial),'.mat'])
    Fs=header.Fs; % loading sampling frequeny
    onset=event_info.sample(1);
    offset=event_info.sample(2);
    signal_tmp=signal(:,onset-4*Fs:offset+4*Fs); % trying to avoid edges
elseif strcmp(ictal_or_inter,'interictal')
    load([Patient_initials,'_SEEG_InterIctal',num2str(trial),'.mat'])
    Fs=header.Fs; % loading sampling frequency
    pre_siz_analysis=1;  %/FS ms; 0s
    post_siz_analysis=60000*(Fs/500); %/FS ms; ~120s
    signal_tmp=signal(:,pre_siz_analysis:post_siz_analysis);
end

%% Downsampling and filtering
signal_ff=[];
fsNew=256; % sampling frequency
c=0;
for ch=1:size(signal_tmp,1)
    c=c+1;
    downsampled=resample(signal_tmp(ch,:),fsNew,Fs);
    signal_ff(c,:) = bandstop(downsampled,[55 65],fsNew);
end

%% Removing extra charecters from channel names: both spike detected and original ones
% From the channels included in spike detection
for ch=1:length(channels_id)
    newStr = erase(channels_id{1,ch},"EEG ");
    channels_id_pruned{1,ch} = erase(newStr,"-Ref");
end
clearvars newStr
% now from the original data labels
for ch=1:size(electrode_loc_info.x0xEF0xBB0xBFname,1)
    newStr = erase((electrode_loc_info.x0xEF0xBB0xBFname(ch,:)),"EEG ");
    rsd_str=erase(newStr,"-Ref");
    electrode_loc_info_pruned.x0xEF0xBB0xBFname(ch,1:length(rsd_str)) = rsd_str;
end
channels_id=channels_id_pruned;
electrode_loc_info.x0xEF0xBB0xBFname=electrode_loc_info_pruned.x0xEF0xBB0xBFname;
clearvars channels_id_pruned electrode_loc_info_pruned

%% First finding channels used in analysis
chans_analysed=[];
for i=1:size(electrode_loc_info.x,1)
    exists=0;
    for j=1:size(channels_id,2)
        if strcmp(strtrim(electrode_loc_info.x0xEF0xBB0xBFname(i,:)),channels_id{1,j})
            exists=1;
        end
    end
    if ~isempty(str2num(electrode_loc_info.x(i,:))) && exists==1
        chans_analysed(i)=1;
    else
        chans_analysed(i)=0;
    end
end
channels=chans_analysed;

%% Finding indices of resected/SOZ channels from the orignal data
patient_data_address=['F:\RESEARCH\Hamid\Multicentre dataset\ds004100\sub-HUP',Patient_initials,'\ses-presurgery\ieeg\'];
folder_address=[patient_data_address,'sub-HUP',Patient_initials,'_ses-presurgery_task-',ictal_or_inter,'_acq-',modality,'_run-',sprintf('%02d',trial)];
channel_info=tdfread([folder_address,'_channels.tsv']);
c=0;
d=0;
channels_soz_inds_tmpp=[];
channels_resctd_inds_tmpp=[];
tmp_chan=find(chans_analysed);
for current_elec=1:size(signal_ff,1)
    chan_currnt=tmp_chan(current_elec);
    if contains(channel_info.status_description(chan_currnt,:),'resect')
        c=c+1;
        channels_resctd_inds_tmpp(c)=current_elec;
    end
    if contains(channel_info.status_description(chan_currnt,:),'soz')
        d=d+1;
        channels_soz_inds_tmpp(d)=current_elec;
    end
end

channels_resctd_inds=zeros(1,size(channels_id,2));
channels_resctd_inds(channels_resctd_inds_tmpp)=1;

channels_soz_inds=zeros(1,size(channels_id,2));
channels_soz_inds(channels_soz_inds_tmpp)=1;

clearvars tmp_chan


%% Detecting white matter contacts
addpath(genpath('F:\RESEARCH\Hamid\Multicentre dataset\Scripts\Coordi_to_atlases'))
addpath(genpath('F:\RESEARCH\Hamid\Multicentre dataset\Scripts\NIfTI_20140122'))

channels_gray_inds=[];

tmp_chan=find(chans_analysed);
for current_elec=1:size(signal_ff,1)
    chan_currnt=tmp_chan(current_elec);
    currnt_elec_coord=[str2num(electrode_loc_info.x(chan_currnt,:)) str2num(electrode_loc_info.y(chan_currnt,:)) str2num(electrode_loc_info.z(chan_currnt,:))];
    try
        [atlas]=mni2atlas(currnt_elec_coord,[2 3],'1mm');
        if ~isempty(atlas(1).label)
            Atlas_location{current_elec}=atlas(1).label{1,1};
        end
        if ~isempty(atlas(2).label)
            GW_matter{current_elec}=atlas(2).label{1,1};
            if contains(atlas(2).label{1,1},'White')
                channels_gray_inds(current_elec)=0;
            else
                channels_gray_inds(current_elec)=1;
            end
        else
            channels_gray_inds(current_elec)=0;
        end
    catch
        channels_gray_inds(current_elec)=0;
    end
end

%% Only keeping channels which have further distances
min_dist=10; % (mm) only channels above this distance are kept

channels_distant_inds_tmpp=1;
tmp_chan=find(chans_analysed);
c=1;
for current_elec1=2:size(signal_ff,1)

    chan_currnt=tmp_chan(current_elec1);
    currnt_elec_coord1=[str2num(electrode_loc_info.x(chan_currnt,:)) str2num(electrode_loc_info.y(chan_currnt,:)) str2num(electrode_loc_info.z(chan_currnt,:))];

    cors=[];
    d=0;
    for current_elec2=channels_distant_inds_tmpp
        d=d+1;
        chan_currnt2=tmp_chan(current_elec2);
        currnt_elec_coord2=[str2num(electrode_loc_info.x(chan_currnt2,:)) str2num(electrode_loc_info.y(chan_currnt2,:)) str2num(electrode_loc_info.z(chan_currnt2,:))];
        dists(d)=pdist2(currnt_elec_coord2,currnt_elec_coord1);
    end
    if sum(dists<min_dist)==0 || channels_soz_inds(current_elec1)==1
        c=c+1;
        channels_distant_inds_tmpp(c)=current_elec1;
    end

end
channels_distant_inds=zeros(1,size(channels_id,2));
channels_distant_inds(channels_distant_inds_tmpp)=1;
clearvars tmp_chan dists channels_distant_inds_tmpp

%% plotting raw
scaling=600;
if strcmp(ictal_or_inter,'interictal')
    span=[size(signal_tmp,2)*1/3:size(signal_tmp,2)*2/3-1];
else
    span=[1:size(signal_tmp,2)];
end
gap=2;
for ch=1:size(signal_tmp,1)
    y=signal_tmp(ch,span);
    x=mean(y);
    if channels_soz_inds(ch)==1
        plot([1:length(y)]+gap*Fs,(y-x)+ch*scaling,'color',[0.9 0.3 0.3],'linewidth',2);
    elseif channels_gray_inds(ch)==0
        plot([1:length(y)]+gap*Fs,(y-x)+ch*scaling,'color',[0.3 0.3 0.9],'linewidth',2);
    else
        plot([1:length(y)]+gap*Fs,(y-x)+ch*scaling,'color',[0.3 0.3 0.3],'linewidth',2);
    end
    h=text(Fs/2.5,ch*scaling,channels_id{ch});
    if channels_soz_inds(ch)==1
        h.Color=[0.9 0.3 0.3];
    elseif channels_gray_inds(ch)==0
        h.Color=[0.3 0.3 0.9];
    else
        h.Color=[0.3 0.3 0.3];
    end
    set(h,'Rotation',10);
    hold on;
end
set(gca,'TickDir','out','Fontsize',16)
xlim([1 20000])
ylim([1 (size(signal_tmp,1)+1)*scaling])
box off
xticks([1 5000 10000 15000 20000])
xticklabels([0 5 60 90 120])
yticklabels([])
xlabel('Time (s)')
ylabel('\it Amplitude (a.u.)')

%% plotting filtered using electrodes use in analysis
trial=1;
epoch=1;
if strcmp(ictal_or_inter,'ictal')
    load([num2str(Patient),'_Project2_data_for_PyConnectivity_long20s_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.mat'],'channels_for_connectiv_inds');
elseif strcmp(ictal_or_inter,'interictal')
    load([num2str(Patient),'_Project2_data_for_PyConnectivity_long20s_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.mat'],'channels_for_connectiv_inds');
end
% channels_to_plot=sort([randsample(find(channels_gray_inds & channels_soz_inds==0),30-sum(channels_soz_inds==1)),find(channels_soz_inds)]);
channels_to_plot=sort(find(channels_for_connectiv_inds));
close all
scaling=600;
span=[size(signal_ff,2)*1/3:size(signal_ff,2)*2/3-1];
gap=2;
c=0;
for ch=channels_to_plot
    c=c+1;
    y=signal_ff(ch,span);
    x=mean(y);
    if channels_soz_inds(ch)==1
        plot([1:length(y)]+gap*fsNew,(y-x)+c*scaling,'color',[0.9 0.3 0.3],'linewidth',2);
    else
        plot([1:length(y)]+gap*fsNew,(y-x)+c*scaling,'color',[0.3 0.3 0.75],'linewidth',2);
    end
    h=text(round(Fs/2.5),c*scaling,channels_id{ch});
    h.Color=[0.3 0.3 0.3];
    set(h,'Rotation',10);
    hold on;
end
set(gca,'TickDir','out','Fontsize',16)
xlim([1 length(y)+gap*fsNew/Fs])
ylim([1 (length(channels_to_plot)+1)*scaling])
box off
xticks(round([1 5000 10000 15000 20000]*fsNew/Fs))
xticklabels([0 5 60 90 120])
yticklabels([])
xlabel('Time (s)')
ylabel('\it Amplitude (a.u.)')
ccc
%% Plotting connectivity matrix of first epoch for the same channels
clc
ictal_or_inter='interictal';

Patient=23

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


% without corr, te, gpfit and ccm
connectivities={'cce','di',...
    'dtf','dcoh','pdcoh',...
    'sgc','gd','psi','lmfit',...
    'anm','cds','reci',...
    'igci'};
trial=1
epoch=1
conect='pdcoh';
if strcmp(ictal_or_inter,'ictal')
    load([num2str(Patient),'_Project2_data_for_PyConnectivity_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.mat']);
elseif strcmp(ictal_or_inter,'interictal')
    load([num2str(Patient),'_Project2_data_for_PyConnectivity_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.mat'],...
        'signal_ff','Data','channels_for_connectiv_inds','channels_resctd_inds','channels_soz_inds','channels','channels_id','channels_distant_inds','channels_uncorltd_inds','channels_gray_inds','span','channel_info');
end
chans_used=find((channels_for_connectiv_inds));
chan_labels=channels_soz_inds(chans_used);


if strcmp(ictal_or_inter,'ictal')
    filePath=[num2str(Patient),'_',conect,'_Project2_PyConnectivity_Seizure',num2str(trial),'_epoch_',num2str(epoch),'.csv'];

elseif strcmp(ictal_or_inter,'interictal')
    filePath=[num2str(Patient),'_',conect,'_Project2_PyConnectivity_Interictal',num2str(trial),'_epoch_',num2str(epoch),'.csv'];
end

opts = detectImportOptions(filePath);
opts.EmptyLineRule = 'read'; % Specify how empty values should be treated (NaN in this case)
for i=1:size(opts.VariableTypes,2)
    opts.VariableTypes{1,i}='double';
end
opts.DataLines=[3,inf];
opts.VariableNamesLine=2;
T = readtable(filePath, opts);
TT=table2array(T);

% Preparing the data for connectivity analysis
if strcmp(conect,'corr')  ||  strcmp(conect,'te')  ||  strcmp(conect,'cce')  ||  strcmp(conect,'gd')  ||  strcmp(conect,'psi')
    TT = abs(TT);
elseif strcmp(conect,'igci')  ||  strcmp(conect,'lmfit')  ||  strcmp(conect,'gpfit')  ||  strcmp(conect,'cds')  ||  strcmp(conect,'reci')
    TT = weight_conversion(abs(TT), 'lengths');
end
% TT = weight_conversion(TT, 'normalize');

TT(isnan(TT))=0;
clearvars opts
% Removing all-zero columns
TT=TT(1:sum(channels_for_connectiv_inds),1:sum(channels_for_connectiv_inds));

TT_rearranged=nan(size(TT,1));
for i=1:size(TT,1)
    for j=1:size(TT,2)
        if i==j
            TT_rearranged(i,j)=nan;
        else
            TT_rearranged(size(TT,1)+1-i,size(TT,1)+1-j)=TT(i,j);
        end
    end
end

h=imagesc(TT_rearranged);
caxis([min(TT_rearranged(:)), max(TT_rearranged(~isnan(TT_rearranged)))]);
alphaData = ~isnan(TT_rearranged); % Alpha data: 1 (opaque) for non-NaN, 0 (transparent) for NaN
set(h, 'AlphaData', alphaData);
set(gca, 'Color', [1 1 1],'Fontsize',14);
hcb=colorbar;
hcb.Title.String = "Connectivity";

xticks([1:30])
yticks([1:30])
xticklabels([channels_id(sort(chans_used,'descend'))])
yticklabels([channels_id(sort(chans_used,'descend'))])
xlabel('\it From')
ylabel('\it To')
title('PDCOH')
