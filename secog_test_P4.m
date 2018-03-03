%% get the behavioral data
subjname = {'P2' , 'P4'};
subjNum = 2;
Dall=secog_analyze ('sing_subj' , 'subjCode' , 'P4');
mainDir = ['/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/' subjname{subjNum} ,'/'];

save([mainDir , 'AllData_Behav.mat'] , 'Dall')

%% Visualize behavior
out  = secog_visualize(Dall , subjNum, 'IPI_MT_seq', [], [] , [] , [1 2], [1]);
out  = secog_visualize(Dall , subjNum, 'IPI_MT_chunk', [], [] , [] , [1 2], [1]);

%% Pack EEG
secog_packEEG('PackEEG' , 2);
%% Calculate the PSDs
subjNum = 2;

% Good Channels for P2
% chans = [3:15 , 99:109 , 121:129 ,36];

% Good Channels for P4
chans = [15:23 25 80:88 90:99 101:110 112 28];


load('AllData_Behav.mat')

load('AllData_Events.mat')
Dout  = secog_parseEEG_PSD('ParseEEG-freqDecomp' , Dall, subjNum , 'Channels' , chans);
Dout  = secog_parseEEG_PSD('ParseEEG-calc_norm_PSD' , Dall, subjNum , 'Channels' , chans);
% time warp
Dout  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned' , Events, subjNum);
Dout  = secog_parseEEG_PSD('TimeWarpPSD_Raw_Binned_seqType' , Events, subjNum);
Dout  = secog_parseEEG_PSD('AlignEvents_SeqType' , Events, subjNum);
Dout  = secog_parseEEG_PSD('AlignEvents_SeqType_Warped' , Events, subjNum);


cd('/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/P4')
load('AllData_PSD_Warped.mat')
secog_BlockGroup(Pall,subjNum,'raw_BlockGroup', 'Channels' , chans)
%% average pattern calculation
load('AllData_PSD_Warped.mat')
secog_BlockGroup(Pall,subjNum,'binned_BlockGroup', 'Channels' , chans)
secog_BlockGroup(Pall,subjNum,'raw_BlockGroup', 'Channels' , chans)

load('AllData_PSD_Warped_SeqType.mat')
secog_BlockGroup(Pall,subjNum,'raw_SeqType', 'Channels' , chans)  
secog_BlockGroup(Pall,subjNum,'binned_SeqType', 'Channels' , chans)  



secog_BlockGroup([],subjNum,'raw_AvgPower_BlockGroup', 'Channels' , chans)
secog_BlockGroup([],subjNum,'raw_AvgPower_SeqType', 'Channels' , chans)
secog_BlockGroup([],subjNum,'raw_Power_SeqType', 'Channels' , chans)
secog_BlockGroup([],subjNum,'binned_Power_SeqType', 'Channels' , chans)
secog_BlockGroup([],subjNum,'Aligned_SeqType', 'Channels' , chans)
secog_BlockGroup([],subjNum,'AlignedWarped_SeqType', 'Channels' , chans)




%% plotting case summary
%      'binned_SingleTrial'
%      'binned_SingleTrial_AvgChann'
%      'binned_BlockGroup_AvgChann'
%      'binned_BlockGroup'
%      'raw_SingleTrial'
%      'raw_BlockGroup'
%      'raw_BlockGroup_AvgChann'

BG(1,:) = {[ ] , [2 8], [14 20 26], [29 38], [], [1 7],[13 19 25], [28 37], [] , [3:5] , [9:11] , [6 12] , [15:17] , [21:23] , [],...
    [18 24] , [30:32] , [34:36] , [], [27 33],[]}';
BG(2,:) = {'[]' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'[]',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', '[]' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , '[]', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , '[]', 'ChunkDay3', '[]'}';
BG = BG';



secog_visualizePSD(Pall,subjNum,'raw_BlockGroup' , 'BlockGroup' , 'SingleFingSlow1','Chan2Plot' , [12, 14],'Channels' , chans) % PMd
secog_visualizePSD(Pall,subjNum,'raw_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [27:30],'Channels' , chans) % SMA
load('AllData_PSD_Warped.mat')
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'Intermixed5','Chan2Plot' , [20 :22, 24, 12, 14],'Channels' , chans); % PMd
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow1','Chan2Plot' , [27:30],'Channels' , chans); % SMA

secog_visualizePSD(Pall,subjNum,'raw_BlockGroup_SeqType' , 'BlockGroup' , 'SingleFingSlow3','Chan2Plot' , [34:37],'Channels' , chans); % SMA
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup_SeqType' , 'BlockGroup' , 'SingleFingSlow3','Chan2Plot' , [34 :37],'Channels' , chans); % PMd




secog_visualizePSD(Pall,subjNum,'binned_BlockGroup_SeqType' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [27:30],'Channels' , chans); % SMA
secog_visualizePSD([],subjNum,'AvgPower_SeqType' , 'BlockGroup' , 'Intermixed5','Chan2Plot' , [27:30],'Channels' , chans); % SMA

secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed1','Intermixed4','Intermixed7'},'Chan2Plot' , [1:40],'Channels' , chans); % SMA\



secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},'Chan2Plot' , [1:3],'Channels' , chans); % SMA\
secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},'Chan2Plot' , [10 12],'Channels' , chans) % PMd

secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'SingleFingFast1','SingleFingFast2','SingleFingFast3'},'Chan2Plot' , [1:40],'Channels' , chans); % SMA\
secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'SingleFingSlow1','SingleFingSlow2','SingleFingSlow3'},'Chan2Plot' , [10 12],'Channels' , chans) % PMd


secog_visualizePSD([],subjNum,'Compare_SeqType_TemporalPattern' ,'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' , [27:30],'Channels' , chans , 'bandofInterest' , 6); % SMA\

secog_visualizePSD([],subjNum,'Compare_SeqType_TemporalPattern' ,'BlockGroup' , {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' , [10 12],'Channels' , chans , 'bandofInterest' , 4); % SMA\

secog_visualizePSD([],subjNum,'Compare_SeqType_TemporalPattern' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'},...
    'Chan2Plot' , [10 12],'Channels' , chans , 'bandofInterest' , 6); % SMA\

secog_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed4','Intermixed7'},...
    'Chan2Plot' , [1:3],'Channels' , chans , 'bandofInterest' , 4); % SMA\

secog_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3'},...
    'Chan2Plot' , [27:30],'Channels' , chans , 'bandofInterest' , 4); % SMA\

secog_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' },...
    'Chan2Plot' , [27:30],'Channels' , chans , 'bandofInterest' , 6); % SMA\

secog_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},...
    'Chan2Plot' , [27:30],'Channels' , chans , 'bandofInterest' , 6); % SMA\

secog_visualizePSD([],subjNum,'ChunkAligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed9'},...
    'Chan2Plot' , [27:30],'Channels' , chans , 'bandofInterest' , 6); % SMA\


secog_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},...
    'Chan2Plot' , [27:30],'Channels' , chans , 'bandofInterest' , 6); % SMA\
secog_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4'},...
    'Chan2Plot' , [27:30],'Channels' , chans , 'bandofInterest' , 6); % SMA\
secog_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3'},...
    'Chan2Plot' , [1:40],'Channels' , chans , 'bandofInterest' , 6); % SMA\
secog_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' , [27:30],'Channels' , chans , 'bandofInterest' , 6); % SMA\



secog_visualizePSD(Pall,subjNum,'raw_SingleTrial' , 'TBNum', [10 1], 'Chan2Plot' , [23 8],'Channels' , chans)

load('AllData_PSD_StimNorm.mat')
secog_visualizePSD(Pall,subjNum,'binned_SingleTrial' , 'TBNum', [10 4], 'Chan2Plot' , [23 8],'Channels' , chans)
secog_visualizePSD(Pall,subjNum,'binned_SingleTrial_AvgChann' , 'TBNum', [10 40], 'Chan2Plot' , [104:107 , 109])
secog_visualizePSD(Pall,subjNum,'raw_SingleTrial' , 'TBNum',  [10 4], 'Chan2Plot' , [33 38],'Channels' , chans)



%%
load('AllData_PSD_Warped.mat')
% PMD layer 2
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup_AvgChann' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [104:107 , 109])
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [104:107 , 109])
% PMD layer 1
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [12 14])
% PMv
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [6])
% SMA
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , [123:126])
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , blockGroupNames{3},'Chan2Plot' , [123:126])
%% sig test

secog_sigTest('seqs_across_days' , 1, 'Groups' , [1 2 3 4] , 'Chan' , [27:30]);% SMA
secog_sigTest('seqs_across_days' , 1, 'Groups' , [1 2 3 4] , 'Chan' , [12, 10])% PMD 
secog_sigTest('chunks_across_days' , 1, 'Groups' , [1 2 3] , 'Chan' , [27:30]) 
secog_sigTest('chunks_across_days' , 1, 'Groups' , [1 2 3] , 'Chan' , [12, 10])% PMD 
secog_sigTest('Alined_Sequences' , 1, 'Groups' , [1 2 3 4] , 'Chan' , [27:30])% PMD 
secog_sigTest('Alined_SingleFinger' , 1, 'Groups' , [1 2 3 4] , 'Chan' , [27:30])% PMD 