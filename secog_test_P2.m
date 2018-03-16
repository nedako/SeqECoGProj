%% Make all the data structures
% % subjNum = 1;
% % 
% % % Good Channels for P2
% % chans = [4:12 ,14, 122:129 ,36] ;
% % secog_makeData(subjNum , chans)

%% Visualize Behavior
subjNum = 1;

% Good Channels for P2
chans = [4:12 ,14, 122:129 ,36] ;

cd('/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/P2')
load('AllData_Behav.mat')

load('AllData_Events.mat')

out  = secog_visualize(Dall , subjNum, 'IPI_MT_seq' ,[1 2]);
out  = secog_visualize(Dall , subjNum, 'IPI_MT_chunk' ,[1 2]);
out  = secog_visualize(Dall , subjNum, 'RT_seq' ,[1 2]);
out  = secog_visualize(Dall , subjNum, 'RT_chunk' ,[1 2]);
out  = secog_visualize(Dall , subjNum, 'RT_SF' ,[1 2]);
out  = secog_visualize(Dall , subjNum, 'MT_SF' ,[1 2]);


%% plotting case summary
%      'binned_SingleTrial'
%      'binned_SingleTrial_AvgChann'
%      'binned_BlockGroup_AvgChann'
%      'binned_BlockGroup'
%      'raw_SingleTrial'
%      'raw_BlockGroup'
%      'raw_BlockGroup_AvgChann'


% LEFT SMA  --> CH = [1 , 11:18];
% LEFT PMd  --> CH = [2:5];
% LEFT PMv  --> CH = [6:10];

% ALL LEFT  --> CH = [1:18];

% ALL MOTOR --> CH = [1:18];

{'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';


secog_visualizePSD(Pall,subjNum,'raw_BlockGroup' , 'BlockGroup' , 'SingleFingSlow1','Chan2Plot' , CH,'Channels' , chans) 
load('AllData_PSD_Warped.mat')
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'Intermixed5','Chan2Plot' , CH,'Channels' , chans); 
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow1','Chan2Plot' , CH,'Channels' , chans); 

secog_visualizePSD(Pall,subjNum,'raw_BlockGroup_SeqType' , 'BlockGroup' , 'SingleFingSlow1','Chan2Plot' , CH,'Channels' , chans); 
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup_SeqType' , 'BlockGroup' , 'Intermixed5','Chan2Plot' , CH,'Channels' , chans); 




secog_visualizePSD(Pall,subjNum,'binned_BlockGroup_SeqType' , 'BlockGroup' , 'SingleFingFast1','Chan2Plot' , CH,'Channels' , chans);
secog_visualizePSD([],subjNum,'AvgPower_SeqType' , 'BlockGroup' , 'Intermixed5','Chan2Plot' , CH,'Channels' , chans);

secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'},'Chan2Plot' , CH,'Channels' , chans); 
secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'},'Chan2Plot' , CH,'Channels' , chans);

secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},'Chan2Plot' , CH,'Channels' , chans); 
secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},'Chan2Plot' , CH,'Channels' , chans) 



secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},'Chan2Plot' , CH,'Channels' , chans); 

secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'SingleFingSlow1','SingleFingSlow2','SingleFingSlow3','SingleFingSlow4'},'Chan2Plot' , CH,'Channels' , chans); 
secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' },'Chan2Plot' , CH,'Channels' , chans); 

secog_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'},'Chan2Plot' ,CH,'Channels' , chans); 

secog_visualizePSD([],subjNum,'Compare_SeqType_TemporalPattern' ,'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 


secog_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed9'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 

secog_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 4); 

secog_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' },...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 

secog_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 

secog_visualizePSD([],subjNum,'ChunkAligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed9'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 4); 


secog_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 
secog_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 
secog_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'},...
    'Chan2Plot' ,CH,'Channels' , chans , 'bandofInterest' , 6); 
secog_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' ,CH,'Channels' , chans , 'bandofInterest' , 6);



secog_visualizePSD(Pall,subjNum,'raw_SingleTrial' , 'TBNum', [10 1], 'Chan2Plot' , CH,'Channels' , chans)

load('AllData_PSD_StimNorm.mat')
secog_visualizePSD(Pall,subjNum,'binned_SingleTrial' , 'TBNum', [10 4], 'Chan2Plot' , CH)
secog_visualizePSD(Pall,subjNum,'binned_SingleTrial_AvgChann' , 'TBNum', [10 40], 'Chan2Plot' , CH)
secog_visualizePSD(Pall,subjNum,'raw_SingleTrial' , 'TBNum', [1 1], 'Chan2Plot' , CH)
%%
load('AllData_PSD_Warped.mat')
% PMD layer 2
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup_AvgChann' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH)
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH)
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH)
secog_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , blockGroupNames{3},'Chan2Plot' , CH)
%% sig test
secog_sigTest('seqs_across_days' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH);
secog_sigTest('chunks_across_days' , subjNum, 'Groups' , [1 2 3] , 'Chan' , CH) 
secog_sigTest('Alined_Sequences' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
secog_sigTest('Alined_SingleFinger' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
secog_sigTest('ChunkPlacement' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
secog_sigTest('singleFing_across_days' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)

