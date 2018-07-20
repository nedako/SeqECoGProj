%% Make all the data structures
subjNum = 1;

% Good Channels for P2
chans = [4:12 ,14, 122:129 ,36] ;
%%
seqeeg_makeData(subjNum , chans)

%% Visualize Behavior
subjNum = 1;

% Good Channels for P2
chans = [4:12 ,14, 122:129 ,36] ;

cd('/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/P2')
load('AllData_Behav.mat')

load('AllData_Events.mat')

out  = seqeeg_visualize(Dall , subjNum, 'IPI_MT_seq' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'IPI_MT_chunk' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'RT_seq' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'RT_chunk' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'RT_SF' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'MT_SF' ,[1 2]);


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
CH = 19;
{'SingleFingNat' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , 'Intermixed5', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , 'Intermixed8', 'ChunkDay3', 'Intermixed9'}';


seqeeg_visualizePSD(Pall,subjNum,'raw_BlockGroup' , 'BlockGroup' , 'SingleFingSlow1','Chan2Plot' , CH,'Channels' , chans) 
load('AllData_PSD_Warped.mat')
seqeeg_visualizePSD(Pall,subjNum,'raw_BlockGroup' , 'BlockGroup' , 'Intermixed1','Chan2Plot' , CH,'Channels' , chans , 'Rep2Plot' , [1]) 
seqeeg_visualizePSD(Pall,subjNum,'raw_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH,'Channels' , chans, 'Rep2Plot' , [1,2]) 

seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'Intermixed7','Chan2Plot' , CH,'Channels' , chans,'separateReps',1); 
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow1','Chan2Plot' , CH,'Channels' , chans); 

load('AllData_PSD_Warped_SeqType.mat')
seqeeg_visualizePSD(Pall,subjNum,'raw_BlockGroup_SeqType' , 'BlockGroup' , 'Intermixed9','Chan2Plot' ,CH,'Channels' , chans); 
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup_SeqType' , 'BlockGroup' , 'Intermixed9','Chan2Plot' , CH,'Channels' , chans,'separateReps',1); 




seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup_SeqType' , 'BlockGroup' , 'SingleFingFast1','Chan2Plot' , CH,'Channels' , chans);
seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_Repetition' , 'BlockGroup' , {'Intermixed9'},'Chan2Plot' , CH,'Channels' , chans);

seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'},'Chan2Plot' , CH,'Channels' , chans); 




seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},'Chan2Plot' , CH,'Channels' , chans); 

seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'SingleFingSlow1','SingleFingSlow2','SingleFingSlow3','SingleFingSlow4'},'Chan2Plot' , CH,'Channels' , chans); 
seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' },'Chan2Plot' , CH,'Channels' , chans); 

seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed1','Intermixed5','Intermixed8','Intermixed9'},'Chan2Plot' ,CH,'Channels' , chans); 

seqeeg_visualizePSD([],subjNum,'Compare_SeqType_TemporalPattern' ,'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 


seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed9'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 4); 

seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 4); 

seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' },...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 

seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 

seqeeg_visualizePSD([],subjNum,'ChunkAligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed9'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 


seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 
seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 
seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'SingleFingSlow4'},...
    'Chan2Plot' ,CH,'Channels' , chans , 'bandofInterest' , 6); 
seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' ,CH,'Channels' , chans , 'bandofInterest' , 6);
seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' ,CH,'Channels' , chans , 'bandofInterest' , 6);





C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'SingleFingSlow1'});
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4'});
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'SingleFingFast2','SingleFingFast3'});
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'Intermixed3','Intermixed4','Intermixed6','Intermixed7'})
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'Intermixed1'})
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'Intermixed3'})
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'Intermixed7'})


seqeeg_visualizePSD([],subjNum,'Aligned_SeqType_average' ,'Chan2Plot' ,CH,'Channels' , chans );

seqeeg_visualizePSD([],subjNum,'Aligned_singleFing_average' ,'Chan2Plot' ,CH,'Channels' , chans);



seqeeg_visualizePSD([],subjNum,'plot_BlocksPower' , 'Chan2Plot' ,CH,'Channels' , chans );

seqeeg_visualizePSD(Pall,subjNum,'raw_SingleTrial' , 'TBNum', [10 1], 'Chan2Plot' , CH,'Channels' , chans)

load('AllData_PSD_StimNorm.mat')
seqeeg_visualizePSD(Pall,subjNum,'binned_SingleTrial' , 'TBNum', [10 4], 'Chan2Plot' , CH)
seqeeg_visualizePSD(Pall,subjNum,'binned_SingleTrial_AvgChann' , 'TBNum', [10 40], 'Chan2Plot' , CH)
seqeeg_visualizePSD(Pall,subjNum,'raw_SingleTrial' , 'TBNum', [1 1], 'Chan2Plot' , CH)
%%
load('AllData_PSD_Warped.mat')
% PMD layer 2
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup_AvgChann' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH)
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH)
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH)
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , blockGroupNames{3},'Chan2Plot' , CH)
%% sig test
seqeeg_sigTest('seqs_across_days' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH);
seqeeg_sigTest('chunks_across_days' , subjNum, 'Groups' , [1 2 3] , 'Chan' , CH) 
seqeeg_sigTest('Alined_Sequences' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
seqeeg_sigTest('Alined_SingleFinger' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
seqeeg_sigTest('ChunkPlacement' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
seqeeg_sigTest('singleFing_across_days' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
seqeeg_sigTest('SeqFingRepetition' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
