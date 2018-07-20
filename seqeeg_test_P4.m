subjNum = 2;
% Good Channels for P4
chans = [15:21 93:99 101:110 112 28];
%%
seqeeg_makeData(subjNum , chans)


%% Visualize Behavior
subjNum = 2;
% Good Channels for P4
chans = [15:21 93:99 101:110 112 28];

cd('/Volumes/MotorControl/data/SeqECoG/ecog1/iEEG data/P4')
load('AllData_Behav.mat')

load('AllData_Events.mat')
out  = seqeeg_visualize(Dall , subjNum, 'IPI_MT_seq' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'IPI_MT_chunk' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'RT_seq' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'RT_chunk' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'RT_SF' ,[1 2]);
out  = seqeeg_visualize(Dall , subjNum, 'MT_SF' ,[1 2]);


%% plotting case summary
% 'binned_SingleTrial'
% 'binned_SingleTrial_AvgChann'
% 'binned_BlockGroup_AvgChann'
% 'binned_BlockGroup'
% 'raw_SingleTrial'
% 'raw_BlockGroup'
%      'raw_BlockGroup_AvgChann'

% RIGHT SMA --> CH = [1:3];
% RIGHT PMd --> CH = [4:6];
% RIGHT PMv --> CH = [7];
% ALL RIGHT --> CH = [1:7];


% LEFT SMA  --> CH = [8 9];
% LEFT PMd  --> CH = [11:15];
% LEFT PMd+SMA --> CH = [8 9 11:15];
% LEFT M1   --> CH = [16 : 20];
% LEFT M1 HandKnob  --> CH = [18:20];
% LEFT S1   --> CH = [21:25];
% ALL LEFT  --> CH = [8:25];

% ALl PMd   --> CH = [4:6 , 11:15;]
% ALL SMA   --> CH = [1:3 ,8 9];
% ALL MOTOR --> CH = [1:25];

BG(1,:) = {[ ] , [2 8], [14 20 26], [29 38], [], [1 7],[13 19 25], [28 37], [] , [3:5] , [9:11] , [6 12] , [15:17] , [21:23] , [],...
    [18 24] , [30:32] , [34:36] , [], [27 33],[]}';
BG(2,:) = {'[]' , 'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3' ,'[]',...
    'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', '[]' , 'Intermixed1' , 'Intermixed2' , ...
    'ChunkDay1' , 'Intermixed3' , 'Intermixed4' , '[]', 'ChunkDay2' , 'Intermixed6' , ...
    'Intermixed7' , '[]', 'ChunkDay3', '[]'}';
BG = BG';



load('AllData_PSD_Warped.mat')
seqeeg_visualizePSD(Pall,subjNum,'raw_BlockGroup' , 'BlockGroup' , 'Intermixed1','Chan2Plot' , CH,'Channels' , chans , 'Rep2Plot' , [1]) 
seqeeg_visualizePSD(Pall,subjNum,'raw_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH,'Channels' , chans, 'Rep2Plot' , [1,2]) 

seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'Intermixed7','Chan2Plot' , CH,'Channels' , chans,'separateReps',1); 
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow1','Chan2Plot' , CH,'Channels' , chans); 

load('AllData_PSD_Warped_SeqType.mat')
seqeeg_visualizePSD(Pall,subjNum,'raw_BlockGroup_SeqType' , 'BlockGroup' , 'Intermixed7','Chan2Plot' ,CH,'Channels' , chans); 
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup_SeqType' , 'BlockGroup' , 'Intermixed7','Chan2Plot' , CH,'Channels' , chans,'separateReps',1); 

seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup_SeqType' , 'BlockGroup' , 'SingleFingFast1','Chan2Plot' , CH,'Channels' , chans); 
seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_Repetition' , 'BlockGroup' , {'Intermixed1'},'Chan2Plot' ,CH,'Channels' , chans); 
seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_Repetition' , 'BlockGroup' , {'Intermixed5'},'Chan2Plot' ,CH,'Channels' , chans); 
seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_Repetition' , 'BlockGroup' , {'Intermixed7'},'Chan2Plot' ,CH,'Channels' , chans); 

seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'Intermixed1','Intermixed4','Intermixed7'},'Chan2Plot' , CH,'Channels' , chans); 
seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},'Chan2Plot' , CH,'Channels' , chans); 
seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'SingleFingFast1','SingleFingFast2','SingleFingFast3'},'Chan2Plot' , CH,'Channels' , chans); 
seqeeg_visualizePSD([],subjNum,'AvgPower_SeqType_comp' , 'BlockGroup' , {'SingleFingSlow1','SingleFingSlow2','SingleFingSlow3'},'Chan2Plot' , CH,'Channels' , chans) 


seqeeg_visualizePSD([],subjNum,'Compare_SeqType_TemporalPattern' ,'BlockGroup' , {'Intermixed1','Intermixed4','Intermixed7'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6);
seqeeg_visualizePSD([],subjNum,'Compare_SeqType_TemporalPattern' ,'BlockGroup' , {'Intermixed1','Intermixed7'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 4);


seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed4','Intermixed7'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 
seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed7'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 

seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'Intermixed1','SingleFingSlow1'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6);

seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 4);

seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4' },...
    'Chan2Plot' ,CH,'Channels' , chans , 'bandofInterest' , 4);

seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6);

seqeeg_visualizePSD([],subjNum,'Aligned_SeqType' ,'BlockGroup' , {'Intermixed1','Intermixed4','Intermixed7'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6); 


seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'ChunkDay1','ChunkDay2','ChunkDay3'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6);
seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'SingleFingFast1' , 'SingleFingFast2' , 'SingleFingFast3', 'SingleFingFast4'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6);
seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'SingleFingSlow1' , 'SingleFingSlow2'  , 'SingleFingSlow3'},...
    'Chan2Plot' , CH,'Channels' , chans , 'bandofInterest' , 6);
seqeeg_visualizePSD([],subjNum,'AlignedWarped_SeqType' ,'BlockGroup' , {'Intermixed2','Intermixed5','Intermixed8','Intermixed9'},...
    'Chan2Plot' ,CH,'Channels' , chans , 'bandofInterest' , 6);


seqeeg_visualizePSD([],subjNum,'plot_BlocksPower' , 'Chan2Plot' ,CH,'Channels' , chans );


seqeeg_visualizePSD([],subjNum,'Aligned_SeqType_average' ,'Chan2Plot' ,CH,'Channels' , chans );

seqeeg_visualizePSD([],subjNum,'Aligned_singleFing_average' ,'Chan2Plot' ,CH,'Channels' , chans);

seqeeg_visualizePSD(Pall,subjNum,'raw_SingleTrial' , 'TBNum', [10 1], 'Chan2Plot' , CH,'Channels' , chans)

load('AllData_PSD_StimNorm.mat')
seqeeg_visualizePSD(Pall,subjNum,'binned_SingleTrial' , 'TBNum', [10 4], 'Chan2Plot' , CH,'Channels' , chans)
seqeeg_visualizePSD(Pall,subjNum,'binned_SingleTrial_AvgChann' , 'TBNum', [10 40], 'Chan2Plot' , CH)
seqeeg_visualizePSD(Pall,subjNum,'raw_SingleTrial' , 'TBNum',  [10 4], 'Chan2Plot' , CH,'Channels' , chans)

C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'SingleFingSlow1'});
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'SingleFingSlow2','SingleFingSlow3'});
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'SingleFingFast2','SingleFingFast3'});
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'Intermixed3','Intermixed4','Intermixed6','Intermixed7'})
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'Intermixed1'})
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'Intermixed3'})
C = seqeeg_ForcePower('binned_AlignTrials' ,Pall,subjNum, 'Channels' , chans,'BlockGroup',{'Intermixed7'})

%%
load('AllData_PSD_Warped.mat')
% PMD layer 2
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup_AvgChann' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH)
seqeeg_visualizePSD(Pall,subjNum,'binned_BlockGroup' , 'BlockGroup' , 'SingleFingSlow','Chan2Plot' , CH)

%% sig test

seqeeg_sigTest('seqs_across_days' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH);
seqeeg_sigTest('chunks_across_days' , subjNum, 'Groups' , [1 2 3] , 'Chan' , CH) 
seqeeg_sigTest('Alined_Sequences' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
seqeeg_sigTest('Alined_SingleFinger' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
seqeeg_sigTest('ChunkPlacement' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
seqeeg_sigTest('singleFing_across_days' , subjNum, 'Groups' , [1 2 3 4] , 'Chan' , CH)
